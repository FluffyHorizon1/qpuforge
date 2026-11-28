#!/usr/bin/env python3
"""Independent check of flux_coupling (include/qpuforge/layout.hpp): the same thin superconducting film, solved by a
second implementation.

Usage: flux_check.py <window.json> [--h UM] [--lk-ph PH] [--margin F] [--zone UM] [--tol T] [--maxit N] [--floating]

<window.json> comes from `layout_gates flux <chip> <cells.json> <owner> <line> --dump <window.json>`: the apertures of
the solve (the etch with straps and pad cuts taken out), their kinds (bounded holes and the line's two gaps; dumps of
the earlier model also hold apertures at the window edge), the SQUID loop, the owner's etch, the reference of the angle
that weights M_cm, and flux_coupling's result (by its grid solver, or by its method of moments with --moments).

The sheet current is K = curl(g z) for a stream function g on a uniform grid of square cells of side h (aligned to the
SQUID loop's corner, so that the launcher's edges, at whole and half micrometres from the loop, fall on cell edges). g
is constant in each aperture: one unknown in each bounded hole (zero net flux), zero in any held one, and 1 A apart in
the line's two gaps (with one unknown offset when both float). The film fills the window and a margin around it
(--margin, a fraction of the window's larger side, default 0.1 as flux_coupling's grid), on a periodic grid. The energy
is the magnetic energy of a sheet in free space, (mu0/4) sum |k| |g(k)|^2, evaluated with FFTs, plus with --lk-ph the
kinetic energy (L_k/2) sum |grad g|^2 of a sheet inductance L_k per square (default 0: a perfect thin film, as
flux_coupling assumes). Conjugate gradients with a Jacobi preconditioner minimise it (flux_coupling's grid solver uses
1/|k|; the two share no code). The normal field B_z = (mu0/2) |k| g vanishes in the film (L_k = 0) and threads the
apertures: the flux in the loop gives M; the flux in the owner's etch weighted by (theta / 2 pi - 1/2) gives M_cm,
theta as flux_coupling defines it (the polar angle, about the island's centre, of the ground metal nearest each point,
along 16 rays in 0.5 um steps up to 80 um; found here on a 0.1 um raster), with the branch rule (--zone).
When the dump holds the owner's charge distribution (island_charge, from the control window's extraction), M_cm is
also given by reciprocity with that charge: the displacement current leaving the island at each angle in proportion to
its charge instead of evenly in angle. The tool also prints the line's net flux through the owner's etch (zero for a
bounded hole), which sets how much the split of the return between the two ways round would matter if that etch were
held.

Needs numpy only. Prints both results for M and M_cm and their difference; exit status 0, or 2 on a usage error.
"""
import argparse
import json
import math
import sys
import time

import numpy as np

MU0 = 4e-7 * math.pi


def good_size(n):
    """smallest m >= n whose only prime factors are 2, 3 and 5 (fast FFT sizes)"""
    m = max(2, n)
    while True:
        k = m
        for p in (2, 3, 5):
            while k % p == 0:
                k //= p
        if k == 1:
            return m
        m += 1


def cell_range(lo, hi, origin, h, n):
    """indices of the cells (of side h from `origin`) whose centres lie in [lo, hi)"""
    i0 = int(math.ceil((lo - origin) / h - 0.5 - 1e-9))
    i1 = int(math.ceil((hi - origin) / h - 0.5 - 1e-9))
    return max(0, i0), min(n, i1)


def ring_theta(d, pts):
    """theta of flux_coupling's M_cm weight at the points pts (N x 2): the polar angle about the island's centre of the
    nearest ground metal (16 rays, 0.5 um steps, up to 80 um), less the loop's angle, in [0, 2 pi)"""
    x0, y0, x1, y1 = d['window']
    qx, qy = d['island_centre']
    lo = pts.min(axis=0) - 81.0
    hi = pts.max(axis=0) + 81.0
    lo = np.maximum(lo, [x0, y0])
    hi = np.minimum(hi, [x1, y1])
    r = 0.1
    nx = int(math.ceil((hi[0] - lo[0]) / r)) + 1
    ny = int(math.ceil((hi[1] - lo[1]) / r)) + 1
    ground = np.ones((ny, nx), dtype=bool)                  # ground: in the window, not etch, not the owner's metal
    for e in d['etch'] + d['island']:
        i0, i1 = cell_range(e[0], e[2], lo[0], r, nx)
        j0, j1 = cell_range(e[1], e[3], lo[1], r, ny)
        if i1 > i0 and j1 > j0:
            ground[j0:j1, i0:i1] = False
    def is_ground(x, y):
        i = np.floor((x - lo[0]) / r).astype(np.int64)
        j = np.floor((y - lo[1]) / r).astype(np.int64)
        inside = (x >= x0) & (x <= x1) & (y >= y0) & (y <= y1) & (i >= 0) & (i < nx) & (j >= 0) & (j < ny)
        out = np.zeros(x.shape, dtype=bool)
        out[inside] = ground[j[inside], i[inside]]
        return out
    n = len(pts)
    best = np.full(n, 1e300)
    bx = pts[:, 0].copy()
    by = pts[:, 1].copy()
    for dd in range(16):
        ux, uy = math.cos(2 * math.pi * dd / 16), math.sin(2 * math.pi * dd / 16)
        found = np.zeros(n, dtype=bool)
        t = 0.5
        while t < 80.0:
            todo = (~found) & (t < best)
            if not todo.any():
                break
            idx = np.nonzero(todo)[0]
            hit = is_ground(pts[idx, 0] + t * ux, pts[idx, 1] + t * uy)
            h_idx = idx[hit]
            best[h_idx] = t
            bx[h_idx] = pts[h_idx, 0] + t * ux
            by[h_idx] = pts[h_idx, 1] + t * uy
            found[h_idx] = True
            t += 0.5
    th = np.arctan2(by - qy, bx - qx) - d['theta0']
    th = np.mod(th, 2 * math.pi)
    th[th >= 2 * math.pi] -= 2 * math.pi
    return th


def solve(d, h, lk_ph, margin, tol, maxit, log, floating=False, zone=0.0):
    x0, y0, x1, y1 = d['window']
    m = margin * max(x1 - x0, y1 - y0)
    lp = d['loop']                                           # grid aligned to the loop's corner: the launcher's edges (at whole and half
    gx0 = lp[0] - h * math.ceil((lp[0] - (x0 - m)) / h)      # micrometres from the SQUID's axis) fall on cell edges
    gy0 = lp[1] - h * math.ceil((lp[1] - (y0 - m)) / h)
    nx = good_size(int(math.ceil((x1 + m - gx0) / h)))
    ny = good_size(int(math.ceil((y1 + m - gy0) / h)))
    lab = np.full((ny, nx), -1, dtype=np.int64)              # -1: film; else the aperture of the cell
    for r, k in zip(d['pieces'], d['aperture']):
        i0, i1 = cell_range(r[0], r[2], gx0, h, nx)
        j0, j1 = cell_range(r[1], r[3], gy0, h, ny)
        if i1 > i0 and j1 > j0:
            lab[j0:j1, i0:i1] = int(k)
    kind = [int(k) for k in d['kind']]
    if floating:                                             # every aperture a bounded hole: no aperture held at the window edge
        kind = [0 if k == 1 else k for k in kind]
    ga, gb = int(d['gap_a']), int(d['gap_b'])
    fa, fb = kind[ga] == 2, kind[gb] == 2
    # fixed part g0 (amperes) and the unknowns: every film cell, each floating aperture, the gap pair's offset
    g0 = np.zeros((ny, nx))
    dof_of_ap = {}
    nap = len(kind)
    nd = 0
    present = set(np.unique(lab[lab >= 0]).tolist())         # an aperture with no cell (a sliver the window cuts off) has no unknown
    for a in range(nap):
        if a in (ga, gb):
            continue
        if kind[a] == 0 and a in present:
            dof_of_ap[a] = nd
            nd += 1
    if fa and fb:
        g0[lab == ga] = 0.5
        g0[lab == gb] = -0.5
        dof_of_ap[ga] = dof_of_ap[gb] = nd
        nd += 1
    else:
        g0[lab == ga] = 1.0 if fa else 0.0
        g0[lab == gb] = -1.0 if fb else 0.0
    film = (lab < 0).ravel()
    film_idx = np.nonzero(film)[0]
    nf = len(film_idx)
    labf = lab.ravel()
    ap_cell = np.nonzero(np.isin(labf, list(dof_of_ap.keys())))[0]
    ap_dof = np.array([dof_of_ap[int(a)] for a in labf[ap_cell]], dtype=np.int64)
    hm = h * 1e-6
    kx = 2 * math.pi * np.fft.rfftfreq(nx, d=h)
    ky = 2 * math.pi * np.fft.fftfreq(ny, d=h)
    K = np.sqrt(kx[None, :] ** 2 + ky[:, None] ** 2)         # 1/um
    cmag = hm * hm * (MU0 / 2) * 1e6                         # (J / A^2 per cell pair): h^2 B_z per unit g
    lk = lk_ph * 1e-12

    def bz(g):                                              # normal field in tesla for g in amperes
        return (MU0 / 2) * 1e6 * np.fft.irfft2(K * np.fft.rfft2(g), s=g.shape)

    def op(g):                                              # gradient of the energy: h^2 B_z + L_k (-laplacian g)
        out = cmag * np.fft.irfft2(K * np.fft.rfft2(g), s=g.shape)
        if lk > 0:
            out += lk * (4 * g - np.roll(g, 1, 0) - np.roll(g, -1, 0) - np.roll(g, 1, 1) - np.roll(g, -1, 1))
        return out

    def expand(u):
        g = np.zeros(ny * nx)
        g[film_idx] = u[:nf]
        g[ap_cell] = u[nf + ap_dof]
        return g.reshape(ny, nx)

    def gather(v):
        v = v.ravel()
        out = np.empty(nf + nd)
        out[:nf] = v[film_idx]
        out[nf:] = np.bincount(ap_dof, weights=v[ap_cell], minlength=nd)
        return out

    # Jacobi preconditioner: the operator's diagonal (film cells: its mean over the spectrum; apertures: chi^T A chi)
    kfx = 2 * math.pi * np.fft.fftfreq(nx, d=h)
    kfull_mean = float(np.mean(np.sqrt(kfx[None, :] ** 2 + ky[:, None] ** 2)))
    diag = np.empty(nf + nd)
    diag[:nf] = cmag * kfull_mean + 4 * lk
    for j in range(nd):
        chi = np.zeros(ny * nx)
        chi[ap_cell[ap_dof == j]] = 1.0
        chi = chi.reshape(ny, nx)
        diag[nf + j] = float(np.sum(chi * op(chi)))
    b = -gather(op(g0))
    u = np.zeros(nf + nd)
    r = b.copy()
    z = r / diag
    p = z.copy()
    rz = float(r @ z)
    bnorm = float(np.linalg.norm(b))
    t0 = time.time()
    it = 0
    res = 1.0
    for it in range(1, maxit + 1):
        Ap = gather(op(expand(p)))
        alpha = rz / float(p @ Ap)
        u += alpha * p
        r -= alpha * Ap
        res = float(np.linalg.norm(r)) / bnorm
        if res < tol:
            break
        z = r / diag
        rz_new = float(r @ z)
        p = z + (rz_new / rz) * p
        rz = rz_new
        if log and it % 200 == 0:
            print(f"    iteration {it}: relative residual {res:.2e} ({time.time() - t0:.0f} s)", flush=True)
    g = expand(u) + g0
    B = bz(g)
    xc = gx0 + (np.arange(nx) + 0.5) * h
    yc = gy0 + (np.arange(ny) + 0.5) * h
    li0, li1 = cell_range(lp[0], lp[2], gx0, h, nx)
    lj0, lj1 = cell_range(lp[1], lp[3], gy0, h, ny)
    in_loop = np.zeros((ny, nx), dtype=bool)
    in_loop[lj0:lj1, li0:li1] = True
    phi_loop = float(B[in_loop].sum()) * hm * hm             # Wb for 1 A
    M = abs(phi_loop) * 1e12
    moat = (lab == int(d['own_moat'])) & ~in_loop
    jj, ii = np.nonzero(moat)
    pts = np.stack([xc[ii], yc[jj]], axis=1)
    th = ring_theta(d, pts)
    wgt = th / (2 * math.pi) - 0.5
    branch = None                                            # (branch cell, zone-edge cell whose weight it takes)
    if zone > 0:
        # branches of the owner's etch beyond its cut-out (the readout line's gaps, joined to it at the claw) carry no
        # crossing current: there the qubit's stream function keeps the value it has where the branch leaves the zone
        isl = np.array(d['island'], dtype=float)
        zx0, zy0 = isl[:, 0].min() - zone, isl[:, 1].min() - zone
        zx1, zy1 = isl[:, 2].max() + zone, isl[:, 3].max() + zone
        inz = (pts[:, 0] >= zx0) & (pts[:, 0] <= zx1) & (pts[:, 1] >= zy0) & (pts[:, 1] <= zy1)
        if (~inz).any() and inz.any():
            edge = inz & ((pts[:, 0] - zx0 < 3 * h) | (zx1 - pts[:, 0] < 3 * h) | (pts[:, 1] - zy0 < 3 * h) | (zy1 - pts[:, 1] < 3 * h))
            ep = pts[edge]
            ew = wgt[edge]
            eidx = np.nonzero(edge)[0]
            branch = []
            for k in np.nonzero(~inz)[0]:                    # nearest zone-edge cell of the owner's etch
                j = int(np.argmin((ep[:, 0] - pts[k, 0]) ** 2 + (ep[:, 1] - pts[k, 1]) ** 2))
                wgt[k] = ew[j]
                branch.append((int(k), int(eidx[j])))
    cm = float(np.sum(B[jj, ii] * wgt)) * hm * hm
    M_cm = abs(cm) * 1e12
    sgn = -1 if phi_loop < 0 else 1
    M_cm_signed = cm * 1e12 * sgn
    net = (float(np.sum(B[jj, ii])) * hm * hm + phi_loop) * 1e12 * sgn   # the line's net flux through the owner's etch, loop included (pH)
    recip = None
    if 'island_charge' in d:
        # reciprocity with the island's own charge: the displacement current leaving the island at angle theta (the
        # same angle convention, measured at the island's panels) replaces the even spread of the ansatz; F(theta) is
        # the fraction of the island's charge between the SQUID and theta
        ic = np.array(d['island_charge']['panels'], dtype=float)
        thi = ring_theta(d, ic[:, :2])
        order = np.argsort(thi)
        qs = ic[order, 2]
        cum = np.cumsum(qs) / qs.sum()
        F = np.interp(th, thi[order], cum, left=0.0, right=1.0)
        if branch:
            for k, j in branch:
                F[k] = F[j]
        cmr = float(np.sum(B[jj, ii] * (F - 0.5))) * hm * hm
        quad = [float(qs[(thi[order] >= q * math.pi / 2) & (thi[order] < (q + 1) * math.pi / 2)].sum() / qs.sum()) for q in range(4)]
        recip = dict(M_cm_signed=cmr * 1e12 * sgn, quadrants=quad, ground=d['island_charge']['ground_fF'], other=d['island_charge']['other_conductors_fF'],
                     owner=d['island_charge']['owner_fF'])
    film_bz = float(np.max(np.abs(B.ravel()[film_idx]))) if lk == 0 else float('nan')
    return dict(nx=nx, ny=ny, h=h, unknowns=nf + nd, iterations=it, residual=res, seconds=time.time() - t0, M=M, M_cm=M_cm,
                M_cm_signed=M_cm_signed, film_bz=film_bz, max_bz=float(np.max(np.abs(B))), net=net, recip=recip)


def main():
    ap = argparse.ArgumentParser(description="independent check of flux_coupling: the thin film solved for its stream function on a grid")
    ap.add_argument('window', help="window dump from layout_gates flux ... --dump")
    ap.add_argument('--h', type=float, default=1.0, help="cell size in um (default 1)")
    ap.add_argument('--lk-ph', type=float, default=0.0, help="kinetic sheet inductance in pH per square (default 0)")
    ap.add_argument('--margin', type=float, default=0.1, help="film margin around the window, a fraction of its larger side (default 0.1, as flux_coupling)")
    ap.add_argument('--tol', type=float, default=1e-7, help="relative residual of the conjugate gradients (default 1e-7)")
    ap.add_argument('--maxit', type=int, default=6000)
    ap.add_argument('--floating', action='store_true', help="treat every aperture as a bounded hole (zero net flux), none held at the window edge (for dumps "
                    "of the earlier model, which held the apertures the window clips)")
    ap.add_argument('--zone', type=float, default=40.0, help="the branch rule: weight the owner's etch beyond the island's bounding box grown by this (um) with "
                    "the value where it leaves (default 40, as flux_coupling; 0: no branch rule)")
    a = ap.parse_args()
    if not (0.05 <= a.h <= 5) or a.lk_ph < 0 or not (0 <= a.margin <= 2):
        print("flux_check.py: --h must be in [0.05, 5] um, --lk-ph >= 0, --margin in [0, 2]", file=sys.stderr)
        return 2
    try:
        d = json.load(open(a.window))
    except (OSError, ValueError) as e:
        print(f"flux_check.py: {a.window}: {e}", file=sys.stderr)
        return 2
    mom = d['method_of_moments']
    print(f"{d['chip']}: {d['owner']} <- {d['line']} ({d['cells']}); window {d['window'][2] - d['window'][0]:.0f} x {d['window'][3] - d['window'][1]:.0f} um, "
          f"{len(d['kind'])} apertures")
    method = mom.get('method', 'method of moments')           # the dump's own solve: flux_coupling's grid or its method of moments
    unit = 'cells' if method == 'grid' else 'edge panels'
    print(f"  flux_coupling, {method} ({unit} {mom['mesh_um']:.2f} um, {mom['panels']} {'cells' if method == 'grid' else 'panels'}): M {mom['M_pH']:.4f} pH, "
          f"M_cm {mom['M_cm_pH']:.4f} pH (signed {mom['M_cm_signed_pH']:+.4f})")
    s = solve(d, a.h, a.lk_ph, a.margin, a.tol, a.maxit, True, a.floating, a.zone)
    print(f"  stream function on a {s['nx']} x {s['ny']} grid of {s['h']:g} um cells ({s['unknowns']} unknowns, L_k {a.lk_ph:g} pH/sq): "
          f"M {s['M']:.4f} pH, M_cm {s['M_cm']:.4f} pH (signed {s['M_cm_signed']:+.4f}); {s['iterations']} iterations, residual {s['residual']:.1e}, {s['seconds']:.0f} s")
    if a.lk_ph == 0:
        print(f"    largest |B_z| in the film {s['film_bz']:.3g} T against {s['max_bz']:.3g} T in the apertures (per A)")
    print(f"  difference (this tool - flux_coupling): M {s['M'] - mom['M_pH']:+.4f} pH ({100 * (s['M'] / mom['M_pH'] - 1):+.1f}%), M_cm {s['M_cm_signed'] - mom['M_cm_signed_pH']:+.4f} pH")
    kind = 0 if a.floating else int(d['kind'][int(d['own_moat'])])
    print(f"  the owner's etch is {'held (it reaches the window edge)' if kind == 1 else 'a bounded hole'}; the line's net flux through it (loop included) {s['net']:+.4f} pH: "
          f"moving the split of the return between the two ways round from 1/2 by delta (the loop at the matching share) changes M_cm by {-s['net']:+.4f} delta pH")
    r = s['recip']
    if r:
        q = r['quadrants']
        print(f"  reciprocity with the island's own charge (control window, 1.5 um panels): charge to ground {r['ground']:.1f} fF, to other conductors {r['other']:.1f} fF "
              f"of {r['owner']:.1f} fF; by quarter turns from the SQUID {100 * q[0]:.0f} / {100 * q[1]:.0f} / {100 * q[2]:.0f} / {100 * q[3]:.0f} %")
        print(f"    M_cm with that charge in place of the even spread: {r['M_cm_signed']:+.4f} pH (ansatz {s['M_cm_signed']:+.4f} pH, {100 * (r['M_cm_signed'] / s['M_cm_signed'] - 1):+.1f}%)")
    return 0


if __name__ == '__main__':
    sys.exit(main())
