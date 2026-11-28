#!/usr/bin/env python3
"""Order-of-magnitude estimate of what a flip-chip top chip does to a qubit's capacitance, which the extraction in
include/qpuforge/layout.hpp does not model (it puts vacuum above the chip and no metal).

Usage: topchip_estimate.py [--w UM] [--s UM] [--keepout UM] [--gaps UM,UM,...] [--h UM]

A two-dimensional cross-section of one arm of the qubit island: a strip of width w (16 um) between two ground planes
at s (30 um) from its edges, on silicon (eps_r 11.45, semi-infinite below). Above it: (a) vacuum, as the extraction
assumes; (b) a top chip of silicon whose face is at height g (the bump height); (c) the same with a ground plane on the
top chip's face except within `keepout` (40 um) of the strip's edges, as the bump lattice of assemble_chip implies
(the ground bumps keep 40 um from qubits, couplers and buses); (d) a continuous ground plane on the top chip's face.
Finite volumes on a uniform grid of cell size h (1 um) over +-300 um, grounded outer boundary, conjugate gradients.
Prints the capacitance per unit length in each case and its change against (a). The change of an island's E_C is
of the order of minus the change of its arms' capacitance to ground; coupling capacitances across 4-5 um gaps change
less. This is an estimate for one cross-section, not an extraction of the chip with its top chip.
Needs numpy only. Exit status 0, or 2 on a usage error.
"""
import argparse
import sys

import numpy as np

EPS0 = 8.8541878128e-12


def capacitance(w, s, h, eps_si, top, gap, keepout, half=300.0, tol=1e-9, maxit=20000):
    """capacitance per unit length (fF/um) of the strip; top: 'none', 'si', 'window', 'plane'"""
    n = int(round(half / h))
    x = np.arange(-n, n + 1) * h                             # node coordinates
    z = np.arange(-n, n + 1) * h
    nz, nx = len(z), len(x)
    kz0 = n                                                  # z index of the chip's face
    kg = n + int(round(gap / h))                             # z index of the top chip's face
    # cell permittivity: cell (k, i) lies between nodes k, k+1 and i, i+1
    zc = (z[:-1] + z[1:]) / 2
    eps_cell = np.ones((nz - 1, nx - 1))
    eps_cell[zc < 0, :] = eps_si
    if top != 'none':
        eps_cell[zc > gap, :] = eps_si
    # link permittivities (finite volumes): x-links between (k, i) and (k, i+1), z-links between (k, i) and (k+1, i)
    pad = np.pad(eps_cell, 1, mode='edge')
    ex = 0.5 * (pad[:-1, 1:-1] + pad[1:, 1:-1])              # shape (nz, nx-1)
    ez = 0.5 * (pad[1:-1, :-1] + pad[1:-1, 1:])              # shape (nz-1, nx)
    fixed = np.zeros((nz, nx), dtype=bool)
    V0 = np.zeros((nz, nx))
    fixed[0, :] = fixed[-1, :] = fixed[:, 0] = fixed[:, -1] = True
    strip = (np.abs(x) <= w / 2 + 1e-9)
    ground = (np.abs(x) >= w / 2 + s - 1e-9)
    fixed[kz0, strip | ground] = True
    V0[kz0, strip] = 1.0
    if top == 'window':
        fixed[kg, np.abs(x) >= w / 2 + keepout - 1e-9] = True
    elif top == 'plane':
        fixed[kg, :] = True

    def apply(V):                                            # (A V)_node = sum over links eps (V_node - V_neighbour)
        out = np.zeros_like(V)
        dx = V[:, 1:] - V[:, :-1]
        out[:, :-1] -= ex * dx
        out[:, 1:] += ex * dx
        dz = V[1:, :] - V[:-1, :]
        out[:-1, :] -= ez * dz
        out[1:, :] += ez * dz
        return out

    free = ~fixed
    diag = np.zeros((nz, nx))
    diag[:, :-1] += ex
    diag[:, 1:] += ex
    diag[:-1, :] += ez
    diag[1:, :] += ez
    b = -apply(V0)
    b[fixed] = 0
    u = np.zeros((nz, nx))
    r = b.copy()
    zv = np.where(free, r / diag, 0)
    p = zv.copy()
    rz = float(np.sum(r * zv))
    bn = float(np.sqrt(np.sum(b * b)))
    for it in range(maxit):
        Ap = apply(p)
        Ap[fixed] = 0
        a = rz / float(np.sum(p * Ap))
        u += a * p
        r -= a * Ap
        if float(np.sqrt(np.sum(r * r))) / bn < tol:
            break
        zv = np.where(free, r / diag, 0)
        rz2 = float(np.sum(r * zv))
        p = zv + (rz2 / rz) * p
        rz = rz2
    V = u + V0
    Q = float(np.sum(apply(V)[kz0, strip]))                  # eps0 * (charge per unit length) / eps0, in units of eps0 * V
    return Q * EPS0 * 1e15 * 1e-6, it + 1                    # F/m -> fF/um


def main():
    ap = argparse.ArgumentParser(description="flip-chip top chip: change of a qubit arm's capacitance, 2D estimate")
    ap.add_argument('--w', type=float, default=16.0)
    ap.add_argument('--s', type=float, default=30.0)
    ap.add_argument('--keepout', type=float, default=40.0)
    ap.add_argument('--gaps', default='5,10,20')
    ap.add_argument('--h', type=float, default=1.0)
    a = ap.parse_args()
    try:
        gaps = [float(v) for v in a.gaps.split(',')]
    except ValueError:
        print("topchip_estimate.py: --gaps takes numbers in um, e.g. 5,10,20", file=sys.stderr)
        return 2
    if a.h <= 0 or a.h > 2 or a.w <= 0 or a.s <= 0 or a.keepout < 0 or any(g <= a.h for g in gaps):
        print("topchip_estimate.py: sizes must be positive, h <= 2 um and every gap larger than h", file=sys.stderr)
        return 2
    c0, it = capacitance(a.w, a.s, a.h, 11.45, 'none', 10.0, a.keepout)
    print(f"arm cross-section: strip {a.w:g} um, gaps {a.s:g} um, on silicon; cells {a.h:g} um")
    print(f"  (a) vacuum above (the extraction's model): {c0:.5f} fF/um")
    for g in gaps:
        row = []
        for top in ('si', 'window', 'plane'):
            c, _ = capacitance(a.w, a.s, a.h, 11.45, top, g, a.keepout)
            row.append(c)
        print(f"  top chip at {g:g} um: (b) silicon {row[0]:.5f} fF/um ({100 * (row[0] / c0 - 1):+.1f}%), (c) and a ground plane {a.keepout:g} um from the "
              f"strip's edges {row[1]:.5f} ({100 * (row[1] / c0 - 1):+.1f}%), (d) a continuous ground plane {row[2]:.5f} ({100 * (row[2] / c0 - 1):+.1f}%)")
    return 0


if __name__ == '__main__':
    sys.exit(main())
