# Lumped circuit of a grounded transmon whose junction lands on a pocket P that reaches the ground G only through a
# strap L_s, with a 50 ohm line (centre conductor shorted to P, matched far end) across the strap:
#   I -[L_J || C_qp]- P ,  I -[C_ret]- G ,  P -[L_s || Z0]- G
# Natural frequencies from the nodal admittance determinant, solved numerically (mpmath-free: numpy eig of the
# state-space matrix), not from the closed-form formula. Compare with
#   (a) the finding's formula Gamma = w^4 L_s^2 C_ret^2 / (Z0 C_sigma)
#   (b) the capacitive-only model (pocket at the end of the line, current through the strap = C_qp displacement only)
import numpy as np

def eig_t1(f01_GHz, Csig_fF, Cqp_fF, Ls_pH, Z0=50.0):
    w0 = 2 * np.pi * f01_GHz * 1e9
    Cs, Cqp, Ls = Csig_fF * 1e-15, Cqp_fF * 1e-15, Ls_pH * 1e-12
    Cret = Cs - Cqp
    LJ = 1 / (w0 ** 2 * Cs)
    # state: x = [V_I, V_P, i_J (I->P through L_J), i_s (P->G through L_s)]
    # KCL at I: Cret dV_I/dt + Cqp d(V_I - V_P)/dt + i_J = 0
    # KCL at P: -Cqp d(V_I - V_P)/dt - i_J + i_s + V_P / Z0 = 0
    # L_J di_J/dt = V_I - V_P ;  L_s di_s/dt = V_P
    # Mass matrix form: Mx' = Ax
    M = np.array([[Cret + Cqp, -Cqp, 0, 0],
                  [-Cqp, Cqp, 0, 0],
                  [0, 0, LJ, 0],
                  [0, 0, 0, Ls]])
    A = np.array([[0, 0, -1, 0],
                  [0, -1 / Z0, 1, -1],
                  [1, -1, 0, 0],
                  [0, 1, 0, 0]])
    ev = np.linalg.eigvals(np.linalg.solve(M, A))
    ev = ev[np.abs(ev.imag) > 0]
    r = ev[np.argmin(np.abs(ev.imag - w0))]
    return -1 / (2 * r.real) * 1e6, r.imag / (2 * np.pi * 1e9)

def formula_t1(f01_GHz, Csig_fF, Cqp_fF, Ls_pH, Z0=50.0, C=None):
    w = 2 * np.pi * f01_GHz * 1e9
    Cc = (Csig_fF - Cqp_fF) * 1e-15 if C is None else C * 1e-15
    return 1 / (w ** 4 * (Ls_pH * 1e-12) ** 2 * Cc ** 2 / (Z0 * Csig_fF * 1e-15)) * 1e6

cases = [("AF-17 Q16 (measure cell)", 5.99, 95.90, 4.899, 65.5), ("BF-49 Q10 (measure cell)", 5.99, 96.05, 4.374, 136.4)]
for name, f, Cs, Cqp, T1d in cases:
    print(name, f"f01 {f} GHz, C_sigma {Cs} fF, C_qp {Cqp} fF, design T1 {T1d} us")
    for Ls in [0.3, 0.5, 0.7, 1, 3, 10, 30, 70, 88]:
        t, fr = eig_t1(f, Cs, Cqp, Ls)
        print(f"   L_s {Ls:5.1f} pH: eigen T1 {t:10.3f} us = {t / T1d:8.3f}x design | formula(C_ret) {formula_t1(f, Cs, Cqp, Ls):10.3f} us"
              f" | capacitive-only(C_qp) {formula_t1(f, Cs, Cqp, Ls, C=Cqp):12.1f} us | f01 {fr:.4f} GHz ({(fr / f - 1) * 100:+.3f} %)")
    # L_s at which the 10x budget is just met (eigen model)
    lo, hi = 0.01, 10.0
    for _ in range(60):
        mid = (lo * hi) ** 0.5
        (lo, hi) = (mid, hi) if eig_t1(f, Cs, Cqp, mid)[0] / T1d > 10 else (lo, mid)
    print(f"   10x budget needs L_s < {lo:.3f} pH")
# strap partial self-inductance (thin flat strip): mu0 l/(2 pi) [ln(2l/w) + 1/2 + 0.2235 w/l]
mu0 = 4e-7 * np.pi
for l in [88, 90, 98, 100, 108, 110, 118, 120]:
    w = 10.0
    print(f"strip {l:3d} x {w:.0f} um: partial self-inductance {mu0 * l * 1e-6 / (2 * np.pi) * (np.log(2 * l / w) + 0.5 + 0.2235 * w / l) * 1e12:5.1f} pH")
