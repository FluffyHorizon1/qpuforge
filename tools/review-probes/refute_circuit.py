# Lumped check of the galvanic channel: island -[L_J || C_qp]- pocket -[L_s || Z0 (flux line)]- ground, island -C_ret- ground.
# Energy decay rate = -2 Re(s) of the qubit root of det Y(s) = 0.
import numpy as np
P = np.polynomial.polynomial
def t1_us(f01_GHz, Csig_fF, Cqp_fF, Ls_pH, Z0=50.0):
    w0 = 2*np.pi*f01_GHz*1e9
    Cs, Cqp, Ls = Csig_fF*1e-15, Cqp_fF*1e-15, Ls_pH*1e-12
    Cret = Cs - Cqp
    LJ = 1/(w0**2*Cs)            # junction inductance giving f01 with an ideal ground
    # (1 + s^2 LJ Cqp)(s^2 Cret Ls Z0 + s Ls + Z0) + s^2 LJ Cret (Z0 + s Ls)
    a = [1, 0, LJ*Cqp]
    b = [Z0, Ls, Cret*Ls*Z0]
    c = P.polymul([0, 0, LJ*Cret], [Z0, Ls])
    poly = P.polyadd(P.polymul(a, b), c)
    r = P.polyroots(poly)
    r = r[np.argmin(np.abs(r.imag - w0))]
    return 1/(-2*r.real)*1e6, r.imag/(2*np.pi*1e9)
def formula_us(f01_GHz, Csig_fF, Cqp_fF, Ls_pH, Z0=50.0):
    w = 2*np.pi*f01_GHz*1e9; Cs = Csig_fF*1e-15; Cret = Cs - Cqp_fF*1e-15; Ls = Ls_pH*1e-12
    return 1/(w**4*Ls**2*Cret**2/(Z0*Cs))*1e6
for name, f, Cs, Cqp, T1d in [("AF-17 Q16", 5.99, 95.90, 4.899, 65.5), ("BF-49 Q10", 5.99, 96.05, 4.374, 136.4)]:
    for Ls in [0.3, 0.66, 1, 3, 10, 30, 70]:
        t, fr = t1_us(f, Cs, Cqp, Ls)
        print(f"{name}: L_s {Ls:5.2f} pH: T1 {t:10.3f} us ({t/T1d:7.3f}x design {T1d} us), formula {formula_us(f, Cs, Cqp, Ls):10.3f} us, f01 {fr:.4f} GHz ({(fr/f-1)*100:+.2f} %)")
# strap partial self-inductance (thin flat strip, length l, width w): mu0 l / 2pi [ln(2l/w) + 1/2 + 0.2235 w/l]
mu0 = 4e-7*np.pi
for l in [88, 98, 100, 108, 118]:
    w = 10.0
    print(f"strip {l} x {w} um: partial self-inductance {mu0*l*1e-6/(2*np.pi)*(np.log(2*l/w)+0.5+0.2235*w/l)*1e12:.1f} pH")
