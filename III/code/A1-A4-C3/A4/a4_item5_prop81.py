"""a4_item5_prop81.py -- Item 5: DJX Proposition 8.1 (positive conformal-vector identity) verified exactly in V_L:
parafermion conformal vectors (7.1) of Jiao's triples, transport by Phi^{-1} (8.3), Casimir projection Pi_R (4.18),
the recovered Virasoro vectors e, f of (4.19)-(4.20), the identities (8.4), (8.1), (8.2), and the scalar tables (A.11)-(A.12)."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from a4_common import *
from a4_voa import sugawara, eta_k, virasoro_checks, casimir

T = Timer()
rep = Report("A4 item 5 : DJX Proposition 8.1 and (4.19)-(4.22) -- exact Casimir projections in V_L")
L, voa = make_voa()
mats = matrices()
g = huang_generators(voa)
j = jiao_generators(voa)
kH = D.kappa_H
kJ = D.kappa_J
EL, FL, HL, ER, FR, HR, u, om = (g[k] for k in ("EL", "FL", "HL", "ER", "FR", "HR", "u", "omega"))
PE, PF, PH, BE, BF, BH = (j[k] for k in ("PE", "PF", "PH", "BE", "BF", "BH"))


def hH(v, w):
    return voa.herm(v, w, kH)


def hJ(v, w):
    return voa.herm(v, w, kJ)


# ---- parafermion conformal vectors in Jiao's realization
eta3J = eta_k(voa, PE, PF, PH, 3)
eta13J = eta_k(voa, BE, BF, BH, 13)
wP = sugawara(voa, PE, PF, PH, 3)
wB = sugawara(voa, BE, BF, BH, 13)
rep.check("eta_3^J = omega_P - (1/12) (PH)_{-1}PH  ((7.1) agrees with Omega_k - H_{-1}H/(4k))", eta3J == wP - voa.normal(PH, PH).scale(mpq(1, 12)))
rep.check("eta_13^J = omega_B - (1/52) (BH)_{-1}BH", eta13J == wB - voa.normal(BH, BH).scale(mpq(1, 52)))
for name, eta, cexp, H in (("eta_3^J", eta3J, mpq(4, 5), PH), ("eta_13^J", eta13J, mpq(8, 5), BH)):
    res, c = virasoro_checks(voa, eta, name + ": ")
    rep.check("%s is a Virasoro vector of central charge %s" % (name, cexp), all(r[1] for r in res) and c == cexp, "c = %s" % c)
    rep.check("%s is annihilated by H_0, H_1, H_2 of its own Cartan (commutant of the Heisenberg algebra)" % name,
              all(voa.mode(H, n, eta).is_zero() for n in (0, 1, 2)))
    rep.check("%s has zero charge under both Jiao Cartans and is sigma-fixed (momenta in L_0)" % name,
              voa.mode(PH, 0, eta).is_zero() and voa.mode(BH, 0, eta).is_zero() and all(in_L0(x) for x in eta.momenta()))
rep.check("eta_3^J and eta_13^J commute: (eta_3)_n eta_13 = 0 for n = 0..3", all(voa.mode(eta3J, n, eta13J).is_zero() for n in range(4)))
rep.check("the Heisenberg Virasoro (PH)_{-1}PH/12 commutes with eta_3^J (K(sl2,3) x Heisenberg decomposition)", all(voa.mode(voa.normal(PH, PH).scale(mpq(1, 12)), n, eta3J).is_zero() for n in range(4)))
rep.info("number of Fock monomials of eta_3^J, eta_13^J", (len(eta3J), len(eta13J)))
print(T("eta^J"))

# ---- transport (8.3)
eta3 = Phi_inv(voa, mats, eta3J)
eta13 = Phi_inv(voa, mats, eta13J)
rep.check("Phi(Phi^{-1}(eta_k^J)) = eta_k^J (the inverse lift (8.3) is consistent with (6.13))", Phi(voa, mats, eta3) == eta3J and Phi(voa, mats, eta13) == eta13J)
for name, eta, cexp in (("eta_3", eta3, mpq(4, 5)), ("eta_13", eta13, mpq(8, 5))):
    res, c = virasoro_checks(voa, eta, name + ": ")
    rep.check("%s = Phi^{-1}(%s^J) in V_L (Huang copy) is Virasoro with c = %s" % (name, name, cexp), all(r[1] for r in res) and c == cexp, "c = %s" % c)
    rep.check("%s has zero HR(0)- and HL(0)-weight and is tau-fixed (momenta in R_*^{-1}L_0)" % name,
              voa.mode(HR, 0, eta).is_zero() and voa.mode(HL, 0, eta).is_zero() and
              all(sum(mpq(3 * D.ELL[k] + D.HVEC[k]) * L.ipb(x, k) for k in range(8)).__class__ is mpq and
                  (sum(mpq(3 * D.ELL[k] + D.HVEC[k], 4) * L.ipb(x, k) for k in range(8))).denominator == 1 for x in eta.momenta()))
    rep.check("%s is NOT annihilated by (ER)_0 (so the average E_R is a genuine projection)" % name, not voa.mode(ER, 0, eta).is_zero())
rep.check("eta_3 and eta_13 commute in V_L", all(voa.mode(eta3, n, eta13).is_zero() for n in range(4)))
rep.check("theta_H(eta_k) = eta_k  (unitary parafermion conformal vectors, Phi theta_H = theta_J Phi)", voa.theta(eta3, kH) == eta3 and voa.theta(eta13, kH) == eta13)
print(T("transport"))

# ---- Casimir projections (4.18)
def CR(v):
    return casimir(voa, ER, FR, HR, v)


def CL(v):
    return casimir(voa, EL, FL, HL, v)


def PiL(v):
    return v - CL(v).scale(mpq(1, 8))


def Dj_sequence(v):
    """return dict j -> D_j v for j = 0, 2, ..., 16 with D_j = (1 - C_R/(j(j+2))) D_{j-2}"""
    out = {0: v}
    cur = v
    for jj in range(2, 17, 2):
        cur = cur - CR(cur).scale(mpq(1, jj * (jj + 2)))
        out[jj] = cur
    return out


# self-adjointness of the Casimirs (sanity, on the vectors at hand)
rep.check("C_R is self-adjoint for (.,.)_kappaH on {eta_3, eta_13, u_1 theta_H u}: (C v, w) = (v, C w)",
          all(hH(CR(v), w) == hH(v, CR(w)) for v in (eta3, eta13) for w in (eta3, eta13, voa.mode(u, 1, voa.theta(u, kH)))))
# a0 = Pi_L(u_1 theta_H u)
u1thu = voa.mode(u, 1, voa.theta(u, kH))
rep.check("u_1 theta_H(u) has weight 2 and zero Cartan weights for both triples", voa.is_homogeneous(u1thu, 2) and voa.mode(HL, 0, u1thu).is_zero() and voa.mode(HR, 0, u1thu).is_zero())
a0 = PiL(u1thu)
rep.check("a0 = Pi_L(u_1 theta_H u) is annihilated by (EL)_0, (FL)_0 (left invariant)", voa.mode(EL, 0, a0).is_zero() and voa.mode(FL, 0, a0).is_zero())
rep.check("a0 is a 'highest weight' check: C_L a0 = 0", CL(a0).is_zero())
print(T("a0"))
D0 = Dj_sequence(a0)
D3 = Dj_sequence(eta3)
D13 = Dj_sequence(eta13)
print(T("Casimir sequences"))
utilde = D0[16]
Pi3 = D3[16]
Pi13 = D13[16]
rep.check("Pi_R a0, Pi_R eta_3, Pi_R eta_13 are annihilated by (ER)_0, (FR)_0, (HR)_0 (right invariant)",
          all(voa.mode(X_, 0, v).is_zero() for v in (utilde, Pi3, Pi13) for X_ in (ER, FR, HR)))
rep.check("Pi_R is idempotent on these vectors: Pi_R(Pi_R v) = Pi_R v", all(Dj_sequence(v)[16] == v for v in (utilde, Pi3, Pi13)))
# (A.12) scalar moments
A12 = {0: {0: mpq(672), 2: mpq(-1044), 4: mpq(2700), 6: mpq(-4164), 8: mpq(20988, 5), 10: mpq(-7092, 5), 12: mpq(6516, 35), 14: mpq(6516, 35), 16: mpq(6516, 35)},
       3: {0: mpq(2, 5), 2: mpq(-152, 75), 4: mpq(1429, 225), 6: mpq(-2432, 225), 8: mpq(11357, 1125), 10: mpq(-4438, 1125), 12: mpq(524, 7875), 14: mpq(524, 7875), 16: mpq(524, 7875)},
       13: {0: mpq(4, 5), 2: mpq(-958, 325), 4: mpq(2727, 325), 6: mpq(-1836, 325), 8: mpq(2481, 1625), 10: mpq(-54, 1625), 12: mpq(4692, 11375), 14: mpq(4692, 11375), 16: mpq(4692, 11375)}}
for idx, (a_, Dseq) in {0: (a0, D0), 3: (eta3, D3), 13: (eta13, D13)}.items():
    vals = {jj: hH(a_, Dseq[jj]) for jj in range(0, 17, 2)}
    ok = all(vals[jj] == K(A12[idx][jj]) for jj in range(0, 17, 2))
    rep.check("(A.12) column a_%s: (a, D_j a) for j = 0,2,...,16" % idx, ok, {jj: str(vals[jj]) for jj in range(0, 17, 2)} if not ok else "all nine moments agree")
rep.check("D_12 a = D_14 a = D_16 a for a = a0, eta_3, eta_13 (no components of highest weight 14, 16: max HR-weight 16 but weights 14,16 absent)",
          D0[12] == D0[16] and D3[12] == D3[16] and D13[12] == D13[16])
rep.check("(Pi_R a, Pi_R a) = (a, Pi_R a) for a = a0, eta_3, eta_13 (orthogonal projection)", all(hH(Dseq[16], Dseq[16]) == hH(a_, Dseq[16]) for a_, Dseq in ((a0, D0), (eta3, D3), (eta13, D13))))
print(T("(A.12)"))

# ---- (4.19)-(4.20): e, f
wL = sugawara(voa, EL, FL, HL, 1)
wR = sugawara(voa, ER, FR, HR, 39)
wc = om - wL - wR
v0 = utilde - wL.scale(6) - wR.scale(10)
e = (v0.scale(35) - wc.scale(63)).scale(mpq(1, 93))
f = wc - e
rep.check("(A.11): (a0, omega_L) = 3, (a0, omega_R) = 585/41", hH(a0, wL) == K(3) and hH(a0, wR) == K(mpq(585, 41)), (hH(a0, wL), hH(a0, wR)))
rep.check("(A.11): (eta_3, omega_L) = 0, (eta_3, omega_R) = 13/205, (eta_13, omega_L) = 0, (eta_13, omega_R) = 21/205",
          hH(eta3, wL).is_zero() and hH(eta3, wR) == K(mpq(13, 205)) and hH(eta13, wL).is_zero() and hH(eta13, wR) == K(mpq(21, 205)),
          (hH(eta3, wR), hH(eta13, wR)))
rep.check("Prop 4.1 table: (u~, omega_L)/(omega_L,omega_L) = 6 and (u~, omega_R)/(omega_R,omega_R) = 10", hH(utilde, wL) / hH(wL, wL) == K(6) and hH(utilde, wR) / hH(wR, wR) == K(10))
rep.check("(u~, u~) = 6516/35", hH(utilde, utilde) == K(mpq(6516, 35)), hH(utilde, utilde))
for name, w, cexp in (("e", e, mpq(70, 31)), ("f", f, mpq(2400, 1271))):
    res, c = virasoro_checks(voa, w, name + ": ")
    rep.check("(4.22) %s is a Virasoro vector with central charge %s" % (name, cexp), all(r[1] for r in res) and c == cexp, "c = %s" % c)
    rep.check("(4.22) a_n %s = 0 for the six affine generators a and n = 0,1,2" % name, all(voa.mode(a_, n, w).is_zero() for a_ in (EL, FL, HL, ER, FR, HR) for n in (0, 1, 2)))
rep.check("(4.22) e_n f = 0 and f_n e = 0 for n = 0..3", all(voa.mode(e, n, f).is_zero() and voa.mode(f, n, e).is_zero() for n in range(4)))
rep.check("(4.22) e_1 u = (13/31) u and f_1 u = (180/1271) u", voa.mode(e, 1, u) == u.scale(mpq(13, 31)) and voa.mode(f, 1, u) == u.scale(mpq(180, 1271)))
rep.check("(4.22) e + f = omega_c, (e,e) = 35/31, (f,f) = 1200/1271, (e,f) = 0", e + f == wc and hH(e, e) == K(mpq(35, 31)) and hH(f, f) == K(mpq(1200, 1271)) and hH(e, f).is_zero())
rep.check("e_3 e = (35/31) 1 and f_3 f = (1200/1271) 1", voa.mode(e, 3, e) == voa.vacuum().scale(mpq(35, 31)) and voa.mode(f, 3, f) == voa.vacuum().scale(mpq(1200, 1271)))
rep.check("theta_H(e) = e, theta_H(f) = f and e, f are tau-fixed", voa.theta(e, kH) == e and voa.theta(f, kH) == f)
rep.check("omega = omega_L + omega_R + e + f is a decomposition into four mutually commuting Virasoro vectors",
          all(voa.mode(x_, n, y_).is_zero() for x_ in (wL, wR, e, f) for y_ in (wL, wR, e, f) if x_ is not y_ for n in range(4)))
rep.check("(omega_c, e) = (e,e) and (omega_c, f) = (f,f)", hH(wc, e) == hH(e, e) and hH(wc, f) == hH(f, f))
rep.check("the four vectors omega_L, omega_R, e, f are nonzero and (Prop 4.1) u~ = 6 omega_L + 10 omega_R + (156/35) e + (9/5) f",
          utilde == wL.scale(6) + wR.scale(10) + e.scale(mpq(156, 35)) + f.scale(mpq(9, 5)))
print(T("e, f"))

# ---- (8.4), (8.1), (8.2)
T3 = wR.scale(mpq(2, 45)) + e.scale(mpq(39, 175)) + f.scale(mpq(9, 100))
T13 = wR.scale(mpq(14, 195)) + e.scale(mpq(13, 175)) + f.scale(mpq(13, 20))
rep.check("(8.4): Pi_R eta_3 = (2/45) omega_R + (39/175) e + (9/100) f  (exact vector identity in V_L)", Pi3 == T3, "" if Pi3 == T3 else "difference has %d terms, norm %s" % (len(Pi3 - T3), hH(Pi3 - T3, Pi3 - T3)))
rep.check("(8.4): Pi_R eta_13 = (14/195) omega_R + (13/175) e + (13/20) f", Pi13 == T13, "" if Pi13 == T13 else "difference has %d terms" % len(Pi13 - T13))
q0 = eta3.scale(mpq(1, 3)) + eta13.scale(mpq(1, 13))
Piq0 = Pi3.scale(mpq(1, 3)) + Pi13.scale(mpq(1, 13))
rep.check("(8.1): Pi_R q0 = (464/22815) omega_R + (2/25) omega_c", Piq0 == wR.scale(mpq(464, 22815)) + wc.scale(mpq(2, 25)))
rep.check("coefficient arithmetic: (2/45)/3 + (14/195)/13 = 464/22815 ; (39/175)/3 + (13/175)/13 = 2/25 ; (9/100)/3 + (13/20)/13 = 2/25",
          mpq(2, 45) / 3 + mpq(14, 195) / 13 == mpq(464, 22815) and mpq(39, 175) / 3 + mpq(13, 175) / 13 == mpq(2, 25) and mpq(9, 100) / 3 + mpq(13, 20) / 13 == mpq(2, 25))
lhs82 = om
rhs82 = wL + wR.scale(mpq(3403, 4563)) + Pi3.scale(mpq(25, 6)) + Pi13.scale(mpq(25, 26))
rep.check("(8.2): omega = omega_L + (3403/4563) omega_R + (25/6) Pi_R eta_3 + (25/26) Pi_R eta_13  (exact vector identity)", lhs82 == rhs82,
          "" if lhs82 == rhs82 else "difference has %d terms" % len(lhs82 - rhs82))
rep.check("(8.2) coefficients: 3403/4563 + (25/6)(2/45) + (25/26)(14/195) = 1 ; all positive",
          mpq(3403, 4563) + mpq(25, 6) * mpq(2, 45) + mpq(25, 26) * mpq(14, 195) == 1)
rep.check("scalar checks: (Pi_R q0, omega_R) = 232/7995, (Pi_R q0, omega_c) = 34/205, (q0, omega) = 38/195, (omega_R,omega_R) = 117/82, (omega_c,omega_c) = 85/41",
          hH(Piq0, wR) == K(mpq(232, 7995)) and hH(Piq0, wc) == K(mpq(34, 205)) and hH(q0, om) == K(mpq(38, 195)) and hH(wR, wR) == K(mpq(117, 82)) and hH(wc, wc) == K(mpq(85, 41)),
          (hH(Piq0, wR), hH(Piq0, wc), hH(q0, om)))
rep.check("(eta_3, omega) = (1/2)(4/5) = 2/5 and (eta_13, omega) = 4/5 (central charge / 2)", hH(eta3, om) == K(mpq(2, 5)) and hH(eta13, om) == K(mpq(4, 5)))
rep.check("Appendix A.3: c_i^T G_B c_i = c_i^T s_i : (T_3,T_3) = (eta_3, T_3) = 524/7875 and (T_13,T_13) = (eta_13,T_13) = 4692/11375",
          hH(T3, T3) == K(mpq(524, 7875)) and hH(eta3, T3) == K(mpq(524, 7875)) and hH(T13, T13) == K(mpq(4692, 11375)) and hH(eta13, T13) == K(mpq(4692, 11375)))
rep.check("the averaged vectors Pi_R eta_k are not Virasoro vectors (as DJX remark: 'the average itself need not be a Virasoro vector')",
          not all(r[1] for r in virasoro_checks(voa, Pi3)[0]) and not all(r[1] for r in virasoro_checks(voa, Pi13)[0]))
print(T("(8.1), (8.2)"))
rep.summary()
print(T("done"))
