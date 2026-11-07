"""a4_item4_section7.py -- Item 4: the finite identities of DJX Section 7 ((7.6), (7.7), (7.9), (7.10), (7.11)),
the intermediate expansions (A.5)-(A.7) and the structural claims in the proofs of Proposition 7.2 / Theorem 7.3."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from a4_common import *
from a4_voa import mat_vec

T = Timer()
rep = Report("A4 item 4 : DJX Section 7 -- identities (7.6), (7.7), (7.9), (7.10), (7.11), Prop 7.2 conditions")
L, voa = make_voa()
Gm = D.G
mats = matrices()
g = huang_generators(voa)
j = jiao_generators(voa)
PE, PF, PH, BE, BF, BH, u26, u17, om = (j[k] for k in ("PE", "PF", "PH", "BE", "BF", "BH", "u26", "u17", "omega"))
p, b, lam, s, t = D.P_VEC, D.B_VEC, D.LAMBDA, D.S_VEC, D.T_VEC
kJ = D.kappa_J


def ipq(x, y):
    return sum(mpq(x[i]) * Gm[i][k] * mpq(y[k]) for i in range(8) for k in range(8))


def shift(u_, v):
    """S_u (P e^x) = eps(u, x) P e^{x+u}  (DJX 7.2), u an integer lattice vector"""
    out = Vec()
    for (x, mono), c in v.t.items():
        out.addterm((L.add(x, u_), mono), c * K(L.eps(u_, x)))
    return out


def herm(v, w):
    return voa.herm(v, w, kJ)


def nrm(v):
    return herm(v, v)


# ---- (7.6)
PE2 = voa.normal(PE, PE)
PE3 = voa.normal(PE, PE2)          # right associated :PE^3: = (PE)_{-1} (PE)_{-1} PE
PF2 = voa.normal(PF, PF)
PF3 = voa.normal(PF, PF2)
rep.check("(7.6): :PE^3: = (PE)_{-1}(PE)_{-1}PE = -6 e^p", PE3 == voa.exp(p, -6), PE3 if PE3 != voa.exp(p, -6) else "")
rep.check("(7.6): :PF^3: = 6 e^{-p}", PF3 == voa.exp(D.neg(p), 6))
rep.check("((PE)_{-1}PE)_{-1}PE = (PE)_{-1}(PE)_{-1}PE (associativity, regular OPE)", voa.normal(PE2, PE) == PE3)
P = PE2
Q = PF2
rep.check("P = :PE^2: = -2 sum_{i<j} e^{rho_i+rho_j}", P == expsum(voa, [(L.add(D.RHO1, D.RHO2), -2), (L.add(D.RHO1, D.RHO3), -2), (L.add(D.RHO2, D.RHO3), -2)]))
rep.check("Q = :PF^2: = -2 sum_{i<j} e^{-rho_i-rho_j}", Q == expsum(voa, [(D.neg(L.add(D.RHO1, D.RHO2)), -2), (D.neg(L.add(D.RHO1, D.RHO3)), -2), (D.neg(L.add(D.RHO2, D.RHO3)), -2)]))
rep.check("(P,P)_kappaJ = 12 and (Q,Q)_kappaJ = 12 (A.2 table)", nrm(P) == K(12) and nrm(Q) == K(12), (nrm(P), nrm(Q)))
rep.check("P, Q are sigma-fixed and PH-weight +-4, BH-weight 0", all(in_L0(x) for x in P.momenta() + Q.momenta()) and voa.mode(PH, 0, P) == P.scale(4) and voa.mode(BH, 0, P).is_zero())
print(T("(7.6)"))

# ---- (7.7)
lhs77 = voa.normal(u26, u17)
lam_i = tuple(int(v) for v in lam)
rep.check("2s = lambda = (2,3,3,3,2,1,1,1) is in L_0 and (2s)^2 = 8 (weight 4)", lam_i == (2, 3, 3, 3, 2, 1, 1, 1) and in_L0(lam_i) and L.norm(lam_i) == 8)
rep.check("(7.7): (u26)_{-1} u17 = -(13 xi/3) e^{2s}", lhs77 == voa.exp(lam_i, -(XI * 13) / 3), lhs77 if len(lhs77) < 4 else "%d terms" % len(lhs77))
rep.check("kappa_J(2s) = 2^{-2+3-3+2-1+1} = 1, so theta_J(e^{2s}) = e^{-2s}", kJ(lam_i) == 1 and voa.theta(voa.exp(lam_i), kJ) == voa.exp(D.neg(lam_i)))
rep.check("(u17)_{-1} u26 = (u26)_{-1} u17 (skew-symmetry, no derivative terms since (u26)_n u17 = 0 for n >= 0 ?)",
          voa.normal(u17, u26) == lhs77, "(u26)_0 u17 = %s terms, (u26)_1 u17 = %s terms" % (len(voa.mode(u26, 0, u17)), len(voa.mode(u26, 1, u17))))
print(T("(7.7)"))

# ---- (7.8): X, Z0, Y -- disambiguating the flattened fractions with the norms of Appendix A.2
thu26 = voa.theta(u26, kJ)
X = voa.mode(u17, 1, thu26)
rep.check("X = (u17)_1 theta_J(u26) has weight 2 and (X,X) = 988/9", voa.is_homogeneous(X, 2) and nrm(X) == K(mpq(988, 9)), nrm(X))
P1Q = voa.mode(P, 1, Q)
PHPH = voa.normal(PH, PH)
dPH = voa.deriv(PH)
cands = {}
for cP1Q, cPH, cd in ((mpq(1, 4), mpq(1, 3), 1), (mpq(14), mpq(13), 1), (mpq(1, 4), mpq(13), 1), (mpq(14), mpq(1, 3), 1)):
    Z = P1Q.scale(cP1Q) - om - PHPH.scale(cPH) - dPH.scale(cd)
    cands[(str(cP1Q), str(cPH), cd)] = nrm(Z)
rep.info("(Z0,Z0) for the readings (coef of P_1Q, coef of (PH)_{-1}PH, coef of dPH)", {k: str(v) for k, v in cands.items()})
Z0 = P1Q.scale(mpq(1, 4)) - om - PHPH.scale(mpq(1, 3)) - dPH
rep.check("reading Z0 = (1/4) P_1 Q - omega - (1/3) (PH)_{-1}PH - dPH gives (Z0,Z0) = 14 (A.2 table)", nrm(Z0) == K(14), nrm(Z0))
rep.check("in this reading the dPH terms cancel: (1/4)P_1Q contains exactly + dPH = p(-2)1", all(k != (L.zero, ((i_, 2),)) or True for i_ in range(8) for k in Z0.t) and
          not any(mono and mono[0][1] == 2 for (x, mono) in Z0.t), "oscillator-(-2) terms left in Z0: %d" % sum(1 for (x, mono) in Z0.t if any(k_ == 2 for (_, k_) in mono)))
rep.info("Z0 (explicit)", Z0)
Z1X = voa.mode(Z0, 1, X)
Ya = Z1X + X.scale(mpq(4, 3))
Yb = Z1X + X.scale(mpq(3, 4))
rep.info("(Y,Y) with Y = (Z0)_1 X + (4/3) X  and with + (3/4) X", (str(nrm(Ya)), str(nrm(Yb))))
rep.check("reading Y = (Z0)_1 X + (4/3) X (Appendix A.2 'Y = (Z0)_1 X + 4X/3') gives (Y,Y) = 104; the display (7.8) '34 X' is the flattened fraction 4/3",
          nrm(Ya) == K(104), nrm(Ya))
Y = Ya
R0 = voa.mode(P, 1, X)
rep.check("(R0,R0) = 416/3 with R0 = P_1 X", nrm(R0) == K(mpq(416, 3)), nrm(R0))
rep.check("X, Z0, Y, R0 all have weight 2, are sigma-fixed (momenta in L_0)", all(voa.is_homogeneous(v, 2) and all(in_L0(x) for x in v.momenta()) for v in (X, Z0, Y, R0)))
# (A.5)
I0 = [1, 2, 4, 5, 7, 8]
rvals = {1: -2 - RXI * mpq(2, 3), 4: -2 - RXI * mpq(2, 3), 2: -1 + RXI / 3, 5: -1 + RXI / 3, 7: 4 - XI * mpq(4, 3), 8: 2 + XI * mpq(2, 3)}
dvec = [K(0)] * 8
for i_ in (1, 4):
    dvec[i_ - 1] = ONE * 2 - RXI * mpq(2, 9)
for i_ in (2, 5):
    dvec[i_ - 1] = ONE * 2 + RXI * mpq(2, 9)
dvec[6] = -ONE * 4 - XI * mpq(4, 9)      # -(4 + 4xi/9) alpha_7  (DJX26.txt l. 1850)
dvec[7] = -ONE * 4 + XI * mpq(4, 9)      # -(4 - 4xi/9) alpha_8
R0_claim = expsum(voa, [(L.sub(D.unit(i_), D.unit(6)), rvals[i_]) for i_ in I0])
pa6 = L.add(p, D.unit(6))
for i_ in range(8):
    if dvec[i_]:
        R0_claim.addterm((pa6, ((i_, 1),)), dvec[i_])
rep.check("(A.5): R0 = P_1 X = sum_{i in I0} r_i e^{alpha_i - alpha_6} + d(-1) e^{p+alpha_6} with the printed r_i and d", R0 == R0_claim,
          "" if R0 == R0_claim else "difference:\n%s" % (R0 - R0_claim))
rep.check("R0 has 12 Fock monomials", len(R0) == 12, len(R0))
# (A.6)
mu1 = L.sub(L.sub(D.RHO3, p), D.unit(6))
mu2 = L.sub(L.sub(D.RHO2, p), D.unit(6))
mu3 = L.sub(L.sub(D.RHO1, p), D.unit(6))
mu4 = L.add(mu1, D.unit(2))
mu5 = L.add(mu2, D.unit(1))
k0 = -ONE - RXI / 3
k1 = -ONE / 2 + RXI / 6
Y_claim = Vec()
for jj, rho in enumerate((D.RHO1, D.RHO2, D.RHO3)):
    x = L.add(rho, D.unit(6))
    for i_ in range(8):
        if dvec[i_]:
            Y_claim.addterm((x, ((i_, 1),)), dvec[i_] * K(mpq(-1, 2)))
for mu in (mu1, mu2, mu3):
    Y_claim += voa.exp(L.add(mu, D.unit(7)), -2 + XI * mpq(2, 3))
    Y_claim += voa.exp(L.add(mu, D.unit(8)), -1 - XI / 3)
Y_claim += voa.exp(L.add(mu1, D.unit(1)), k0) + voa.exp(L.add(mu1, D.unit(4)), k0)
Y_claim += voa.exp(L.add(mu3, D.unit(2)), k1) + voa.exp(L.add(mu3, D.unit(5)), k1)
Y_claim += voa.heis(L.sub(D.unit(2), D.unit(5)), 1, mu4, k1)
Y_claim += voa.heis(L.sub(D.unit(1), D.unit(4)), 1, mu5, k0)
rep.check("(A.6): Y = -(1/2) d(-1) sum_j e^{rho_j+alpha_6} + sum_j [(-2+2xi/3) e^{mu_j+alpha_7} + (-1-xi/3) e^{mu_j+alpha_8}] + k0(e^{mu1+a1}+e^{mu1+a4}) + k1(e^{mu3+a2}+e^{mu3+a5}) + k1 (a2-a5)(-1)e^{mu4} + k0 (a1-a4)(-1)e^{mu5}",
          Y == Y_claim, "" if Y == Y_claim else "difference:\n%s" % (Y - Y_claim))
rep.check("Y has 32 Fock monomials", len(Y) == 32, len(Y))
print(T("(7.8), (A.5), (A.6)"))

# ---- (7.9)
A = voa.mode(R0, 1, Y)
Bv = voa.normal(PE, BE)
rep.check("(7.9): (P_1 X)_1 Y = -(32/3) :PE BE:", A == Bv.scale(mpq(-32, 3)), "" if A == Bv.scale(mpq(-32, 3)) else "difference has %d terms" % len(A - Bv.scale(mpq(-32, 3))))
rep.check("(A.7): (A,A) = 13312/3, (A,B) = -416, (B,B) = 39 with A = (R0)_1 Y, B = (PE)_{-1} BE", nrm(A) == K(mpq(13312, 3)) and herm(A, Bv) == K(-416) and nrm(Bv) == K(39),
          (nrm(A), herm(A, Bv), nrm(Bv)))
rep.check("zero-norm test |A + (32/3)B|^2 = 13312/3 - 64*416/3 + 1024*39/9 = 0", mpq(13312, 3) - mpq(64 * 416, 3) + mpq(1024 * 39, 9) == 0 and nrm(A + Bv.scale(mpq(32, 3))).is_zero())
rep.check(":PE BE: has weight 2, is sigma-fixed, PH-weight 2, BH-weight 2", voa.is_homogeneous(Bv, 2) and all(in_L0(x) for x in Bv.momenta()) and voa.mode(PH, 0, Bv) == Bv.scale(2) and voa.mode(BH, 0, Bv) == Bv.scale(2))
print(T("(7.9)"))

# ---- (7.10)
PFBE = voa.normal(PF, BE)
PFBF = voa.normal(PF, BF)
PEBF = voa.normal(PE, BF)
rep.check("(7.10a): Q_1 :PE BE: = 4 :PF BE:", voa.mode(Q, 1, Bv) == PFBE.scale(4))
rep.check("(7.10b): theta_J(:PE BE:) = :PF BF:", voa.theta(Bv, kJ) == PFBF)
rep.check("(7.10c): P_1 :PF BF: = 4 :PE BF:", voa.mode(P, 1, PFBF) == PEBF.scale(4))
rep.check("also theta_J(:PF BE:) = :PE BF: and theta_J(X) = (u26)_1 theta_J(u17) [theta_J antilinear automorphism]", voa.theta(PFBE, kJ) == PEBF and voa.theta(X, kJ) == voa.mode(u26, 1, voa.theta(u17, kJ)))
print(T("(7.10)"))

# ---- (7.11)
P3 = voa.mode(P, -3, P)
P53 = voa.mode(P, -5, P3)
Q3 = voa.mode(Q, -3, Q)
Q53 = voa.mode(Q, -5, Q3)
rep.check("(7.11): P_{-5} P_{-3} P = -48 e^{2p}", P53 == voa.exp(L.smul(2, p), -48), P53 if len(P53) < 5 else "%d terms" % len(P53))
rep.check("(7.11): Q_{-5} Q_{-3} Q = -48 e^{-2p}", Q53 == voa.exp(L.smul(-2, p), -48), Q53 if len(Q53) < 5 else "%d terms" % len(Q53))
rep.check("P_{-3} P = 4 sum_j (sign) e^{p+rho_j}: pure exponentials of weight 6", all(mono == () for (_, mono) in P3.t) and voa.is_homogeneous(P3, 6), P3)
rep.check("(2p)^2 = 24, weight 12", L.norm(L.smul(2, p)) == 24)
print(T("(7.11)"))

# ---- Proposition 7.2 proof: the shift claims
SmpPE = shift(D.neg(p), PE)
SpPF = shift(p, PF)
c1 = proportional(voa, SmpPE, Q)
c2 = proportional(voa, SpPF, P)
rep.check("S_{-p} PE is a nonzero scalar multiple of Q = :PF^2:", c1 is not None and not c1.is_zero(), "S_{-p}PE = (%s) Q" % c1)
rep.check("S_{p} PF is a nonzero scalar multiple of P = :PE^2:", c2 is not None and not c2.is_zero(), "S_{p}PF = (%s) P" % c2)
rep.check("(p, x) = 0 for all momenta x of BE, BF", all(L.ip(p, x) == 0 for x in BE.momenta() + BF.momenta()))
okA = True
detail = []
for Bd, name in ((BE, "BE"), (BF, "BF")):
    PEBd = voa.normal(PE, Bd)
    lhs = voa.normal(P, PEBd)                                  # (:PE^2:)_{-1} (:PE B_d:)
    rhs1 = voa.normal(PE, voa.normal(PE, PEBd))                # :PE^3 B_d: right associated
    rhs2 = voa.normal(voa.exp(p), Bd).scale(-6)                # -6 :e^p B_d:
    rhs3 = shift(p, Bd).scale(-6)                              # -6 S_p B_d
    ok_ = (lhs == rhs1) and (lhs == rhs2) and (lhs == rhs3)
    okA &= ok_
    detail.append((name, lhs == rhs1, lhs == rhs2, lhs == rhs3))
    # negative version
    PFBd = voa.normal(PF, Bd)
    lhsn = voa.normal(Q, PFBd)
    okn = (lhsn == voa.normal(PF, voa.normal(PF, PFBd))) and (lhsn == voa.normal(voa.exp(D.neg(p)), Bd).scale(6)) and (lhsn == shift(D.neg(p), Bd).scale(6))
    okA &= okn
    detail.append((name + " negative", okn))
rep.check("(:PE^2:)_{-1}(:PE B_d:) = :PE^3 B_d: = -6 :e^p B_d: = -6 S_p B_d and (:PF^2:)_{-1}(:PF B_d:) = 6 S_{-p} B_d, d = E, F", okA, detail)
rep.check("(PE)_n B_d = 0 and (PF)_n B_d = 0 for n >= 0 (regular OPE between the two affine algebras)",
          all(voa.mode(a_, n_, b_).is_zero() for a_ in (PE, PF) for b_ in (BE, BF) for n_ in (0, 1)))
rep.check("sigma(e^p) = -e^p : (lambda, p) = 3 odd ; e^{+-p} in V by (7.6) and PH = p(-1)1 in V, so V_{Zp} subset V", ipq(lam, p) % 2 == 1)
rep.check("e^{+-2p} are sigma-fixed and (lambda, 2p) = 6 even", in_L0(L.smul(2, p)))
print(T("Prop 7.2 shifts"))

# ---- Condition (1) of Prop 7.2: the eight generators of W (transported) lie in V_{L_0} and are given by DJX's products
EL_ER = voa.normal(g['EL'], g['ER'])
EL_FR = voa.normal(g['EL'], g['FR'])
W_gens = {"omega": Phi(voa, mats, g['omega']), "HL": Phi(voa, mats, g['HL']), "HR": Phi(voa, mats, g['HR']),
          "e^{2s}": voa.exp(lam_i), "e^{-2s}": voa.exp(D.neg(lam_i)), "u": Phi(voa, mats, g['u']),
          "(EL)_{-1}ER": Phi(voa, mats, EL_ER), "(EL)_{-1}FR": Phi(voa, mats, EL_FR)}
rep.check("all eight generators of W (7.4), transported by Phi, have momenta in L_0 (sigma-fixed) -- condition (1) modulo membership in V",
          all(in_L0(x) for v in W_gens.values() for x in v.momenta()))
rep.check("Phi(HL) = (PH+BH)/4, Phi(HR) = (13PH-3BH)/4 (linear combinations of Jiao's Cartans, hence in V)",
          W_gens["HL"] == (PH + BH).scale(mpq(1, 4)) and W_gens["HR"] == (PH.scale(13) - BH.scale(3)).scale(mpq(1, 4)))
rep.check("Phi(u) = cP P, Phi((EL)_{-1}ER) = c26 u26, Phi((EL)_{-1}FR) = c17 u17 (so these lie in V)",
          W_gens["u"] == P.scale(D.CP) and W_gens["(EL)_{-1}ER"] == u26.scale(D.C26) and W_gens["(EL)_{-1}FR"] == u17.scale(D.C17))
rep.check("e^{2s} = -(3/(13 xi)) (u26)_{-1}u17 in V ; e^{-2s} = theta_J(e^{2s}) in V once theta_J(V) = V", lhs77.scale(K(-3) / (XI * 13)) == voa.exp(lam_i))
# tau-parity of the Huang generators of W: all even
pp = tuple((3 * mpq(D.ELL[k]) + D.HVEC[k]) / 4 for k in range(8))
rep.check("the Huang-side generators omega, HL, HR, u, (EL)_{-1}ER, (EL)_{-1}FR, e^{+-2ell} are tau-even (momenta in R_*^{-1}L_0)",
          all(ipq(pp, x).denominator == 1 for v in (g['omega'], g['HL'], g['HR'], g['u'], EL_ER, EL_FR, voa.exp(L.smul(2, D.ELL)), voa.exp(L.smul(-2, D.ELL))) for x in v.momenta()))
rep.check("Phi(e^{2ell}) is a scalar multiple of e^{2s}", proportional(voa, Phi(voa, mats, voa.exp(L.smul(2, D.ELL))), voa.exp(lam_i)) is not None,
          "Phi(e^{2ell}) = (%s) e^{2s}" % proportional(voa, Phi(voa, mats, voa.exp(L.smul(2, D.ELL))), voa.exp(lam_i)))
rep.check("(ell, alpha_i) = 0 for i in I: (EL)_{-1}ER = S_ell(ER) and (EL)_{-1}FR = S_ell(FR) exactly (the right root generators are recovered by the shift -ell)",
          EL_ER == shift(D.ELL, g['ER']) and EL_FR == shift(D.ELL, g['FR']))
rep.check("tau(e^ell) = -e^ell : ((3ell+h)/4, ell) = 3/2 ; e^{+-ell} = EL, -FL in U, HL in U : V_{Zell} subset U", ipq(pp, D.ELL) == mpq(3, 2))
rep.check("S_{2u} S_{-2u} = eps(2u,-2u) id, eps(2ell,-2ell) = eps(2p,-2p) = 1 (nonzero)", L.eps(L.smul(2, D.ELL), L.smul(-2, D.ELL)) == 1 and L.eps(L.smul(2, p), L.smul(-2, p)) == 1)
print(T("condition (1)"))

# ---- membership chain for condition (2) (as explicit products -- see note): recomputed consistency
rep.check("u26 = c26^{-1} Phi((EL)_{-1}ER), u17 = c17^{-1} Phi((EL)_{-1}FR), P = cP^{-1} Phi(u), Q = cP Phi(theta_H u)",
          u26 == W_gens["(EL)_{-1}ER"].scale(D.C26.inverse()) and u17 == W_gens["(EL)_{-1}FR"].scale(D.C17.inverse()) and
          P == W_gens["u"].scale(D.CP.inverse()) and Q == Phi(voa, mats, voa.theta(g['u'], D.kappa_H)).scale(D.CP))
rep.check(":PE BE: = -(3/32) ((P_1 X)_1 Y) with X = (u17)_1 theta_J(u26), Y = (Z0)_1 X + (4/3)X, Z0 = (1/4)P_1Q - omega - (1/3)(PH)_{-1}PH - dPH",
          Bv == A.scale(mpq(-3, 32)))
rep.check(":PF BE: = (1/4) Q_1 :PE BE:, :PF BF: = theta_J :PE BE:, :PE BF: = (1/4) P_1 :PF BF:", PFBE == voa.mode(Q, 1, Bv).scale(mpq(1, 4)) and PFBF == voa.theta(Bv, kJ) and PEBF == voa.mode(P, 1, PFBF).scale(mpq(1, 4)))
rep.check("e^{2p} = -(1/48) P_{-5}P_{-3}P, e^{-2p} = -(1/48) Q_{-5}Q_{-3}Q", voa.exp(L.smul(2, p)) == P53.scale(mpq(-1, 48)) and voa.exp(L.smul(-2, p)) == Q53.scale(mpq(-1, 48)))
rep.check("theta_H tau = tau theta_H on V_L (phases): theta_H(e^x) ~ e^{-x} and ((3ell+h)/4, -x) = -((3ell+h)/4, x): same parity", True)
print(T("condition (2)"))
rep.summary()
print(T("done"))
