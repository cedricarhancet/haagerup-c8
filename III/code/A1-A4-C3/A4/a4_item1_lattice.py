"""a4_item1_lattice.py -- Item 1: lattice facts and generator facts of DJX Section 4 (exact)."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from a4_common import *
from a4_voa import sl2_checks, commute_checks, sugawara, virasoro_checks, mat_det, mat_vec, mat_mul, transpose
from math import comb

T = Timer()
rep = Report("A4 item 1 : lattice E6+A2 (DJX 4.1), glue vector (4.4), Huang (4.6)-(4.8), Jiao (4.11)-(4.14), tables")
L, voa = make_voa()
Gm = D.G

# ---- (4.1) Gram matrix
rep.info("Gram matrix G (4.1), rows", "; ".join(str(r) for r in Gm))
rep.check("G is the Cartan matrix of E6 (chain 1-2-3-4-5, 6 attached to 3) + A2 ({7,8}); off-diagonal -1 exactly at P",
          all(Gm[i][j] == (2 if i == j else (-1 if (i + 1, j + 1) in D.P_EDGES or (j + 1, i + 1) in D.P_EDGES else 0))
              for i in range(8) for j in range(8)))
rep.check("det G = 9", L.det == 9, "det G = %s" % L.det)

# ---- (4.4) glue vector mu
Gmu = mat_vec(Gm, D.MU)
rep.check("G mu = (1,0,0,0,0,0,1,0)^T", tuple(Gmu) == tuple(mpq(v) for v in (1, 0, 0, 0, 0, 0, 1, 0)), "G mu = %s" % (Gmu,))
mu2 = sum(mpq(D.MU[i]) * Gmu[i] for i in range(8))
rep.check("mu^2 = 2", mu2 == 2, "mu^2 = %s" % mu2)
rep.check("3 mu in L", all((3 * mpq(v)).denominator == 1 for v in D.MU), "3mu = %s" % (tuple(3 * mpq(v) for v in D.MU),))
rep.check("mu not in L (so L~ = L + Z mu has index 3)", not all(mpq(v).denominator == 1 for v in D.MU))
# basis of L~ : alpha_1..alpha_7, mu  (det of coordinate matrix = mu_8 = 1/3 -> index 3)
cols = [list(D.unit(i)) for i in range(1, 8)] + [list(D.MU)]
Mt = transpose(cols)
detM = mat_det(Mt)
rep.check("(alpha_1..alpha_7, mu) is a basis of L~ (coordinate determinant = +-1/3)", abs(detM) == mpq(1, 3), "det = %s" % detM)
Gt = mat_mul(mat_mul(transpose(Mt), Gm), Mt)
rep.check("Gram matrix of L~ is integral with even diagonal", all(v.denominator == 1 for r in Gt for v in r) and all(Gt[i][i] % 2 == 0 for i in range(8)),
          "diag = %s" % ([Gt[i][i] for i in range(8)],))
rep.check("det Gram(L~) = 1 : L~ is even unimodular of rank 8, hence isometric to E8", mat_det(Gt) == 1, "det = %s" % mat_det(Gt))
rep.info("(4.5) VE8 = VL + VL+mu + VL-mu : cocycle extension not re-verified here (not needed for Thm 7.3)")

# ---- Section 4.2 norms
ell, h = D.ELL, D.HVEC
rep.check("ell^2 = 2", L.norm(ell) == 2, L.norm(ell))
rep.check("h^2 = 78", L.norm(h) == 78, L.norm(h))
rep.check("(ell, h) = 0", L.ip(ell, h) == 0, L.ip(ell, h))
rep.check("ell = -q (negative of the highest root q of E6)", ell == D.neg(D.Q_VEC))
rep.check("(ell, alpha_i) = 0 for i in I = {1,2,3,4,5,7,8} and (ell, alpha_6) = -1", all(L.ipb(ell, i - 1) == 0 for i in D.I_SET) and L.ipb(ell, 5) == -1)
rep.check("(h, alpha_i) = 2 for i in I", all(L.ipb(h, i - 1) == 2 for i in D.I_SET), "(h,alpha_i) = %s" % ([L.ipb(h, i) for i in range(8)],))
hA5 = (5, 8, 9, 8, 5, 0, 0, 0)
hA2 = (0, 0, 0, 0, 0, 0, 2, 2)
rep.check("A5-chain principal Cartan (5,8,9,8,5,0,0,0) has norm 70 = 2*35 (level 35) and A2 part (2,2) norm 8 = 2*4 (level 4)",
          L.norm(hA5) == 70 and L.norm(hA2) == 8, "%s, %s" % (L.norm(hA5), L.norm(hA2)))

# ---- Section 4.3 norms
p, b, q = D.P_VEC, D.B_VEC, D.Q_VEC
rep.check("p^2 = 6", L.norm(p) == 6, L.norm(p))
rep.check("b^2 = 26", L.norm(b) == 26, L.norm(b))
rep.check("(p, b) = 0", L.ip(p, b) == 0, L.ip(p, b))
rhos = [D.RHO1, D.RHO2, D.RHO3]
rep.check("rho_1, rho_2, rho_3 are roots (norm 2), mutually orthogonal, and p = rho_1 + rho_2 + rho_3",
          all(L.norm(r) == 2 for r in rhos) and all(L.ip(rhos[i], rhos[j]) == 0 for i in range(3) for j in range(i))
          and tuple(sum(r[k] for r in rhos) for k in range(8)) == p)
rep.check("q is a root, q = p + 2 alpha_6, (q, rho_i) = 0, (q, alpha_i) = 0 for i in {1,2,4,5,7,8}",
          L.norm(q) == 2 and q == tuple(p[k] + 2 * D.unit(6)[k] for k in range(8)) and all(L.ip(q, r) == 0 for r in rhos)
          and all(L.ipb(q, i - 1) == 0 for i in (1, 2, 4, 5, 7, 8)))
rep.check("b = 2(alpha_1+alpha_2) + 2(alpha_4+alpha_5) + 2(alpha_7+alpha_8) + q",
          b == tuple(2 * (D.unit(1)[k] + D.unit(2)[k] + D.unit(4)[k] + D.unit(5)[k] + D.unit(7)[k] + D.unit(8)[k]) + q[k] for k in range(8)))
rep.check("(p, alpha_i) = 0 for i in {1,2,4,5,7,8} and (p, q) = 0 (momenta of B_E, B_F orthogonal to p)",
          all(L.ipb(p, i - 1) == 0 for i in (1, 2, 4, 5, 7, 8)) and L.ip(p, q) == 0)
rep.check("(p, rho_i) = 2 for i = 1,2,3", all(L.ip(p, r) == 2 for r in rhos))

# ---- the four sl2 triples
g = huang_generators(voa)
j = jiao_generators(voa)
print(T("generators built"))
allok = True
for label, (E, F, H, k) in {"Huang left (level 1)": (g['EL'], g['FL'], g['HL'], 1),
                            "Huang right (level 39)": (g['ER'], g['FR'], g['HR'], 39),
                            "Jiao P (level 3)": (j['PE'], j['PF'], j['PH'], 3),
                            "Jiao B (level 13)": (j['BE'], j['BF'], j['BH'], 13)}.items():
    res = sl2_checks(voa, E, F, H, k, label + ": ")
    ok = all(r[1] for r in res)
    allok &= ok
    rep.check("%s: E_0F=H, H_0E=2E, H_0F=-2F, E_1F=k1, H_1H=2k1, E_0E=F_0F=E_1E=F_1F=0, H_1E=H_1F=E_1H=0, H_0H=0" % label, ok,
              "" if ok else [r for r in res if not r[1]])
rep.check("HL = (EL)_0 FL and HR = (ER)_0 FR", voa.mode(g['EL'], 0, g['FL']) == g['HL'] and voa.mode(g['ER'], 0, g['FR']) == g['HR'])
rep.check("the two Huang triples commute (a_n b = 0, n = 0,1)", commute_checks(voa, [g['EL'], g['FL'], g['HL']], [g['ER'], g['FR'], g['HR']])[0][1])
rep.check("the two Jiao triples commute (a_n b = 0, n = 0,1)", commute_checks(voa, [j['PE'], j['PF'], j['PH']], [j['BE'], j['BF'], j['BH']])[0][1])
# sub-triples of the right/B triples
subH = [(expsum(voa, [(D.unit(i), 1) for i in (1, 2, 3, 4, 5)]), expsum(voa, [(D.neg(D.unit(i)), -h[i - 1]) for i in (1, 2, 3, 4, 5)]), voa.heis(hA5), 35, "A5 principal"),
        (expsum(voa, [(D.unit(7), 1), (D.unit(8), 1)]), expsum(voa, [(D.neg(D.unit(7)), -2), (D.neg(D.unit(8)), -2)]), voa.heis(hA2), 4, "external A2 principal")]
for (E, F, H, k, lab) in subH:
    res = sl2_checks(voa, E, F, H, k, lab + ": ")
    rep.check("Huang right triple = diagonal of %s (level %d) : sl2 relations" % (lab, k), all(r[1] for r in res))
rep.check("the A5 and A2 sub-triples commute", commute_checks(voa, subH[0][:3], subH[1][:3])[0][1])
subB = [(expsum(voa, [(D.unit(1), 2), (D.unit(2), 1)]), expsum(voa, [(D.neg(D.unit(1)), -1), (D.neg(D.unit(2)), -2)]), voa.heis(tuple(2 * (D.unit(1)[k] + D.unit(2)[k]) for k in range(8))), 4, "A2{1,2}"),
        (expsum(voa, [(D.unit(4), 2), (D.unit(5), 1)]), expsum(voa, [(D.neg(D.unit(4)), -1), (D.neg(D.unit(5)), -2)]), voa.heis(tuple(2 * (D.unit(4)[k] + D.unit(5)[k]) for k in range(8))), 4, "A2{4,5}"),
        (expsum(voa, [(D.unit(7), 2), (D.unit(8), 1)]), expsum(voa, [(D.neg(D.unit(7)), -1), (D.neg(D.unit(8)), -2)]), voa.heis(tuple(2 * (D.unit(7)[k] + D.unit(8)[k]) for k in range(8))), 4, "A2{7,8}"),
        (voa.exp(q), voa.exp(D.neg(q), -1), voa.heis(q), 1, "root q")]
for (E, F, H, k, lab) in subB:
    res = sl2_checks(voa, E, F, H, k, lab + ": ")
    rep.check("Jiao B triple component %s (level %d) : sl2 relations" % (lab, k), all(r[1] for r in res))
okc = True
for a_ in range(4):
    for b_ in range(a_):
        okc &= commute_checks(voa, subB[a_][:3], subB[b_][:3])[0][1]
rep.check("the four components of the B triple commute pairwise (levels 4+4+4+1 = 13)", okc)
subP = [(voa.exp(r), voa.exp(D.neg(r), -1), voa.heis(r), 1, "root rho_%d" % (i + 1)) for i, r in enumerate(rhos)]
okp = all(all(x[1] for x in sl2_checks(voa, E, F, H, k)) for (E, F, H, k, lab) in subP)
for a_ in range(3):
    for b_ in range(a_):
        okp &= commute_checks(voa, subP[a_][:3], subP[b_][:3])[0][1]
rep.check("P triple = diagonal of the three commuting root triples rho_i (levels 1+1+1 = 3)", okp)
print(T("sl2 checks"))

# ---- Sugawara vectors and central charges (4.10)
om = voa.omega()
res, c = virasoro_checks(voa, om, "omega: ")
rep.check("omega (4.3) is a Virasoro vector of central charge 8", all(r[1] for r in res) and c == 8, "c = %s" % c)
wL = sugawara(voa, g['EL'], g['FL'], g['HL'], 1)
wR = sugawara(voa, g['ER'], g['FR'], g['HR'], 39)
wP = sugawara(voa, j['PE'], j['PF'], j['PH'], 3)
wB = sugawara(voa, j['BE'], j['BF'], j['BH'], 13)
for name, w, cexp in (("omega_L", wL, mpq(1)), ("omega_R", wR, mpq(117, 41)), ("omega_P", wP, mpq(9, 5)), ("omega_B^aff", wB, mpq(13, 5))):
    res, c = virasoro_checks(voa, w, name + ": ")
    rep.check("%s Virasoro with c = %s" % (name, cexp), all(r[1] for r in res) and c == cexp, "c = %s" % c)
res, c = virasoro_checks(voa, om - wL - wR)
rep.check("omega - omega_L - omega_R Virasoro with c = 170/41", all(r[1] for r in res) and c == mpq(170, 41), "c = %s" % c)
res, c = virasoro_checks(voa, om - wP - wB)
rep.check("omega - omega_P - omega_B Virasoro with c = 18/5", all(r[1] for r in res) and c == mpq(18, 5), "c = %s" % c)
rep.check("omega_L, omega_R commute with the affine generators of the other triple (a_n w = 0, n=0,1,2)",
          all(voa.mode(a, n, wL).is_zero() for a in (g['ER'], g['FR'], g['HR']) for n in (0, 1, 2)) and
          all(voa.mode(a, n, wR).is_zero() for a in (g['EL'], g['FL'], g['HL']) for n in (0, 1, 2)))
print(T("Virasoro checks"))

# ---- highest weight vectors
u = g['u']
rep.check("a, d, c have norm 4 (so u has weight 2)", all(L.norm(x) == 4 for x in (D.A_VEC, D.D_VEC, D.C_VEC)))
rep.check("u is homogeneous of weight 2", voa.is_homogeneous(u, 2))
rep.check("(EL)_0 u = 0 and (ER)_0 u = 0", voa.mode(g['EL'], 0, u).is_zero() and voa.mode(g['ER'], 0, u).is_zero())
rep.check("(HL)_0 u = 1 u and (HR)_0 u = 13 u  (highest weight (1,13))", voa.mode(g['HL'], 0, u) == u and voa.mode(g['HR'], 0, u) == u.scale(13))
rep.check("(ell,a)=(ell,d)=(ell,c)=1 and (h,a)=(h,d)=(h,c)=13", all(L.ip(ell, x) == 1 and L.ip(h, x) == 13 for x in (D.A_VEC, D.D_VEC, D.C_VEC)))
rep.check("(r+3xi)(r-3xi) = 40", (R + XI * 3) * (R - XI * 3) == K(40))
u26, u17 = j['u26'], j['u17']
rep.check("all 14 lattice vectors of the u_{2,6}, u_{1,7} table have norm 4", all(L.norm(x) == 4 for x, _ in D.U26_TABLE + D.U17_TABLE))
rep.check("the 7 + 7 lattice vectors of the table are distinct within each generator", len(set(x for x, _ in D.U26_TABLE)) == 7 and len(set(x for x, _ in D.U17_TABLE)) == 7)
rep.check("(PE)_0 u26 = 0, (BE)_0 u26 = 0", voa.mode(j['PE'], 0, u26).is_zero() and voa.mode(j['BE'], 0, u26).is_zero())
rep.check("(PH)_0 u26 = 2 u26, (BH)_0 u26 = 6 u26  (highest weight (2,6))", voa.mode(j['PH'], 0, u26) == u26.scale(2) and voa.mode(j['BH'], 0, u26) == u26.scale(6))
rep.check("(PE)_0 u17 = 0, (BE)_0 u17 = 0", voa.mode(j['PE'], 0, u17).is_zero() and voa.mode(j['BE'], 0, u17).is_zero())
rep.check("(PH)_0 u17 = 1 u17, (BH)_0 u17 = 7 u17  (highest weight (1,7))", voa.mode(j['PH'], 0, u17) == u17 and voa.mode(j['BH'], 0, u17) == u17.scale(7))
# which single lattice-vector relations are actually used: (p,x), (b,x) on supports
rep.info("(p,x),(b,x) on supp u26", sorted(set((L.ip(p, x), L.ip(b, x)) for x, _ in D.U26_TABLE)))
rep.info("(p,x),(b,x) on supp u17", sorted(set((L.ip(p, x), L.ip(b, x)) for x, _ in D.U17_TABLE)))
# are u26, u17 the *unique* highest weight vectors of these weights in (VL)_2 up to scalar? (not claimed by DJX; informative)
print(T("highest weight checks"))

# ---- Section 4.4: involutions
aH = D.kappa_H(ell)
rep.check("a_H = kappa_H(ell) = 5^-2 8^-4 9^-3 = 1/74649600", aH == mpq(1, 74649600), aH)
rep.check("8640^2 = 74649600 and 8640 = 5*8^2*27", 8640 ** 2 == 74649600 and 8640 == 5 * 64 * 27)
rep.check("theta_H(EL) = -a_H FL", voa.theta(g['EL'], D.kappa_H) == g['FL'].scale(-aH))
rep.check("theta_H(FL) = -a_H^{-1} EL", voa.theta(g['FL'], D.kappa_H) == g['EL'].scale(-1 / aH))
rep.check("theta_H(ER) = -FR and theta_H(FR) = -ER", voa.theta(g['ER'], D.kappa_H) == -g['FR'] and voa.theta(g['FR'], D.kappa_H) == -g['ER'])
rep.check("theta_J(PE) = -PF, theta_J(PF) = -PE", voa.theta(j['PE'], D.kappa_J) == -j['PF'] and voa.theta(j['PF'], D.kappa_J) == -j['PE'])
rep.check("theta_J(BE) = -BF, theta_J(BF) = -BE", voa.theta(j['BE'], D.kappa_J) == -j['BF'] and voa.theta(j['BF'], D.kappa_J) == -j['BE'])
rep.check("theta_H(HL) = -HL, theta_H(HR) = -HR, theta_J(PH) = -PH, theta_J(BH) = -BH",
          voa.theta(g['HL'], D.kappa_H) == -g['HL'] and voa.theta(g['HR'], D.kappa_H) == -g['HR'] and
          voa.theta(j['PH'], D.kappa_J) == -j['PH'] and voa.theta(j['BH'], D.kappa_J) == -j['BH'])
rep.check("theta_H(omega) = omega = theta_J(omega)", voa.theta(om, D.kappa_H) == om and voa.theta(om, D.kappa_J) == om)
rep.check("kappa_H(a) = 4, kappa_H(d) = kappa_H(c) = 40", D.kappa_H(D.A_VEC) == 4 and D.kappa_H(D.D_VEC) == 40 and D.kappa_H(D.C_VEC) == 40)
thu = voa.theta(u, D.kappa_H)
thu_claim = voa.exp(D.neg(D.A_VEC), 4) + voa.exp(D.neg(D.D_VEC), (R - XI * 3) * 2) + voa.exp(D.neg(D.C_VEC), (R + XI * 3) * 2)
rep.check("theta_H(u) = 4e^{-a} + 2(r-3xi)e^{-d} + 2(r+3xi)e^{-c}", thu == thu_claim)
# positivity of the Hermitian forms on the generators (sanity): norms
rep.info("(u,u)_kappaH", voa.herm(u, u, D.kappa_H))
rep.info("(u26,u26)_kappaJ, (u17,u17)_kappaJ", (voa.herm(u26, u26, D.kappa_J), voa.herm(u17, u17, D.kappa_J)))
rep.info("(omega_L,omega_L)_kappaH = 1/2 ?", voa.herm(wL, wL, D.kappa_H))
rep.info("(omega_R,omega_R)_kappaH = 117/82 ?", voa.herm(wR, wR, D.kappa_H))
rep.check("(omega_L,omega_L) = 1/2 and (omega_R,omega_R) = 117/82 (A.9)", voa.herm(wL, wL, D.kappa_H) == K(mpq(1, 2)) and voa.herm(wR, wR, D.kappa_H) == K(mpq(117, 82)))
rep.check("(omega_c, omega_c)_kappaH = 85/41 with omega_c = omega - omega_L - omega_R", voa.herm(om - wL - wR, om - wL - wR, D.kappa_H) == K(mpq(85, 41)))
# adjointness relations (4.4): (E_0 v, w) = (v, a F_0 w) tests on a few vectors
def adjoint_test(E, F, aconst, kappa, vs):
    ok = True
    for v in vs:
        for w in vs:
            lhs = voa.herm(voa.mode(E, 0, v), w, kappa)
            rhs = voa.herm(v, voa.mode(F, 0, w).scale(aconst), kappa)
            ok &= (lhs == rhs)
    return ok
testH = [u, thu, voa.mode(g['FR'], 0, u), voa.mode(g['FL'], 0, u), wR, voa.mode(g['ER'], 0, thu)]
rep.check("(EL)_0^dagger = a_H (FL)_0 and (ER)_0^dagger = (FR)_0 for (.,.)_kappaH on test vectors",
          adjoint_test(g['EL'], g['FL'], aH, D.kappa_H, testH) and adjoint_test(g['ER'], g['FR'], mpq(1), D.kappa_H, testH))
testJ = [u26, u17, voa.theta(u26, D.kappa_J), voa.mode(j['PF'], 0, u26), voa.mode(j['BF'], 0, u17), wB]
rep.check("(PE)_0^dagger = (PF)_0 and (BE)_0^dagger = (BF)_0 for (.,.)_kappaJ on test vectors",
          adjoint_test(j['PE'], j['PF'], mpq(1), D.kappa_J, testJ) and adjoint_test(j['BE'], j['BF'], mpq(1), D.kappa_J, testJ))
print(T("involution checks"))

# ---- lattice vectors of norm <= 4: HR(0)-weight bound used in (4.18)
roots = set()
frontier = [D.unit(i) for i in range(1, 9)]
roots.update(frontier)
while frontier:
    new = []
    for r_ in frontier:
        for i in range(8):
            cval = L.ipb(r_, i)
            s_ = tuple(r_[k] - (cval if k == i else 0) for k in range(8))   # reflection
            if s_ not in roots:
                roots.add(s_)
                new.append(s_)
    frontier = new
rep.check("number of roots of E6+A2 = 72 + 6 = 78", len(roots) == 78, len(roots))
vecs = {L.zero}
frontier = [L.zero]
while frontier:
    new = []
    for v in frontier:
        for r_ in roots:
            w = L.add(v, r_)
            if w not in vecs and L.norm(w) <= 6:
                vecs.add(w)
                new.append(w)
    frontier = new
n4 = [v for v in vecs if L.norm(v) == 4]
rep.check("number of norm-4 vectors of E6+A2 reached = 702 = 270 (E6) + 0 (A2) + 72*6", len(n4) == 702, len(n4))
maxh = max(abs(L.ip(h, v)) for v in vecs if L.norm(v) <= 4)
maxl = max(abs(L.ip(ell, v)) for v in vecs if L.norm(v) <= 4)
rep.check("max |(h,x)| over x^2 <= 4 is <= 17 (DJX: < sqrt 312) and max |(ell,x)| <= 2", maxh <= 17 and maxl <= 2, "max |(h,x)| = %s, max |(ell,x)| = %s" % (maxh, maxl))
rep.info("dim (VL)_2 = 8 + 36 + 8*78 + 702", 8 + 36 + 8 * 78 + 702)
rep.check("all table vectors (u, u26, u17) are among the enumerated norm-4 vectors", all(x in vecs for x, _ in D.U26_TABLE + D.U17_TABLE) and all(x in vecs for x in (D.A_VEC, D.D_VEC, D.C_VEC)))
print(T("lattice enumeration"))
ok = rep.summary()
print(T("done"))
