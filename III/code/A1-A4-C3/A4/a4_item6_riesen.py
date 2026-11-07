"""a4_item6_riesen.py -- complement (not in items 1-7 of the task): DJX Section 4.5 and Appendix A.3, Riesen-Huang
identification: lattice map P, sign function F (4.17), transported vectors B_k = Psi(b_k) from Tables 2-3 (A.8),
Gram matrix (A.9), (A.10), e = e_* (4.21) with e from (4.19)-(4.20), pairings (A.11)."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from a4_common import *
from a4_voa import Lattice, LatticeVOA, sugawara, virasoro_checks, casimir, mat_det, mat_mul, mat_vec, transpose

T = Timer()
rep = Report("A4 item 6 (complement) : DJX Section 4.5 / Appendix A.3 -- Riesen's generators transported to Huang's coordinates")
L, voa = make_voa()
g = huang_generators(voa)
kH = D.kappa_H
EL, FL, HL, ER, FR, HR, u, om = (g[k] for k in ("EL", "FL", "HL", "ER", "FR", "HR", "u", "omega"))
hH = lambda v, w: voa.herm(v, w, kH)

# ---- Riesen's lattice: Gram G_Rie (A2 + E6, pairs (1,2),(3,5),(4,6),(5,6),(6,7),(7,8))
GR_ = [[0] * 8 for _ in range(8)]
for i in range(8):
    GR_[i][i] = 2
for (i, k) in [(1, 2), (3, 5), (4, 6), (5, 6), (6, 7), (7, 8)]:
    GR_[i - 1][k - 1] = GR_[k - 1][i - 1] = -1
LR = Lattice(GR_)
rep.check("det G_Rie = 9", LR.det == 9)
# zeta_i and m_i, h_Rie
def dl(i):
    return tuple(1 if k == i - 1 else 0 for k in range(8))
zeta = [dl(1), dl(2), dl(4), dl(6), dl(7), dl(8), tuple(dl(3)[k] + dl(4)[k] + 2 * dl(5)[k] + 2 * dl(6)[k] + dl(7)[k] for k in range(8))]
m = [2, 2, 5, 8, 9, 8, 5]
hRie = tuple(sum(m[i] * zeta[i][k] for i in range(7)) for k in range(8))
rep.check("h_Rie = sum m_i zeta_i = (2,2,5,10,10,18,14,8)", hRie == (2, 2, 5, 10, 10, 18, 14, 8), hRie)
rep.check("all zeta_i are roots and their Gram matrix equals that of (alpha7, alpha8, alpha1, ..., alpha5) (zeta_1..zeta_7 = A2 + A5 simple roots)",
          all(LR.norm(z) == 2 for z in zeta) and all(LR.ip(zeta[i], zeta[k]) == L.ip(D.unit(idx_i), D.unit(idx_k)) for i, idx_i in enumerate((7, 8, 1, 2, 3, 4, 5)) for k, idx_k in enumerate((7, 8, 1, 2, 3, 4, 5))))
rep.check("delta_3^2 = 2, h_Rie^2 = 78, (delta_3, h_Rie) = 0", LR.norm(dl(3)) == 2 and LR.norm(hRie) == 78 and LR.ip(dl(3), hRie) == 0)
# P matrix (table): delta1,2 -> alpha7, alpha8 ; delta3 -> -(1,2,3,2,1,2,0,0) ; delta4 -> alpha1 ; delta5 -> alpha3+alpha4+alpha5+alpha6 ; delta6,7,8 -> alpha2, alpha3, alpha4
Pcols = [D.unit(7), D.unit(8), (-1, -2, -3, -2, -1, -2, 0, 0), D.unit(1), (0, 0, 1, 1, 1, 1, 0, 0), D.unit(2), D.unit(3), D.unit(4)]
Pm = transpose([list(c) for c in Pcols])
rep.check("det P = -1", mat_det(Pm) == -1, mat_det(Pm))
rep.check("P^T G P = G_Rie", mat_mul(mat_mul(transpose(Pm), D.G), Pm) == [[mpq(v) for v in r] for r in GR_])
Pv = lambda x: tuple(int(v) for v in mat_vec(Pm, x))
rep.check("P delta_3 = ell, P h_Rie = h", Pv(dl(3)) == D.ELL and Pv(hRie) == D.HVEC)
rep.check("(P zeta_1..7) = (alpha7, alpha8, alpha1, alpha2, alpha3, alpha4, alpha5)", [Pv(z) for z in zeta] == [D.unit(i) for i in (7, 8, 1, 2, 3, 4, 5)])
rep.check("kappa_H(P zeta_i) = m_i", all(kH(Pv(zeta[i])) == m[i] for i in range(7)))


def F(x):
    x1, x2, x3, x4, x5, x6, x7, x8 = x
    return (x5 + x3 * (x4 + x5 + x6 + x7 + x8) + x5 * (x6 + x7)) % 2


ok = True
for i in range(8):
    for k in range(8):
        x, y = dl(i + 1), dl(k + 1)
        xy = tuple(x[t] + y[t] for t in range(8))
        ok &= LR.eps(x, y) * (-1) ** F(xy) == (-1) ** (F(x) + F(y)) * L.eps(Pv(x), Pv(y))
rep.check("(4.17): eps_Rie(x,y)(-1)^{F(x+y)} = (-1)^{F(x)+F(y)} eps(Px,Py) on all 64 basis pairs (both sides bimultiplicative in (x,y))", ok)
import random
random.seed(3)
ok = True
for _ in range(2000):
    x = tuple(random.randint(-3, 3) for _ in range(8))
    y = tuple(random.randint(-3, 3) for _ in range(8))
    xy = tuple(x[t] + y[t] for t in range(8))
    ok &= LR.eps(x, y) * (-1) ** F(xy) == (-1) ** (F(x) + F(y)) * L.eps(Pv(x), Pv(y))
rep.check("(4.17) on 2000 random pairs", ok)
rep.check("F(+-delta_3) = F(+-zeta_i) = 0 mod 2", all(F(x) == 0 for x in [dl(3), tuple(-v for v in dl(3))] + zeta + [tuple(-v for v in z) for z in zeta]))
beta = [(0, 0, 1, 2, 1, 3, 2, 1), (0, 0, 2, 2, 3, 4, 3, 2), (1, 1, 1, 1, 1, 2, 2, 1)]
rep.check("beta_i have norm 4; P beta_1 = c, P beta_2 = d, P beta_3 = a", all(LR.norm(b_) == 4 for b_ in beta) and Pv(beta[0]) == D.C_VEC and Pv(beta[1]) == D.D_VEC and Pv(beta[2]) == D.A_VEC)
rep.check("((-1)^F(beta_1), (-1)^F(beta_2), (-1)^F(beta_3)) = (-1, 1, 1)", [(-1) ** F(b_) for b_ in beta] == [-1, 1, 1])
# Psi(v) = u: v = ((-r+3xi)/20) e^{beta1} + ((r+3xi)/20) e^{beta2} + e^{beta3}
Psi_v = voa.exp(D.C_VEC, (-R + XI * 3) / 20 * (-1)) + voa.exp(D.D_VEC, (R + XI * 3) / 20) + voa.exp(D.A_VEC, 1)
rep.check("Psi(v) = u (term by term with the signs (-1)^F)", Psi_v == u)
ok = True
for _ in range(2000):
    n = tuple(random.randint(-3, 3) for _ in range(7))
    z = tuple(sum(n[i] * zeta[i][k] for i in range(7)) for k in range(8))
    ok &= F(z) == (n[6] * (n[2] + n[3] + n[4] + n[5])) % 2
rep.check("F(sum n_i zeta_i) = n7 (n3+n4+n5+n6) mod 2 (2000 random n)", ok)
print(T("lattice map P"))

# ---- B_k from Tables 2-3 (A.8)
NU = D.NU
Bk = [Vec() for _ in range(4)]
for pairs, qs in D.TABLE2:
    for (i, k) in pairs:
        for kk in range(4):
            if qs[kk]:
                # nu_i(-1) nu_j(-1) 1 : product of two Heisenberg vectors
                vv = voa.heis_create(NU[i - 1], 1, voa.heis(NU[k - 1]))
                Bk[kk] += vv.scale(mpq(qs[kk]))
for vecs, gs in D.TABLE3:
    for x in vecs:
        Sx = voa.exp(x) + voa.exp(D.neg(x), kH(x))
        for kk in range(4):
            if gs[kk]:
                Bk[kk] += Sx.scale(mpq(gs[kk]))
rep.check("all lattice vectors of Table 3 have norm 4", all(L.norm(x) == 4 for vecs, _ in D.TABLE3 for x in vecs))
rep.check("Table 3 lists 21 distinct lattice vectors (2+2+2+2+4+2+2+1+1+1+1+1)", len(set(x for vecs, _ in D.TABLE3 for x in vecs)) == 21, len(set(x for vecs, _ in D.TABLE3 for x in vecs)))
rep.check("B_0..B_3 are homogeneous of weight 2 and theta_H-fixed", all(voa.is_homogeneous(B_, 2) and voa.theta(B_, kH) == B_ for B_ in Bk))
GB_claim = [[mpq(10), mpq(-2, 3), mpq(8), mpq(-4)], [mpq(-2, 3), mpq(104, 9), mpq(32, 3), mpq(-16, 3)], [mpq(8), mpq(32, 3), mpq(1965), mpq(-2245, 2)], [mpq(-4), mpq(-16, 3), mpq(-2245, 2), mpq(2805, 4)]]
GB = [[hH(Bk[i], Bk[k]) for k in range(4)] for i in range(4)]
rep.check("(A.9): Gram matrix (B_i, B_j) = G_B as printed", all(GB[i][k] == K(GB_claim[i][k]) for i in range(4) for k in range(4)), [[str(v) for v in r] for r in GB])
wL = sugawara(voa, EL, FL, HL, 1)
wR = sugawara(voa, ER, FR, HR, 39)
wc = om - wL - wR
rep.check("(omega_L, B_i) = (omega_R, B_i) = 0", all(hH(wL, B_).is_zero() and hH(wR, B_).is_zero() for B_ in Bk))
GBfull = [[mpq(1, 2), 0, 0, 0, 0, 0], [0, mpq(117, 82), 0, 0, 0, 0]] + [[0, 0] + GB_claim[i] for i in range(4)]
rep.check("det G_B(6x6) = 9586850", mat_det(GBfull) == 9586850, mat_det(GBfull))
rep.check("the B_i lie in the affine commutant: a_n B_i = 0 for the six affine generators, n = 0,1,2",
          all(voa.mode(a_, n, B_).is_zero() for B_ in Bk for a_ in (EL, FL, HL, ER, FR, HR) for n in (0, 1, 2)))
rep.check("(omega_c, B_i) = (2, 8/3, 45, -45/2)", [hH(wc, B_) for B_ in Bk] == [K(2), K(mpq(8, 3)), K(45), K(mpq(-45, 2))], [str(hH(wc, B_)) for B_ in Bk])
rep.check("(A.10): omega_c = (8 B0 + 9 B1 + 2 B2 + 2 B3)/41 (exact vector identity)", wc == (Bk[0].scale(8) + Bk[1].scale(9) + Bk[2].scale(2) + Bk[3].scale(2)).scale(mpq(1, 41)))
print(T("B_k"))

# ---- e_*, f_* and e, f
estar = (Bk[0].scale(2) + Bk[1].scale(6)).scale(mpq(1, 31)) + Bk[2].scale((ONE * 2 - R) / 93) + Bk[3].scale((ONE * 2 - R * 2) / 93)
fstar = wc - estar
def PiL(v):
    return v - casimir(voa, EL, FL, HL, v).scale(mpq(1, 8))
def PiR(v):
    cur = v
    for jj in range(2, 17, 2):
        cur = cur - casimir(voa, ER, FR, HR, cur).scale(mpq(1, jj * (jj + 2)))
    return cur
a0 = PiL(voa.mode(u, 1, voa.theta(u, kH)))
utilde = PiR(a0)
v0 = utilde - wL.scale(6) - wR.scale(10)
e = (v0.scale(35) - wc.scale(63)).scale(mpq(1, 93))
f = wc - e
rep.check("(4.21) / Prop 4.1: e = e_* = (2B0+6B1)/31 + ((2-r)/93) B2 + ((2-2r)/93) B3  (exact vector identity)", e == estar, "" if e == estar else "difference has %d terms" % len(e - estar))
rep.check("f = f_* = omega_c - e_*", f == fstar)
for name, w, cexp in (("e_*", estar, mpq(70, 31)), ("f_*", fstar, mpq(2400, 1271))):
    res, c = virasoro_checks(voa, w, name + ": ")
    rep.check("%s is Virasoro with c = %s (Riesen's two Virasoro vectors e_R, f_R transported)" % (name, cexp), all(r[1] for r in res) and c == cexp, "c = %s" % c)
rep.check("Prop 4.1 table: (e_*,e_*) = 35/31, (f_*,f_*) = 1200/1271, (u~,e_*) = 156/31, (u~,f_*) = 2160/1271, (e_*,f_*) = 0",
          hH(estar, estar) == K(mpq(35, 31)) and hH(fstar, fstar) == K(mpq(1200, 1271)) and hH(utilde, estar) == K(mpq(156, 31)) and hH(utilde, fstar) == K(mpq(2160, 1271)) and hH(estar, fstar).is_zero())
T0 = wL.scale(6) + wR.scale(10) + estar.scale(mpq(156, 35)) + fstar.scale(mpq(9, 5))
rep.check("(T0,T0) = 36/2 + 100*117/82 + (156/35)^2 (35/31) + (9/5)^2 (1200/1271) = 6516/35 and u~ = T0",
          hH(T0, T0) == K(mpq(6516, 35)) and utilde == T0)
# (A.11) pairings with B_k
# (eta_3, eta_13 transported as in item 5)
from a4_voa import eta_k
jj_ = jiao_generators(voa)
mats = matrices()
eta3 = Phi_inv(voa, mats, eta_k(voa, jj_['PE'], jj_['PF'], jj_['PH'], 3))
eta13 = Phi_inv(voa, mats, eta_k(voa, jj_['BE'], jj_['BF'], jj_['BH'], 13))
A11 = {"a0": (a0, [K(mpq(26, 5)), K(mpq(164, 15)), ONE * 136 + R * 8, -ONE * 68 - R * 8]),
       "a3": (eta3, [K(mpq(13, 50)), K(mpq(41, 75)), (ONE * 34 + R * 2) / 5, (-ONE * 17 - R * 2) / 5]),
       "a13": (eta13, [K(mpq(143, 150)), K(mpq(91, 225)), (ONE * 260 - R * 26) / 15, (-ONE * 130 + R * 26) / 15])}
for name, (a_, vals) in A11.items():
    got = [hH(a_, B_) for B_ in Bk]
    rep.check("(A.11) column %s: pairings with B_0..B_3" % name, got == vals, [str(v) for v in got] if got != vals else "agree")
# the two averages in the B-basis (A.3): T_3, T_13 columns
T3 = wR.scale(mpq(2, 45)) + estar.scale(mpq(39, 175)) + fstar.scale(mpq(9, 100))
T13 = wR.scale(mpq(14, 195)) + estar.scale(mpq(13, 175)) + fstar.scale(mpq(13, 20))
rep.check("E_R(eta_3) = T_3 and E_R(eta_13) = T_13 with e_*, f_* from Riesen's tables (Appendix A.3 conclusion)", PiR(eta3) == T3 and PiR(eta13) == T13)
print(T("e_*, f_*"))
rep.summary()
print(T("done"))
