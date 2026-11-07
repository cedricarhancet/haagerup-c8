"""a4_item2_section5.py -- Item 2: DJX Section 5 (lattice neighbor L_H, s, t, basis (5.4), roots (5.5), involutions)."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from a4_common import *
from a4_voa import mat_det, mat_vec, mat_mul, transpose, mat_inverse
import itertools

T = Timer()
rep = Report("A4 item 2 : DJX Section 5 (L_0, L_H, s, t, T, sigma, tau)")
L, voa = make_voa()
Gm = D.G
p, b = D.P_VEC, D.B_VEC
lam, s, t = D.LAMBDA, D.S_VEC, D.T_VEC


def ipq(x, y):
    return sum(mpq(x[i]) * Gm[i][j] * mpq(y[j]) for i in range(8) for j in range(8))


rep.info("lambda = (p+b)/2", lam)
rep.info("s = (p+b)/4", s)
rep.info("t = (13p-3b)/4", t)
rep.check("lambda^2 = 8", ipq(lam, lam) == 8, ipq(lam, lam))
rep.check("s^2 = 2", ipq(s, s) == 2, ipq(s, s))
rep.check("t^2 = 78", ipq(t, t) == 78, ipq(t, t))
rep.check("(s, t) = 0", ipq(s, t) == 0, ipq(s, t))
rep.check("p = (3s+t)/4 and b = (13s-t)/4  (5.2)", tuple((3 * s[k] + t[k]) / 4 for k in range(8)) == tuple(mpq(v) for v in p)
          and tuple((13 * s[k] - t[k]) / 4 for k in range(8)) == tuple(mpq(v) for v in b))
rep.check("(s, p) = 3/2", ipq(s, p) == mpq(3, 2), ipq(s, p))
rep.check("t has HALF-INTEGER alpha^J coordinates: t is NOT in L_J (Table 1 'Here t is a lattice vector' must be read as t in L_H)",
          not all(v.denominator == 1 for v in t), t)
rep.check("t in s + L_0 subset L_H (t - s has integer coordinates with even sum), lambda = 2s in L_0, s not in L_J",
          in_L0(tuple(t[k] - s[k] for k in range(8))) and all(v.denominator == 1 for v in lam) and not all(v.denominator == 1 for v in s),
          "t - s = %s" % (tuple(t[k] - s[k] for k in range(8)),))
rep.check("(s, b) = 13/2, (t, p) = 39/2, (t, b) = -39/2 (consequences)", ipq(s, b) == mpq(13, 2) and ipq(t, p) == mpq(39, 2) and ipq(t, b) == mpq(-39, 2))
Glam = mat_vec(Gm, lam)
rep.check("G lambda = (1,1,-1,1,1,-1,1,1)^T", tuple(Glam) == tuple(mpq(v) for v in (1, 1, -1, 1, 1, -1, 1, 1)), Glam)
rep.check("hence L_0 = {x in L_J : (lambda,x) in 2Z} = integer vectors with even coordinate sum (G lambda = all-ones mod 2)", all(v % 2 == 1 for v in Glam))
rep.check("p has odd coordinate sum (represents the odd coset), L_J = L_0 u (p + L_0)", sum(p) % 2 == 1, sum(p))
rep.check("coordinate sum of lambda = 2s is even: 2s in L_0", sum(lam) % 2 == 0, sum(lam))

# basis B of L_0 (Section 6.1) and B_H of L_H (5.4)
mats = matrices()
Bm = mats["B"]
detB = mat_det(Bm)
rep.check("columns of B = (2e1, e2-e1, e3+e1, e4-e1, e5-e1, 2s, e7-e1, e8-e1) lie in L_0 and det B = +-2 : B is a basis of L_0",
          all(in_L0(col) for col in D.B_BASIS) and abs(detB) == 2, "det B = %s" % detB)
MBH = mats["MBH"]
detBH = mat_det(MBH)
rep.check("det of the coordinate matrix of B_H (5.4) is +-1 (covolume of L_H = covolume of L_J)", abs(detBH) == 1, "det = %s" % detBH)
GBH = mat_mul(mat_mul(transpose(MBH), Gm), MBH)
rep.check("Gram matrix of B_H is integral, even on the diagonal, det 9 : L_H is an even lattice of determinant 9",
          all(v.denominator == 1 for r in GBH for v in r) and all(GBH[i][i] % 2 == 0 for i in range(8)) and mat_det(GBH) == 9,
          "det = %s, diag = %s" % (mat_det(GBH), [GBH[i][i] for i in range(8)]))
rep.check("(s, b_i) in Z for the basis b_i of L_0, i.e. (s, L_0) in Z", all(ipq(s, col).denominator == 1 for col in D.B_BASIS))
rep.check("L_H = L_0 + Zs = L_0 u (s + L_0): 2s in L_0 and s not in L_0 (index 2)", in_L0(lam) and not all(v.denominator == 1 for v in s))

# (5.5): the r_i in B_H coordinates; their alpha^J coordinates, Gram matrix and determinant
Rm = mats["R"]
ri = [tuple(Rm[r][c] for r in range(8)) for c in range(8)]
rep.info("r_i in alpha^J coordinates (columns of R)", "; ".join(str(tuple(int(v) if v.denominator == 1 else v for v in col)) for col in ri))
rep.check("all r_i lie in L_H (= L_0 or s + L_0): r_i - (s if odd) has even integer coordinates",
          all(in_L0(col) or in_L0(tuple(col[k] - s[k] for k in range(8))) for col in ri))
rep.info("which r_i lie in the odd coset s + L_0", [i + 1 for i, col in enumerate(ri) if not in_L0(col)])
GR = mat_mul(mat_mul(transpose(Rm), Gm), Rm)
rep.check("Gram matrix of (r_1..r_8) equals G : R^T G R = G", all(GR[i][j] == Gm[i][j] for i in range(8) for j in range(8)))
detRBH = mat_det(transpose([list(c) for c in D.R_IN_BH]))
rep.check("determinant of the r_i relative to B_H is +-1 : (r_i) is a basis of L_H", abs(detRBH) == 1, "det = %s" % detRBH)
rep.check("det R = +-1 (R maps L (det 9) onto L_H (det 9))", abs(mat_det(Rm)) == 1, mat_det(Rm))
Rl = mat_vec(Rm, D.ELL)
Rh = mat_vec(Rm, D.HVEC)
rep.check("T ell = R ell = s", tuple(Rl) == tuple(mpq(v) for v in s), Rl)
rep.check("T h = R h = t", tuple(Rh) == tuple(mpq(v) for v in t), Rh)
rep.check("hence T is a lattice isometry L -> L_H with T ell = s, T h = t (Proposition 5.1)", True)
# since det R = +-1 and r_i in L_H: R(L) = L_H
Rinv = mats["R_inv"]
rep.check("R^{-1} has integer entries on L_H: R^{-1} b_i in Z^8 for the basis of L_0 and R^{-1} s in Z^8",
          all(all(v.denominator == 1 for v in mat_vec(Rinv, col)) for col in D.B_BASIS) and all(v.denominator == 1 for v in mat_vec(Rinv, s)))
rep.check("R^{-1} p = (3 ell + h)/4 (so tau = exp(2 pi i p(0)) reads exp(pi i (3HL + HR)_0 / 2) on the Huang copy)",
          tuple(mat_vec(Rinv, p)) == tuple((3 * mpq(D.ELL[k]) + D.HVEC[k]) / 4 for k in range(8)), mat_vec(Rinv, p))
rep.check("R^{-1} lambda = 2 ell (lambda = 2s = 2 T ell); note R^{-1} b = (13 ell - h)/4",
          tuple(mat_vec(Rinv, lam)) == tuple(2 * mpq(v) for v in D.ELL) and tuple(mat_vec(Rinv, b)) == tuple((13 * mpq(D.ELL[k]) - D.HVEC[k]) / 4 for k in range(8)))

# (5.6)-(5.7): the involutions on the lattice cosets
rep.check("tau = exp(2 pi i p(0)) is +1 on L_0 : (p, b_i) in Z for all basis vectors of L_0", all(ipq(p, col).denominator == 1 for col in D.B_BASIS))
rep.check("tau is -1 on s + L_0 : (p, s) = 3/2 (half-integer)", ipq(p, s).denominator == 2)
rep.check("sigma = exp(pi i lambda(0)) is +1 on L_0 : (lambda, b_i) even for all basis vectors", all(ipq(lam, col) % 2 == 0 for col in D.B_BASIS))
rep.check("sigma is -1 on p + L_0 : (lambda, p) = 3 odd", ipq(lam, p) == 3, ipq(lam, p))
rep.check("(V_LH)^tau = V_L0 = (V_LJ)^sigma as graded spaces (5.7) : both even parts are the L_0 part", True)

# tau on the Huang generators (pull-back (3ell+h)/4): phases on momenta
pp = tuple((3 * mpq(D.ELL[k]) + D.HVEC[k]) / 4 for k in range(8))
g = huang_generators(voa)


def phase_set(v, vec):
    return sorted(set(ipq(vec, x) for x in v.momenta()))


rep.check("tau fixes HL, HR, omega (momentum 0) and u: ((3ell+h)/4, x) in Z for x in {a, d, c}",
          all(val.denominator == 1 for val in phase_set(g['u'], pp)), phase_set(g['u'], pp))
rep.check("tau negates EL, FL: ((3ell+h)/4, +-ell) = +-3/2", phase_set(g['EL'], pp) == [mpq(3, 2)] and phase_set(g['FL'], pp) == [mpq(-3, 2)])
rep.check("tau negates ER, FR: ((3ell+h)/4, +-alpha_i) = +-1/2", phase_set(g['ER'], pp) == [mpq(1, 2)] and phase_set(g['FR'], pp) == [mpq(-1, 2)])
j = jiao_generators(voa)
rep.check("sigma fixes PH, BH, omega, u26, u17: (lambda, x) even on all momenta of u26, u17",
          all(val % 2 == 0 for val in phase_set(j['u26'], lam) + phase_set(j['u17'], lam)), (phase_set(j['u26'], lam), phase_set(j['u17'], lam)))
rep.check("sigma negates PE, PF, BE, BF: (lambda, x) odd on all their momenta",
          all(val % 2 == 1 for v in (j['PE'], j['PF'], j['BE'], j['BF']) for val in phase_set(v, lam)))
rep.check("sigma = exp(pi i (PH+BH)_0 / 2) : lambda = (p+b)/2", tuple(lam) == tuple((mpq(p[k]) + b[k]) / 2 for k in range(8)))
# Remark 5.2 : the comparison is not the identity embedding; R is not the identity
rep.check("R is not a permutation of the root basis (Remark 5.2: the neighbor comparison is a genuine isometry)", any(abs(Rm[i][j]) > 1 for i in range(8) for j in range(8)))
rep.summary()
print(T("done"))
