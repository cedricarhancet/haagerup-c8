#!/usr/bin/env python3
# Exact checks of the arithmetic data of Section 6 of
# "Haagerup symmetry in a rank-six lattice conformal field theory" (v06).
# Requires Python 3 and cypari2 (PARI/GP), used for qfminim, qfisom and matkerint.
# Every assertion is checked in exact rational arithmetic. Run: python3 verify_c6_section6_v01.py
from fractions import Fraction as F
from itertools import product
from collections import defaultdict
import cypari2

pari = cypari2.Pari()
out = []


def report(msg):
    out.append(msg)
    print(msg)


def pmat(G):
    return pari.matrix(len(G), len(G), [x for r in G for x in r])


def ip(x, y, G):
    return sum(F(x[i]) * G[i][j] * F(y[j]) for i in range(len(x)) for j in range(len(y)))


def nroots(G):
    return int(pari.qfminim(pmat(G), 2, 0)[0])


# ---------------------------------------------------------------------------
# 1. The lattice Q = E6 + A2, the overlattices E8^+-, the Cartan columns
GQ = [[0] * 8 for _ in range(8)]
for i in range(8):
    GQ[i][i] = 2
for a, b in [(1, 3), (2, 4), (3, 4), (4, 5), (5, 6), (7, 8)]:
    GQ[a - 1][b - 1] = GQ[b - 1][a - 1] = -1
g = [F(x, 3) for x in (4, 3, 5, 6, 4, 2, 2, 1)]
sg = g[:6] + [g[7], g[6]]
cases = {'I': ([0, 0, 0, 1, 0, 0, 0, 0], [8, 10, 14, 19, 14, 8, 2, 2]),
         'II': ([0, -1, 1, 0, 1, 0, 0, 0], [4, 5, 7, 10, 7, 4, 2, 2])}
GK = [[2, -1, 0, 0, 0, -1], [-1, 2, 0, 0, 0, 0], [0, 0, 2, -1, 0, -1],
      [0, 0, -1, 2, 0, 0], [0, 0, 0, 0, 2, -1], [-1, 0, -1, 0, -1, 4]]
GKf = [[2, -1, -1, 1, -1, -1], [-1, 2, 1, 0, 0, 1], [-1, 1, 2, 0, 0, 0],
       [1, 0, 0, 2, -1, -1], [-1, 0, 0, -1, 2, 0], [-1, 1, 0, -1, 0, 8]]
GNdisp = [[4, -1, -2, -2, -1, -1], [-1, 4, -1, -1, 2, 2], [-2, -1, 4, 1, -1, -1],
          [-2, -1, 1, 4, 1, 1], [-1, 2, -1, 1, 4, 1], [-1, 2, -1, 1, 1, 4]]
BI = [[1, 1, -2, 1, 1, 0], [0, 0, 0, -3, 0, 0], [2, -1, -1, -1, -1, 3], [0, 0, 0, -3, 0, 3],
      [-2, 1, 1, -2, 1, 3], [-1, 2, 2, -1, -1, 0], [-1, -1, -1, 2, -1, 0], [1, -2, 1, 1, 1, -3]]
BII = [[1, 1, -1, 2, -2, 1], [0, 0, 0, 3, 0, -3], [-1, 2, -2, 4, -1, -1], [0, 3, 0, 3, 0, -6],
       [1, 1, 2, 2, 1, -5], [-1, 2, 1, 1, 2, -1], [-1, -1, 1, -2, -1, 2], [1, -2, -1, -1, 1, 1]]


def in_E8(x, gg):
    return any(all((x[i] - t * gg[i]).denominator == 1 for i in range(8)) for t in range(3))


for c, (h1, h2) in cases.items():
    Lb = [h1, [F(a + b, 2) for a, b in zip(h1, h2)]]
    assert all(t.denominator == 1 for t in Lb[1])
    GL = [[ip(x, y, GQ) for y in Lb] for x in Lb]
    B = BI if c == 'I' else BII
    cols = [[F(B[r][k], 3) for r in range(8)] for k in range(6)]
    assert all(in_E8(x, g) for x in cols)
    assert all(ip(x, h1, GQ) == 0 and ip(x, h2, GQ) == 0 for x in cols)
    assert [[ip(x, y, GQ) for y in cols] for x in cols] == [[F(v) for v in r] for r in GK]
    report('case %s: Gram of L = %s, columns of B_%s span a copy of K in E8^+' % (c, GL, c))
assert pari.matdet(pmat(GK)) == 39
GKi = pari.matsolve(pmat(GK), pari.matid(6))
gK = [F(str(GKi[i, 0])) for i in range(6)]
assert [x * 39 for x in gK] == [34, 17, 8, 4, 6, 12]
assert (ip(gK, gK, GK) / 2) % 1 == F(17, 39)
assert nroots(GK) == 14 and nroots(GKf) == 30 and pari.matdet(pmat(GKf)) == 39
th = pari.qfrep(pmat(GK), 3, 1)
assert [2 * int(th[k]) for k in range(3)] == [14, 96, 156]
report('K: det 39, q(g_K) = 17/39, shells 14, 96, 156, K-flat has 30 roots')

# ---------------------------------------------------------------------------
# 2. N = Q cap h^perp, its rootless Gram matrix and its four isotropic lines
for c, (h1, h2) in cases.items():
    rows = [[sum(h[i] * GQ[i][j] for i in range(8)) for j in range(8)] for h in (h1, h2)]
    ker = pari.matkerint(pari.matrix(2, 8, [x for r in rows for x in r]))
    BN = [[int(ker[i, j]) for i in range(8)] for j in range(6)]
    GN = [[int(ip(x, y, GQ)) for y in BN] for x in BN]
    assert pari.matdet(pmat(GN)) == 351 and nroots(GN) == 0
    assert pari.qfisom(pmat(GN), pmat(GNdisp)) != 0
report('N = Q cap h^perp: det 351, no roots, isometric to (eq:rootless-N-gram) in cases I and II')

GN = GNdisp
iso = [v for v in product(range(3), repeat=6) if any(v)
       and all(sum(GN[i][j] * v[j] for j in range(6)) % 3 == 0 for i in range(6))
       and sum(v[i] * GN[i][j] * v[j] for i in range(6) for j in range(6)) % 18 == 0]
lines = sorted({min(v, tuple((-t) % 3 for t in v)) for v in iso})
assert lines == [(1, 0, 2, 0, 1, 2), (1, 0, 2, 0, 2, 1), (1, 1, 2, 1, 1, 2), (1, 1, 2, 1, 2, 1)]
o3 = [v for v in product(range(3), repeat=6) if any(v)
      and all(sum(GN[i][j] * v[j] for j in range(6)) % 3 == 0 for i in range(6))]
assert len(o3) == 26
for v in lines:
    gens = [[3 * int(i == j) for j in range(6)] for i in range(6)] + [list(v)]
    M = pari.matrix(6, 7, [gens[cc][r] for r in range(6) for cc in range(7)])
    H = pari.mathnf(M)
    Bs = [[F(int(H[r, cc]), 3) for r in range(6)] for cc in range(6)]
    G = [[ip(b1, b2, GN) for b2 in Bs] for b1 in Bs]
    assert all(x.denominator == 1 for r in G for x in r)
    G = [[int(x) for x in r] for r in G]
    assert all(G[i][i] % 2 == 0 for i in range(6))
    assert pari.matdet(pmat(G)) == 39 and nroots(G) == 14
    assert pari.qfisom(pmat(G), pmat(GK)) != 0 and pari.qfisom(pmat(G), pmat(GKf)) == 0
report('N\'/N has 26 elements of order 3, exactly four isotropic lines, as displayed')
report('each of the four overlattices is even, of determinant 39, with 14 roots, isometric to K')

# genus argument: the 3-primary form (xy+2z^2)/3 and the discriminant of l^perp/l
q3 = lambda x, y, z: F(x * y + 2 * z * z, 3) % 1
b3 = lambda u, v: (q3(*[(a + b) % 3 for a, b in zip(u, v)]) - q3(*u) - q3(*v)) % 1
vec = [u for u in product(range(3), repeat=3) if any(u)]
iso3 = sorted({min(u, tuple((-t) % 3 for t in u)) for u in vec if q3(*u) == 0})
assert iso3 == [(0, 1, 0), (1, 0, 0), (1, 1, 1), (1, 1, 2)]
for l in iso3:
    perp = [u for u in vec if b3(u, l) == 0]
    span = {tuple((k * t) % 3 for t in l) for k in range(3)}
    vals = {q3(*u) for u in perp if u not in span}
    assert vals == {F(2, 3)}
report('3-primary part: four isotropic lines, each l^perp/l has quadratic value 2/3 = q_K(13)')
w = lambda u: (u[0], u[1], (-u[2]) % 3)
assert {min(w(l), tuple((-t) % 3 for t in w(l))) for l in iso3[2:]} == set(iso3[2:])
assert all(min(w(l), tuple((-t) % 3 for t in w(l))) == l for l in iso3[:2])
assert min(w(iso3[2]), tuple((-t) % 3 for t in w(iso3[2]))) == iso3[3]
report('(x,y,z) -> (x,y,-z) fixes the two pure lines and exchanges the two mixed lines')

# ---------------------------------------------------------------------------
# 3. Pointed labels, minimal energies and the relevant invariant spectrum (case (3,13))
e = lambda m, n: F(m * m, 12) + F(n * n, 52)
lab = lambda m, n: (5 * (m + 24 * ((m + n) // 2))) % 39
qK = lambda p: F(17 * p * p, 39) % 1
best = defaultdict(lambda: (F(99), []))
for m in range(-12, 13):
    for n in range(-40, 41):
        if (m - n) % 2:
            continue
        p, en = lab(m, n), e(m, n)
        assert p == (26 * m + 21 * n) % 39
        assert (qK(p) + en) % 1 == 0
        b0, l0 = best[p]
        if en < b0:
            best[p] = (en, [(m, n)])
        elif en == b0:
            l0.append((m, n))
assert len(best) == 39 and max(v[0] for v in best.values()) == F(16, 13)
high = {p: v for p, v in best.items() if v[0] >= 1}
assert sorted(high) == [4, 12, 17, 22, 27, 35]
assert all(high[p][0] == F(40, 39) for p in (4, 17, 22, 35))
assert high[12][0] == high[27][0] == F(16, 13)
report('39 classes, all minima <= 16/13, classes of minimum >= 1: +-4, +-17 (40/39), +-12 (16/13)')
for p in sorted(high):
    report('  p = %2d, min energy %s, all minimizing momenta %s' % (p, high[p][0], sorted(high[p][1])))
# U_1 weights in case (3,13) and U_2 affine-primary types (ArhI14, Theorem 10.5 and Table 2)
U1 = [(2, 0), (-2, 0), (0, 2), (0, -2)]
prim = [(0, 4), (1, 1), (1, 3), (1, 5), (1, 7), (2, 2), (2, 6), (3, 3)]
desc = [(4, 0), (2, 0), (0, 0), (0, 4), (0, 2), (2, 2), (0, 0), (0, 0), (0, 0)]


def weights(a, b):
    return [(x, y) for x in range(-a, a + 1, 2) for y in range(-b, b + 1, 2)]


def mult_U2(mn):
    return sum(weights(a, b).count(mn) for (a, b) in prim + desc)


low = {}
for p, (en, mins) in best.items():
    if p == 0:
        continue
    N = en + qK(p)
    assert N.denominator == 1 and N in (1, 2)
    for mn in mins:
        d = U1.count(mn) if N == 1 else mult_U2(mn)
        low.setdefault(p, set()).add(d)
    assert len(low[p]) == 1
rel = sorted(p for p in low if low[p] == {1})
assert rel == [3, 4, 13, 17, 22, 26, 35, 36]
assert all(low[p] == {0} for p in low if p not in rel)
report('relevant labels {+-3, +-4, +-13, +-17}, one-dimensional lowest spaces, all minimizers agree')
assert {qK(p) for p in (13, 26)} == {F(2, 3)} and {qK(p) for p in (3, 36)} == {F(12, 13)}
assert {qK(p) for p in (4, 35, 17, 22)} == {F(38, 39)}
units = [u for u in range(39) if (u * u) % 39 == 1]
assert units == [1, 14, 25, 38]
for u in units:
    for S in ({13, 26}, {3, 36}, {4, 35, 17, 22}):
        assert {(u * s) % 39 for s in S} == S
assert (14 * 26 - (-26)) % 39 == 0 and (14 * 21 - 21) % 39 == 0
assert (25 * 26 - 26) % 39 == 0 and (25 * 21 + 21) % 39 == 0
report('isometries of (Z/39,q): u in {1,14,25,38}, they preserve the three relevant sets')

# cubic support among the eight relevant fields
S = [3, 36, 4, 35, 13, 26, 17, 22]
trip = sorted({tuple(sorted(t)) for t in product(S, repeat=3) if sum(t) % 39 == 0})
assert trip == [(4, 13, 22), (13, 13, 13), (17, 26, 35), (26, 26, 26)]
report('nonzero cubic sector triples: (13,13,13), (-13,-13,-13), (4,-17,13), (-4,17,-13)')

# ---------------------------------------------------------------------------
# 4. Normalized constants, stress weights and the I14 refinement of the coset stress
k = 3
assert 2 * k * (k - 1) == 12 and F(12, 9) == F(4, 3) and F(1, k) == F(1, 3)
sug = lambda a, kk: F(a * (a + 2), 4 * (kk + 2))
assert sug(1, 3) + sug(7, 13) == F(6, 5) and sug(2, 3) + sug(6, 13) == F(6, 5)
assert sug(1, 3) - F(1, 12) == F(1, 15) and sug(7, 13) - F(49, 52) == F(7, 65)
assert sug(2, 3) - F(4, 12) == F(1, 15) and sug(6, 13) - F(36, 52) == F(7, 65)
assert F(4, 5) + F(1, 15) + F(7, 65) == F(38, 39)
assert F(9, 5) + F(13, 5) + F(18, 5) == 8 and F(182, 95) + F(32, 19) == F(18, 5)
assert F(56, 95) + F(4, 19) == F(4, 5)
assert F(4, 3) - 2 * F(76, 39) == F(-100, 39) and 1 - F(100, 39) == F(-61, 39)
report('|C| = 2(k-1)/k = 4/3 and 1/k = 1/3 at k = 3, stress weights 4/5 = 56/95 + 4/19')

# ---------------------------------------------------------------------------
# 5. Trivial-twist objects of Z(H3) x P and the dimension bookkeeping
tw = {'0': F(0), 'b': F(0), 'a': F(0), 'c0': F(0), 'c1': F(1, 3), 'c2': F(2, 3)}
for l in range(1, 7):
    tw['d%d' % l] = F(6 * l * l, 13) % 1
triv = sorted((z, p) for z in tw for p in range(39) if (tw[z] + qK(p)) % 1 == 0)
assert triv == sorted([('0', 0), ('b', 0), ('a', 0), ('c0', 0), ('c1', 13), ('c1', 26)])
sols = [(r, s) for r in range(6) for s in range(6) if 1 + r + 2 * s == 6 and 3 * r + 3 * s == 9]
assert sols == [(1, 2)]
report('trivial-twist simple objects: (0,0), (b,0), (a,0), (c0,0), (c1,13), (c1,26), r = 1, s = 2')

# ---------------------------------------------------------------------------
# 6. Weak-field generation table
assert (4 * 10) % 39 == 1 and (17 * 23) % 39 == 1
ell = lambda p, r: min(abs(kk) for kk in range(-19, 20) if (kk * p - r) % 39 == 0)
assert [ell(4, 3), ell(4, 13), ell(4, 17)] == [9, 13, 14]
assert [ell(17, 3), ell(17, 13), ell(17, 4)] == [9, 13, 14]
for p in (4, 17):
    assert len({(kk * p) % 39 for kk in range(-19, 20)}) == 39
    assert len({(kk * p) % 39 for kk in range(-18, 19)}) == 37
report('generation table: 9, 13, 14 for both initial pairs, F^[19] = F, F^[18] proper')
# ---------------------------------------------------------------------------
# 7. The realization of levels (1,39) (Proposition prop:second-spectator-stress)
eI = lambda m, n: F(m * m, 4) + F(n * n, 156)
# second dual generator of [[2,1],[1,20]] is 37 times the first one modulo L = span{(2,0),(1,39)}
# in (m,n) coordinates g1 = (1,-1), g2 = (0,2) and L = span{(2,0),(1,39)}
g1, g2 = (1, -1), (0, 2)
dvec = (g2[0] - 37 * g1[0], g2[1] - 37 * g1[1])
assert dvec == (-37, 39) and (dvec[0] - 1) % 2 == 0 and dvec[1] == 39
labI = lambda m, n: (2 * (m + 37 * ((m + n) // 2))) % 39
bestI = defaultdict(lambda: (F(99), []))
for m in range(-8, 9):
    for n in range(-60, 61):
        if (m - n) % 2:
            continue
        p, en = labI(m, n), eI(m, n)
        assert p == (-2 * n) % 39 and (qK(p) + en) % 1 == 0
        b0, l0 = bestI[p]
        if en < b0:
            bestI[p] = (en, [(m, n)])
        elif en == b0:
            l0.append((m, n))
assert len(bestI) == 39
primI = [(0, 4), (0, 8), (0, 8), (0, 12), (1, 1), (1, 5), (1, 7), (1, 13)]
descI = [(2, 0), (0, 2), (0, 4), (2, 2)]
U1I = [(2, 0), (-2, 0), (0, 2), (0, -2)]
multI = lambda mn: sum(weights(a, b).count(mn) for (a, b) in primI + descI)
relI = []
for p, (en, mins) in bestI.items():
    if p == 0:
        continue
    N = en + qK(p)
    if N == 1:
        d = {U1I.count(x) for x in mins}
    elif N == 2:
        d = {multI(x) for x in mins}
    else:
        continue
    assert len(d) == 1
    if d == {1}:
        relI.append(p)
assert sorted(relI) == [4, 13, 17, 22, 26, 35]
assert bestI[35] == (F(1, 39), [(0, 2)]) and bestI[4] == (F(1, 39), [(0, -2)])
assert bestI[13][0] == F(4, 3) and sorted(bestI[13][1]) == [(-1, 13), (1, 13)]
assert bestI[17][0] == F(40, 39) and sorted(bestI[17][1]) == [(-1, 11), (1, 11)]
assert bestI[3] == (F(27, 13), [(0, 18)]) and bestI[36] == (F(27, 13), [(0, -18)])
assert 2 - F(1 * 3, 12) - F(13 * 15, 164) == F(23, 41)
assert 8 - 1 - F(117, 41) == F(170, 41) and F(18, 5) < F(170, 41)
assert all(b * (b + 2) != 492 for b in range(0, 40))
assert {(14 * x) % 39 for x in (4, 35)} == {17, 22}
report('levels (1,39): labels -2n, relevant classes of grade <= 2 are +-4, +-13, +-17, eigenvalue 23/41')
report('levels (1,39): class +-3 at (0,+-18) in grade 3, coset weight 3 - b(b+2)/164 is never zero')
report('ALL CHECKS PASSED')
