#!/usr/bin/env python3
# Exact complements to lattice/verify_c6_section6_v01.py for
# "Haagerup symmetry in a rank-six lattice conformal field theory" (source v12).
# It checks the rational mass certificate of Proposition 6.2 with the local factors of
# Conway--Sloane IV, the data of the realization of levels (1,39) used in Proposition 7.14,
# the signs and values of the rescaling-invariant cubic quantities of Section 7.3, the
# cocycle parity identities of Section 7.3 and the pointed labels of Proposition 7.16.
# Requires Python 3 and cypari2 (PARI/GP). Run: python3 -I verify_c6_complements_v12.py
from fractions import Fraction as F
from itertools import product
from math import gcd
import cypari2

pari = cypari2.Pari()


def report(msg):
    print(msg)


def pmat(G):
    return pari.matrix(len(G), len(G), [x for r in G for x in r])


def jacobi(a, n):
    return int(pari.kronecker(a, n))


GK = [[2, -1, 0, 0, 0, -1], [-1, 2, 0, 0, 0, 0], [0, 0, 2, -1, 0, -1],
      [0, 0, -1, 2, 0, 0], [0, 0, 0, 0, 2, -1], [-1, 0, -1, 0, -1, 4]]
GKf = [[2, -1, -1, 1, -1, -1], [-1, 2, 1, 0, 0, 1], [-1, 1, 2, 0, 0, 0],
       [1, 0, 0, 2, -1, -1], [-1, 0, 0, -1, 2, 0], [-1, 1, 0, -1, 0, 8]]

# ---------------------------------------------------------------------------
# 1. Mass of the genus of K (Proposition 6.2), Conway--Sloane IV
chi = lambda a: 0 if gcd(a, 39) > 1 else jacobi(a, 39)
B3 = lambda x: x ** 3 - F(3, 2) * x ** 2 + F(1, 2) * x
B3chi = 39 ** 2 * sum(chi(a) * B3(F(a, 39)) for a in range(1, 40))
assert B3chi == 528
# (-39/m) = (m/39) for odd m, and chi(2) = 1, so zeta_{-39}(3) = (1 - 2^-3) L(3,chi)
assert all(jacobi(-39, m) == chi(m) for m in range(1, 2000, 2))
assert chi(2) == 1
# 2-adic: type II, dimension 6, det 39 = -1 mod 8, octane 0, free, species 6+
assert 39 % 8 == 7
M2 = 1 / (2 * (1 - F(1, 4)) * (1 - F(1, 16)) * (1 - F(1, 8)))
std2 = 1 / (2 * (1 - F(1, 4)) * (1 - F(1, 16)))
ratio2 = M2 * F(1, 2 ** 6) / std2
assert ratio2 == F(1, 56)
# p in {3,13}: constituents of dimensions 5 and 1, cross term p^(5/2): ratio p^(5/2)/2.
# L(3,chi) = 2 pi^3 B3chi / (3 * 39^(5/2)), so the product of the two odd ratios, 39^(5/2)/4,
# cancels the factor 39^(-5/2) and the mass is rational:
mass = F(7, 8) * F(2 * B3chi, 3) / 360 * ratio2 * F(1, 4)
assert mass == F(B3chi, 138240) == F(11, 2880)
oK = int(pari.qfauto(pmat(GK))[0])
oKf = int(pari.qfauto(pmat(GKf))[0])
assert (oK, oKf) == (288, 2880)
assert F(1, oK) + F(1, oKf) == mass
assert int(pari.qfisom(pmat(GK), pmat(GKf)) == 0) == 1
report('mass: B_{3,chi} = 528, 2-adic ratio 1/56, mass 11/2880 = 1/288 + 1/2880 with |O(K)| = 288, |O(K-flat)| = 2880')

# ---------------------------------------------------------------------------
# 2. Realization of levels (1,39) (Proposition 7.14)
qK = lambda p: F(17 * p * p, 39) % 1
eI = lambda m, n: F(m * m, 4) + F(n * n, 156)
lab = lambda n: (-2 * n) % 39
sug = lambda a, k: F(a * (a + 2), 4 * (k + 2))
cc = lambda k: F(3 * k, k + 2)
assert 8 - cc(1) - cc(39) == F(170, 41)
assert (lab(13), lab(11), lab(18), lab(2), lab(-2)) == (13, 17, 3, 35, 4)
assert eI(1, 13) == F(4, 3) and eI(1, 11) == F(40, 39) and eI(0, 18) == F(27, 13) and eI(0, 2) == F(1, 39)
# grade-two relevant vectors of momenta (+-1,+-13) and (+-1,+-11) lie in modules of type V(1,13)
w2 = 2 - sug(1, 1) - sug(13, 39)
assert w2 == F(23, 41)
# grade three, momentum (0,18): b is even, b >= 18, and b = 20 would put the momentum (0,20),
# of label -1, at Omega-weight 3 - 100/39 = 17/39 = q(-1) < 1, which is excluded since -1 is not relevant
assert lab(20) == 38 and 3 - eI(0, 20) == F(17, 39) == qK(-1)
w3 = 3 - sug(18, 39)
assert w3 == F(33, 41)
assert all(3 - sug(b, 39) != 0 for b in range(0, 60))
# root currents of level 39: momenta (0,+-2), labels -+4, Omega-weight 1 - 1/39 = 38/39
assert 1 - eI(0, 2) == F(38, 39)
assert sorted({(14 * x) % 39 for x in (4, 35)}) == [17, 22]
report('levels (1,39): central charge 170/41, eigenvalues 23/41 on +-13 and +-17, 33/41 on +-3 (b = 18 forced)')

# ---------------------------------------------------------------------------
# 3. Rescaling-invariant cubic quantities (Section 7.3 and Proposition 7.8)
k = 3
c13sq = F(2 * k * (k - 1), k * k)       # |e_{(-1)}^2 1|^2 / 3^2 with |e|^2 = k
c4sq = F(1, k)                          # |a_{(2)} b|^2 on the unit vector f/sqrt(3)
assert (c13sq, c4sq) == (F(4, 3), F(1, 3))
sign = lambda h1, h2, h3: (-1) ** (h1 + h2 + h3)
rhoU13 = sign(1, 1, 2) * c13sq           # weights of x_13, x_13, x_26 in U
rhoU4 = sign(2, 2, 1) * c4sq             # weights of x_4, x_{-17}, x_{-13} in U
assert (rhoU13, rhoU4) == (F(4, 3), F(-1, 3))
assert (rhoU13 ** 2, rhoU4 ** 2) == (F(16, 9), F(1, 9))
report('cubic quantities: rho^U = 4/3 and -1/3, rho = 16/9 and 1/9, |C| = 4/3 and 1/3')

# ---------------------------------------------------------------------------
# 4. Cocycle parities (Section 7.3): 39^2 is odd
for a, b in product(range(-78, 79), repeat=2):
    s, t = F(a, 39), F(b, 39)
    if (s + t).denominator == 1:
        assert (1521 * (s + t) - (s + t)) % 2 == 0
    if (s - t).denominator == 1:
        assert (1521 * (s + t) - (s - t)) % 2 == 0
report('cocycles: commutator parities on E8 and on Lambda_K checked for all pairs in (1/39)Z, |numerators| <= 78')

# ---------------------------------------------------------------------------
# 5. Mixed pointed bulk (Proposition 7.16)
mono = lambda p, s: (qK(p + s) - qK(p) - qK(s)) % 1
local = [p for p in range(39) if mono(p, 13) == mono(p, 26) == 0]
assert local == list(range(0, 39, 3))
assert all(F(34 * p * 13, 39) % 1 == mono(p, 13) for p in range(39))
assert all(qK(3 * r) == (-F(r * r, 13)) % 1 for r in range(13))
relevant = {3, 36, 13, 26, 4, 35, 17, 22}
assert sorted(p for p in relevant if p in local) == [3, 36]
report('mixed bulk: local labels 3Z/39Z, monodromy exp(2 pi i 34p/3), pointed relevant labels +-3')
report('ALL CHECKS PASSED')
