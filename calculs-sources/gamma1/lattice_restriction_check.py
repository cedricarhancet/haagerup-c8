"""Characters of the nine modules V_{(i,j)} of V(E6+A2) versus the Evans--Gannon characters X^(gamma).
Checks (exact integers, to order q^N):
  Theta_{(0,0)}/eta^8 = X_0 + X_b   (for both gamma),
  Theta_{(1,1)}/eta^8 = X_a = X_c0 = Theta_{(1,2)}/eta^8,
  Theta_{(0,1)}/eta^8 = X_c1,  Theta_{(1,0)}/eta^8 = X_c2.
E6 cosets are read off from E8 = union of (E6 + i mu) x (A2 + i nu): fixing the A2-component of an E8
vector fixes the E6 coset.  Coordinates: E8 in doubled coordinates, A2 spanned by e1-e2, e2-e3.
Usage: python3 lattice_restriction_check.py
"""
import math, itertools
from fractions import Fraction as Fr
import xi_exact as X
N = 4                                     # number of q-powers checked after the leading one
def gen(maxnorm2):                        # doubled coords, |x|^2 = sum(x_i^2)/4 <= maxnorm2/2 ... (norm = |x|^2)
    lim = 4*maxnorm2; out = []
    for parity in (0, 1):
        stack = [((), 0)]
        while stack:
            cur, s2 = stack.pop()
            if len(cur) == 8:
                if sum(cur) % 4 == 0: out.append(cur)
                continue
            r = math.isqrt(lim - s2)
            for x in range(-r, r + 1):
                if (x - parity) % 2 == 0: stack.append((cur + (x,), s2 + x*x))
    return out
V = gen(2*N + 4)
# A2 component: projection onto span(e1-e2, e2-e3) = vectors in first three coordinates with zero sum
def a2_part(v):            # v doubled coords; return doubled-coordinate projection (first 3 coords minus mean)
    m = Fr(v[0] + v[1] + v[2], 3)
    return (Fr(v[0]) - m, Fr(v[1]) - m, Fr(v[2]) - m)
def norm(w):               # actual norm of doubled-coordinate vector
    return Fr(sum(x*x for x in w)) / 4
from collections import defaultdict
theta_e6 = {}              # A2-component -> theta series of the E6 coset (coeff dict by norm/2)
groups = defaultdict(lambda: defaultdict(int))
for v in V:
    p = a2_part(v)
    e6n = norm(v) - norm(p)
    groups[p][e6n/2] += 1
zero = (Fr(0), Fr(0), Fr(0))
# representatives of A2*/A2 classes: omega1 = (2/3,-1/3,-1/3) in real coords -> doubled (4/3,-2/3,-2/3)
w1 = (Fr(4, 3), Fr(-2, 3), Fr(-2, 3)); w2 = (Fr(-4, 3), Fr(2, 3), Fr(2, 3))
E6 = {0: groups[zero], 1: groups[w1], 2: groups[w2]}
# A2 cosets: theta series of A2 + j*omega
def a2_theta(shift, maxh):
    out = defaultdict(int)
    for a, b in itertools.product(range(-12, 13), repeat=2):
        # A2 vector a*(e1-e2) + b*(e2-e3) in doubled coords
        x = (2*a + shift[0], 2*(b - a) + shift[1], -2*b + shift[2])
        h = norm(x)/2
        if h <= maxh: out[h] += 1
    return out
A2 = {0: a2_theta(zero, N + 2), 1: a2_theta(w1, N + 2), 2: a2_theta(w2, N + 2)}
def coset_theta(i, j, maxh):
    out = defaultdict(int)
    for h1, c1 in E6[i].items():
        for h2, c2 in A2[j].items():
            if h1 + h2 <= maxh: out[h1 + h2] += c1*c2
    return out
# 1/eta^8 = q^{-1/3} prod (1-q^n)^{-8}
P = [1] + [0]*(N + 3)
for n in range(1, N + 4):
    for _ in range(8):
        for k in range(n, N + 4): P[k] += P[k - n]
def char(i, j):
    th = coset_theta(i, j, N + 2)
    hmin = min(th)
    ser = [0]*(N + 1)
    for h, c in th.items():
        k = h - hmin
        if k.denominator != 1: raise ValueError('non-integral spacing')
        for m in range(N + 1 - int(k)):
            ser[int(k) + m] += c*P[m]
    return hmin - Fr(1, 3), ser
def Xchar(z, gamma, nmax):
    # X_z = sum_n (Xi[n][z][0] + gamma Xi[n][z][d3]) q^{Lam_z + n}
    d3 = X.names.index('d3')
    ser = [X.Xi[n][z][0] + gamma*X.Xi[n][z][d3] for n in range(-1, nmax)]
    while ser and ser[0] == 0: ser = ser[1:]; 
    return ser
names = X.names
ok = True
for gamma in (0, 1):
    for (i, j), zs in {(0, 0): ['0', 'b'], (1, 1): ['a'], (2, 2): ['a'], (1, 2): ['c0'], (2, 1): ['c0'],
                       (0, 1): ['c1'], (0, 2): ['c1'], (1, 0): ['c2'], (2, 0): ['c2']}.items():
        e, ser = char(i, j)
        tot = None
        for z in zs:
            zi = names.index(z)
            s = [X.Xi[n][zi][0] + gamma*X.Xi[n][zi][names.index('d3')] for n in range(-1, N + 2)]
            tot = s if tot is None else [a + b for a, b in zip(tot, s)]
        while tot and tot[0] == 0: tot = tot[1:]
        match = [int(x) for x in tot[:N + 1]] == ser[:N + 1]
        ok &= match
        print('gamma=%d (%d,%d) leading exponent %s : theta/eta^8 %s  vs  X_%s %s  %s' % (
            gamma, i, j, e, ser[:N + 1], '+'.join(zs), [int(x) for x in tot[:N + 1]], 'OK' if match else 'MISMATCH'))
print('ALL OK' if ok else 'SOME MISMATCH')
