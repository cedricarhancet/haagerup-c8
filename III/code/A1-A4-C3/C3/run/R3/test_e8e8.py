"""Tests of e8e8_voa: dimensions, Virasoro relations, central charge, skew-symmetry, bilinear form, phi_0."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import time, random
from gmpy2 import mpq
import e8e8_voa as V
t0 = time.time()
B0, B1, B2 = V.basis(0), V.basis(1), V.basis(2)
print('dims', len(B0), len(B1), len(B2)); assert (len(B0), len(B1), len(B2)) == (1, 496, 69752)
om = V.omega()
print('omega_(3) omega =', V.mode_vec(om, 3, om))   # c/2 = 8
print('L1 omega =', V.Lmode(1, om), ' L0 omega = 2 omega:', V.Lmode(0, om) == V.scale(om, 2))
random.seed(1)
fails = 0
for k in random.sample(B1, 20) + random.sample(B2, 40):
    v = {k: mpq(1)}
    if V.Lmode(0, v) != V.scale(v, V.weight(k)): fails += 1
print('L0 = weight on 60 random states, failures:', fails)
# skew symmetry u_(n) v = sum_j (-1)^{n+j+1} / j! L_{-1}^j (v_(n+j) u)   for u, v in W_1 (n >= 0 into W_1, W_0)
def Lm1(v): return V.Lmode(-1, v)
fails = 0
for _ in range(30):
    u = {random.choice(B1): mpq(1)}; v = {random.choice(B1): mpq(1)}
    for n in (0, 1):
        lhs = V.mode_vec(u, n, v)
        rhs = {}; x = V.mode_vec(v, n, u); fact = 1; j = 0
        terms = []
        for j in range(0, 3):
            x = V.mode_vec(v, n + j, u)
            for _ in range(j): x = Lm1(x) if x else x
            from math import factorial
            V.vadd(rhs, x, mpq((-1)**(n + j + 1), factorial(j)))
        if V.clean(lhs) != V.clean(rhs): fails += 1
print('skew-symmetry W1 x W1 (n=0,1): failures', fails)
# Jacobi / commutator on weight-1: [u_(0), v_(0)] w = (u_(0) v)_(0) w on W_1 and W_2 samples
fails = 0
for _ in range(30):
    u = {random.choice(B1): mpq(1)}; v = {random.choice(B1): mpq(1)}; w = {random.choice(B1 + random.sample(B2, 50)): mpq(1)}
    lhs = V.vadd(V.mode_vec(u, 0, V.mode_vec(v, 0, w)), V.mode_vec(v, 0, V.mode_vec(u, 0, w)), -1)
    rhs = V.mode_vec(V.mode_vec(u, 0, v), 0, w)
    if V.clean(lhs) != V.clean(rhs): fails += 1
print('commutator formula, failures:', fails)
# commutator [u_(1), v_(0)] for u in W_2, v in W_1 on W_1: [v_(0), u_(1)] = (v_(0) u)_(1)
fails = 0
for _ in range(20):
    u = {random.choice(B2): mpq(1)}; v = {random.choice(B1): mpq(1)}; w = {random.choice(B1): mpq(1)}
    lhs = V.vadd(V.mode_vec(v, 0, V.mode_vec(u, 1, w)), V.mode_vec(u, 1, V.mode_vec(v, 0, w)), -1)
    rhs = V.mode_vec(V.mode_vec(v, 0, u), 1, w)
    if V.clean(lhs) != V.clean(rhs): fails += 1
print('[v_(0), u_(1)] = (v_(0)u)_(1) on W_1 (u in W_2), failures:', fails)
# phi_0 automorphism and positivity of (phi_0 u, u)
def phi0(v): return {(V.NEGI[b], m): c*(-1)**len(m) for (b, m), c in v.items()}
fails = 0; tests = 0
for _ in range(25):
    u = {random.choice(B1 + random.sample(B2, 30)): mpq(1)}; v = {random.choice(B1 + random.sample(B2, 30)): mpq(1)}
    wu, wv = V.weight(next(iter(u))), V.weight(next(iter(v)))
    for n in range(wu + wv - 3, wu + wv):
        if wu + wv - n - 1 > 2 or wu + wv - n - 1 < 0: continue
        a = V.mode_vec(u, n, v); b = V.mode_vec(phi0(u), n, phi0(v)); tests += 1
        if phi0(a) != V.clean(b): fails += 1
print('phi_0 automorphism tests', tests, 'failures', fails)
bad = 0
for k in random.sample(B1, 40) + random.sample(B2, 40):
    f = V.form(phi0({k: mpq(1)}), {k: mpq(1)})
    if f <= 0: bad += 1
print('(phi_0 u, u) > 0 on 80 random basis states: failures', bad)
print('time %.1fs' % (time.time() - t0))
