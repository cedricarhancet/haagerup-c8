"""Tests of lattice_voa.py (Axis A0): weights, Virasoro relations, Sugawara, translation, Borcherds commutator."""
import random, time
from gmpy2 import mpq
from math import comb
import lattice_voa as V
random.seed(1)
t0 = time.time()
ok = True
def eq(a, b): return V.clean(a) == V.clean(b)
def sub(a, b):
    out = dict(a); V.vadd(out, b, -1); return out
# 1. L_0 = weight on random basis states
for w in range(4):
    B = V.basis(w)
    for key in random.sample(B, min(40, len(B))):
        if not eq(V.Lmode(0, {key: mpq(1)}), {key: mpq(w)}): ok = False; print('L0 fail', key)
print('L0 = weight: ok' if ok else 'L0 FAIL')
# 2. omega_3 omega = (c/2) 1, c = 8 ; L_1 omega = 0, L_2 omega = (c/2) 1
c3 = V.mode_vec(V.OMEGA, 3, V.OMEGA)
print('omega_3 omega =', dict(c3), ' (expected 4 * vacuum)')
assert c3 == {V.VAC: mpq(4)}
assert not V.clean(V.mode_vec(V.OMEGA, 2, V.OMEGA))
# 3. Virasoro commutators on random states of weight <= 2
def L(n, v): return V.Lmode(n, v)
cnt = 0
for w in (1, 2):
    for key in random.sample(V.basis(w), 15):
        v = {key: mpq(1)}
        for m, n in ((1, -1), (2, -1), (1, -2), (2, -2), (0, -1), (3, -2), (2, 1), (-1, 1)):
            lhs = sub(L(m, L(n, v)), L(n, L(m, v)))
            rhs = {k: x*(m - n) for k, x in L(m + n, v).items()}
            if m + n == 0:
                V.vadd(rhs, v, mpq(8*(m**3 - m), 12))
            if not eq(lhs, rhs): ok = False; print('Vir fail', key, m, n)
            cnt += 1
print('Virasoro commutators:', cnt, 'checks', 'ok' if ok else 'FAIL')
# 4. translation: (L_{-1} u)_k = -k u_{k-1}
for key in random.sample(V.basis(2), 10):
    u = {key: mpq(1)}; Lu = L(-1, u)
    for vk in random.sample(V.basis(1), 5):
        v = {vk: mpq(1)}
        for k in range(0, 4):
            if not eq(V.mode_vec(Lu, k, v), {kk: -k*x for kk, x in V.mode_vec(u, k - 1, v).items()}):
                ok = False; print('translation fail')
print('translation property', 'ok' if ok else 'FAIL')
# 5. Borcherds commutator [u_m, v_n] w = sum_j C(m,j) (u_j v)_{m+n-j} w on random triples
cnt = 0
for _ in range(60):
    wu = random.choice((1, 1, 2)); wv = random.choice((1, 2)); ww = random.choice((0, 1, 2))
    u = {random.choice(V.basis(wu)): mpq(1)}; v = {random.choice(V.basis(wv)): mpq(random.randint(1, 3))}
    w = {random.choice(V.basis(ww)): mpq(1)}
    m = random.randint(-1, 3); n = random.randint(-1, 3)
    lhs = sub(V.mode_vec(u, m, V.mode_vec(v, n, w)), V.mode_vec(v, n, V.mode_vec(u, m, w)))
    rhs = {}
    for j in range(0, wu + wv + 2):
        if m >= 0 and j > m: break
        x = V.mode_vec(u, j, v)
        if x: V.vadd(rhs, V.mode_vec(x, m + n - j, w), comb(m, j) if m >= 0 else (-1)**j*comb(-m + j - 1, j))
    if not eq(lhs, rhs): ok = False; print('Borcherds fail', u, v, w, m, n)
    cnt += 1
print('Borcherds commutator:', cnt, 'random triples', 'ok' if ok else 'FAIL')
# 6. Sugawara of e6 + sl3 at level one equals omega
def sugawara():
    out = {}
    # Cartan parts: (1/(2(1+h))) sum G^{-1}_ij alpha_i(-1) alpha_j(-1) on each factor
    for (idx, hv) in ((range(6), 12), (range(6, 8), 3)):
        for i in idx:
            for j in idx:
                if V.GINV[i][j]:
                    k = (V.ZERO, tuple(sorted(((1, i), (1, j)))))
                    out[k] = out.get(k, 0) + V.GINV[i][j]/(2*(1 + hv))
    for a in V.ROOTS:
        hv = 12 if all(x == 0 for x in a[6:]) else 3
        dual = {(V.neg(a), ()): mpq(1, V.eps(a, V.neg(a)))}
        x = V.mode_vec({(a, ()): mpq(1)}, -1, dual)
        V.vadd(out, x, mpq(1, 2*(1 + hv)))
    return V.clean(out)
S = sugawara()
print('Sugawara(e6 + sl3, level 1) == omega:', S == V.OMEGA)
ok = ok and S == V.OMEGA
# 7. invariant form: <omega, omega> = c/2, <alpha_i(-1)1, alpha_j(-1)1> = -G_ij
print('<omega,omega> =', V.form(V.OMEGA, V.OMEGA), '; <a1(-1), a3(-1)> =', V.form(V.state_h(0), V.state_h(2)))
print('ALL TESTS', 'PASSED' if ok else 'FAILED', ' time %.1fs' % (time.time() - t0))
