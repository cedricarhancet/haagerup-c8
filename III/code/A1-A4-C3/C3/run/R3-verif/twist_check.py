# Extra: the field of definition of u_k depends on the rational structure. Twisting the standard Q-form of
# V_L by sign characters t_a, t_b : L -> {+-1} that centralize sl2^4 (e'^beta = sqrt(-2)^{a(beta)} sqrt(3)^{b(beta)} e^beta)
# gives a Q-form in which A x A is still rational; we test in which twist u_0 becomes rational / Q(sqrt(-6))-rational.
import sys, pickle, itertools
sys.path.insert(0, '.')
from voa import *
from lat import *
from q4 import Q4
S = pickle.load(open('sl2_4.pkl', 'rb')); G = pickle.load(open('G.pkl', 'rb'))
E, F, H = S['E'], S['F'], S['H']
sols = pickle.load(open('solutions.pkl', 'rb'))['sols']
def tostate(xc):
    st = {}
    for c, g in zip(xc, G['Gb']):
        if c:
            for k, v in g.items():
                nv = st.get(k, Q4(0)) + c*v
                if nv: st[k] = nv
                else: st.pop(k, None)
    return st
U = [tostate(x) for x in sols]
# sign characters trivial on the support roots (as in phi_chi.py)
supp = [a[1] for X in range(4) for a in E[X]]
def gf2_null(rows, n):
    rows = [r[:] for r in rows]; piv = []; r = 0
    for c in range(n):
        p = next((i for i in range(r, len(rows)) if rows[i][c]), None)
        if p is None: continue
        rows[r], rows[p] = rows[p], rows[r]
        for i in range(len(rows)):
            if i != r and rows[i][c]:
                rows[i] = [(x+y) % 2 for x, y in zip(rows[i], rows[r])]
        piv.append(c); r += 1
    out = []
    for f in [c for c in range(16) if c not in piv]:
        v = [0]*n; v[f] = 1
        for i, c in enumerate(piv): v[c] = rows[i][f]
        out.append(v)
    return out
NB = gf2_null([[x % 2 for x in coords16(a)] for a in supp], 16)
chars = []
for bits in itertools.product((0, 1), repeat=len(NB)):
    chars.append([sum(b*v[i] for b, v in zip(bits, NB)) % 2 for i in range(16)])
def ev(s, beta): return sum(x*y for x, y in zip(s, coords16(beta))) % 2
inv2 = Q4(0, Fr(-1, 2), 0, 0)   # 1/sqrt(-2)
inv3 = Q4(0, 0, Fr(1, 3), 0)    # 1/sqrt(3)
def field_of(st):
    used = [False]*4
    for k, c in st.items():
        for j in range(4):
            if c.a[j]: used[j] = True
    return ['1', 'sqrt(-2)', 'sqrt(3)', 'sqrt(-6)'], used
best = []
for sa in chars:
    for sb in chars:
        # twisted coordinates of u_0
        st = {}
        for (m, beta), c in U[0].items():
            f = c
            if ev(sa, beta): f = f*inv2
            if ev(sb, beta): f = f*inv3
            st[(m, beta)] = f
        _, used = field_of(st)
        if used == [True, False, False, False]:
            best.append(('Q', sa, sb))
        elif used[1] is False and used[2] is False:
            best.append(('Q(sqrt(-6))', sa, sb))
kinds = {}
for kd, sa, sb in best: kinds[kd] = kinds.get(kd, 0) + 1
print('twists (t_a, t_b) among 16 x 16 for which u_0 has coordinates in Q:', kinds.get('Q', 0), '; in Q(sqrt(-6)) but not Q:', kinds.get('Q(sqrt(-6))', 0))
# show one rational twist, and verify E,F,H stay rational (they are: chars trivial on support roots)
for kd, sa, sb in best:
    if kd == 'Q':
        print('example rational twist: a =', sa, ' b =', sb)
        ok = all(ev(sa, k[1]) == 0 and ev(sb, k[1]) == 0 for x in E + F for k in x)
        print('  sl2^4 generators untouched by the twist:', ok)
        # coordinates of all 4 solutions in this twisted Q-form
        for j in range(4):
            st = {}
            for (m, beta), c in U[j].items():
                f = c
                if ev(sa, beta): f = f*inv2
                if ev(sb, beta): f = f*inv3
                st[(m, beta)] = f
            print('  u%d field generators used:' % j, [n for n, u in zip(*field_of(st)) if u])
        break
# closure of the twisted Q-form under products: e'^b_(n) e'^g must have rational coefficients in the e' basis
import random
random.seed(3)
sa, sb = [kd for kd in best if kd[0] == 'Q'][0][1:]
def twist(beta):
    f = Q4(1)
    if ev(sa, beta): f = f*Q4(0, 1, 0, 0)
    if ev(sb, beta): f = f*Q4(0, 0, 1, 0)
    return f
bad = 0; ntest = 0
pool = ROOTS + NORM4
for _ in range(400):
    b = random.choice(ROOTS); g = random.choice([x for x in pool if ip(x, b) < 0])
    for n in range(-1, 2):
        r = mode_vec({((), b): twist(b)}, n, {((), g): twist(g)})
        for (m, lat), c in r.items():
            # coefficient in the e' basis: c / twist(lat)
            t = twist(lat)
            # divide by t (t is 1, s2, s3 or s6, invertible): multiply by conj-like inverse
            inv = {Q4(1): Q4(1), Q4(0, 1, 0, 0): Q4(0, Fr(-1, 2), 0, 0), Q4(0, 0, 1, 0): Q4(0, 0, Fr(1, 3), 0), Q4(0, 0, 0, 1): Q4(0, 0, 0, Fr(-1, 6))}[t]
            q = Q4.lift(c)*inv
            ntest += 1
            if any(q.a[1:]): bad += 1
print('twisted Q-form closed under products (random e\'^b_(n) e\'^g, %d coefficients): failures %d' % (ntest, bad))
