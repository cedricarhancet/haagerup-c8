# Step 1: explicit sl2^4 in e8+e8 = W_1, verification of relations, levels, W_1 decomposition
import sys, pickle, json, itertools, time, os
sys.path.insert(0, '.')
from voa import *
from lat import *
t0 = time.time()
CO = {(3, 4): Fr(-2), (3, 6): Fr(-2)}   # from solve_signs.py (torus gauge: all others = 1)
E = []; F = []; H = []
for X in range(4):
    e = {}; f = {}
    for j, (a, k) in enumerate(support(X)):
        c = CO.get((X, j), Fr(1))
        e[((), a)] = c
        f[((), vneg(a))] = Fr(-k)/c
    E.append(e); F.append(f)
    Hv = Hvec(X)
    H.append({(((1, i),), ZERO): Fr(Hv[i], 2) for i in range(16) if Hv[i]})
def br(x, y): return mode_vec(x, 0, y)
def scal(c, x): return {k: c*v for k, v in x.items()}
ok = True
for X in range(4):
    ok &= br(E[X], F[X]) == H[X]
    ok &= br(H[X], E[X]) == scal(2, E[X])
    ok &= br(H[X], F[X]) == scal(-2, F[X])
    lev = mode_vec(E[X], 1, F[X])
    print('sl2 #%d: E_(1)F =' % X, lev, ' H_(1)H =', mode_vec(H[X], 1, H[X]))
    ok &= lev == {((), ZERO): Fr(LEVELS[X])}
    ok &= mode_vec(H[X], 1, H[X]) == {((), ZERO): Fr(2*LEVELS[X])}
    ok &= mode_vec(E[X], 1, E[X]) == {} and mode_vec(F[X], 1, F[X]) == {}
    for Y in range(4):
        if X != Y:
            for a in (E[X], F[X], H[X]):
                for b in (E[Y], F[Y], H[Y]):
                    ok &= br(a, b) == {}
                    ok &= mode_vec(a, 1, b) == {}
print('Chevalley relations + commutation + levels OK:', ok)
# indices in each copy
for X in range(4):
    Hv = Hvec(X)
    print('sl2 #%d index copy1 = %s, copy2 = %s' % (X, Fr(sum(x*x for x in Hv[:8]), 8), Fr(sum(x*x for x in Hv[8:]), 8)))
# W_1 decomposition by highest weight vectors: ker of all E_X(0) on each weight space
from collections import defaultdict
wsp = defaultdict(list)
for a in ROOTS:
    wsp[tuple(ip(a, Hvec(X)) for X in range(4))].append(a)
import sympy as sp
dec = {}
HW = {}
for tau, roots in wsp.items():
    if not all(x >= 0 for x in tau) or tau == (0, 0, 0, 0):
        continue
    # matrix of E_X(0) on span{e^a : a in roots}
    rows = {}
    cols = []
    for a in roots:
        col = {}
        for X in range(4):
            r = br(E[X], e_state(a))
            for k, v in r.items():
                col[(X, k)] = v
        cols.append(col)
    keys = sorted(set(k for c in cols for k in c), key=str)
    M = sp.Matrix(len(keys), len(roots), lambda i, j: cols[j].get(keys[i], 0))
    ns = M.nullspace()
    if ns:
        dec[tau] = len(ns)
        HW[tau] = [{((), roots[j]): Fr(int(sp.fraction(v[j])[0]), int(sp.fraction(v[j])[1])) for j in range(len(roots)) if v[j] != 0} for v in ns]
exp = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 't12_W1.json')))['W1']
exp = {tuple(int(x) for x in k.split(',')): v for k, v in exp.items()}
print('hw-vector decomposition:', sorted(dec.items()))
print('W1 decomposition matches expected (20 types):', dec == exp, len(dec))
print('dim check', sum(m*(t[0]+1)*(t[1]+1)*(t[2]+1)*(t[3]+1) for t, m in dec.items()))
pickle.dump({'E': E, 'F': F, 'H': H, 'HW': HW}, open('sl2_4.pkl', 'wb'))
print('time %.1fs' % (time.time()-t0))
