# Step 3d/4: build the 4 solutions exactly over K = Q(sqrt(-2), sqrt(3)), verify them directly in W,
# intrinsic invariants, the antilinear involution phi_chi.
import sys, pickle, time, itertools, random
sys.path.insert(0, '.')
from voa import *
from lat import *
from q4 import Q4
t0 = time.time()
S = pickle.load(open('sl2_4.pkl', 'rb'))
G = pickle.load(open('G.pkl', 'rb'))
E, F, H, HW = S['E'], S['F'], S['H'], S['HW']
Gb, pivots, omCC = G['Gb'], G['pivots'], G['omCC']
s2 = Q4(0, 1, 0, 0); s3 = Q4(0, 0, 1, 0)
sols = []
for e13, e14 in itertools.product((1, -1), repeat=2):
    x13 = s2*Fr(e13, 39); x14 = s3*Fr(e14, 39)
    x = {0: Fr(7, 195), 1: Fr(1, 15), 2: Fr(-2, 65), 3: Fr(0), 5: Fr(-4, 39), 6: Fr(-2, 65), 7: Fr(7, 39), 8: Fr(-4, 39), 9: Fr(1, 10), 15: Fr(0)}
    x[4] = -x[5]/2 - x[6]/2 - Fr(1, 5)
    x = {k: Q4(v) for k, v in x.items()}
    x[13] = x13; x[14] = x14; x[10] = -x13; x[11] = 39*x13*x14; x[12] = 39*x13*x14
    sols.append([x[i] for i in range(16)])
def tostate(xc):
    st = {}
    for c, g in zip(xc, Gb):
        if c:
            for k, v in g.items():
                nv = st.get(k, Q4(0)) + c*v
                if nv: st[k] = nv
                else: st.pop(k, None)
    return st
def prod(u, v, k):
    out = {}
    for s, a in u.items():
        for t, b in v.items():
            if not s[0] and not t[0] and ip(s[1], t[1]) > -(k+1):
                continue
            r = mode(s, k, {t: Fr(1)})
            for kk, vv in r.items():
                nv = out.get(kk, Q4(0)) + a*b*vv
                if nv: out[kk] = nv
                else: out.pop(kk, None)
    return out
def eqst(a, b):
    keys = set(a) | set(b)
    return all(Q4.lift(a.get(k, 0)) == Q4.lift(b.get(k, 0)) for k in keys)
def scal(c, a): return {k: c*v for k, v in a.items()}
gens = E + F + H
U = [tostate(x) for x in sols]
omCCq = {k: Q4(v) for k, v in omCC.items()}
print('solutions (G-coordinates x0..x15):')
for i, x in enumerate(sols):
    print('  u%d:' % i, x)
hz_ok = True
def hA(a, b): return Fr(a*(a+2), 20) + Fr(b*(b+2), 60)
import json
expected = pickle.load(open('system.pkl', 'rb'))['expected']
import sympy as sp
sq2 = sp.sqrt(-2); sq3 = sp.sqrt(3)
def tosp(q):
    a, b, c, d = q.a
    return sp.Rational(a.numerator, a.denominator) + sp.Rational(b.numerator, b.denominator)*sq2 + sp.Rational(c.numerator, c.denominator)*sq3 + sp.Rational(d.numerator, d.denominator)*sq2*sq3
for k, u in enumerate(U):
    inG = all(not prod(scal(Q4(1), x), u, n) for x in gens for n in (0, 1, 2))
    uu1 = prod(u, u, 1)
    idem = eqst(uu1, scal(2, u))
    uu3 = prod(u, u, 3)
    # o(u) on every HW_tau : o(u)v = v_(1)u
    eig_ok = True
    for tau, vecs in HW.items():
        m = len(vecs)
        keys = sorted(set(kk for v in vecs for kk in v), key=str)
        Amat = sp.Matrix(len(keys), m, lambda r, c: sp.Rational(vecs[c].get(keys[r], Fr(0)).numerator, vecs[c].get(keys[r], Fr(0)).denominator))
        cols = []
        for v in vecs:
            img = prod(scal(Q4(1), v), u, 1)
            assert set(img) <= set(keys)
            b = sp.Matrix([tosp(img.get(kk, Q4(0))) for kk in keys])
            sol = (Amat.T*Amat).LUsolve(Amat.T*b)
            assert sp.simplify(Amat*sol - b) == sp.zeros(len(keys), 1)
            cols.append(sol)
        M = sp.Matrix.hstack(*cols)
        tq = sp.Symbol('t')
        cp = sp.expand((tq*sp.eye(m) - M).det())
        tgt = sp.expand(sp.prod([tq - sp.Rational(l.numerator, l.denominator) for l in expected[tau]]))
        good = sp.simplify(cp - tgt) == 0
        # diagonalizable: eigenvalues are distinct within each tau except possibly repeated ones
        if len(set(expected[tau])) < m:
            good &= M.is_diagonalizable()
        eig_ok &= good
    up = sadd(dict(omCCq), u, -1)
    ortho = not prod(u, up, 1) and not prod(u, up, 3) if True else None
    print('u%d: in Com(AxA): %s | u_(1)u = 2u: %s | u_(3)u = %s | o(u) spectrum on W_1 as required: %s | u_(n)(omega_CC - u)=0 n=1,3: %s' % (
        k, inG, idem, uu3, eig_ok, ortho))
# intrinsic invariants: Gram matrix of the 4 solutions, (u, omega_CC)
print('Gram matrix (u_k, u_l) = (u_k)_(3) u_l:')
for k in range(4):
    print('   ', [prod(U[k], U[l], 3).get(((), ZERO), Q4(0)) for l in range(4)])
print('(u_k, omega_CC):', [prod(U[k], omCCq, 3).get(((), ZERO), Q4(0)) for k in range(4)])
pickle.dump({'sols': sols}, open('solutions.pkl', 'wb'))
print('time %.1fs' % (time.time()-t0))
