# Step 3a: algebra structure on G: Gram form, products u_(1)u', o(u) on highest-weight spaces of W_1
import sys, pickle, time
sys.path.insert(0, '.')
from voa import *
from lat import *
import flint
t0 = time.time()
S = pickle.load(open('sl2_4.pkl', 'rb'))
G = pickle.load(open('G.pkl', 'rb'))
Gb, pivots = G['Gb'], G['pivots']
HW = S['HW']
n = len(Gb)
cache = {}
def bprod(s, t, k):
    """s_(k) t for basis keys"""
    key = (s, t, k)
    r = cache.get(key)
    if r is None:
        # quick zero test for pure lattice states: e^b_(k) e^g = 0 unless <b,g> <= -(k+1) (z^{<b,g>} z^{d}, d>=0)
        if not s[0] and not t[0] and ip(s[1], t[1]) > -(k+1):
            r = {}
        else:
            r = mode(s, k, {t: Fr(1)})
        cache[key] = r
    return r
def prod(u, v, k):
    out = {}
    for s, a in u.items():
        for t, b in v.items():
            r = bprod(s, t, k)
            if r:
                sadd(out, r, a*b)
    return out
def coords_in_G(v):
    co = [v.get(p, Fr(0)) for p in pivots]
    rec = {}
    for c, g in zip(co, Gb):
        if c: sadd(rec, g, c)
    return co, clean(rec) == clean(v)
gram = [[None]*n for _ in range(n)]
P = [[None]*n for _ in range(n)]
allin = True
for i in range(n):
    for j in range(i, n):
        r3 = prod(Gb[i], Gb[j], 3)
        assert set(r3) <= {((), ZERO)}
        gram[i][j] = gram[j][i] = r3.get(((), ZERO), Fr(0))
        r1 = prod(Gb[i], Gb[j], 1)
        co, ok = coords_in_G(r1)
        allin &= ok
        P[i][j] = P[j][i] = co
    print('row', i, 'done, %.1fs' % (time.time()-t0), flush=True)
print('all products u_i(1)u_j lie in G:', allin)
# symmetric check u_(1)v = v_(1)u (+ L(-1) term which vanishes in weight 2 for k=1? check numerically a few)
r_a = prod(Gb[0], Gb[5], 1); r_b = prod(Gb[5], Gb[0], 1)
print('u0_(1)u5 == u5_(1)u0 :', clean(r_a) == clean(r_b))
# o(u_i) on HW_tau: o(u)v = v_(1)u  (skew-symmetry, weights 1 and 2)
Mt = {}
for tau, vecs in sorted(HW.items()):
    m = len(vecs)
    # coordinates in HW basis: use pivot keys of vecs
    keys = sorted(set(k for v in vecs for k in v), key=str)
    A = flint.fmpq_mat(len(keys), m)
    for c, v in enumerate(vecs):
        for r, k in enumerate(keys):
            A[r, c] = flint.fmpq(v.get(k, Fr(0)).numerator, v.get(k, Fr(0)).denominator)
    mats = []
    for i in range(n):
        Mi = [[Fr(0)]*m for _ in range(m)]
        for c, v in enumerate(vecs):
            img = prod(v, Gb[i], 1)
            # solve A x = img
            bvec = flint.fmpq_mat(len(keys), 1)
            assert set(img) <= set(keys), 'image leaves HW_tau weight space?'
            for r, k in enumerate(keys):
                x = img.get(k, Fr(0)); bvec[r, 0] = flint.fmpq(x.numerator, x.denominator)
            # least squares exact: solve via normal eq on full-rank A
            At = A.transpose()
            sol = (At*A).solve(At*bvec)
            assert A*sol == bvec, 'o(u) does not preserve HW_tau'
            for r in range(m):
                Mi[r][c] = Fr(int(sol[r, 0].p), int(sol[r, 0].q))
        mats.append(Mi)
    Mt[tau] = mats
print('o(u_i) matrices on HW spaces computed; preserved: True')
# direct cross-check of skew symmetry for one tau, i: u_(1)v computed with composite modes
tau0 = (0, 0, 0, 2); v = HW[tau0][0]
print('cross-check o(u)v = v_(1)u :', clean(prod(Gb[3], v, 1)) == clean(prod(v, Gb[3], 1)))
pickle.dump({'gram': gram, 'P': P, 'Mt': Mt}, open('griess.pkl', 'wb'))
print('time %.1fs' % (time.time()-t0))
