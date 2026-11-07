# Step 5 (optional claim): S_2 x 1 = {v in W_2 : (1 x x)_(n) v = 0 (x in g'_2+g''_2, n>=0),
#                                     (omega_CC - u)_(n) v = 0, n = 1,2,3}   for u = u_k
# Linear algebra over K = Q(sqrt(-2), sqrt(3)) done via the regular representation over Q.
import sys, pickle, time, itertools
sys.path.insert(0, '.')
from voa import *
from lat import *
from q4 import Q4
import flint
from collections import defaultdict
t0 = time.time()
KSOL = int(sys.argv[1]) if len(sys.argv) > 1 else 0
ONLY = sys.argv[2] if len(sys.argv) > 2 else 'all'   # 'all' or 'dominant'
S = pickle.load(open('sl2_4.pkl', 'rb'))
G = pickle.load(open('G.pkl', 'rb'))
E, F, H = S['E'], S['F'], S['H']
Gb, omCC = G['Gb'], G['omCC']
sols = pickle.load(open('solutions.pkl', 'rb'))['sols']
xs = sols[KSOL]
# rational components of u and of omega' = omega_CC - u
U = [dict() for _ in range(4)]
for c, g in zip(xs, Gb):
    for j in range(4):
        if c.a[j]:
            sadd(U[j], g, c.a[j])
Wc = [sadd(dict(omCC), U[0], -1), {k: -v for k, v in U[1].items()}, {k: -v for k, v in U[2].items()}, {k: -v for k, v in U[3].items()}]
Wc = [clean(w) for w in Wc]
allW = sorted(set(k for w in Wc for k in w), key=str)
print('omega\' support size', len(allW), ' components sizes', [len(w) for w in Wc])
# multiplication table of K basis (1, s2, s3, s6): kappa_j kappa_l = sum_m MT[j][l][m] kappa_m
MT = [[[0]*4 for _ in range(4)] for _ in range(4)]
for j in range(4):
    for l in range(4):
        a = Q4(*[int(i == j) for i in range(4)]); b = Q4(*[int(i == l) for i in range(4)])
        p = a*b
        MT[j][l] = [x for x in p.a]
Hs = [Hvec(X) for X in range(4)]
R = [a for a in ROOTS if ip(a, Hs[2]) == 0 and ip(a, Hs[3]) == 0]
N = [a for a in NORM4 if ip(a, Hs[2]) == 0 and ip(a, Hs[3]) == 0]
blocks = defaultdict(list)
for i in range(16):
    for j in range(i, 16):
        blocks[(0, 0)].append((((1, i), (1, j)), ZERO))
for i in range(16):
    blocks[(0, 0)].append((((2, i),), ZERO))
for a in R:
    w = (ip(a, Hs[0]), ip(a, Hs[1]))
    for i in range(16):
        blocks[w].append((((1, i),), a))
for a in N:
    blocks[(ip(a, Hs[0]), ip(a, Hs[1]))].append(((), a))
gens = [E[2], F[2], H[2], E[3], F[3], H[3]]
def kernel_rational(cols, nrows_keys):
    """cols: list of dict(rowkey -> Fraction). returns list of kernel vectors (list of Fractions)"""
    import math
    keys = sorted(set(k for c in cols for k in c), key=str)
    if not keys:
        return [[Fr(int(i == j)) for i in range(len(cols))] for j in range(len(cols))]
    ridx = {k: i for i, k in enumerate(keys)}
    den = 1
    for c in cols:
        for v in c.values():
            den = den*v.denominator//math.gcd(den, v.denominator)
    M = flint.fmpz_mat(len(keys), len(cols))
    for j, c in enumerate(cols):
        for k, v in c.items():
            M[ridx[k], j] = int(v*den)
    NS, nul = M.nullspace()
    return [[Fr(int(NS[i, a])) for i in range(len(cols))] for a in range(nul)]
cache = {}
def bmode(s, n, t):
    key = (s, n, t)
    r = cache.get(key)
    if r is None:
        # necessary condition: output weight >= norm of output lattice / 2
        lat = vadd(s[1], t[1])
        if state_wt(s) + state_wt(t) - n - 1 < Fr(sum(x*x for x in lat), 8):
            r = {}
        else:
            r = mode(s, n, {t: Fr(1)})
        cache[key] = r
    return r
total = 0
dims = {}
Ybases = {}
Sbases = {}
wlist = sorted(blocks)
if ONLY == 'dominant':
    wlist = [w for w in wlist if w[0] >= 0 and w[1] >= 0]
for w in wlist:
    basis = blocks[w]
    cols = []
    for bk in basis:
        col = {}
        for gi, x in enumerate(gens):
            for n in (0, 1, 2):
                for k, v in mode_vec(x, n, {bk: Fr(1)}).items():
                    col[(gi, n, k)] = v
        cols.append(col)
    ker = kernel_rational(cols, None)
    Y = [clean({basis[i]: v[i] for i in range(len(basis)) if v[i]}) for v in ker]
    Ybases[w] = Y
    # omega' constraints, regular representation
    big = []   # columns indexed by (c, j)
    Tn = []
    for y in Y:
        T = [dict() for _ in range(4)]   # T[l] = dict((n,key) -> Fraction) : (W_l)_(n) y
        for t, b in y.items():
            for s in allW:
                for n in (1, 2, 3):
                    r = bmode(s, n, t)
                    if not r: continue
                    for l in range(4):
                        a = Wc[l].get(s)
                        if a:
                            for kk, vv in r.items():
                                key = (n, kk)
                                T[l][key] = T[l].get(key, 0) + a*b*vv
        Tn.append(T)
    for c, T in enumerate(Tn):
        for j in range(4):
            col = {}
            for l in range(4):
                for key, val in T[l].items():
                    if not val: continue
                    for m in range(4):
                        f = MT[j][l][m]
                        if f:
                            kk = (key, m)
                            col[kk] = col.get(kk, 0) + f*val
            big.append(clean(col))
    kerK = kernel_rational(big, None)
    assert len(kerK) % 4 == 0
    d = len(kerK)//4
    dims[w] = (len(basis), len(Y), d)
    # K-vectors: lambda_c = sum_j v[4c+j] kappa_j
    Sv = []
    for v in kerK:
        st = {}
        for c, y in enumerate(Y):
            lam = Q4(*[v[4*c+j] for j in range(4)])
            if lam:
                for k, val in y.items():
                    nv = st.get(k, Q4(0)) + lam*val
                    if nv: st[k] = nv
                    else: st.pop(k, None)
        Sv.append(st)
    Sbases[w] = Sv
    total += d
    print('block', w, ' |basis| =', len(basis), ' dim Y =', len(Y), ' dim_K S_2-block =', d, '  (%.1fs)' % (time.time()-t0), flush=True)
print('TOTAL dim_K (S x 1) cap W_2 over computed blocks =', total)
pickle.dump({'dims': dims, 'Sbases': Sbases, 'Ybases': Ybases, 'KSOL': KSOL}, open('s2_space_%d_%s.pkl' % (KSOL, ONLY), 'wb'))
print('time %.1fs' % (time.time()-t0))
