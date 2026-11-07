# Step 5b: checks on S_2 x 1 : total dim, weight multiplicities vs expected decomposition U2 (t12_table.json),
# phi_chi-stability ; positivity of (phi_chi x, y) on the 5-dim slice R_{d4} x w0 of W_1.
import sys, pickle, time, json, itertools, os
sys.path.insert(0, '.')
from voa import *
from lat import *
from q4 import Q4
import flint, math
t0 = time.time()
KSOL = int(sys.argv[1]) if len(sys.argv) > 1 else 0
D = pickle.load(open('s2_space_%d_all.pkl' % KSOL, 'rb'))
dims, Sb = D['dims'], D['Sbases']
tot = sum(d for (_, _, d) in dims.values())
print('u%d : dim_K (S x 1) cap W_2 =' % KSOL, tot)
T = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'data', 't12_table.json')))
U2 = {tuple(map(int, k.split(','))): v for k, v in T['U2'].items()}
expw = {}
for (a, b), m in U2.items():
    for wa in range(-a, a+1, 2):
        for wb in range(-b, b+1, 2):
            expw[(wa, wb)] = expw.get((wa, wb), 0) + m
got = {w: d for w, (_, _, d) in dims.items() if d}
print('weight multiplicities equal those of U2 =', sorted(U2.items()), ':', got == expw)
# phi_chi
import sympy as sp
mu = [0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 1, 0, -1, 0, 0, 0]   # chi(beta) = 2^{<mu,beta>} (from phi_chi.out)
def chi(b):
    e = sum(m*x for m, x in zip(mu, b))
    assert e % 2 == 0
    return Fr(2)**(e//2)
def phi(st):
    return {(m, vneg(b)): (c.conj() if isinstance(c, Q4) else c)*chi(b)*(-1)**len(m) for (m, b), c in st.items()}
def expand(st, keys):
    v = []
    for k in keys:
        c = Q4.lift(st.get(k, 0))
        v.extend(c.a)
    return v
allok = True
for w, vecs in Sb.items():
    mw = (-w[0], -w[1])
    tgt = Sb.get(mw, [])
    imgs = [phi(v) for v in vecs]
    keys = sorted(set(k for v in tgt for k in v) | set(k for v in imgs for k in v), key=str)
    if not vecs:
        continue
    cols = [expand(v, keys) for v in tgt]
    A = flint.fmpq_mat(len(keys)*4, len(cols))
    for j, c in enumerate(cols):
        for i, x in enumerate(c):
            if x: A[i, j] = flint.fmpq(x.numerator, x.denominator)
    r0 = A.rank()
    cols2 = cols + [expand(v, keys) for v in imgs]
    B = flint.fmpq_mat(len(keys)*4, len(cols2))
    for j, c in enumerate(cols2):
        for i, x in enumerate(c):
            if x: B[i, j] = flint.fmpq(x.numerator, x.denominator)
    r1 = B.rank()
    allok &= (r0 == r1)
print('phi_chi maps (S x 1) cap W_2 into itself (all blocks):', allok)
# ---- positivity on the slice R_{d4} x w0 ----
S = pickle.load(open('sl2_4.pkl', 'rb'))
E, F, H, HW = S['E'], S['F'], S['H'], S['HW']
G = pickle.load(open('G.pkl', 'rb'))
sols = pickle.load(open('solutions.pkl', 'rb'))['sols']
u = {}
for c, g in zip(sols[KSOL], G['Gb']):
    for k, v in g.items():
        nv = u.get(k, Q4(0)) + c*v
        if nv: u[k] = nv
        else: u.pop(k, None)
def o(v):  # o(u) v = v_(1) u
    return clean(mode_vec(v, 1, u))
def lin(vs, cs):
    out = {}
    for v, c in zip(vs, cs):
        for k, x in v.items():
            nv = out.get(k, Q4(0)) + c*x
            if nv: out[k] = nv
            else: out.pop(k, None)
    return out
# (0,0,2,0): multiplicity 2, eigenvalues {0, 5/13}; find eigenvector for 5/13
v1, v2 = HW[(0, 0, 2, 0)]
# solve o(u) (a v1 + b v2) = 5/13 (a v1 + b v2): 2x2 system over K
o1, o2 = o(v1), o(v2)
keys = sorted(set(v1) | set(v2) | set(o1) | set(o2), key=str)
# find coordinates of o1, o2 in basis v1, v2 using two independent keys
def coord(vec):
    for k1, k2 in itertools.combinations(keys, 2):
        a11, a12 = Q4.lift(v1.get(k1, 0)), Q4.lift(v2.get(k1, 0))
        a21, a22 = Q4.lift(v1.get(k2, 0)), Q4.lift(v2.get(k2, 0))
        det = a11*a22 - a12*a21
        if det:
            # only rational det handled (v1, v2 rational)
            dq = det.a[0]; assert det == Q4(dq)
            b1, b2 = Q4.lift(vec.get(k1, 0)), Q4.lift(vec.get(k2, 0))
            x = (b1*a22 - a12*b2)/dq; y = (a11*b2 - a21*b1)/dq
            assert lin([v1, v2], [x, y]) == {k: Q4.lift(c) for k, c in vec.items() if c} or all(not c for c in sadd(lin([v1, v2], [x, y]), vec, -1).values())
            return x, y
M = [coord(o1), coord(o2)]   # columns
# M = [[m11, m12],[m21, m22]] with o(v1) = m11 v1 + m21 v2 etc.
m11, m21 = M[0]; m12, m22 = M[1]
lam = Fr(5, 13)
# eigenvector (a, b): (m11 - lam) a + m12 b = 0
if m12:
    a, b = m12, lam - m11
else:
    a, b = lam - m22, m21
x0 = lin([v1, v2], [a, b])
assert all(not c for c in sadd(o(x0), {k: lam*c for k, c in x0.items()}, -1).values()), 'eigenvector check'
y = HW[(1, 1, 2, 0)][0]
oy = o(y)
assert all(not c for c in sadd(oy, {k: Fr(12, 65)*c for k, c in y.items()}, -1).values()), 'eigenvalue 12/65 check on (1,1,2,0)'
slice_ = [x0, y, clean(mode_vec(F[0], 0, y)), clean(mode_vec(F[1], 0, y)), clean(mode_vec(F[0], 0, mode_vec(F[1], 0, y)))]
slice_ = [{k: Q4.lift(c) for k, c in v.items()} for v in slice_]
# check slice vectors are annihilated by E_2, E_3 and have (H_2, H_3) weight (2, 0), o(u)-eigen
for v in slice_:
    assert not clean(mode_vec(E[2], 0, v)) and not clean(mode_vec(E[3], 0, v))
def herm(x, yv):
    r = mode_vec(phi(x), 1, yv)
    return -Q4.lift(r.get(((), ZERO), 0))
Gm = [[herm(a, b) for b in slice_] for a in slice_]
print('Hermitian Gram matrix on the slice R_{d4} x w0 (5x5):')
for row in Gm: print('   ', row)
herm_ok = all(Gm[i][j] == Gm[j][i].conj() for i in range(5) for j in range(5))
print('Hermitian:', herm_ok)
# leading principal minors (exact in K, must be real and > 0)
import sympy as sp
def tosp(q):
    a, b, c, d = q.a
    return sp.Rational(a.numerator, a.denominator) + sp.Rational(b.numerator, b.denominator)*sp.sqrt(-2) + sp.Rational(c.numerator, c.denominator)*sp.sqrt(3) + sp.Rational(d.numerator, d.denominator)*sp.sqrt(-6)
Ms = sp.Matrix(5, 5, lambda i, j: tosp(Gm[i][j]))
minors = [sp.nsimplify(sp.simplify(Ms[:k, :k].det())) for k in range(1, 6)]
print('leading principal minors:', minors)
print('positive definite:', all(sp.re(m) > 0 and sp.simplify(sp.im(m)) == 0 for m in minors))
print('time %.1fs' % (time.time()-t0))
