"""Axis A1 (case II, modulo p3 = 8388571, coherent embedding): bi-weights of the commutant C = Com_S(A) in weight 4.
C_4 = common kernel of y_(k) (y in U_1, k >= 0) on the (0,0)-block of S_4 (checkpoint of closure_w4_modp.py, choice 0).
The two Virasoro vectors of C_2 (Griess subalgebra 9) are omega_i = 2 f_i, f_i the two orthogonal idempotents of
span(e, x) (exact over Q(sqrt13)).  We compute L_i(0) = (omega_i)_(1) on C_4 and its spectrum.
Prediction of the conjecture C = K(sl2,17)^{Z2} (x) GKO(4,13;17): on C_4 (dim 7) the bi-weights are
(4,0) x 3 [d^2 f, f_(-1) f, primary W], (0,4) x 3, (2,2) x 1."""
import sys, os, json, pickle, time
HERE = os.path.abspath(os.path.dirname(__file__)); sys.path.insert(0, os.path.abspath(HERE + '/../generation-forte'))
os.chdir(HERE + '/../generation-forte')
import numpy as np, scipy.sparse as sps, sympy as sp
from flint import nmod_mat, nmod_poly, nmod
from gmpy2 import mpq
import lattice_voa as V, exact_engine as E
t0 = time.time()
def log(*a): print(*a, '(%.0fs)' % (time.time() - t0), flush=True)
S = E.setup(0); T = S['T']; P = E.P; F = E.F
assert P == 8388571
z = np.load('../seconde-realisation/ckpt_S4_0_p8388571.npz')
B4 = E.Blocks(T, 4); mu0 = (0, 0); keys = B4.blk[mu0]; kidx = {k: i for i, k in enumerate(keys)}
rows4 = np.array(z['0_0'], dtype=np.int64) % P            # 280 x 3078
log('S_4 (0,0) block:', rows4.shape)
# --- C_4 basis: kernel of y_(k) ---
outkeys = {}; R_, C_, D_ = [], [], []
for j, kk in enumerate(keys):
    for a, y in enumerate(S['U1']):
        for k in range(0, 5):
            for k2, c in V.mode_vec(y, k, {kk: mpq(1)}).items():
                if c:
                    key = (a, k, k2)
                    if key not in outkeys: outkeys[key] = len(outkeys)
                    R_.append(outkeys[key]); C_.append(j); D_.append(F.q(c))
M = sps.csr_matrix((np.array(D_, dtype=np.int64), (np.array(R_), np.array(C_))), shape=(len(outkeys), len(keys)))
Z = np.asarray((M @ rows4.T)) % P                         # (#out) x 280
Z = Z[np.any(Z, axis=1)]
Zm = nmod_mat(Z.shape[0], Z.shape[1], [int(x) for x in Z.flatten()], P)
Rr, rank = Zm.rref()
piv = []; r = 0
for i in range(rank):
    j = next(j for j in range(Z.shape[1]) if int(Rr[i, j]) != 0); piv.append(j)
free = [j for j in range(Z.shape[1]) if j not in piv]
lam = np.zeros((Z.shape[1], len(free)), dtype=np.int64)
for c, fj in enumerate(free):
    lam[fj, c] = 1
    for i, pj in enumerate(piv): lam[pj, c] = (-int(Rr[i, fj])) % P
C4 = (rows4.T @ lam) % P                                  # 3078 x dimC4
log('dim C_4 =', C4.shape[1])
# --- idempotents of the subalgebra 9 (exact) ---
gb = json.load(open('griess_B_II.json'))
n = len(gb['gram']); prod = {tuple(map(int, k.split(','))): [sp.sympify(c) for c in v] for k, v in gb['prod'].items()}
gram = sp.Matrix([[sp.sympify(c) for c in r] for r in gb['gram']])
G9 = json.load(open('griess_subalgebras.json'))['II']['subalgebras'][9]
xv = sp.Matrix([0] + [sp.sympify(c) for c in G9['x']]); ev = sp.Matrix([1] + [0]*(n - 1))
def mult(a, b):
    out = sp.zeros(n, 1)
    for i in range(n):
        for j in range(n):
            if a[i] != 0 and b[j] != 0: out += a[i]*b[j]*sp.Matrix(prod[(i, j)])
    return out.applyfunc(lambda t: sp.expand(sp.radsimp(t)))
# products in span(e, x): e.e = e (e = omega/2), e.x = x, x.x = alpha e + beta x  (checked exactly)
r13 = sp.sqrt(13)
def nf(t): return sp.expand(sp.radsimp(sp.expand(t)))
ee, ex, xx = mult(ev, ev), mult(ev, xv), mult(xv, xv)
assert (ee - ev).applyfunc(nf) == sp.zeros(n, 1) and (ex - xv).applyfunc(nf) == sp.zeros(n, 1)
i_ = next(i for i in range(1, n) if xv[i] != 0)
beta = nf(xx[i_]/xv[i_]); alpha = nf(xx[0] - beta*xv[0])
assert (xx - alpha*ev - beta*xv).applyfunc(nf) == sp.zeros(n, 1)
# w = a e + b x, b != 0:  a = (1 - b beta)/2 ,  a^2 + alpha b^2 = a  ->  b^2 (beta^2/4 + alpha) = 1/4
Dq = nf(sp.Rational(1, 4)/(beta**2/4 + alpha))
def sqrt13(D):   # exact square root in Q(sqrt13)
    D = nf(D); p_ = D.coeff(r13, 0); q_ = D.coeff(r13, 1)
    u, v = sp.symbols('u v'); sol = sp.solve([u**2 + 13*v**2 - p_, 2*u*v - q_], [u, v], dict=True)
    sol = [s_ for s_ in sol if s_[u].is_Rational and s_[v].is_Rational]
    assert sol, ('not a square in Q(sqrt13)', D); return sol[0][u] + sol[0][v]*r13
b0 = sqrt13(Dq)
idem = []
for b in (b0, -b0):
    a = nf((1 - b*beta)/2); v = (a*ev + b*xv).applyfunc(nf)
    assert (mult(v, v) - v).applyfunc(nf) == sp.zeros(n, 1); idem.append(v)
assert ((idem[0] + idem[1]) - ev).applyfunc(nf) == sp.zeros(n, 1)
log('alpha, beta =', alpha, beta, ' b =', b0)
cs = [sp.nsimplify(sp.radsimp(8*(v.T*gram*v)[0])) for v in idem]      # c = 2 (omega, omega) = 8 (f, f)
log('idempotents of span(e,x):', [list(v) for v in idem], 'central charges', cs)
# --- matrices of b_(1) on the block, then on C_4 ---
Bb = pickle.load(open('griess_B_II_basis.pkl', 'rb'))
def kred(x):  # sympy element of Q(sqrt13) -> F_p
    x = sp.expand(sp.radsimp(x)); q = x.coeff(sp.sqrt(13)) if x.has(sp.sqrt(13)) else 0; p_ = sp.expand(x - q*sp.sqrt(13))
    return (F.q(mpq(int(sp.Rational(p_).p), int(sp.Rational(p_).q))) + F.q(mpq(int(sp.Rational(q).p), int(sp.Rational(q).q)))*F.s13) % P
Xb = []
for b in Bb:
    R_, C_, D_ = [], [], []
    for j, kk in enumerate(keys):
        for k2, c in V.mode_vec(b, 1, {kk: mpq(1)}).items():
            if c:
                assert k2 in kidx
                R_.append(kidx[k2]); C_.append(j); D_.append(F.q(c))
    Xb.append(sps.csr_matrix((np.array(D_, dtype=np.int64), (np.array(R_), np.array(C_))), shape=(len(keys), len(keys))))
log('b_(1) matrices built')
# coordinates in C_4 : solve C4 * X = Y  (pivot rows)
C4m = nmod_mat(C4.shape[0], C4.shape[1], [int(x) for x in C4.flatten()], P)
def coords(Y):
    A = np.hstack([C4, Y]) % P
    Am = nmod_mat(A.shape[0], A.shape[1], [int(x) for x in A.flatten()], P)
    Rr, rk = Am.rref()
    d = C4.shape[1]; assert rk == d, ('image not in C_4', rk)
    return np.array([[int(Rr[i, d + j]) for j in range(Y.shape[1])] for i in range(d)], dtype=np.int64)
Y0 = np.asarray(Xb[0] @ C4) % P
log('DIAG e_(1) C4 == 2 C4 ?', bool(not ((Y0 - 2*C4) % P).any()))
for i in range(1, n):
    Yi = np.asarray(Xb[i] @ C4) % P
    A = np.hstack([C4, Yi]) % P; Am = nmod_mat(A.shape[0], A.shape[1], [int(x) for x in A.flatten()], P)
    log('DIAG rank [C4 | b_%d(1) C4] =' % i, Am.rank())
for sg in (1, -1):
    co = []
    for t in xv:
        t = sp.expand(t); q_ = t.coeff(sp.sqrt(13), 1); p_ = t.coeff(sp.sqrt(13), 0)
        co.append((F.q(mpq(int(sp.Rational(p_).p), int(sp.Rational(p_).q))) + sg*F.q(mpq(int(sp.Rational(q_).p), int(sp.Rational(q_).q)))*F.s13) % P)
    Yx = sum((co[i]*(np.asarray(Xb[i] @ C4) % P) for i in range(n)), np.zeros(C4.shape, dtype=np.int64)) % P
    A = np.hstack([C4, Yx]) % P; Am = nmod_mat(A.shape[0], A.shape[1], [int(x) for x in A.flatten()], P)
    log('DIAG sign', sg, 'rank [C4 | x_(1) C4] =', Am.rank())
res = {'dim_C4': int(C4.shape[1]), 'central_charges': [str(c) for c in cs]}
Lmats = []
for v, c in zip(idem, cs):
    co = [2*kred(t) % P for t in v]                       # omega_i = 2 f_i
    Y = sum((co[i]*(np.asarray(Xb[i] @ C4) % P) for i in range(n)), np.zeros(C4.shape, dtype=np.int64)) % P
    L0 = coords(Y); Lmats.append(L0)
    Lm = nmod_mat(L0.shape[0], L0.shape[1], [int(x) for x in L0.flatten()], P)
    cp = Lm.charpoly()
    roots = [(int(r), m) for r, m in cp.roots()] if hasattr(cp, 'roots') else None
    def recon(a):
        for den in range(1, 2000):
            num = a*den % P
            if num < 10**5: return '%d/%d' % (num, den)
            if P - num < 10**5: return '-%d/%d' % (P - num, den)
        return str(a)
    res['L0 spectrum (c=%s)' % c] = [(recon(r), m) for r, m in roots] if roots is not None else str(cp)
    log('c =', c, 'L(0) on C_4: char poly roots', res['L0 spectrum (c=%s)' % c])
# commutation and joint spectrum
L1, L2 = Lmats
comm = (L1 @ L2 - L2 @ L1) % P
res['[L1(0), L2(0)] = 0'] = bool(not comm.any())
S_ = (L1 + L2) % P
res['trace(L1+L2) mod p'] = int(np.trace(S_) % P)
log(res)
json.dump(res, open(HERE + '/c1_biweights.json', 'w'), indent=1)
np.savez_compressed(HERE + '/c1_C4_basis_p8388571.npz', C4=C4, L1=L1, L2=L2)
