"""Axis B (independent check of Riesen's algebra U^R = <A, e, f, v>, case (I)), modulo a prime p.
Riesen's generators e, f (Griess subalgebra 5) and v (type (1,13), isotropic line 1) transported into our lattice
conventions (Calculs/cas-I/r1_riesen_transported.pkl, coefficients in Q(sqrt13, sqrt3, i)); A = currents of case (I)
(identified exactly with Riesen's, r1).  Reduction: p < 2^23 with p = 1 mod 12 and (13/p) = 1; i, sqrt3, sqrt13 in F_p,
sqrt39 = sqrt13 sqrt3 (a ring homomorphism Q(sqrt13, sqrt3, i) -> F_p; dimensions do not depend on the choice).
 W_2 : smallest g-stable subspace of V'_2 containing A_2, e, f, v and closed under the Griess products b_(1) w
       (b highest weight vectors of W_2 suffice, by g-equivariance).  Expected 120; C_2 = Com(A) in W_2 expected = span(e, f).
 S_3 : y(-1)W_2, y(-2)U_1, L(-1)W_2, u_0 W_2, u_{-1}U_1 (u highest weight vectors of W_2 outside A_2), g-closed (956?).
 S_4 : y(-1)S_3, y(-2)W_2, y(-3)U_1, L(-1)S_3, u_{-1}W_2, u_{-2}U_1, u_0 S_3, g-closed (6010?); tests y_1 S_4 c S_3,
       y_2 S_4 c W_2.  Same families as closure_w4_modp.py (case II).  Ranks modulo p are lower bounds over the field.
Usage (from this folder): PSTART=8388608 [SAVEALL=dir] [W3P=1] [NOW4=1] [VCHOICE=conj|rand] python3 b1_UR_closure_modp.py"""
import sys, os, json, time, gc, pickle
HERE = os.path.abspath(os.path.dirname(__file__)); GF = os.path.abspath(HERE + '/../cas-II')
sys.path.insert(0, GF); os.chdir(GF)
import numpy as np
import scipy.sparse as sps
import sympy as sp
from flint import nmod_mat
from gmpy2 import mpq
import lattice_voa as V, cases as C, commutant_griess as CG
t0 = time.time()
def log(*a):
    print(*a, '(%.0fs)' % (time.time() - t0), flush=True)
def find_prime(start):
    p = sp.prevprime(start)
    while not (p % 12 == 1 and sp.legendre_symbol(13, p) == 1): p = sp.prevprime(p)
    return p
P = find_prime(int(os.environ.get('PSTART', 2**23))); assert P < 2**23
SGN = int(os.environ.get('SGN', 1))
I_ = sp.sqrt_mod(P - 1, P); S3_ = sp.sqrt_mod(3, P); S13_ = sp.sqrt_mod(13, P)
if SGN < 0: S13_ = P - S13_          # other embedding (control)
class _F:
    p = P
    def q(self, x):
        n, d = (int(x.numerator), int(x.denominator)) if hasattr(x, 'numerator') else (int(x.p), int(x.q))
        return n*pow(d, -1, P) % P
F = _F()
R13, R3, II = sp.symbols('R13 R3 II')
def red(c):
    e = sp.expand(sp.sympify(c)).xreplace({sp.sqrt(39): R13*R3, sp.sqrt(13): R13, sp.sqrt(3): R3, sp.I: II})
    e = sp.expand(e); Pl = sp.Poly(e, R13, R3, II); tot = 0
    for (a, b, c_), co in Pl.terms():
        co = sp.Rational(co)
        tot += F.q(co)*pow(S13_, a, P)*pow(S3_, b, P)*pow(I_, c_, P)
    return tot % P
T = C.build('I'); U1 = CG.U1_basis(T)
A2 = []
for x in U1:
    A2.append(V.clean(V.mode_vec(x, -2, {V.VAC: mpq(1)})))
    for y in U1: A2.append(V.clean(V.mode_vec(x, -1, y)))
FULL = pickle.load(open(HERE + '/../cas-I/r1_riesen_transported.pkl', 'rb'))
GEN = {n: [(red(c), {k: mpq(1)}) for k, c in FULL[n].items()] for n in ('e', 'f', 'v')}
GEN = {n: [(c, d) for c, d in g if c] for n, g in GEN.items()}
VCH = os.environ.get('VCHOICE', '')      # controls: 'conj' = the other isotropic line (complex conjugate of v);
if VCH:                                  # 'rand' = v + 3 conj(v), a non-isotropic vector of the admissible plane E
    conjv = {k: sp.conjugate(sp.sympify(c)) for k, c in FULL['v'].items()}
    comb = {}
    for k, c in FULL['v'].items():
        comb[k] = (0 if VCH == 'conj' else red(c))
    for k, c in conjv.items():
        comb[k] = (comb.get(k, 0) + (1 if VCH == 'conj' else 3)*red(c)) % P
    GEN['v'] = [(c, {k: mpq(1)}) for k, c in comb.items() if c]
log('p =', P, ' sqrt13 =', S13_, ' sqrt3 =', S3_, ' i =', I_, '; terms of e, f, v:', {n: len(g) for n, g in GEN.items()})
KEYS = {w: V.basis(w) for w in range(5)}
IDX = {w: {k: i for i, k in enumerate(KEYS[w])} for w in range(5)}
WT = {w: [C.weight_of_key(T, k) for k in KEYS[w]] for w in range(5)}
BLK = {}
for w in range(5):
    b = {}
    for i, mu in enumerate(WT[w]): b.setdefault(mu, []).append(i)
    BLK[w] = {mu: np.array(c) for mu, c in b.items()}
POS = {w: {} for w in range(5)}          # key index -> position inside its block
for w in range(5):
    for mu, c in BLK[w].items():
        for j, i in enumerate(c): POS[w][int(i)] = j
def comb_to_blocks(comb, w):
    """comb = [(coef, rational dict)] (a weight vector) -> (mu, array over the block)"""
    z = {}
    for c, v in comb:
        for k, x in v.items():
            i = IDX[w][k]; z[i] = (z.get(i, 0) + c*F.q(x)) % P
    z = {i: x for i, x in z.items() if x}
    if not z: return None
    mu = WT[w][next(iter(z))]
    a = np.zeros(len(BLK[w][mu]), dtype=np.int64)
    for i, x in z.items():
        assert WT[w][i] == mu
        a[POS[w][i]] = x
    return mu, a
def rref_rows(rows, ncols):
    """rows: list of int64 arrays -> reduced basis (2D int64 array)"""
    if not rows: return np.zeros((0, ncols), dtype=np.int64)
    M = nmod_mat(len(rows), ncols, [int(x) for r in rows for x in r], P)
    R, rk = M.rref()
    return np.array([[int(R[i, j]) for j in range(ncols)] for i in range(rk)], dtype=np.int64).reshape(rk, ncols)
class Space:
    def __init__(self, w): self.w = w; self.B = {}
    def dim(self): return sum(b.shape[0] for b in self.B.values())
    def add(self, new):            # new: dict mu -> list of arrays
        grown = False
        for mu, rows in new.items():
            n = len(BLK[self.w][mu]); old = self.B.get(mu, np.zeros((0, n), dtype=np.int64))
            nb = rref_rows(list(old) + rows, n)
            if nb.shape[0] > old.shape[0]: grown = True
            self.B[mu] = nb
        return grown
    def contains(self, mu, rows):
        """number of rows (dict-free) not in the span"""
        n = len(BLK[self.w][mu]); old = self.B.get(mu, np.zeros((0, n), dtype=np.int64))
        r0 = old.shape[0]
        return rref_rows(list(old) + rows, n).shape[0] - r0
def op_matrix(comb, n, ws, wt):
    """sparse (rows keys_wt, cols keys_ws) matrix of sum c v_n, comb = [(c, rational dict)]"""
    R_, C_, D_ = [], [], []
    for c, v in comb:
        for key, x0 in v.items():
            rows, cols, vals = [], [], []
            cx = c*F.q(x0) % P
            for j, k in enumerate(KEYS[ws]):
                r = V.mode_vec({key: mpq(1)}, n, {k: mpq(1)})
                for kk, x in r.items():
                    if x: rows.append(IDX[wt][kk]); cols.append(j); vals.append(cx*F.q(x) % P)
            V._MODE.clear(); V._EMODE.clear()
            R_.append(np.array(rows, dtype=np.int64)); C_.append(np.array(cols, dtype=np.int64)); D_.append(np.array(vals, dtype=np.int64))
    M = sps.csc_matrix((np.concatenate(D_) if D_ else np.zeros(0, dtype=np.int64),
                        (np.concatenate(R_) if R_ else np.zeros(0, dtype=np.int64), np.concatenate(C_) if C_ else np.zeros(0, dtype=np.int64))),
                       shape=(len(KEYS[wt]), len(KEYS[ws])))
    M.sum_duplicates(); M.data %= P; M.eliminate_zeros()
    return M
def apply_op(M, space, wt, src=None):
    """images of all basis vectors of `space` (or of src: dict mu -> 2D array) under M -> dict mu -> list of arrays"""
    out = {}
    S = space.B if src is None else src
    for mu, Bm in S.items():
        if Bm.shape[0] == 0: continue
        sub = M[:, BLK[space.w][mu]]
        Z = (sub @ Bm.T) % P                      # (n_wt, r)
        nz = np.nonzero(np.any(Z, axis=1))[0]
        if len(nz) == 0: continue
        for nu in set(WT[wt][i] for i in nz):
            cols = BLK[wt][nu]
            Zb = Z[cols, :]
            for j in range(Zb.shape[1]):
                if np.any(Zb[:, j]): out.setdefault(nu, []).append(Zb[:, j].copy())
    return out
def merge(a, b):
    for k, v in b.items(): a.setdefault(k, []).extend(v)
    return a
G_OPS = {}
def g_ops(w):
    if w not in G_OPS: G_OPS[w] = [op_matrix([(1, T[k])], 0, w, w) for k in ('f1', 'f2', 'e1', 'e2')]
    return G_OPS[w]
def g_close(space):
    for it in range(30):
        new = {}
        for M in g_ops(space.w): merge(new, apply_op(M, space, space.w))
        if not space.add(new): break
    return space
def mult_decomp(space):
    per = {mu: b.shape[0] for mu, b in space.B.items()}
    rr = lambda a, b: per.get((a, b), 0)
    m = {f'{a},{b}': rr(a, b) - rr(a + 2, b) - rr(a, b + 2) + rr(a + 2, b + 2) for (a, b) in per if a >= 0 and b >= 0}
    return {k: v for k, v in m.items() if v}
res = {'p': P, 'sqrt13': int(S13_), 'sqrt3': int(S3_), 'i': int(I_)}
S1 = Space(1); tmp1 = {}
for u in U1:
    mu, a = comb_to_blocks([(1, u)], 1); tmp1.setdefault(mu, []).append(a)
S1.add(tmp1)
def gvec_to_comb(mu, a, w):
    return [(int(a[j]), {KEYS[w][int(i)]: mpq(1)}) for j, i in enumerate(BLK[w][mu]) if a[j]]
def hw_space(space):
    """highest weight vectors (killed by e1_(0), e2_(0)) of a g-stable space, by block: dict mu -> rows"""
    Me = [op_matrix([(1, T[k])], 0, space.w, space.w) for k in ('e1', 'e2')]
    out = {}
    for mu, Bm in space.B.items():
        if Bm.shape[0] == 0 or mu[0] < 0 or mu[1] < 0: continue
        cols = BLK[space.w][mu]
        imgs = [((M[:, cols] @ Bm.T) % P) for M in Me]
        Z = np.vstack(imgs)                   # (2 n_w, r)
        Zm = nmod_mat(Z.shape[0], Z.shape[1], [int(x) for x in Z.flatten()], P)
        X, nul = Zm.nullspace()
        if nul:
            lam = np.array([[int(X[i, j]) for j in range(nul)] for i in range(Z.shape[1])], dtype=np.int64)
            out[mu] = (lam.T @ Bm) % P
    return out
# ---------- weight 2: closure
A2s = Space(2); tmp = {}
for a in A2:
    r = comb_to_blocks([(1, a)], 2)
    if r: tmp.setdefault(r[0], []).append(r[1])
A2s.add(tmp); log('dim A_2 =', A2s.dim(), '(expected 22: level one factor, e\'(-1)e\' = 0)'); res['dimA2'] = A2s.dim()
W2 = Space(2); W2.add({mu: list(b) for mu, b in A2s.B.items()})
tmp = {}
for n_, g in GEN.items():
    r = comb_to_blocks(g, 2); tmp.setdefault(r[0], []).append(r[1])
W2.add(tmp); g_close(W2); log('A_2 + g.{e, f, v}:', W2.dim())
KM = {}
def key_op(i):
    if i not in KM: KM[i] = op_matrix([(1, {KEYS[2][i]: mpq(1)})], 1, 2, 2).tocsr()
    return KM[i]
def M_of(mu, a):
    M = None
    for j, i in enumerate(BLK[2][mu]):
        if a[j]:
            X = key_op(int(i)) * int(a[j])
            M = X if M is None else M + X
    M = M.tocsc(); M.data %= P; M.eliminate_zeros()
    return M
for it in range(10):
    H = hw_space(W2); nh = sum(b.shape[0] for b in H.values())
    new = {}
    for mu, rows in H.items():
        for a in rows: merge(new, apply_op(M_of(mu, a), W2, 2))
    grown = W2.add(new); g_close(W2)
    log('  Griess closure, iteration %d: %d highest weight vectors, dim W_2 = %d' % (it, nh, W2.dim()))
    if not grown: break
res['dimW2'] = W2.dim(); res['W2_decomposition'] = mult_decomp(W2)
log('dim W_2 =', W2.dim(), ' (expected 120); g-decomposition', res['W2_decomposition'])
# C_2 = Com(A) in W_2: block (0,0) vectors killed by y_(0), y_(1)
mu0 = (0, 0); B0 = W2.B[mu0]
Zs = []
for y in U1:
    for k in (0, 1):
        M = op_matrix([(1, y)], k, 2, 2 - k); Zs.append((M[:, BLK[2][mu0]] @ B0.T) % P)
Z = np.vstack(Zs); Zm = nmod_mat(Z.shape[0], Z.shape[1], [int(x) for x in Z.flatten()], P)
X, nul = Zm.nullspace(); res['dimC2'] = nul
# is span(e, f) = C_2 ?  (rank of [C_2 ; e ; f])
C2rows = [(np.array([int(X[i, j]) for i in range(Z.shape[1])], dtype=np.int64) @ B0) % P for j in range(nul)]
ef = [comb_to_blocks(GEN[n_], 2)[1] for n_ in ('e', 'f')]
rk = nmod_mat(nul + 2, len(BLK[2][mu0]), [int(x) for r in C2rows + ef for x in r], P).rank()
res['C2_equals_span_ef'] = bool(rk == nul == 2)
log('dim C_2 =', nul, '; C_2 = span(e, f):', res['C2_equals_span_ef'])
HW2 = hw_space(W2); HA = hw_space(A2s)
hw = []
for mu, rows in HW2.items():
    base = list(HA.get(mu, np.zeros((0, len(BLK[2][mu])), dtype=np.int64)))
    r0 = len(rref_rows(base, len(BLK[2][mu]))) if base else 0
    for a in rows:
        r1 = rref_rows(base + [a], len(BLK[2][mu])).shape[0]
        if r1 > r0: base.append(a); r0 = r1; hw.append(('%s#%d' % (mu, len(hw)), gvec_to_comb(mu, a, 2)))
log('highest weight vectors of W_2 outside A_2:', len(hw), [h[0] for h in hw]); res['hw_outside_A2'] = [h[0] for h in hw]
json.dump(res, open(HERE + '/b1_UR_closure_p%d%s.json' % (P, ('_' + VCH) if VCH else ''), 'w'), indent=1)
KM.clear(); gc.collect()
S2 = W2
# ---------- weight 3
S3 = Space(3); gens = {}
for x in U1:
    merge(gens, apply_op(op_matrix([(1, x)], -1, 2, 3), S2, 3))
    merge(gens, apply_op(op_matrix([(1, x)], -2, 1, 3), S1, 3))
merge(gens, apply_op(op_matrix([(1, V.OMEGA)], 0, 2, 3), S2, 3))
for name, u in hw:
    merge(gens, apply_op(op_matrix(u, 0, 2, 3), S2, 3))
    merge(gens, apply_op(op_matrix(u, -1, 1, 3), S1, 3))
S3.add(gens); del gens; gc.collect()
g_close(S3); res['dimS3'] = S3.dim(); res['S3_decomposition'] = mult_decomp(S3)
log('dim S_3 (generated, g-closed) =', S3.dim(), ' (expected 956)')
json.dump(res, open(HERE + '/b1_UR_closure_p%d%s.json' % (P, ('_' + VCH) if VCH else ''), 'w'), indent=1)
if os.environ.get('W3P'):              # coherence test (W3'): u_(1) S_3 c S_3 for u = e, f, v
    w3p = {}
    for nm in ('e', 'f', 'v'):
        img = apply_op(op_matrix(GEN[nm], 1, 3, 3), S3, 3)
        w3p[nm] = sum(S3.contains(mu, rows) for mu, rows in img.items())
        log("(W3') extra rank of %s_(1) S_3 over S_3:" % nm, w3p[nm])
    res['W3prime_u1S3_extra'] = w3p
    json.dump(res, open(HERE + '/b1_UR_closure_p%d_W3prime.json' % P, 'w'), indent=1)
if os.environ.get('NOW4'): sys.exit(0)
# ---------- weight 4
S4 = Space(4)
def add_gens(g):
    S4.add(g); gc.collect()
for x in U1:
    add_gens(apply_op(op_matrix([(1, x)], -1, 3, 4), S3, 4))
    add_gens(apply_op(op_matrix([(1, x)], -2, 2, 4), S2, 4))
    add_gens(apply_op(op_matrix([(1, x)], -3, 1, 4), S1, 4))
add_gens(apply_op(op_matrix([(1, V.OMEGA)], 0, 3, 4), S3, 4))
log('S_4: affine and L(-1) part:', S4.dim())
for name, u in hw:
    add_gens(apply_op(op_matrix(u, -1, 2, 4), S2, 4))
    add_gens(apply_op(op_matrix(u, -2, 1, 4), S1, 4))
    add_gens(apply_op(op_matrix(u, 0, 3, 4), S3, 4))
    log('  after u =', name, ':', S4.dim())
g_close(S4)
res['dimS4'] = S4.dim(); res['S4_decomposition'] = mult_decomp(S4)
log('dim S_4 (generated, g-closed) =', S4.dim(), ' (expected 6010)')
json.dump(res, open(HERE + '/b1_UR_closure_p%d%s.json' % (P, ('_' + VCH) if VCH else ''), 'w'), indent=1)
w4 = {}
for i, y in enumerate(U1):
    a = sum(S3.contains(mu, r) for mu, r in apply_op(op_matrix([(1, y)], 1, 4, 3), S4, 3).items())
    b = sum(S2.contains(mu, r) for mu, r in apply_op(op_matrix([(1, y)], 2, 4, 2), S4, 2).items())
    w4['U1_%d' % i] = (a, b)
res['W4_tests_extra_ranks'] = w4
log('(W4) y_1 S_4 c S_3, y_2 S_4 c W_2 (extra ranks):', w4)
np.savez_compressed(HERE + '/b1_S4_00_p%d.npz' % P, S4_00=S4.B[(0, 0)])
if os.environ.get('SAVEALL'):          # all blocks of W_2, S_3, S_4 (large: scratch only)
    od = os.environ['SAVEALL']; os.makedirs(od, exist_ok=True)
    for w, Sp in ((2, S2), (3, S3), (4, S4)):
        np.savez_compressed(od + '/b1_S%d_allblocks_p%d.npz' % (w, P), **{'%d_%d' % mu: b for mu, b in Sp.B.items()})
json.dump(res, open(HERE + '/b1_UR_closure_p%d%s.json' % (P, ('_' + VCH) if VCH else ''), 'w'), indent=1)
log('done')
