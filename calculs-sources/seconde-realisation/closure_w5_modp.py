"""Weight-five obstruction test modulo p (analogue of [Rie26, (4.2)]): with S_2 = U_2, S_3, S_4 the constructed spaces of
closure_w4_modp.py, put N_5 = T S_4 + H(-1) S_4 + (S_2)(-1) S_3 (closed under g' + g''). We test u_(0) S_4 c N_5 for u in
{e, x, highest weight vectors of the eight types}; by zero-mode covariance and the normal-product expansion of the
modes of A_2 this gives (S_2)(0) S_4 c N_5.  We also compare dim N_5 with dim U_5 = 29512 (Haagerup vacuum character).
Usage: PSTART=8388608 python3 closure_w5_modp.py CHOICE"""
import sys, os, json, time, gc
import numpy as np
import scipy.sparse as sps
import sympy as sp
from flint import nmod_mat
from gmpy2 import mpq
import lattice_voa as V, cases as C, closure_explore as CE, griess_subalgebras as GS, exact_isotropy_s13 as ES, modp
_R, _VE = V.lattice_vectors(10)
V.VECS[:] = _VE
for _v in _VE: V.NORM[_v] = V.ip(_v, _v)
F = modp.Fp(int(os.environ.get('PSTART', 2**23))); P = F.p
assert P < 2**23
t0 = time.time()
def log(*a):
    print(*a, '(%.0fs)' % (time.time() - t0), flush=True)
chi = int(sys.argv[1]); TEST4 = sys.argv[2].split(';') if len(sys.argv) > 2 else ['x']
T, D, out, types, N1, N2, s13f, U1, A2 = CE.setup('II', 9)
info = json.load(open('closure_explore_II_9.json')); Lj = json.load(open('isotropic_lines_II_9.json'))
G = json.load(open('griess_subalgebras.json'))['II']['subalgebras'][9]
_, basis, _, _ = GS.griess_data('II')
tl = sorted(types)
choice = info['combos'][chi] if chi >= 0 else [0]*len(tl)
xs = [ES.split(c) for c in G['x']]
x_p = {}; x_q = {}
for i, (p_, q_) in enumerate(xs):
    if p_: V.vadd(x_p, basis[i + 1], mpq(int(p_.p), int(p_.q)))
    if q_: V.vadd(x_q, basis[i + 1], mpq(int(q_.p), int(q_.q)))
x_p, x_q = V.clean(x_p), V.clean(x_q)
lines = {tau: [[F.alg(c) for c in sp.sympify(l)] for l in Lj[str(tau)][0]['lines']] for tau in tl}
U2 = [[(1, a)] for a in A2 if a] + [[(1, basis[0])], [(1, x_p), (F.s13, x_q)]]
hw = [('e', U2[-2]), ('x', U2[-1])]
for tau, l in zip(tl, choice):
    ty = types[tau]; lv = lines[tau][l]; ml = []
    for w in range(len(ty['mods'][0][0])):
        comb = []
        for k in range(ty['d']):
            if ty['mods'][k][0][w]: comb.append((lv[k], ty['mods'][k][0][w]))
            if ty['mods'][k][1][w]: comb.append((lv[k]*F.s13 % P, ty['mods'][k][1][w]))
        ml.append(comb); U2.append(comb)
    hw.append((str(tau), ml[0]))
HW = dict(hw)
log('p =', P, 'choice', choice, '; terms of highest weight vectors:', {n: sum(len(v) for _, v in u) for n, u in hw})
KEYS = {w: V.basis(w) for w in range(6)}
assert len(KEYS[5]) == 352480, len(KEYS[5])
IDX = {w: {k: i for i, k in enumerate(KEYS[w])} for w in range(6)}
WT = {w: [C.weight_of_key(T, k) for k in KEYS[w]] for w in range(6)}
BLK = {}
for w in range(6):
    b = {}
    for i, mu in enumerate(WT[w]): b.setdefault(mu, []).append(i)
    BLK[w] = {mu: np.array(c) for mu, c in b.items()}
POS = {w: {} for w in range(6)}          # key index -> position inside its block
for w in range(6):
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
    def add(self, new):            # new: dict mu -> list of arrays (processed in chunks to bound memory)
        grown = False
        for mu in list(new.keys()):
            rows = new.pop(mu)
            n = len(BLK[self.w][mu])
            for c0 in range(0, len(rows), 300):
                old = self.B.get(mu, np.zeros((0, n), dtype=np.int64))
                nb = rref_rows(list(old) + rows[c0:c0 + 300], n)
                if nb.shape[0] > old.shape[0]: grown = True
                self.B[mu] = nb
            del rows
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
    """images of all basis vectors of `space` (or of src: dict mu -> 2D array) under M -> dict mu -> list of arrays.
    Memory-light: only the nonzero rows of each column slice are multiplied."""
    out = {}
    S = space.B if src is None else src
    for mu, Bm in S.items():
        if Bm.shape[0] == 0: continue
        sub = M[:, BLK[space.w][mu]]
        rows_nz = np.unique(sub.indices)
        if len(rows_nz) == 0: continue
        Z = (sub[rows_nz, :] @ Bm.T) % P              # (len(rows_nz), r), dense, small
        if not isinstance(Z, np.ndarray): Z = Z.toarray() if hasattr(Z, 'toarray') else np.asarray(Z)
        Z = np.asarray(Z) % P
        byblock = {}
        for idx_, i in enumerate(rows_nz): byblock.setdefault(WT[wt][i], []).append((POS[wt][int(i)], idx_))
        for nu, lst in byblock.items():
            n = len(BLK[wt][nu]); Zb = np.zeros((n, Z.shape[1]), dtype=np.int64)
            for pos, idx_ in lst: Zb[pos, :] = Z[idx_, :]
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
res = {'p': P, 'choice': choice}
# ---------- weight 1, 2
S1 = Space(1); tmp1 = {}
for u in U1:
    mu, a = comb_to_blocks([(1, u)], 1); tmp1.setdefault(mu, []).append(a)
S1.add(tmp1)
S2 = Space(2); tmp = {}
for c in U2:
    r = comb_to_blocks(c, 2)
    if r: tmp.setdefault(r[0], []).append(r[1])
S2.add(tmp)
log('dim U_1 =', S1.dim(), '; dim U_2 =', S2.dim())
# ---------- weight 3
S3 = Space(3); gens = {}
for x in U1:
    merge(gens, apply_op(op_matrix([(1, x)], -1, 2, 3), S2, 3))
    merge(gens, apply_op(op_matrix([(1, x)], -2, 1, 3), S1, 3))
merge(gens, apply_op(op_matrix([(1, V.OMEGA)], 0, 2, 3), S2, 3))
OPS0_23 = {}
for name, u in hw:
    OPS0_23[name] = op_matrix(u, 0, 2, 3); merge(gens, apply_op(OPS0_23[name], S2, 3))
    merge(gens, apply_op(op_matrix(u, -1, 1, 3), S1, 3))
S3.add(gens); del gens; gc.collect()
log('S_3 before g-closure:', S3.dim())
g_close(S3); res['dimS3'] = S3.dim(); res['S3_decomposition'] = mult_decomp(S3)
log('dim S_3 (generated, g-closed) =', S3.dim(), ' (expected 956)')
res = {'p': P, 'choice': choice, 'dimS3': S3.dim()}
# ---------- weight 4
S4 = Space(4)
_ck = f'ckpt_S4_{chi}_p{P}.npz'
if os.path.exists(_ck):
    _z = np.load(_ck)
    for _k in _z.files:
        _a, _b = map(int, _k.split('_')); S4.B[(_a, _b)] = _z[_k]
    log('S_4 loaded from checkpoint:', S4.dim())
def add_gens(g):
    S4.add(g); gc.collect()
if S4.dim() == 0:
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
np.savez_compressed(f'ckpt_S4_{chi}_p{P}.npz', **{f'{a}_{b}': B_ for (a, b), B_ in S4.B.items()})
log('dim S_4 (generated, g-closed) =', S4.dim(), ' (expected 6010)')
json.dump(res, open(f'closure_w5_modp_{chi}_p{P}.json', 'w'), indent=1)
# ---------- weight 5
del G_OPS[3]; gc.collect()
N5 = Space(5)
N5.add(apply_op(op_matrix([(1, V.OMEGA)], 0, 4, 5), S4, 5)); gc.collect()
log('N_5: T S_4 :', N5.dim())
for x in U1:
    N5.add(apply_op(op_matrix([(1, x)], -1, 4, 5), S4, 5)); gc.collect()
log('N_5: + H(-1) S_4 :', N5.dim())
for name, u in hw:
    N5.add(apply_op(op_matrix(u, -1, 3, 5), S3, 5)); gc.collect()
    log('  N_5: + u(-1) S_3, u =', name, ':', N5.dim())
del G_OPS[4]; gc.collect()
np.savez_compressed(f'ckpt_N5pre_{chi}_p{P}.npz', **{f'{a}_{b}': B_ for (a, b), B_ in N5.B.items()})
g_close(N5)
res['dimN5'] = N5.dim(); res['N5_decomposition'] = mult_decomp(N5)
log('dim N_5 (g-closed) =', N5.dim(), ' (dim U_5 = 29512)')
json.dump(res, open(f'closure_w5_modp_{chi}_p{P}.json', 'w'), indent=1)
del G_OPS[5]; gc.collect()
# test u_(0) S_4 c N_5, key by key in chunks
t5 = {}
for name, u in hw:
    terms = [(c, {k: x}) for c, v in u for k, x in v.items()]
    acc = {}
    for i in range(0, len(terms), 10):
        M = op_matrix(terms[i:i + 10], 0, 4, 5)
        for mu, Bm in S4.B.items():
            if Bm.shape[0] == 0: continue
            sub = M[:, BLK[4][mu]]
            rows_nz = np.unique(sub.indices)
            if len(rows_nz) == 0: continue
            Z = np.asarray((sub[rows_nz, :] @ Bm.T)) % P
            byblock = {}
            for idx_, i in enumerate(rows_nz): byblock.setdefault(WT[5][i], []).append((POS[5][int(i)], idx_))
            for nu, lst in byblock.items():
                key_ = (mu, nu)
                if key_ not in acc: acc[key_] = np.zeros((len(BLK[5][nu]), Bm.shape[0]), dtype=np.int64)
                A_ = acc[key_]
                for pos, idx_ in lst: A_[pos, :] = (A_[pos, :] + Z[idx_, :]) % P
        del M; gc.collect()
    extra = 0
    for (mu, nu), Zb in acc.items():
        rows = [Zb[:, j].copy() for j in range(Zb.shape[1]) if np.any(Zb[:, j])]
        if rows: extra += N5.contains(nu, rows)
    t5[name] = int(extra)
    log('(4.2) u =', name, ': extra rank of u_(0) S_4 over N_5 =', extra)
    res['W5_test'] = t5
    json.dump(res, open(f'closure_w5_modp_{chi}_p{P}.json', 'w'), indent=1)
log('done')
