"""Exact engine over K = Q(sqrt13, sqrt-3) for the case (II) algebra S (subalgebra 9, choice 0).
K-vectors: dict key -> K element (quadruple of mpq in the basis 1, r, q, rq; r = sqrt13, q = sqrt-3).
Exact operators: dict source key -> K-vector (image of the basis state).  Membership u in span_K(B) is decided
blockwise (Cartan weights of g' + g'') through the restriction of scalars K = Q^4: pivot rows found modulo p
(nonzero minor => independence over Q), square system solved exactly (fmpz_mat.solve -> fmpq_mat), and the
identity B C = X then VERIFIED EXACTLY on all rows."""
import json, time, gc, sys
from gmpy2 import mpq
import sympy as sp
from flint import fmpz_mat, fmpq_mat, nmod_mat, fmpq
import lattice_voa as V, cases as C, closure_explore as CE, griess_subalgebras as GS, Kfield as K
import exact_isotropy_s13 as ES, modp
BAS = [K.ONE, K.R, K.Q, K.RQ]
F = modp.Fp(2**23); P = F.p
def kred(a):            # K -> F_p, coherent embedding r -> s13, q -> s3
    s13, s3 = F.s13, F.s3
    return (F.q(a[0]) + F.q(a[1])*s13 + F.q(a[2])*s3 + F.q(a[3])*s13 % P*s3) % P
def kv_terms(kv):       # K-vector -> list of (K coef, rational dict), one per Weil component
    comps = [{}, {}, {}, {}]
    for k, a in kv.items():
        for i in range(4):
            if a[i]: comps[i][k] = a[i]
    return [(BAS[i], comps[i]) for i in range(4) if comps[i]]
def kv_from_terms(terms):
    out = {}
    for c, v in terms:
        if K.iszero(c) or not v: continue
        for k, x in v.items():
            out[k] = K.add(out.get(k, K.Z), K.scal(x, c))
    return {k: a for k, a in out.items() if not K.iszero(a)}
def kv_scale(alpha, v): return {k: K.mul(alpha, a) for k, a in v.items()}
def kv_add(acc, v, alpha=None):
    for k, a in v.items():
        b = a if alpha is None else K.mul(alpha, a)
        acc[k] = K.add(acc.get(k, K.Z), b)
    return acc
def kv_clean(v): return {k: a for k, a in v.items() if not K.iszero(a)}
def build_op(terms, n, ws):
    """exact operator of the field u = sum c v (terms) in mode n, on the basis states of weight ws"""
    O = {}
    for k in V.basis(ws):
        acc = {}
        for c, v in terms:
            r = V.mode_vec(v, n, {k: mpq(1)})
            for kk, x in r.items():
                if x: acc[kk] = K.add(acc.get(kk, K.Z), K.scal(x, c))
        acc = kv_clean(acc)
        if acc: O[k] = acc
    V._MODE.clear(); V._EMODE.clear()
    return O
def apply_op(O, kv):
    out = {}
    for k, a in kv.items():
        img = O.get(k)
        if img is None: continue
        if a == K.ONE:
            for kk, b in img.items(): out[kk] = K.add(out.get(kk, K.Z), b)
        else:
            for kk, b in img.items(): out[kk] = K.add(out.get(kk, K.Z), K.mul(a, b))
    return kv_clean(out)
def kprod(terms_u, n, kv):
    out = {}
    for c, v in terms_u:
        for d, w in kv_terms(kv):
            r = V.mode_vec(v, n, w)
            if not r: continue
            cd = K.mul(c, d)
            for kk, x in r.items():
                if x: out[kk] = K.add(out.get(kk, K.Z), K.scal(x, cd))
    return kv_clean(out)
class Blocks:
    def __init__(self, T, w):
        self.T = T; self.w = w; self.keys = V.basis(w)
        self.blk = {}
        for k in self.keys: self.blk.setdefault(C.weight_of_key(T, k), []).append(k)
        self.pos = {mu: {k: i for i, k in enumerate(ks)} for mu, ks in self.blk.items()}
    def wt(self, kv): return C.weight_of_key(self.T, next(iter(kv)))
def modp_vec(kv, pos, n):
    import numpy as np
    z = np.zeros(n, dtype=np.int64)
    for k, a in kv.items(): z[pos[k]] = kred(a)
    return z
def select_pivots(vecs, pos, n, start=None):
    """indices of a maximal subset (greedy, in order) independent modulo p; start = already chosen vectors"""
    import numpy as np
    chosen = list(start or []); basis_rows = [modp_vec(v, pos, n) for v in chosen]; idx = []
    rk = rank_mod(basis_rows, n)
    for i, v in enumerate(vecs):
        z = modp_vec(v, pos, n)
        if not z.any(): continue
        r2 = rank_mod(basis_rows + [z], n)
        if r2 > rk: basis_rows.append(z); idx.append(i); rk = r2
    return idx
def rank_mod(rows, n):
    if not rows: return 0
    return nmod_mat(len(rows), n, [int(x) for r in rows for x in r], P).rank()
def weil_column(kv, pos, n):
    col = [0]*(4*n); den = 1
    for k, a in kv.items():
        for i in range(4):
            if a[i]: den = den*int(a[i].denominator)//__import__('math').gcd(den, int(a[i].denominator))
    for k, a in kv.items():
        for i in range(4):
            if a[i]: col[4*pos[k] + i] = int(a[i]*den)
    return col, den
def solve_block(B, X, pos, n, verify_chunk=2000):
    """B: list of K-vectors (independent over K expected), X: list of K-vectors in the same block.
    Returns (rankK_ok, list of indices of X NOT in span_K(B)) -- exact."""
    if not X: return True, []
    if not B: return True, [i for i, x in enumerate(X) if x]
    cols = []
    for b in B:
        for z in BAS: cols.append(weil_column(kv_scale(z, b), pos, n)[0])
    m = len(cols); N = 4*n
    # pivot rows modulo p
    Mt = nmod_mat(m, N, [c % P for col in cols for c in col], P)
    R, rk = Mt.rref()
    if rk < m: return False, None
    piv = []
    j = 0
    for i in range(rk):
        while int(R[i, j]) == 0: j += 1
        piv.append(j); j += 1
    Xc = [weil_column(x, pos, n)[0] for x in X]
    BP = fmpz_mat(m, m, [cols[c][r] for r in piv for c in range(m)])
    XP = fmpz_mat(m, len(Xc), [Xc[t][r] for r in piv for t in range(len(Xc))])
    Csol = BP.solve(XP)                        # fmpq_mat, exact
    Num, den = Csol.numer_denom()              # Csol = Num / den, Num integral
    bad = set()
    for r0 in range(0, N, verify_chunk):
        rows = range(r0, min(N, r0 + verify_chunk))
        Bc = fmpz_mat(len(rows), m, [cols[c][r] for r in rows for c in range(m)])
        Xr = fmpz_mat(len(rows), len(Xc), [Xc[t][r] for r in rows for t in range(len(Xc))])
        D = Bc*Num - Xr*den                    # exact integer identity B C = X  <=>  D = 0
        for t in range(len(Xc)):
            if t in bad: continue
            for i in range(len(rows)):
                if D[i, t] != 0: bad.add(t); break
    return True, sorted(bad)
def setup(chi=0):
    T, D, out, types, N1, N2, s13, U1, A2 = CE.setup('II', 9)
    info = json.load(open('closure_explore_II_9.json')); Lj = json.load(open('isotropic_lines_II_9.json'))
    G = json.load(open('griess_subalgebras.json'))['II']['subalgebras'][9]
    _, basis, _, _ = GS.griess_data('II')
    tl = sorted(types); choice = info['combos'][chi]
    lines = {tau: [[K.from_sympy(c) for c in sp.sympify(l)] for l in Lj[str(tau)][0]['lines']] for tau in tl}
    xs = [ES.split(c) for c in G['x']]
    xterms = []
    for i, (p_, q_) in enumerate(xs):
        if p_: xterms.append((K.scal(mpq(int(p_.p), int(p_.q)), K.ONE), basis[i + 1]))
        if q_: xterms.append((K.scal(mpq(int(q_.p), int(q_.q)), K.R), basis[i + 1]))
    hw = {'e': [(K.ONE, basis[0])], 'x': xterms}
    mods = []                    # (name, index, terms) for all module weight vectors
    for tau, l in zip(tl, choice):
        ty = types[tau]; lv = lines[tau][l]
        for w in range(len(ty['mods'][0][0])):
            terms = []
            for k in range(ty['d']):
                if ty['mods'][k][0][w]: terms.append((lv[k], ty['mods'][k][0][w]))
                if ty['mods'][k][1][w]: terms.append((K.mul(lv[k], K.R), ty['mods'][k][1][w]))
            mods.append((str(tau), w, terms))
            if w == 0: hw[str(tau)] = terms
    U2span = [kv_from_terms([(K.ONE, a)]) for a in A2 if a] + [kv_from_terms(hw['e']), kv_from_terms(hw['x'])] + [kv_from_terms(t) for _, _, t in mods]
    U1k = [kv_from_terms([(K.ONE, u)]) for u in U1]
    return dict(T=T, types=types, U1=U1, U1k=U1k, A2=A2, basis=basis, hw=hw, mods=mods, U2span=[u for u in U2span if u], choice=choice)

class Ech:
    """incremental row echelon form modulo p (numpy), for greedy pivot selection"""
    def __init__(self, n):
        import numpy as np
        self.n = n; self.rows = []; self.piv = []
    def add(self, z):
        import numpy as np
        z = z % P
        for pc, row in zip(self.piv, self.rows):
            c = int(z[pc])
            if c: z = (z - c*row) % P
        nz = np.nonzero(z)[0]
        if len(nz) == 0: return False
        pc = int(nz[0]); inv = pow(int(z[pc]), -1, P)
        self.rows.append((z*inv) % P); self.piv.append(pc)
        return True
