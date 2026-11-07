"""EXACT highest weight VOA relations in weights 3 and 4 and their Zhu polynomials over K (case II).
Independent exact rewrite of the modular discovery z2 of Article I.
 Phi_3   : U_2 + g(x)U_2 + Lambda^2 U_2 -> V'_3,  (v, h(x)w, a^b) -> Tv + h(-1)w + B(a,b)
 Phi_sym : Sym^2 U_2 -> V'_4,                      a.b -> N(a,b) = (a(-1)b + b(-1)a)/2
 Phi_mix : S_3 + g(x)S_3 + Sym^2 U_2 -> V'_4,      (u, h(x)v, a.b) -> Tu + h(-1)v + N(a,b)
Highest vectors of the formal domains are computed over Q (integer e-action); kernels over K by restriction of
scalars with exact verification; Zhu polynomials have coefficients in K.  Output: rel_exact.pkl, z2x.json."""
import sys, os, time, pickle, json, glob, gc
from collections import defaultdict
sys.path.insert(0, '/home/claude/hj/work/Calculs/zhu')
sys.path.insert(0, '/home/claude/hj/work/Calculs/zhu-exact')
HERE = '/home/claude/hj/work/Calculs/zhu-exact'
from zhu_common import *          # chdir to generation-forte; E, V, K, mpq, np
import kexact as KX
from flint import fmpz_mat, nmod_mat
t0 = time.time()
def log(*a): print(*a, '(%.0fs)' % (time.time() - t0), flush=True)
S = E.setup(0); T = S['T']
d = pickle.load(open('/home/claude/hj/work/Calculs/zhu/U2_lowering_basis.pkl', 'rb'))
HV = d['hv']; LB = [(a, i, j) for a, i, j, _ in d['basis']]; LBI = {x: n for n, x in enumerate(LB)}
TYPE = [w for _, w in HV]; Z = [b[3] for b in d['basis']]
def wtU2(b): a, i, j = LB[b]; r, l = TYPE[a]; return (r - 2 * i, l - 2 * j)
TO = pickle.load(open(HERE + '/prodx/Tops.pkl', 'rb'))
SMALL = pickle.load(open(HERE + '/prodw/small.pkl', 'rb'))
_cache = {}
def bucket(kind, w):
    key = (kind, w)
    if key not in _cache:
        if len(_cache) > 6: _cache.clear(); gc.collect()
        dct = {}
        for f in glob.glob(HERE + '/prodw/%s_%d_%d/a*.pkl' % (kind, w[0], w[1])): dct.update(pickle.load(open(f, 'rb')))
        _cache[key] = dct
    return _cache[key]
def Bv(a, b):
    w = (wtU2(a)[0] + wtU2(b)[0], wtU2(a)[1] + wtU2(b)[1]); return bucket('B', w)[a, b]
def Nv(a, b):
    w = (wtU2(a)[0] + wtU2(b)[0], wtU2(a)[1] + wtU2(b)[1]); return bucket('N', w)[a, b]
ONE = K.ONE
def kq(c): return K.scal(mpq(c), ONE) if not isinstance(c, tuple) else c
def vadd(acc, v, c):           # acc += c v, c in K
    for k, a in v.items(): acc[k] = K.add(acc.get(k, K.Z), K.mul(c, a))
def vclean(v): return {k: a for k, a in v.items() if not K.iszero(a)}
# ---------------- exact expansion in the lowering basis of U_2 ----------------
POS2 = {k: i for i, k in enumerate(V.basis(2))}
def make_expander(vecs):
    keys = sorted({k for v in vecs for k in v}); idx = {k: i for i, k in enumerate(keys)}; n = len(vecs)
    A = nmod_mat(n, len(keys), [x for v in vecs for x in KX.kv_red(v, idx, len(keys))], KX.P)
    R, r = A.rref(); assert r == n, 'basis not independent mod p'
    J = []; j = 0
    for i in range(r):
        while int(R[i, j]) == 0: j += 1
        J.append(j); j += 1
    M = [[vecs[i].get(keys[c], K.Z) for c in J] for i in range(n)]   # n x n : coefficients c with c M = z_J
    Mt = [[M[i][jj] for i in range(n)] for jj in range(n)]
    def expand(z):
        if not z: return {}
        rhs = [[z.get(keys[c], K.Z)] for c in J] if all(k in idx for k in z) else None
        if rhs is None: raise ValueError('not in span (support)')
        X = KX.solve_K(Mt, rhs)
        c = {i: X[i][0] for i in range(n) if not K.iszero(X[i][0])}
        acc = {}
        for i, a in c.items(): vadd(acc, vecs[i], a)
        for k in set(acc) | set(z):
            if not K.iszero(K.add(acc.get(k, K.Z), K.scal(mpq(-1), z.get(k, K.Z)))): raise ValueError('not in span')
        return c
    return expand
expand_U2_raw = make_expander(Z)
_exp_cache = {}
def expand_U2(z, key):
    if key not in _exp_cache: _exp_cache[key] = {('Z', b): c for b, c in expand_U2_raw(z).items()}
    return _exp_cache[key]
def lc(*pairs):
    out = {}
    for c, v in pairs:
        if v: vadd(out, v, kq(c))
    return vclean(out)
def cminus(a, b): return expand_U2(lc((mpq(1, 2), SMALL[a, b]['Tp2'])), ('-', a, b))
def cplus(a, b): return expand_U2(lc((1, SMALL[a, b]['p1']), (mpq(-1, 2), SMALL[a, b]['Tp2'])), ('+', a, b))
def padd(acc, poly, c):        # c in K
    for m, x in poly.items():
        y = K.add(acc.get(m, K.Z), K.mul(c, x))
        if K.iszero(y): acc.pop(m, None)
        else: acc[m] = y
    return acc
# ---------------- formal domains (exact integer action) ----------------
GW = [(0, 0), (2, 0), (-2, 0), (0, 0), (0, 2), (0, -2)]
ADE = {1: {0: {1: -2}, 2: {0: 1}}, 4: {3: {4: -2}, 5: {3: 1}}}
ADF = {2: {0: {2: 2}, 1: {0: -1}}, 5: {3: {5: 2}, 4: {3: -1}}}
def e_low(types, x, fac):
    a, i, j = x; r, l = types[a]
    if fac == 0: return {} if i == 0 else {(a, i - 1, j): i * (r - i + 1)}
    return {} if j == 0 else {(a, i, j - 1): j * (l - j + 1)}
def f_low(types, x, fac):
    a, i, j = x; r, l = types[a]
    if fac == 0: return {} if i == r else {(a, i + 1, j): 1}
    return {} if j == l else {(a, i, j + 1): 1}
def ad_cur(h, fac, kind):
    tab = (ADE if kind == 'e' else ADF)[(1 if fac == 0 else 4) if kind == 'e' else (2 if fac == 0 else 5)]
    return tab.get(h, {})
class Domain:
    def __init__(self, parts, utypes): self.parts = parts; self.ut = utypes
    def build(self, ulist, uwt, uidx):
        self.ulist, self.uwt, self.uidx = ulist, uwt, uidx
        el = []
        if 'T' in self.parts: el += [('T', u) for u in range(len(ulist))]
        if 'h' in self.parts: el += [('h', h, u) for h in range(6) for u in range(len(ulist))]
        if 'alt' in self.parts: el += [('alt', a, b) for a in range(120) for b in range(a + 1, 120)]
        if 'sym' in self.parts: el += [('sym', a, b) for a in range(120) for b in range(a, 120)]
        self.elts = el; self.idx = {x: n for n, x in enumerate(el)}
        self.wt = [self.weight(x) for x in el]
        self.byw = defaultdict(list)
        for n, w in enumerate(self.wt): self.byw[w].append(n)
    def weight(self, x):
        if x[0] == 'T': return self.uwt(x[1])
        if x[0] == 'h': w = self.uwt(x[2]); return (w[0] + GW[x[1]][0], w[1] + GW[x[1]][1])
        w1, w2 = wtU2(x[1]), wtU2(x[2]); return (w1[0] + w2[0], w1[1] + w2[1])
    def act(self, x, fac, kind):
        out = defaultdict(int)
        lowf = e_low if kind == 'e' else f_low
        if x[0] == 'T':
            for y, c in lowf(self.ut, self.ulist[x[1]], fac).items(): out[('T', self.uidx[y])] += c
        elif x[0] == 'h':
            h, u = x[1], x[2]
            for h2, c in ad_cur(h, fac, kind).items(): out[('h', h2, u)] += c
            for y, c in lowf(self.ut, self.ulist[u], fac).items(): out[('h', h, self.uidx[y])] += c
        else:
            kd, a, b = x
            for y, c in lowf(TYPE, LB[a], fac).items():
                a2 = LBI[y]
                if kd == 'sym': out[('sym', min(a2, b), max(a2, b))] += c
                elif a2 != b: out[('alt', a2, b) if a2 < b else ('alt', b, a2)] += c if a2 < b else -c
            for y, c in lowf(TYPE, LB[b], fac).items():
                b2 = LBI[y]
                if kd == 'sym': out[('sym', min(a, b2), max(a, b2))] += c
                elif a != b2: out[('alt', a, b2) if a < b2 else ('alt', b2, a)] += c if a < b2 else -c
        return {self.idx[y]: c for y, c in out.items() if c}
def highest_rows(dom, eta):
    """integer basis (list of dicts n -> int) of the vectors of weight eta of the formal domain killed by e1, e2"""
    cols = dom.byw[eta]; pos = {n: i for i, n in enumerate(cols)}
    tgts = []
    for fac, sh in ((0, (2, 0)), (1, (0, 2))):
        tgts.append(dom.byw.get((eta[0] + sh[0], eta[1] + sh[1]), []))
    tp = [{n: i for i, n in enumerate(t)} for t in tgts]
    ncol = sum(len(t) for t in tgts)
    if ncol == 0: return [{n: 1} for n in cols]
    A = fmpz_mat(len(cols), ncol)
    off = [0, len(tgts[0])]
    for n in cols:
        for fac in (0, 1):
            for m, c in dom.act(dom.elts[n], fac, 'e').items(): A[pos[n], off[fac] + tp[fac][m]] += c
    X, nul = A.transpose().nullspace()      # columns x with A^T x = 0  <=>  x^T A = 0
    rows = []
    for j in range(nul):
        v = {cols[i]: int(X[i, j]) for i in range(len(cols)) if X[i, j] != 0}
        if v: rows.append(v)
    return rows
# ---------------- currents on U_2 and S_3 ----------------
HOP2 = {h: E.build_op([(ONE, S['U1'][h])], -1, 2) for h in range(6)}
log('h(-1) on V_2 ready')
def img3(x):
    if x[0] == 'T': return E.apply_op(TO[2], Z[x[1]])
    if x[0] == 'h': return E.apply_op(HOP2[x[1]], Z[x[2]])
    return Bv(x[1], x[2])
def image_of(dom, row, img):
    acc = {}
    for n, c in row.items(): vadd(acc, img(dom.elts[n]), kq(c))
    return vclean(acc)
def zeta3(x):
    if x[0] == 'T': return {(('Z', x[1]),): kq(-2)}
    if x[0] == 'h': return {(('Z', x[2]), ('J', x[1])): ONE}
    a, b = x[1], x[2]; out = {(('Z', a), ('Z', b)): ONE}; padd(out, {(('Z', b), ('Z', a)): ONE}, kq(-1))
    for s, c in cminus(a, b).items(): padd(out, {(s,): c}, kq(-1))
    return out
def poly_of(dom, vec, zeta):     # vec: dict domain index -> K
    out = {}
    for n, c in vec.items(): padd(out, zeta(dom.elts[n]), c)
    return out
def comb(rowsH, coeffs):         # sum_k coeffs[k] rowsH[k] (rowsH integer dicts, coeffs K)
    out = {}
    for k, c in coeffs.items():
        for n, x in rowsH[k].items(): out[n] = K.add(out.get(n, K.Z), K.mul(c, kq(x)))
    return {n: x for n, x in out.items() if not K.iszero(x)}
REL = []; STATS = {'bracket': [], 'symmetric': [], 'mixed': []}
def process(kind, dom, img, zeta, collect_S3=False):
    S3sel = []
    for eta in sorted(dom.byw):
        if eta[0] < 0 or eta[1] < 0: continue
        H = highest_rows(dom, eta)
        if not H: continue
        Ims = [image_of(dom, r, img) for r in H]
        ker, r, st = KX.left_kernel(Ims)
        assert st == 'ok', (kind, eta, st)
        STATS[kind].append((eta, len(H), r, len(ker)))
        for c in ker:
            vec = comb(H, c)
            REL.append((kind, eta, poly_of(dom, vec, zeta)))
        if collect_S3 and r:
            # rows of H whose images are independent (pivot rows modulo p)
            keys = sorted({k for v in Ims for k in v}); idx = {k: i for i, k in enumerate(keys)}
            A = nmod_mat(len(keys), len(H), [x for kk in range(len(keys)) for x in [0] * len(H)], KX.P)
            for i, v in enumerate(Ims):
                for k, a in v.items(): A[idx[k], i] = KX.kred(a)
            R, rk = A.rref(); assert rk == r
            j = 0
            for i in range(rk):
                while int(R[i, j]) == 0: j += 1
                S3sel.append((eta, H[j])); j += 1
        log(kind, eta, 'highest', len(H), 'rank', r, 'relations', len(ker))
    return S3sel
D3 = Domain(('T', 'h', 'alt'), TYPE); D3.build(LB, wtU2, LBI)
log('L_3 domain', len(D3.elts))
S3sel = process('bracket', D3, img3, zeta3, collect_S3=True)
S3types = [eta for eta, _ in S3sel]; S3lift = [row for _, row in S3sel]
log('S_3 highest vectors', len(S3types), 'dim', sum((r + 1) * (l + 1) for r, l in S3types))
S3LB = [(c, i, j) for c, (r, l) in enumerate(S3types) for i in range(r + 1) for j in range(l + 1)]
S3LBI = {x: n for n, x in enumerate(S3LB)}
def lower_dom(dom, vec, i, j):
    vec = dict(vec)
    for fac, k in ((0, i), (1, j)):
        for _ in range(k):
            out = defaultdict(int)
            for n, c in vec.items():
                for m, cc in dom.act(dom.elts[n], fac, 'f').items(): out[m] += c * cc
            vec = {m: c for m, c in out.items() if c}
    return vec
S3dom = [lower_dom(D3, S3lift[c], i, j) for c, i, j in S3LB]
# exact vectors of the lowering basis of S_3 (grouped by weight to reuse buckets)
S3vec = [None] * len(S3LB)
order = sorted(range(len(S3LB)), key=lambda u: D3.wt[next(iter(S3dom[u]))])
for u in order:
    S3vec[u] = image_of(D3, S3dom[u], img3)
    assert S3vec[u], 'zero S_3 vector'
log('S_3 lowering basis: exact vectors', len(S3vec))
# independence check modulo p (rank 956)
keys3 = sorted({k for v in S3vec for k in v}); idx3 = {k: i for i, k in enumerate(keys3)}
A3 = nmod_mat(len(S3vec), len(keys3), [x for v in S3vec for x in KX.kv_red(v, idx3, len(keys3))], KX.P)
rkS3 = A3.rank(); log('rank of S_3 lowering basis mod p:', rkS3)
del A3; gc.collect()
_cache.clear(); gc.collect()
# ---------------- Phi_sym ----------------
def zsym(x):
    a, b = x[1], x[2]; h = K.scal(mpq(1, 2), ONE); out = {}
    padd(out, {(('Z', a), ('Z', b)): ONE}, h); padd(out, {(('Z', b), ('Z', a)): ONE}, h)
    for s, c in cplus(a, b).items(): padd(out, {(s,): c}, ONE)
    return out
DS = Domain(('sym',), TYPE); DS.build(LB, wtU2, LBI)
process('symmetric', DS, lambda x: Nv(x[1], x[2]), zsym)
# ---------------- Phi_mix ----------------
TOP3 = TO[3]
HOP3 = {}
def hop3(h):
    if h not in HOP3: HOP3[h] = E.build_op([(ONE, S['U1'][h])], -1, 3); log('h(-1) on V_3 ready', h)
    return HOP3[h]
_img_cache = {}
def imgmix(x):
    if x[0] == 'T':
        k = ('T', x[1])
        if k not in _img_cache: _img_cache[k] = E.apply_op(TOP3, S3vec[x[1]])
        return _img_cache[k]
    if x[0] == 'h':
        k = ('h', x[1], x[2])
        if k not in _img_cache: _img_cache[k] = E.apply_op(hop3(x[1]), S3vec[x[2]])
        return _img_cache[k]
    return Nv(x[1], x[2])
Z3cache = {}
def Z3(u):
    if u not in Z3cache: Z3cache[u] = poly_of(D3, {n: kq(c) for n, c in S3dom[u].items()}, zeta3)
    return Z3cache[u]
def zmix(x):
    if x[0] == 'T': return {m: K.mul(kq(-3), c) for m, c in Z3(x[1]).items()}
    if x[0] == 'h': return {m + (('J', x[1]),): c for m, c in Z3(x[2]).items()}
    return zsym(x)
def S3wt(u): c, i, j = S3LB[u]; r, l = S3types[c]; return (r - 2 * i, l - 2 * j)
DM = Domain(('T', 'h', 'sym'), S3types); DM.build(S3LB, S3wt, S3LBI)
log('L_4 domain', len(DM.elts))
for h in range(6): hop3(h)
process('mixed', DM, imgmix, zmix)
cnt = defaultdict(int)
for k, _, _ in REL: cnt[k] += 1
log('relation polynomials:', dict(cnt))
# integrality at the prime
bad = [(n, m) for n, (_, _, poly) in enumerate(REL) for m, c in poly.items() if not KX.den_ok(c)]
log('coefficients with p in a denominator:', len(bad))
pickle.dump({'REL': REL, 'TYPE': TYPE, 'HV': HV, 'LB': LB, 'S3types': S3types, 'P': KX.P}, open(HERE + '/rel_exact.pkl', 'wb'))
RELp = [(k, eta, {m: KX.kred(c) for m, c in poly.items() if KX.kred(c)}) for k, eta, poly in REL]
pickle.dump({'REL': RELp, 'TYPE': TYPE, 'HV': HV, 'LB': LB, 'S3types': S3types, 'P': KX.P}, open(HERE + '/rel_exact_reduced_modp.pkl', 'wb'))
json.dump({'counts': dict(cnt), 'S3_highest': len(S3types), 'S3_types': [list(t) for t in S3types], 'rank_S3_lowering_modp': rkS3,
           'kernel_blocks': {k: [[list(e), h, r, nk] for e, h, r, nk in v] for k, v in STATS.items()}, 'p_denominators': len(bad), 'prime': KX.P},
          open(HERE + '/z2x.json', 'w'), indent=1)
log('done')
