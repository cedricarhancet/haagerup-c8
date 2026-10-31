"""Modular discovery (p3 = 8388571, coherent embedding of K): highest weight VOA relations in weights 3 and 4 and
their Zhu polynomials (Riesen's three families, adapted to case II).
 Phi_3   : U_2 + g(x)U_2 + Lambda^2 U_2 -> V'_3,  (v, h(x)w, a^b) -> Tv + h(-1)w + B(a,b)
 Phi_sym : Sym^2 U_2 -> V'_4,                      a.b -> N(a,b)
 Phi_mix : S_3 + g(x)S_3 + Sym^2 U_2 -> V'_4,      (u, h(x)v, a.b) -> Tu + h(-1)v + N(a,b)
A highest relation is a vector of the formal domain, killed by e1, e2 (domain action) and by Phi.
Zhu conversion: zeta_3(v) = -2Z(v), zeta_3(h(x)w) = Z(w)J_h, zeta_3(a^b) = Z_aZ_b - Z_bZ_a - Z(c_-(a,b));
zeta_sym(a.b) = (Z_aZ_b + Z_bZ_a)/2 + Z(c_+(a,b)); zeta_mix(u) = -3Z_3(u), zeta_mix(h(x)v) = Z_3(v)J_h,
with Z_3 = zeta_3 o sigma_3 for an equivariant section sigma_3 (lifts of highest vectors of S_3, extended by lowering).
Polynomials: dict {tuple of symbols: coeff mod p}, symbols ('Z', b) (b in the lowering basis of U_2), ('J', h)."""
import pickle, time, json, gc, os
from collections import defaultdict
import numpy as np, scipy.sparse as sp
from flint import nmod_mat
from zhu_common import *
import cases as C
t0 = time.time()
def log(*a): print(*a, '(%.0fs)' % (time.time() - t0), flush=True)
S = E.setup(0); T = S['T']
d = pickle.load(open(ZD + '/U2_lowering_basis.pkl', 'rb'))
HV = d['hv']; LB = [(a, i, j) for a, i, j, _ in d['basis']]; LBI = {x: n for n, x in enumerate(LB)}
TYPE = [w for _, w in HV]
def wtU2(b): a, i, j = LB[b]; r, l = TYPE[a]; return (r - 2 * i, l - 2 * j)
KEYS = {w: V.basis(w) for w in range(5)}
WT = {w: np.array([C.weight_of_key(T, k) for k in KEYS[w]], dtype=np.int64) for w in range(5)}
Tm = pickle.load(open(ZD + '/prod/Tmat.pkl', 'rb'))
def dense(sv, w):
    z = np.zeros(len(KEYS[w]), np.int64)
    if len(sv[0]): z[sv[0]] = sv[1]
    return z
def Tapp(z, w, k=1):
    for _ in range(k): z = (Tm[w] @ z) % P; w += 1
    return z
# products a(n)b, a <= b
PR = {}
for a in range(120): PR.update({(a,) + k: v for k, v in pickle.load(open(ZD + '/prod/a%03d.pkl' % a, 'rb')).items()})
HU = pickle.load(open(ZD + '/prod/hU2.pkl', 'rb'))
inv2, inv4, inv12 = pow(2, -1, P), pow(4, -1, P), pow(12, -1, P)
def pr(a, b, n, w):
    if a > b: raise ValueError
    return dense(PR[a, b, n], w)
def Bvec(a, b): return (pr(a, b, 0, 3) - Tapp(pr(a, b, 1, 2), 2) * inv2 + Tapp(pr(a, b, 2, 1), 1, 2) * inv4) % P
def Nvec(a, b): return (pr(a, b, -1, 4) - Tapp(pr(a, b, 0, 3), 3) * inv2 + Tapp(pr(a, b, 1, 2), 2, 2) * inv4 - Tapp(pr(a, b, 2, 1), 1, 3) * inv12) % P
def cplus(a, b): return (pr(a, b, 1, 2) - Tapp(pr(a, b, 2, 1), 1) * inv2) % P
def cminus(a, b): return (Tapp(pr(a, b, 2, 1), 1) * inv2) % P
# expansion of U_2 vectors (V'_2 coordinates) in the lowering basis
LBM = np.array([dense(modp_sparse(x[3], {k: i for i, k in enumerate(KEYS[2])}), 2) for x in d['basis']], dtype=np.int64)
def solve_rows(M, ncols_hint=None):
    """pivot columns J such that M[:, J] is invertible mod p, and its inverse"""
    A = nmod_mat(M.shape[0], M.shape[1], [int(x) for x in M.flatten()], P)
    R = A.rref()[0]; piv = []
    r = 0
    for c in range(M.shape[1]):
        if r < M.shape[0] and int(R[r, c]) == 1 and all(int(R[rr, c]) == 0 for rr in range(M.shape[0]) if rr != r):
            piv.append(c); r += 1
    sub = nmod_mat(M.shape[0], len(piv), [int(M[i, c]) for i in range(M.shape[0]) for c in piv], P)
    return piv, sub.inv()
LPIV, LINV = solve_rows(LBM)
def expand_U2(z):          # coefficients c with c . LB = z (checked)
    rhs = nmod_mat(1, len(LPIV), [int(z[c]) for c in LPIV], P)
    c = rhs * LINV; c = np.array([int(c[0, i]) for i in range(120)], np.int64)
    assert np.all((c @ LBM - z) % P == 0), 'not in U_2'
    return {('Z', b): int(c[b]) for b in range(120) if c[b]}
def padd(acc, poly, s=1):
    for m, c in poly.items():
        x = (acc.get(m, 0) + s * c) % P
        if x: acc[m] = x
        else: acc.pop(m, None)
    return acc
# ---- g action on the formal domains (lowering bases) ----
GW = [(0, 0), (2, 0), (-2, 0), (0, 0), (0, 2), (0, -2)]          # H1 e1 f1 H2 e2 f2
ADE = {1: {0: {1: -2}, 2: {0: 1}}, 4: {3: {4: -2}, 5: {3: 1}}}     # ad e1, ad e2 on the current basis
ADF = {2: {0: {2: 2}, 1: {0: -1}}, 5: {3: {5: 2}, 4: {3: -1}}}     # ad f1, ad f2
def e_low(types, x, fac):   # e1 (fac=0) or e2 (fac=1) on a lowering basis element (a,i,j)
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
    """formal domain = direct sum of parts: ('U', U-index) / ('gU', h, U-index) / ('alt', b1, b2) / ('sym', b1, b2)"""
    def __init__(self, parts, utypes=None):
        self.elts = []; self.parts = parts; self.ut = utypes
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
    def act(self, x, fac, kind):    # e (kind 'e') or f (kind 'f') acting on a domain element -> {elt: coeff}
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
        return {self.idx[y]: c % P for y, c in out.items() if c % P}
def kernel_rows(M):
    """basis of {c : c M = 0} mod p (rows)"""
    if M.shape[1] == 0: return np.eye(M.shape[0], dtype=np.int64)
    A = nmod_mat(M.shape[1], M.shape[0], [int(x) for x in M.T.flatten()], P)
    X, nul = A.nullspace()
    return np.array([[int(X[i, j]) for i in range(M.shape[0])] for j in range(nul)], dtype=np.int64).reshape(nul, M.shape[0])
def highest_kernel(dom, image, wimg, use_phi=True):
    """for each dominant weight eta: (domain indices, highest vectors of the domain (rows), kernel of Phi among them)"""
    res = {}
    for eta, cols in sorted(dom.byw.items()):
        if eta[0] < 0 or eta[1] < 0: continue
        pos = {n: i for i, n in enumerate(cols)}
        blocks = []
        for fac, sh in ((0, (2, 0)), (1, (0, 2))):
            tgt = dom.byw.get((eta[0] + sh[0], eta[1] + sh[1]), [])
            tp = {n: i for i, n in enumerate(tgt)}
            Mx = np.zeros((len(cols), len(tgt)), np.int64)
            for n in cols:
                for m, c in dom.act(dom.elts[n], fac, 'e').items(): Mx[pos[n], tp[m]] = c
            blocks.append(Mx)
        Hrows = kernel_rows(np.hstack(blocks)) if sum(b.shape[1] for b in blocks) else np.eye(len(cols), dtype=np.int64)
        if Hrows.shape[0] == 0: continue
        sel = np.nonzero(WT[wimg][:, 0] == eta[0])[0]; sel = sel[WT[wimg][sel, 1] == eta[1]]
        Im = np.array([image(dom.elts[n])[sel] for n in cols], dtype=np.int64).reshape(len(cols), len(sel))
        HI = (Hrows @ Im) % P
        K2 = kernel_rows(HI)
        rel = (K2 @ Hrows) % P if K2.shape[0] else np.zeros((0, len(cols)), np.int64)
        res[eta] = (cols, Hrows, rel, HI)
    return res
# ------------------------------------------------------------------ Phi_3
def img3(x):
    if x[0] == 'T': return dense(HU['T', x[1]], 3)
    if x[0] == 'h': return dense(HU['h', x[1], x[2]], 3)
    return Bvec(x[1], x[2])
D3 = Domain(('T', 'h', 'alt'), TYPE); D3.build(LB, wtU2, LBI)
log('L_3 domain', len(D3.elts))
K3 = highest_kernel(D3, img3, 3)
nrel3 = sum(v[2].shape[0] for v in K3.values()); nhw3 = sum(v[1].shape[0] for v in K3.values())
log('Phi_3: highest domain vectors', nhw3, ' highest relations', nrel3)
def zeta3(x):
    if x[0] == 'T': return {(('Z', x[1]),): (-2) % P}
    if x[0] == 'h': return {(('Z', x[2]), ('J', x[1])): 1}
    a, b = x[1], x[2]; out = {(('Z', a), ('Z', b)): 1}; padd(out, {(('Z', b), ('Z', a)): 1}, -1)
    for s, c in expand_U2(cminus(a, b)).items(): padd(out, {(s,): c}, -1)
    return out
def poly_of(dom, cols, row, zeta):
    out = {}
    for n, c in zip(cols, row):
        if c: padd(out, zeta(dom.elts[n]), int(c))
    return out
REL = []
for eta, (cols, H, rel, HI) in K3.items():
    for r in rel: REL.append(('bracket', eta, poly_of(D3, cols, r, zeta3)))
log('bracket polynomials:', len(REL))
# ---- S_3: highest vectors from HW(L_3) images, lowering basis, section sigma_3
S3types = []; S3lift = []          # per highest vector of S_3: type, lift (dict domain idx -> coeff)
for eta, (cols, H, rel, HI) in K3.items():
    # choose rows of H whose images are independent
    if HI.shape[0] == 0: continue
    A = nmod_mat(HI.shape[1], HI.shape[0], [int(x) for x in HI.T.flatten()], P)
    R = A.rref()[0]; rank = 0; piv = []
    for c in range(HI.shape[0]):
        if rank < HI.shape[1] and int(R[rank, c]) == 1: piv.append(c); rank += 1
    for c in piv:
        S3types.append(eta); S3lift.append({cols[i]: int(H[c, i]) for i in range(len(cols)) if H[c, i]})
log('S_3 highest vectors:', len(S3types), ' dim', sum((r + 1) * (l + 1) for r, l in S3types))
S3LB = [(c, i, j) for c, (r, l) in enumerate(S3types) for i in range(r + 1) for j in range(l + 1)]
S3LBI = {x: n for n, x in enumerate(S3LB)}
def lower_dom(dom, vec, i, j):
    for fac, k in ((0, i), (1, j)):
        for _ in range(k):
            out = defaultdict(int)
            for n, c in vec.items():
                for m, cc in dom.act(dom.elts[n], fac, 'f').items(): out[m] = (out[m] + c * cc) % P
            vec = {m: c for m, c in out.items() if c}
    return vec
S3dom = [lower_dom(D3, S3lift[c], i, j) for c, i, j in S3LB]          # sigma_3 of the lowering basis of S_3
def img_comb(dom, vec, image):
    z = None
    for n, c in vec.items():
        v = image(dom.elts[n]) * c
        z = v if z is None else z + v
    return z % P
S3vec = [img_comb(D3, v, img3) for v in S3dom]
# coordinates in the exact basis of S_3 (mod p), then T u and h(-1)u
SD = pickle.load(open(ZD + '/prod/S3.pkl', 'rb'))
S3M = np.array([dense(s, 3) for s in SD['S3']], np.int64)
SPIV, SINV = solve_rows(S3M)
def spm(lst, w):
    rows, cols, vals = [], [], []
    for i, (ix, vl) in enumerate(lst): rows += [i] * len(ix); cols += list(ix); vals += list(vl)
    return sp.csr_matrix((np.array(vals, np.int64), (rows, cols)), shape=(len(lst), len(KEYS[w])))
def coords_S3(z):
    rhs = nmod_mat(1, len(SPIV), [int(z[c]) for c in SPIV], P); c = rhs * SINV
    c = np.array([int(c[0, i]) for i in range(len(SPIV))], np.int64)
    assert np.all((c @ S3M - z) % P == 0)
    return c
S3C = sp.csr_matrix(np.array([coords_S3(z) for z in S3vec], np.int64))
log('S_3 lowering basis expressed in the exact basis:', S3C.shape)
def mulmod(A, B):   # sparse product mod p without overflow: split B entries in 12-bit pieces
    out = None
    Bc = B.tocsr(); A = A.tocsr()
    for sh in (0, 12):
        Bp = Bc.copy(); Bp.data = (Bp.data >> sh) & 4095
        Cp = (A @ Bp); Cp.data %= P
        Cp = Cp * (1 << sh); Cp.data %= P
        out = Cp if out is None else out + Cp
    out.data %= P; out.eliminate_zeros(); return out.tocsr()
TU = mulmod(S3C, spm([SD['T', i] for i in range(len(SD['S3']))], 4))
HUs = {h: mulmod(S3C, spm([SD['h', h, i] for i in range(len(SD['S3']))], 4)) for h in range(6)}
del SD; gc.collect()
def Z3(u): return poly_of(D3, list(S3dom[u].keys()), list(S3dom[u].values()), zeta3)
Z3cache = {}
def Z3c(u):
    if u not in Z3cache: Z3cache[u] = Z3(u)
    return Z3cache[u]
def S3wt(u): c, i, j = S3LB[u]; r, l = S3types[c]; return (r - 2 * i, l - 2 * j)
# ------------------------------------------------------------------ Phi_sym and Phi_mix
def zsym(x):
    a, b = x[1], x[2]; out = {}
    padd(out, {(('Z', a), ('Z', b)): inv2}); padd(out, {(('Z', b), ('Z', a)): inv2})
    for s, c in expand_U2(cplus(a, b)).items(): padd(out, {(s,): c})
    return out
DS = Domain(('sym',), TYPE); DS.build(LB, wtU2, LBI)
KS = highest_kernel(DS, lambda x: Nvec(x[1], x[2]), 4)
log('Phi_sym: highest domain vectors', sum(v[1].shape[0] for v in KS.values()), ' highest relations', sum(v[2].shape[0] for v in KS.values()))
for eta, (cols, H, rel, HI) in KS.items():
    for r in rel: REL.append(('symmetric', eta, poly_of(DS, cols, r, zsym)))
del KS; gc.collect()
def imgmix(x):
    if x[0] == 'T': return TU[x[1]].toarray().ravel() % P
    if x[0] == 'h': return HUs[x[1]][x[2]].toarray().ravel() % P
    return Nvec(x[1], x[2])
def zmix(x):
    if x[0] == 'T': return {m: (-3 * c) % P for m, c in Z3c(x[1]).items()}
    if x[0] == 'h': return {m + (('J', x[1]),): c for m, c in Z3c(x[2]).items()}
    return zsym(x)
DM = Domain(('T', 'h', 'sym'), S3types); DM.build(S3LB, S3wt, S3LBI)
log('L_4 domain', len(DM.elts))
KM = highest_kernel(DM, imgmix, 4)
log('Phi_mix: highest domain vectors', sum(v[1].shape[0] for v in KM.values()), ' highest relations', sum(v[2].shape[0] for v in KM.values()))
for eta, (cols, H, rel, HI) in KM.items():
    for r in rel: REL.append(('mixed', eta, poly_of(DM, cols, r, zmix)))
cnt = defaultdict(int)
for k, _, _ in REL: cnt[k] += 1
log('relation polynomials:', dict(cnt))
pickle.dump({'REL': REL, 'TYPE': TYPE, 'HV': HV, 'LB': LB, 'S3types': S3types, 'P': P}, open(ZD + '/relations_modp.pkl', 'wb'))
json.dump({'counts': dict(cnt), 'S3_highest': len(S3types), 'S3_types': [list(t) for t in S3types], 'prime': P}, open(ZD + '/z2_relations_modp.json', 'w'), indent=1)
log('done')
