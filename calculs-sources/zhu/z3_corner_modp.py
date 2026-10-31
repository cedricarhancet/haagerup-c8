"""Modular discovery (p3): projection of the Zhu relations to the full Morita corner of A(S) given by A(A+),
A+ = A + J (14 blocks), followed by linear elimination, quadratic saturation, completion of the rewriting system
(diamond lemma on endpoint corners) and structure of the resulting finite-dimensional algebra (trace form, center,
blocks).  This follows the architecture of Riesen's Section 9 (code of the certificates C11-C15 used as a model)."""
import pickle, json, time, sys, random
from collections import defaultdict
from functools import lru_cache
import numpy as np
from flint import nmod_mat
t0 = time.time()
def log(*a): print(*a, '(%.0fs)' % (time.time() - t0), flush=True)
import os
ZD = '/home/claude/work/Calculs/zhu'
if os.environ.get('HAAG_PDIR'): ZD = os.environ['HAAG_PDIR']
D = pickle.load(open(ZD + '/relations_modp.pkl', 'rb'))
P = D['P']; REL = D['REL']; HV = D['HV']; TYPE = D['TYPE']; LB = D['LB']
inv = lambda x: pow(int(x) % P, -1, P)
# ---- constituents and vertices (blocks of A(A+)) ----
CONS = [(0, 0), (0, 2), (1, 1), (2, 0), (0, 4), (1, 3), (2, 2), (3, 1), (0, 6), (1, 5), (2, 4), (3, 3), (1, 7), (2, 6), (0, 8), (3, 5)]
VERT = list(range(12)) + [12, 12, 13, 13]
NV = 14
DIM = [(r + 1) * (l + 1) for r, l in CONS]
def hA(t): r, l = t; return (r * (r + 2) * inv(20) + l * (l + 2) * inv(60)) % P
# ---- sl2 matrices in the lowering basis v, fv, f^2v, ... ----
def sl2(n):
    e = np.zeros((n + 1, n + 1), np.int64); f = e.copy(); h = e.copy()
    for i in range(n + 1):
        h[i, i] = (n - 2 * i) % P
        if i < n: f[i + 1, i] = 1
        if i > 0: e[i - 1, i] = i * (n - i + 1) % P
    return e, f, h
def rho(t, g):     # g in 0..5 = H1 e1 f1 H2 e2 f2, on V(r) (x) V(l), second factor fastest
    r, l = t; e1, f1, h1 = sl2(r); e2, f2, h2 = sl2(l); I1 = np.eye(r + 1, dtype=np.int64); I2 = np.eye(l + 1, dtype=np.int64)
    m = {0: np.kron(h1, I2), 1: np.kron(e1, I2), 2: np.kron(f1, I2), 3: np.kron(I1, h2), 4: np.kron(I1, e2), 5: np.kron(I1, f2)}[g]
    return m % P
RHO = {(c, g): rho(CONS[c], g) for c in range(len(CONS)) for g in range(6)}
def allow(m, n, l): return abs(m - n) <= l <= m + n and (m + n - l) % 2 == 0
@lru_cache(None)
def hw1(n, m, l):  # normalized highest matrix Hom(V(n), V(m)) of weight l (one factor)
    e_n, f_n, h_n = sl2(n); e_m, f_m, h_m = sl2(m)
    dlt = (l - m + n) // 2
    # unknown F (m+1)x(n+1): e_m F - F e_n = 0, h-weight l, normalized F[0, dlt] = 1
    N = (m + 1) * (n + 1); rows = []
    for i in range(m + 1):
        for j in range(n + 1):
            row = np.zeros(N, np.int64)
            for k in range(m + 1): row[k * (n + 1) + j] += e_m[i, k]
            for k in range(n + 1): row[i * (n + 1) + k] -= e_n[k, j]
            rows.append(row % P)
    for i in range(m + 1):
        for j in range(n + 1):
            if (m - 2 * i) - (n - 2 * j) != l:
                row = np.zeros(N, np.int64); row[i * (n + 1) + j] = 1; rows.append(row)
    A = np.array(rows, np.int64)
    M = nmod_mat(A.shape[0], A.shape[1], [int(x) for x in A.flatten()], P)
    X, nul = M.nullspace(); assert nul == 1, (n, m, l, nul)
    F = np.array([int(X[k, 0]) for k in range(N)], np.int64).reshape(m + 1, n + 1)
    c = int(F[0, dlt]); assert c, (n, m, l)
    return F * inv(c) % P
def hw2(lam, mu, nu): return np.kron(hw1(lam[0], mu[0], nu[0]), hw1(lam[1], mu[1], nu[1])) % P
# ---- arrows ----
FREE = [a for a, (lab, w) in enumerate(HV) if lab.startswith('P:')]
ARROWS = []; AID = {}
for a in FREE:
    nu = tuple(HV[a][1])
    for s, lam in enumerate(CONS):
        for t, mu in enumerate(CONS):
            if allow(lam[0], mu[0], nu[0]) and allow(lam[1], mu[1], nu[1]):
                AID[a, s, t] = len(ARROWS); ARROWS.append({'z': HV[a][0], 'a': a, 'src_c': s, 'tgt_c': t, 'source': VERT[s], 'target': VERT[t]})
log('quiver: %d vertices, %d constituents, %d arrows' % (NV, len(CONS), len(ARROWS)))
ELOOP = {c: AID[[a for a in FREE if HV[a][0] == 'P:e'][0], c, c] for c in range(len(CONS))}
# ---- X-hat of the symbols: dict (tgt_c, src_c) -> list of (path, matrix) ----
def mmul(X, Y):
    out = defaultdict(lambda: defaultdict(lambda: None))
    for (t, m), lx in X.items():
        for (m2, s), ly in Y.items():
            if m != m2: continue
            for px, Mx in lx:
                for py, My in ly:
                    Z = (Mx @ My) % P
                    if not Z.any(): continue
                    cur = out[t, s][px + py]
                    out[t, s][px + py] = Z if cur is None else (cur + Z) % P
    return {k: [(p, M) for p, M in v.items() if M is not None and M.any()] for k, v in out.items()}
def madd(X, Y, c=1):
    out = defaultdict(dict)
    for k, l in X.items():
        for p, M in l: out[k][p] = M % P
    for k, l in Y.items():
        for p, M in l:
            out[k][p] = ((out[k].get(p, 0) + c * M) % P)
    return {k: [(p, M) for p, M in v.items() if np.any(M)] for k, v in out.items() if any(np.any(M) for M in v.values())}
def Jhat(g): return {(c, c): [((), RHO[c, g])] for c in range(len(CONS))}
def ad(g, X): return madd(mmul(Jhat(g), X), mmul(X, Jhat(g)), -1)
AFF = {'A:e1(-2)1': lambda: {k: [(p, (-M) % P) for p, M in l] for k, l in Jhat(1).items()},
       'A:e1(-1)e1': lambda: mmul(Jhat(1), Jhat(1)),
       'A:Cas1': lambda: madd(madd(mmul(Jhat(2), Jhat(1)), mmul(Jhat(1), Jhat(2))), {k: [(p, M * inv(2) % P) for p, M in l] for k, l in mmul(Jhat(0), Jhat(0)).items()}),
       'A:Cas2': lambda: madd(madd(mmul(Jhat(5), Jhat(4)), mmul(Jhat(4), Jhat(5))), {k: [(p, M * inv(2) % P) for p, M in l] for k, l in mmul(Jhat(3), Jhat(3)).items()}),
       'A:e2(-2)1': lambda: {k: [(p, (-M) % P) for p, M in l] for k, l in Jhat(4).items()},
       'A:e2(-1)e2': lambda: mmul(Jhat(4), Jhat(4)),
       'A:e1(-1)e2': lambda: mmul(Jhat(4), Jhat(1))}          # [x(-1)y] = [y]*[x]
ZHAT = {}
for b, (a, i, j) in enumerate(LB):
    lab, nu = HV[a]
    if lab.startswith('A:'): X = AFF[lab]()
    else:
        X = {}
        for (aa, s, t), n in AID.items():
            if aa == a: X[t, s] = [((n,), hw2(CONS[s], CONS[t], tuple(nu)))]
    for _ in range(i): X = ad(2, X)
    for _ in range(j): X = ad(5, X)
    ZHAT[b] = X
SYM = lambda s: ZHAT[s[1]] if s[0] == 'Z' else Jhat(s[1])
log('symbol matrices ready')
# ---- projection of the relations ----
def entry_rows(poly, eta):
    """for each admissible constituent pair: the path polynomial at the normalized entry of X-hat(poly)"""
    out = []
    pairs = [(s, t) for s in range(len(CONS)) for t in range(len(CONS)) if allow(CONS[s][0], CONS[t][0], eta[0]) and allow(CONS[s][1], CONS[t][1], eta[1])]
    for s, t in pairs:
        lam, mu = CONS[s], CONS[t]
        d1 = (eta[0] - mu[0] + lam[0]) // 2; d2 = (eta[1] - mu[1] + lam[1]) // 2
        col = d1 * (lam[1] + 1) + d2
        row = defaultdict(int)
        for mono, c in poly.items():
            state = {(t, ()): np.eye(DIM[t], dtype=np.int64)[0]}
            for sym in mono:
                X = SYM(sym); new = {}
                for (k, p), v in state.items():
                    for (k2, k3), l in X.items():
                        if k2 != k: continue
                        for q, M in l:
                            w = (v @ M) % P
                            if w.any():
                                key = (k3, p + q); new[key] = (new[key] + w) % P if key in new else w
                state = new
                if not state: break
            for (k, p), v in state.items():
                if k == s and v[col]: row[p] = (row[p] + c * int(v[col])) % P
        row = {p: x for p, x in row.items() if x}
        if row: out.append(((VERT[t], VERT[s]), row))
    return out
ROWS = []
for n, (kind, eta, poly) in enumerate(REL):
    for ep, row in entry_rows(poly, eta): ROWS.append((ep, (kind, n), row))
    if n % 200 == 0: log('projected', n, '/', len(REL), 'rows', len(ROWS))
for x, A_ in enumerate(ARROWS):
    s, t = A_['src_c'], A_['tgt_c']
    row = defaultdict(int)
    row[(x,)] += (hA(CONS[t]) - hA(CONS[s])) % P
    row[(ELOOP[t], x)] += 2; row[(x, ELOOP[s])] -= 2
    row = {p: c % P for p, c in row.items() if c % P}
    if row: ROWS.append(((A_['target'], A_['source']), ('central', x), row))
log('rows:', len(ROWS))
pickle.dump({'ROWS': ROWS, 'ARROWS': ARROWS, 'CONS': CONS, 'VERT': VERT, 'P': P}, open(ZD + '/corner_rows_modp.pkl', 'wb'))
log('done')
