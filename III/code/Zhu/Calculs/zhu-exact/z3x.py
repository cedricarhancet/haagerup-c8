"""EXACT projection of the exact Zhu relations (rel_exact.pkl) to the Morita corner of A(S) given by the matrix units of
A(A+), A+ = A + J (16 constituents, 14 vertices), as in z3 of Article I but over K.  The sl_2 intertwiner matrices
are rational (computed over Q); the relation coefficients lie in K.  Output: rows_exact.pkl (exact rows: endpoints,
origin, path polynomial with K coefficients) and rows_exact_reduced_modp.pkl (reduction at the prime of the
coherent embedding, same format as corner_rows_modp.pkl)."""
import pickle, json, time, sys
from collections import defaultdict
from functools import lru_cache
from fractions import Fraction as Fr
sys.path.insert(0, '/home/claude/hj/work/Calculs/zhu-exact')
import kexact as KX
from gmpy2 import mpq
import Kfield as K
from flint import fmpq_mat, fmpq, fmpz_mat
t0 = time.time()
def log(*a): print(*a, '(%.0fs)' % (time.time() - t0), flush=True)
HERE = '/home/claude/hj/work/Calculs/zhu-exact'
D = pickle.load(open(HERE + '/rel_exact.pkl', 'rb'))
REL = D['REL']; HV = D['HV']; TYPE = D['TYPE']; LB = D['LB']
CONS = [(0, 0), (0, 2), (1, 1), (2, 0), (0, 4), (1, 3), (2, 2), (3, 1), (0, 6), (1, 5), (2, 4), (3, 3), (1, 7), (2, 6), (0, 8), (3, 5)]
VERT = list(range(12)) + [12, 12, 13, 13]
NV = 14
DIM = [(r + 1) * (l + 1) for r, l in CONS]
def hA(t): r, l = t; return Fr(r * (r + 2), 20) + Fr(l * (l + 2), 60)
def sl2(n):
    e = [[0] * (n + 1) for _ in range(n + 1)]; f = [[0] * (n + 1) for _ in range(n + 1)]; h = [[0] * (n + 1) for _ in range(n + 1)]
    for i in range(n + 1):
        h[i][i] = n - 2 * i
        if i < n: f[i + 1][i] = 1
        if i > 0: e[i - 1][i] = i * (n - i + 1)
    return e, f, h
def kron(A, B):
    mB, nB = len(B), len(B[0])
    return [[A[i // mB][j // nB] * B[i % mB][j % nB] for j in range(len(A[0]) * nB)] for i in range(len(A) * mB)]
def eye(n): return [[1 if i == j else 0 for j in range(n)] for i in range(n)]
def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(len(B)) if A[i][k]) for j in range(len(B[0]))] for i in range(len(A))]
def rho(t, g):
    r, l = t; e1, f1, h1 = sl2(r); e2, f2, h2 = sl2(l); I1 = eye(r + 1); I2 = eye(l + 1)
    return {0: kron(h1, I2), 1: kron(e1, I2), 2: kron(f1, I2), 3: kron(I1, h2), 4: kron(I1, e2), 5: kron(I1, f2)}[g]
RHO = {(c, g): rho(CONS[c], g) for c in range(len(CONS)) for g in range(6)}
def allow(m, n, l): return abs(m - n) <= l <= m + n and (m + n - l) % 2 == 0
@lru_cache(None)
def hw1(n, m, l):
    """normalized highest matrix F in Hom(V(n), V(m)) of weight l: e_m F = F e_n, weight l, F[0, dlt] = 1 (exact)"""
    e_n, f_n, h_n = sl2(n); e_m, f_m, h_m = sl2(m)
    dlt = (l - m + n) // 2
    N = (m + 1) * (n + 1); rows = []
    for i in range(m + 1):
        for j in range(n + 1):
            row = [0] * N
            for k in range(m + 1): row[k * (n + 1) + j] += e_m[i][k]
            for k in range(n + 1): row[i * (n + 1) + k] -= e_n[k][j]
            rows.append(row)
    for i in range(m + 1):
        for j in range(n + 1):
            if (m - 2 * i) - (n - 2 * j) != l:
                row = [0] * N; row[i * (n + 1) + j] = 1; rows.append(row)
    M = fmpz_mat(len(rows), N, [x for r_ in rows for x in r_])
    X, nul = M.nullspace(); assert nul == 1, (n, m, l, nul)
    F = [[Fr(int(X[k * (n + 1) + j, 0])) for j in range(n + 1)] for k in range(m + 1)]
    c = F[0][dlt]; assert c, (n, m, l)
    return tuple(tuple(x / c for x in row) for row in F)
def hw2(lam, mu, nu): return kron([list(r) for r in hw1(lam[0], mu[0], nu[0])], [list(r) for r in hw1(lam[1], mu[1], nu[1])])
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
def madd(X, Y, c=1):
    out = defaultdict(dict)
    for k, l in X.items():
        for p, M in l: out[k][p] = [list(r) for r in M]
    for k, l in Y.items():
        for p, M in l:
            cur = out[k].get(p)
            out[k][p] = [[(cur[i][j] if cur else 0) + c * M[i][j] for j in range(len(M[0]))] for i in range(len(M))]
    return {k: [(p, M) for p, M in v.items() if any(x for r in M for x in r)] for k, v in out.items()}
def mmul(X, Y):
    out = defaultdict(dict)
    for (t, m), lx in X.items():
        for (m2, s), ly in Y.items():
            if m != m2: continue
            for px, Mx in lx:
                for py, My in ly:
                    Zm = matmul(Mx, My)
                    if not any(x for r in Zm for x in r): continue
                    cur = out[t, s].get(px + py)
                    out[t, s][px + py] = Zm if cur is None else [[cur[i][j] + Zm[i][j] for j in range(len(Zm[0]))] for i in range(len(Zm))]
    return {k: [(p, M) for p, M in v.items() if any(x for r in M for x in r)] for k, v in out.items()}
def Jhat(g): return {(c, c): [((), RHO[c, g])] for c in range(len(CONS))}
def ad(g, X): return madd(mmul(Jhat(g), X), mmul(X, Jhat(g)), -1)
def scal(X, c): return {k: [(p, [[c * x for x in r] for r in M]) for p, M in l] for k, l in X.items()}
AFF = {'A:e1(-2)1': lambda: scal(Jhat(1), -1), 'A:e1(-1)e1': lambda: mmul(Jhat(1), Jhat(1)),
       'A:Cas1': lambda: madd(madd(mmul(Jhat(2), Jhat(1)), mmul(Jhat(1), Jhat(2))), scal(mmul(Jhat(0), Jhat(0)), Fr(1, 2))),
       'A:Cas2': lambda: madd(madd(mmul(Jhat(5), Jhat(4)), mmul(Jhat(4), Jhat(5))), scal(mmul(Jhat(3), Jhat(3)), Fr(1, 2))),
       'A:e2(-2)1': lambda: scal(Jhat(4), -1), 'A:e2(-1)e2': lambda: mmul(Jhat(4), Jhat(4)),
       'A:e1(-1)e2': lambda: mmul(Jhat(4), Jhat(1))}
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
log('symbol matrices ready (exact)')
_mono_cache = {}
def mono_rows(mono, eta):
    """for a monomial: dict (s,t) -> {path: rational} (the normalized entry), cached"""
    key = (mono, eta)
    if key in _mono_cache: return _mono_cache[key]
    res = {}
    pairs = [(s, t) for s in range(len(CONS)) for t in range(len(CONS)) if allow(CONS[s][0], CONS[t][0], eta[0]) and allow(CONS[s][1], CONS[t][1], eta[1])]
    for s, t in pairs:
        lam, mu = CONS[s], CONS[t]
        d1 = (eta[0] - mu[0] + lam[0]) // 2; d2 = (eta[1] - mu[1] + lam[1]) // 2
        col = d1 * (lam[1] + 1) + d2
        state = {(t, ()): [1 if i == 0 else 0 for i in range(DIM[t])]}
        for sym in mono:
            X = SYM(sym); new = {}
            for (k, p), v in state.items():
                for (k2, k3), l in X.items():
                    if k2 != k: continue
                    for q, M in l:
                        w = [sum(v[i] * M[i][j] for i in range(len(v)) if v[i]) for j in range(len(M[0]))]
                        if any(w):
                            key2 = (k3, p + q)
                            if key2 in new: new[key2] = [x + y for x, y in zip(new[key2], w)]
                            else: new[key2] = w
            state = new
            if not state: break
        row = {}
        for (k, p), v in state.items():
            if k == s and v[col]: row[p] = row.get(p, 0) + v[col]
        row = {p: x for p, x in row.items() if x}
        if row: res[s, t] = row
    _mono_cache[key] = res
    return res
def kq(x): x = Fr(x); return K.scal(mpq(x.numerator, x.denominator), K.ONE)
ROWS = []
for n, (kind, eta, poly) in enumerate(REL):
    acc = defaultdict(dict)
    for mono, c in poly.items():
        for (s, t), row in mono_rows(mono, eta).items():
            dd = acc[s, t]
            for p, x in row.items():
                y = K.add(dd.get(p, K.Z), K.mul(c, kq(x)))
                if K.iszero(y): dd.pop(p, None)
                else: dd[p] = y
    for (s, t), row in acc.items():
        if row: ROWS.append(((VERT[t], VERT[s]), (kind, n), row))
    if n % 100 == 0: log('projected', n, '/', len(REL), 'rows', len(ROWS))
for x, A_ in enumerate(ARROWS):
    s, t = A_['src_c'], A_['tgt_c']
    row = defaultdict(lambda: Fr(0))
    row[(x,)] += hA(CONS[t]) - hA(CONS[s])
    row[(ELOOP[t], x)] += 2; row[(x, ELOOP[s])] -= 2
    row = {p: kq(c) for p, c in row.items() if c}
    if row: ROWS.append(((A_['target'], A_['source']), ('central', x), row))
log('rows:', len(ROWS))
bad = sum(1 for _, _, r in ROWS for c in r.values() if not KX.den_ok(c))
log('coefficients with p in a denominator:', bad)
pickle.dump({'ROWS': ROWS, 'ARROWS': ARROWS, 'CONS': CONS, 'VERT': VERT}, open(HERE + '/rows_exact.pkl', 'wb'))
RP = []
for ep, o, r in ROWS:
    rr = {p: KX.kred(c) for p, c in r.items() if KX.kred(c)}
    if rr: RP.append((ep, o, rr))
pickle.dump({'ROWS': RP, 'ARROWS': ARROWS, 'CONS': CONS, 'VERT': VERT, 'P': KX.P}, open(HERE + '/rows_exact_reduced_modp.pkl', 'wb'))
json.dump({'rows': len(ROWS), 'p_denominators': bad, 'arrows': len(ARROWS)}, open(HERE + '/z3x.json', 'w'), indent=1)
log('done')
