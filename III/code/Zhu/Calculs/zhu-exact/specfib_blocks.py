"""Block-wise action of the central class Z modulo p on the special fibre D_p (Proposition 6.1, part 4).
Uses the presentation and the element Z computed by specfib_Z.py (specfib/specfib_Z.pkl).  For each eigenvalue c of
the left multiplication L_Z (Z central, D_p semisimple, so that ker(L_Z - c) is the sum of the blocks on which Z acts
by c), computes dim ker(L_Z - c) and, for each vertex v, dim e_v ker(L_Z - c) e_v = sum over these blocks of the
squares of the ranks of e_v.  When one block carries c, this gives its vertex ranks."""
import pickle, json, time
from functools import lru_cache
from fractions import Fraction as Fr
from flint import nmod_mat
t0 = time.time()
def log(*a): print(*a, '(%.0fs)' % (time.time() - t0), flush=True)
ZD = '/home/claude/hj/work/Calculs/zhu-exact/specfib'
D = pickle.load(open(ZD + '/corner_rows_modp.pkl', 'rb')); P = D['P']; AR = D['ARROWS']
S = pickle.load(open(ZD + '/specfib_Z.pkl', 'rb'))
basis, rules, Zel, roots, verts = S['basis'], S['rules'], S['Zel'], S['roots'], S['verts']
def src(a): return AR[a]['source']
def tgt(a): return AR[a]['target']
def add(out, row, c=1):
    for p, v in row.items():
        z = (out.get(p, 0) + c * v) % P
        if z: out[p] = z
        else: out.pop(p, None)
LENS = sorted({len(q) for q in rules}, reverse=True)
@lru_cache(None)
def nf(p):
    for L in LENS:
        for i in range(len(p) - L + 1):
            q = p[i:i + L]
            if q in rules:
                o = {}
                for w, c in rules[q].items(): add(o, nf(p[:i] + w + p[i + L:]), c)
                return o
    return {p: 1}
N = len(basis)
bidx = {}
for i, b in enumerate(basis): bidx[() if b[0] == 'v' else b, b[1] if b[0] == 'v' else None] = i
def ends(b): return (b[1], b[1]) if b[0] == 'v' else (tgt(b[0]), src(b[-1]))
def mul_elt_basis(el, b):
    out = {}
    eb = ends(basis[b]); pb = () if basis[b][0] == 'v' else basis[b]
    for k, c in el.items():
        if k[0] == 'v':
            if k[1] != eb[0]: continue
            add(out, {b: 1}, c); continue
        if src(k[-1]) != eb[0]: continue
        for w, cc in nf(k + pb).items():
            i = bidx[(w, tgt(k[0])) if w == () else (w, None)]
            out[i] = (out.get(i, 0) + c * cc) % P
    return {i: c for i, c in out.items() if c}
M = [[0] * N for _ in range(N)]
for b in range(N):
    for i, c in mul_elt_basis(Zel, b).items(): M[i][b] = (M[i][b] + c) % P
def ratrec(a, bound=10**3):
    r0, r1, s0, s1 = P, a % P, 0, 1
    while r1 > bound:
        q = r0 // r1; r0, r1 = r1, r0 - q * r1; s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound: return None
    return Fr(r1, s1) if s1 > 0 else Fr(-r1, -s1)
endp = [ends(b) for b in basis]
out = []
for r, e in roots:
    A = nmod_mat(N, N, [(M[i][j] - (r if i == j else 0)) % P for i in range(N) for j in range(N)], P)
    K = A.nullspace()[0]
    kdim = N - A.rank()
    cols = [[int(K[i, j]) for i in range(N)] for j in range(kdim)]     # basis vectors of the kernel
    ranks = {}
    for v in verts:
        idx = [i for i in range(N) if endp[i] == (v, v)]
        if not idx: continue
        sub = nmod_mat(kdim, len(idx), [cols[j][i] for j in range(kdim) for i in idx], P)
        rk = sub.rank()
        if rk: ranks[v] = rk
    out.append({'root': r, 'rational': str(ratrec(r)), 'kernel_dim': kdim, 'multiplicity_charpoly': e, 'corner_dims': ranks})
    log('c =', ratrec(r), 'dim ker', kdim, 'corner dims', ranks)
json.dump(out, open(ZD + '/specfib_blocks.json', 'w'), indent=1)
log('done')
