"""Compute the action of Z in the verified modular normal basis.

The modular completion has already been computed. Its rules are compared
with the reductions of the exact rules before using the exact eliminations.
"""
import os
import pickle
import sys
import time
from functools import lru_cache
from flint import nmod_mat

t0 = time.time()
def log(*a): print(*a, '(%.0fs)' % (time.time() - t0), flush=True)
import os
ZD = '/home/claude/hj/work/Calculs/zhu-exact/specfib'
if os.environ.get('HAAG_PDIR'): ZD = os.environ['HAAG_PDIR']
D = pickle.load(open(ZD + '/corner_rows_modp.pkl', 'rb')); P = D['P']; AR = D['ARROWS']; NVERT = 14
def key(p): return (-len(p), p)
def add(out, row, c=1):
    for p, v in row.items():
        z = (out.get(p, 0) + c * v) % P
        if z: out[p] = z
        else: out.pop(p, None)
sys.path.insert(0, '/home/claude/hj/work/Calculs/zhu-exact')
import kexact as KX

source = pickle.load(open(ZD + '/z4_upper_algebra_modp.pkl', 'rb'))
exact = pickle.load(open(ZD + '/../upper_exact.pkl', 'rb'))
def reduce_coefficients(poly):
    result = {}
    for word, c in poly.items():
        reduced = KX.kred(c)
        if reduced:
            result[word] = reduced
    return result

rules = source['rules']
assert {q: reduce_coefficients(rhs) for q, rhs in exact['rules'].items()} == rules
active = source['active']
verts = source['verts']
assert active == exact['active'] and verts == exact['verts']
alllin = {a: reduce_coefficients(rhs) for a, rhs in exact['alllin'].items()}
killed = set(exact['killed'])
def src(a): return AR[a]['source']
def tgt(a): return AR[a]['target']

LENS = sorted({len(q) for q in rules}, reverse=True)
@lru_cache(None)
def nf(p):
    for length in LENS:
        for i in range(len(p) - length + 1):
            rhs = rules.get(p[i:i + length])
            if rhs is not None:
                out = {}
                for w, c in rhs.items():
                    add(out, nf(p[:i] + w + p[i + length:]), c)
                return out
    return {p: 1}

log('verified exact reductions, using the recomputed modular rules')
# ---- action of the central class Z = [omega] on the corner, modulo p ----
from fractions import Fraction as Fr
LMAX = max(len(q) for q in rules)
def irreducible(w): return all(w[i:i + L] not in rules for L in range(2, LMAX + 1) for i in range(max(0, len(w) - L), len(w) - L + 1))
words = [[(a,) for a in active]]
while words[-1] and len(words) < 12:
    words.append([w + (a,) for w in words[-1] for a in active if src(w[-1]) == tgt(a) and irreducible(w + (a,))])
normalw = [w for ws in words[1:] for w in ws]
basis = [('v', v) for v in verts] + [(a,) for a in active] + normalw
N = len(basis); log('basis', N)
bidx = {}
for i, b in enumerate(basis): bidx[() if b[0] == 'v' else b, b[1] if b[0] == 'v' else None] = i
def ends(b): return (b[1], b[1]) if b[0] == 'v' else (tgt(b[0]), src(b[-1]))
CONS = D['CONS']; VERT = D['VERT']
def hA(t): r, l = t; return Fr(r * (r + 2), 20) + Fr(l * (l + 2), 60)
def fq(x): return x.numerator * pow(x.denominator, -1, P) % P
# loop arrows of the vector e at each constituent
eloop = {}
for x, A_ in enumerate(AR):
    if A_['z'] == 'P:e' and A_['src_c'] == A_['tgt_c']: eloop[A_['src_c']] = x
# full normal form of a path including linear eliminations
@lru_cache(None)
def fullnf(p):
    if any(src(a) in killed or tgt(a) in killed for a in p): return {}
    for i, a in enumerate(p):
        if a in alllin:
            o = {}
            for w, c in alllin[a].items(): add(o, fullnf(p[:i] + w + p[i + 1:]), c)
            return o
    return nf(p)
# Z as a combination of paths: sum over vertices v of h_A e_v plus 2 * loops (constituents of active vertices)
Zel = {}   # key: ('v', v) or path
for v in verts:
    cs = [c for c in range(len(CONS)) if VERT[c] == v]
    Zel[('v', v)] = fq(hA(CONS[cs[0]]))
for c, x in eloop.items():
    if VERT[c] in verts:
        for w, cc in fullnf((x,)).items():
            k = ('v', tgt(x)) if w == () else w
            Zel[k] = (Zel.get(k, 0) + 2 * cc) % P
def mul_elt_basis(el, b):     # el (dict) times basis element b, returns dict basis index -> coeff
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
M = nmod_mat(N, N, [0] * (N * N), P)
for b in range(N):
    for i, c in mul_elt_basis(Zel, b).items(): M[i, b] = (int(M[i, b]) + c) % P
cp = M.charpoly(); fac = cp.factor()
roots = []
for f, e in fac[1]:
    if f.degree() == 1: roots.append(((-int(f.coeffs()[0])) % P, e))
    else: roots.append(('deg%d' % f.degree(), e))
log('charpoly of L_Z factors:', roots)
def ratrec(a, bound=10**3):
    # rational reconstruction of a mod P with |num|, den <= bound
    r0, r1, s0, s1 = P, a % P, 0, 1
    while r1 > bound:
        q = r0 // r1; r0, r1 = r1, r0 - q * r1; s0, s1 = s1, s0 - q * s1
    if s1 == 0 or abs(s1) > bound: return None
    return Fr(r1, s1) if s1 > 0 else Fr(-r1, -s1)
for r, e in roots:
    if isinstance(r, int): log('root', r, 'multiplicity', e, 'rational', ratrec(r))
pickle.dump({'roots': roots, 'Zel': Zel, 'basis': basis, 'alllin': alllin, 'rules': rules, 'killed': sorted(killed), 'verts': verts, 'active': active}, open(ZD + '/specfib_Z.pkl', 'wb'))
log('done')
