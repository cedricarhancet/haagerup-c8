"""Modular discovery (p3): upper presentation of the Zhu corner from the projected rows.
 1. echelon per endpoint corner (leading term = longest, then lexicographically smallest path), linear eliminations
    and vertex exclusions iterated to stability (quadratic saturation, as in Riesen's C13 modular discovery);
 2. completion of the quadratic rewriting system by critical overlaps (a,b),(b,c) until no new relation of length <= 2
    appears and every normal path of length three is a leading term (finite spanning set of length <= 2);
 3. diamond check: all overlaps of the final system resolve; regular trace form; center; blocks and e_v-ranks."""
import pickle, json, time, random
from collections import defaultdict
from functools import lru_cache
import numpy as np
from flint import nmod_mat, nmod_poly
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
def reduce_by(tab, v):
    while v:
        lead = min(v, key=key)
        if lead not in tab: return v, lead
        add(v, tab[lead], -v[lead])
    return v, None
def echelon(rows):
    tabs = defaultdict(dict)
    for ep, origin, raw in rows:
        v, lead = reduce_by(tabs[ep], dict(raw))
        if v:
            iv = pow(v[lead], -1, P); tabs[ep][lead] = {p: c * iv % P for p, c in v.items()}
            # keep table reduced on the new leading term (back substitution not needed for correctness)
    return tabs
def src(a): return AR[a]['source']
def tgt(a): return AR[a]['target']
rows = [(ep, o, r) for ep, o, r in D['ROWS']]
log('input rows', len(rows), 'arrows', len(AR))
killed = set(); alllin = {}; hist = []
for stage in range(40):
    tabs = echelon(rows)
    lin = {(ep, p): r for ep, t in tabs.items() for p, r in t.items() if len(p) <= 1}
    newk = {ep[0] for ep, p in lin if not p} - killed
    killed |= newk
    lr = {p[0]: {w: (-c) % P for w, c in r.items() if w != p} for (ep, p), r in lin.items() if p}
    hist.append({'stage': stage, 'rank': sum(map(len, tabs.values())), 'linear': len(lr), 'killed': sorted(killed)})
    log('saturation', hist[-1])
    if not lr and not newk: break
    alllin.update(lr)
    @lru_cache(None)
    def nf(p):
        if any(src(a) in killed or tgt(a) in killed for a in p): return {}
        for i, a in enumerate(p):
            if a in lr:
                out = {}
                for w, c in lr[a].items(): add(out, nf(p[:i] + w + p[i + 1:]), c)
                return out
        return {p: 1}
    nrows = []
    for ep, tab in tabs.items():
        if set(ep) & killed: continue
        for lead, r in tab.items():
            raw = {}
            for p, c in r.items(): add(raw, nf(p), c)
            if raw: nrows.append((ep, ('sat', stage), raw))
    rows = nrows
active = [a for a in range(len(AR)) if a not in alllin and src(a) not in killed and tgt(a) not in killed]
verts = [v for v in range(NVERT) if v not in killed]
log('after saturation: vertices', verts, 'active arrows', len(active))
# ---- completion ----
for stage in range(40):
    tabs = echelon(rows)
    lr = {p[0]: {w: (-c) % P for w, c in r.items() if w != p} for ep, t in tabs.items() for p, r in t.items() if len(p) == 1}
    k2 = {ep[0] for ep, t in tabs.items() if () in t}
    if lr or k2:
        killed |= k2; verts = [v for v in verts if v not in killed]
        @lru_cache(None)
        def lnf(p):
            if any(src(a) in killed or tgt(a) in killed for a in p): return {}
            for i, a in enumerate(p):
                if a in lr:
                    o = {}
                    for w, c in lr[a].items(): add(o, lnf(p[:i] + w + p[i + 1:]), c)
                    return o
            return {p: 1}
        active = [a for a in active if a not in lr and src(a) not in killed and tgt(a) not in killed]
        alllin.update(lr)
        rows = []
        for ep, t in tabs.items():
            if set(ep) & killed: continue
            for p, r in t.items():
                raw = {}
                for w, c in r.items(): add(raw, lnf(w), c)
                if raw: rows.append((ep, ('lin', stage), raw))
        log('completion: linear stage', stage, 'active', len(active), 'killed', sorted(killed)); continue
    rules = {p: {w: (-c) % P for w, c in r.items() if w != p} for ep, t in tabs.items() for p, r in t.items()}
    assert all(len(p) >= 2 for p in rules)
    @lru_cache(None)
    def nf(p):
        for L in sorted({len(q) for q in rules}, reverse=True):
            for i in range(len(p) - L + 1):
                q = p[i:i + L]
                if q in rules:
                    o = {}
                    for w, c in rules[q].items(): add(o, nf(p[:i] + w + p[i + L:]), c)
                    return o
        return {p: 1}
    starts = defaultdict(list)
    for q in rules: starts[q[0]].append(q)
    new = []; n_ov = 0
    for q1 in list(rules):
        for k in range(1, len(q1)):            # overlap: suffix of q1 = prefix of q2
            suf = q1[k:]
            for q2 in starts[suf[0]]:
                if len(q2) <= len(suf) or q2[:len(suf)] != suf: continue
                n_ov += 1
                w = q1 + q2[len(suf):]
                v = {}
                for x, c in rules[q1].items(): add(v, nf(x + q2[len(suf):]), c)
                for x, c in rules[q2].items(): add(v, nf(q1[:k] + x), -c)
                if v: new.append(((tgt(w[0]), src(w[-1])), ('overlap', stage), v))
    normal2 = [(a, b) for a in active for b in active if src(a) == tgt(b) and (a, b) not in rules]
    normal3 = [(a, b, c) for a, b in normal2 for c in active if src(b) == tgt(c) and (b, c) not in rules and (a, b, c) not in rules]
    log('completion stage', stage, 'rules', len(rules), 'overlaps', n_ov, 'nonresolving', len(new), 'normal2', len(normal2), 'normal3 (not leading)', len(normal3))
    if not new: break
    rows = [(ep, ('rule',), {**{p: 1}, **{w: (-c) % P for w, c in r.items()}}) for p, r in rules.items() for ep in [(tgt(p[0]), src(p[-1]))]] + new
    if stage == 39: log('completion did not stabilize')

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
