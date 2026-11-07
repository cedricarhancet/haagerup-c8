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
ZD = '/home/claude/work/Calculs/zhu'
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
# ---- structure ----
# irreducible words of all lengths (Bergman: a basis once every ambiguity resolves)
LMAX = max(len(q) for q in rules)
def irreducible(w): return all(w[i:i + L] not in rules for L in range(2, LMAX + 1) for i in range(max(0, len(w) - L), len(w) - L + 1))
words = [[(a,) for a in active]]
while words[-1] and len(words) < 12:
    nxt = [w + (a,) for w in words[-1] for a in active if src(w[-1]) == tgt(a) and irreducible(w + (a,))]
    words.append(nxt)
log('irreducible words by length:', [len(w) for w in words])
assert not words[-1], 'infinitely many irreducible words?'
normalw = [w for ws in words[1:] for w in ws]
basis = [('v', v) for v in verts] + [(a,) for a in active] + normalw
bidx = {}
for i, b in enumerate(basis): bidx[() if b[0] == 'v' else b, b[1] if b[0] == 'v' else None] = i
def ends(b):
    if b[0] == 'v': return (b[1], b[1])
    return (tgt(b[0]), src(b[-1]))
def vec(poly, ep):
    z = {}
    for p, c in poly.items():
        k = (p, ep[0]) if p == () else (p, None)
        z[bidx[k]] = (z.get(bidx[k], 0) + c) % P
    return z
def mult(x, y):
    ex, ey = ends(basis[x]), ends(basis[y])
    if ex[1] != ey[0]: return {}
    px = () if basis[x][0] == 'v' else basis[x]; py = () if basis[y][0] == 'v' else basis[y]
    return vec(nf(px + py), (ex[0], ey[1]))
N = len(basis); log('basis of irreducible words: %d = %d vertices + %d arrows + %d longer words' % (N, len(verts), len(active), len(normalw)))
tau = np.zeros(N, np.int64)
for b in range(N):
    for c in range(N):
        tau[b] = (tau[b] + mult(b, c).get(c, 0)) % P
G = np.zeros((N, N), np.int64)
for x in range(N):
    for y in range(N):
        for z, c in mult(x, y).items(): G[x, y] = (G[x, y] + c * tau[z]) % P
rk = nmod_mat(N, N, [int(v) for v in G.flatten()], P).rank()
log('regular trace form: rank %d of %d' % (rk, N))
# center: x with x g = g x for every generator g (vertices and arrows), incremental echelon
gens_idx = [i for i, b in enumerate(basis) if b[0] == 'v' or len(b) == 1]
MULT = {}
def mult_c(x, y):
    k = (x, y)
    if k not in MULT: MULT[k] = mult(x, y)
    return MULT[k]
cur = None
for g in gens_idx:
    M = np.zeros((N, N), np.int64)
    for x in range(N):
        for z, c in mult(x, g).items(): M[z, x] = (M[z, x] + c) % P
        for z, c in mult(g, x).items(): M[z, x] = (M[z, x] - c) % P
    M = M[np.any(M, axis=1)]
    if cur is not None: M = np.vstack([cur, M])
    if M.shape[0] == 0: continue
    R, r = nmod_mat(M.shape[0], N, [int(v) for v in M.flatten()], P).rref()
    cur = np.array([[int(R[i, j]) for j in range(N)] for i in range(r)], np.int64)
X, nul = nmod_mat(cur.shape[0], N, [int(v) for v in cur.flatten()], P).nullspace()
Zb = [np.array([int(X[i, j]) for i in range(N)], np.int64) for j in range(nul)]
log('center dimension', nul)
def mulvec(u, v):
    out = np.zeros(N, np.int64)
    for x in np.nonzero(u)[0]:
        for y in np.nonzero(v)[0]:
            uv = int(u[x]) * int(v[y]) % P
            for z, c in mult(int(x), int(y)).items(): out[z] = (int(out[z]) + uv * c) % P
    return out
random.seed(1)
zr = sum((random.randrange(P) * z for z in Zb), np.zeros(N, np.int64)) % P
Zm = np.array(Zb, np.int64)
R_, rank_ = nmod_mat(len(Zb), N, [int(v) for v in Zm.flatten()], P).rref(); pc = []
for i in range(rank_):
    for j in range(N):
        if int(R_[i, j]): pc.append(j); break
Si = nmod_mat(len(Zb), len(pc), [int(Zm[i, j]) for i in range(len(Zb)) for j in pc], P).inv()
def coordsZ(v): r = nmod_mat(1, len(pc), [int(v[j]) for j in pc], P) * Si; return [int(r[0, i]) for i in range(len(Zb))]
Lz = nmod_mat(len(Zb), len(Zb), sum([coordsZ(mulvec(zr, z)) for z in Zb], []), P).transpose()
fac = Lz.charpoly().factor()
log('charpoly of a random central element factors as', [(f.degree(), e) for f, e in fac[1]])
one = np.zeros(N, np.int64); one[[bidx[((), v)] for v in verts]] = 1
blocks = []
roots = [(-int(f.coeffs()[0])) % P for f, e in fac[1] if f.degree() == 1 and e == 1]
if len(roots) == len(Zb):
    for r0 in roots:
        f = one.copy()
        for r1 in roots:
            if r1 != r0: f = mulvec(f, (zr - r1 * one) % P) * pow((r0 - r1) % P, -1, P) % P
        Lf = np.array([mulvec(f, np.eye(N, dtype=np.int64)[y]) for y in range(N)], np.int64)
        dimb = nmod_mat(N, N, [int(v) for v in Lf.flatten()], P).rank()
        ranks = {}
        for v in verts:
            ws = [y for y, b in enumerate(basis) if ends(b) == (v, v)]
            Fv = np.array([mulvec(f, np.eye(N, dtype=np.int64)[y]) for y in ws], np.int64)
            rr = nmod_mat(Fv.shape[0], N, [int(x) for x in Fv.flatten()], P).rank() if len(ws) else 0
            if rr: ranks[v] = int(round(rr ** 0.5))
        blocks.append({'dim': dimb, 'size': int(round(dimb ** 0.5)), 'vertex_ranks': ranks})
    log('blocks:', sorted(b['size'] for b in blocks))
out = {'status': 'modular discovery (p = %d); not a proof' % P, 'saturation': hist, 'killed_vertices': sorted(killed), 'vertices': verts,
       'active_arrows': len(active), 'normal2': len(normal2), 'dimension': N, 'trace_form_rank': rk, 'center_dim': nul, 'blocks': blocks}
json.dump(out, open(ZD + '/z4_upper_algebra_modp.json', 'w'), indent=1)
pickle.dump({'rules': rules, 'active': active, 'verts': verts, 'basis': basis}, open(ZD + '/z4_upper_algebra_modp.pkl', 'wb'))
log('done')
