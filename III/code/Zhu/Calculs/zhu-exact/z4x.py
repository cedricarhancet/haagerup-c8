"""EXACT upper presentation of the Zhu corner over K from the exact projected rows (rows_exact.pkl): the algorithm of
z4 (echelon per endpoint corner with leading term = longest then lexicographically smallest path, linear eliminations
and vertex exclusions to stability, completion by critical overlaps) replayed over K.  Every division is by a leading
coefficient; each leading coefficient is recorded with its reduction at the prime of the coherent embedding.
Output: upper_exact.pkl (exact eliminations, killed vertices, rules) and z4x.json (integrality and unit checks,
irreducible word counts, exact confluence)."""
import pickle, json, time, sys
from collections import defaultdict
from functools import lru_cache
sys.path.insert(0, '/home/claude/hj/work/Calculs/zhu-exact')
import kfast as KF
import kexact as KX
t0 = time.time()
def log(*a): print(*a, '(%.0fs)' % (time.time() - t0), flush=True)
HERE = '/home/claude/hj/work/Calculs/zhu-exact'
D = pickle.load(open(HERE + '/rows_exact.pkl', 'rb')); AR = D['ARROWS']; NVERT = 14
P = KX.P
STATS = {'lead_nonunit': 0, 'lead_total': 0, 'bad_denominators': 0}
def key(p): return (-len(p), p)
def add(out, row, c=KF.ONE):
    for p, v in row.items():
        z = KF.add(out.get(p, KF.ZERO), KF.mul(c, v))
        if KF.iszero(z): out.pop(p, None)
        else: out[p] = z
def reduce_by(tab, v):
    while v:
        lead = min(v, key=key)
        if lead not in tab: return v, lead
        add(v, tab[lead], KF.neg(v[lead]))
    return v, None
def unit_check(c):
    STATS['lead_total'] += 1
    if not (KX.den_ok(c) and KX.kred(c) != 0): STATS['lead_nonunit'] += 1
def echelon(rows):
    tabs = defaultdict(dict)
    for ep, origin, raw in rows:
        v, lead = reduce_by(tabs[ep], dict(raw))
        if v:
            c = v[lead]; unit_check(c); iv = KF.inv(c)
            tabs[ep][lead] = {p: KF.mul(x, iv) for p, x in v.items()}
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
    lr = {p[0]: {w: KF.neg(c) for w, c in r.items() if w != p} for (ep, p), r in lin.items() if p}
    hist.append({'stage': stage, 'rank': sum(map(len, tabs.values())), 'linear': len(lr), 'killed': sorted(killed)})
    log('saturation', hist[-1], STATS)
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
        return {p: KF.ONE}
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
pickle.dump({'alllin': alllin, 'killed': killed, 'rows': rows, 'active': active, 'verts': verts, 'hist': hist}, open(HERE + '/upper_exact_saturation.pkl', 'wb'))
comp = []
for stage in range(40):
    tabs = echelon(rows)
    lr = {p[0]: {w: KF.neg(c) for w, c in r.items() if w != p} for ep, t in tabs.items() for p, r in t.items() if len(p) == 1}
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
            return {p: KF.ONE}
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
    rules = {p: {w: KF.neg(c) for w, c in r.items() if w != p} for ep, t in tabs.items() for p, r in t.items()}
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
        return {p: KF.ONE}
    starts = defaultdict(list)
    for q in rules: starts[q[0]].append(q)
    new = []; n_ov = 0
    for q1 in list(rules):
        for k in range(1, len(q1)):
            suf = q1[k:]
            for q2 in starts[suf[0]]:
                if len(q2) <= len(suf) or q2[:len(suf)] != suf: continue
                n_ov += 1
                w = q1 + q2[len(suf):]
                v = {}
                for x, c in rules[q1].items(): add(v, nf(x + q2[len(suf):]), c)
                for x, c in rules[q2].items(): add(v, nf(q1[:k] + x), KF.neg(c))
                if v: new.append(((tgt(w[0]), src(w[-1])), ('overlap', stage), v))
    normal2 = [(a, b) for a in active for b in active if src(a) == tgt(b) and (a, b) not in rules]
    normal3 = [(a, b, c) for a, b in normal2 for c in active if src(b) == tgt(c) and (b, c) not in rules and (a, b, c) not in rules]
    comp.append({'stage': stage, 'rules': len(rules), 'overlaps': n_ov, 'nonresolving': len(new), 'normal2': len(normal2), 'normal3': len(normal3)})
    log('completion', comp[-1], STATS)
    if not new: break
    rows = [(ep, ('rule',), {**{p: KF.ONE}, **{w: KF.neg(c) for w, c in r.items()}}) for p, r in rules.items() for ep in [(tgt(p[0]), src(p[-1]))]] + new
# inclusion ambiguities (a leading word inside a longer leading word) must also resolve (Bergman)
incl = 0; incl_bad = 0
for q1 in rules:
    for q2 in rules:
        if len(q2) < len(q1):
            for i in range(len(q1) - len(q2) + 1):
                if q1[i:i + len(q2)] == q2:
                    incl += 1
                    v = {}
                    for x, c in rules[q1].items(): add(v, nf(x), c)
                    for x, c in rules[q2].items(): add(v, nf(q1[:i] + x + q1[i + len(q2):]), KF.neg(c))
                    if v: incl_bad += 1
log('inclusion ambiguities', incl, 'unresolved', incl_bad)
# integrality of all final data at the prime
def all_coeffs():
    for r in rules.values():
        yield from r.values()
    for r in alllin.values():
        yield from r.values()
nonint = sum(1 for c in all_coeffs() if not KX.den_ok(c))
LMAX = max(len(q) for q in rules)
def irreducible(w): return all(w[i:i + L] not in rules for L in range(2, LMAX + 1) for i in range(max(0, len(w) - L), len(w) - L + 1))
words = [[(a,) for a in active]]
while words[-1] and len(words) < 12:
    words.append([w + (a,) for w in words[-1] for a in active if src(w[-1]) == tgt(a) and irreducible(w + (a,))])
counts = [len(w) for w in words]
log('irreducible words by length (exact rules):', counts, 'non-integral coefficients:', nonint, STATS)
# every word of length 4 among active arrows contains a leading word (finiteness over the local ring)
len4_irred = sum(1 for w in words[3]) if len(words) > 3 else 0
pickle.dump({'rules': rules, 'alllin': alllin, 'killed': sorted(killed), 'active': active, 'verts': verts}, open(HERE + '/upper_exact.pkl', 'wb'))
RED = {p: {w: KX.kred(c) for w, c in r.items() if KX.kred(c)} for p, r in rules.items()}
pickle.dump({'rules': RED, 'active': active, 'verts': verts}, open(HERE + '/upper_exact_reduced_modp.pkl', 'wb'))
json.dump({'saturation': hist, 'completion': comp, 'killed': sorted(killed), 'vertices': verts, 'active_arrows': len(active),
           'irreducible_counts': counts, 'length4_irreducible': len4_irred, 'inclusion_ambiguities': incl, 'inclusion_unresolved': incl_bad,
           'non_integral_coefficients': nonint, 'leading_coefficients': STATS, 'max_rule_length': LMAX, 'prime': P}, open(HERE + '/z4x.json', 'w'), indent=1)
log('done')
