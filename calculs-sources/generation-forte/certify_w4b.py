"""Axis A, weight four, EXACT over K: S_4^gen = N'_4.
N'_4 = span{ L(-1)S_3, y(-1)S_3, y(-2)U_2, y(-3)U_1, m(-1)U_2 (m in the eight primary modules) } (g-stable families);
an exact basis is chosen among these words (pivots modulo p => independent over K).  Tests (exact membership):
u(0)s for u = e, x, the eight highest weight vectors, s in the exact basis of S_3; e(-1)w, x(-1)w for w in U_2.
Phase (d): every word of the families lies in the span (exact dimension of N'_4)."""
import time, pickle, json, os, sys, gc
import numpy as np
from gmpy2 import mpq
import lattice_voa as V
import exact_engine as E, Kfield as K
t0 = time.time()
def log(*a): print(*a, '(%.0fs)' % (time.time() - t0), flush=True)
S = E.setup(0); T = S['T']
d = pickle.load(open('exact_S3_basis.pkl', 'rb'))
S3list = [v for vs in d['S3b'].values() for v in vs]; U2list = [v for vs in d['U2b'].values() for v in vs]
B4 = E.Blocks(T, 4)
log('S_3 basis', len(S3list), ' U_2 basis', len(U2list), ' V\'_4 blocks', len(B4.blk))
ech = {mu: E.Ech(len(ks)) for mu, ks in B4.blk.items()}
piv = {mu: 0 for mu in B4.blk}
os.makedirs('piv_w4', exist_ok=True)
bidx = {mu: i for i, mu in enumerate(sorted(B4.blk))}
counts = {}
def offer(v, fam):
    if not v: return
    mu = B4.wt(v)
    if ech[mu].add(E.modp_vec(v, B4.pos[mu], len(B4.blk[mu]))):
        with open('piv_w4/blk_%03d.pkl' % bidx[mu], 'ab') as f: pickle.dump(v, f)
        piv[mu] += 1; counts[fam] = counts.get(fam, 0) + 1
def families(callback):
    O = E.build_op([(K.ONE, V.OMEGA)], 0, 3)
    for s in S3list: callback(E.apply_op(O, s), 'L(-1)S3')
    del O; gc.collect(); log('  family L(-1)S_3 done')
    for y in S['U1']:
        O = E.build_op([(K.ONE, y)], -1, 3)
        for s in S3list: callback(E.apply_op(O, s), 'y(-1)S3')
        del O; gc.collect()
    log('  family y(-1)S_3 done')
    for y in S['U1']:
        O = E.build_op([(K.ONE, y)], -2, 2)
        for w in U2list: callback(E.apply_op(O, w), 'y(-2)U2')
        for y2 in S['U1k']: callback(E.kprod([(K.ONE, y)], -3, y2), 'y(-3)U1')
    log('  families y(-2)U_2, y(-3)U_1 done')
    for name, idx, terms in S['mods']:
        O = E.build_op(terms, -1, 2)
        for w in U2list: callback(E.apply_op(O, w), 'm(-1)U2')
        del O
    gc.collect(); log('  family m(-1)U_2 done')
if os.path.exists('piv_w4/DONE.json'):
    counts = json.load(open('piv_w4/DONE.json')); log('pivots of N\'_4 reloaded from piv_w4 (phase already done)')
    npiv = sum(counts.values())
else:
    families(offer)
    npiv = sum(piv.values())
    json.dump(counts, open('piv_w4/DONE.json', 'w'))
log('N\'_4: %d independent words over K (pivots mod p); by family %s' % (npiv, counts))
# tests: stored per block on disk
os.makedirs('tests_w4', exist_ok=True)
fh = {}
ntests = 0
def store(v, lab):
    global ntests
    if not v: return
    mu = B4.wt(v); i = bidx[mu]
    if i not in fh: fh[i] = open('tests_w4/blk_%03d.pkl' % i, 'ab')
    pickle.dump((lab, v), fh[i]); ntests += 1
for name, terms in S['hw'].items():
    single = [(c, {k: x}) for c, v in terms for k, x in v.items()]
    acc = [dict() for _ in S3list]
    for c0 in range(0, len(single), 12):
        O = E.build_op(single[c0:c0 + 12], 0, 3)
        for j, s in enumerate(S3list): E.kv_add(acc[j], E.apply_op(O, s))
        del O; gc.collect()
    for j in range(len(S3list)): store(E.kv_clean(acc[j]), ('u(0)s', name, j))
    del acc; gc.collect(); log('  tests u(0)S_3, u =', name)
for name in ('e', 'x'):
    O = E.build_op(S['hw'][name], -1, 2)
    for j, w in enumerate(U2list): store(E.apply_op(O, w), ('u(-1)w', name, j))
for f in fh.values(): f.close()
log('tests stored:', ntests)
def load(i, d='tests_w4'):
    out = []
    fn = d + '/blk_%03d.pkl' % i
    if not os.path.exists(fn): return out
    with open(fn, 'rb') as f:
        while True:
            try: out.append(pickle.load(f))
            except EOFError: break
    return out
res = {'N4_pivots': npiv, 'N4_by_family': counts, 'tests': ntests}
fails = []
for mu in sorted(B4.blk, key=lambda m: len(B4.blk[m])):
    tests = load(bidx[mu])
    if not tests: continue
    ok, bad = E.solve_block(load(bidx[mu], 'piv_w4'), [v for _, v in tests], B4.pos[mu], len(B4.blk[mu]))
    assert ok, mu
    fails += [tests[i][0] for i in bad]
log('EXACT: tests u(0)S_3 and e(-1)U_2, x(-1)U_2 in N\'_4 :', ntests, 'failures', len(fails))
res['failures'] = [str(f) for f in fails[:50]]; res['n_failures'] = len(fails)
json.dump(res, open('certify_w4.json', 'w'), indent=1)
if os.environ.get('PHASE_D'):
    buf = {}
    def collect(v, fam):
        if v: buf.setdefault(B4.wt(v), []).append(v)
        if sum(len(x) for x in buf.values()) > 1500: flush()
    nd = [0, 0]
    def flush():
        for mu, X in buf.items():
            ok, bad = E.solve_block(load(bidx[mu], 'piv_w4'), X, B4.pos[mu], len(B4.blk[mu]))
            nd[0] += len(X); nd[1] += len(bad)
        buf.clear()
    families(collect); flush()
    log('EXACT phase (d): family words in span of the %d pivots:' % npiv, nd[0], 'failures', nd[1])
    res['phase_d'] = nd
    json.dump(res, open('certify_w4.json', 'w'), indent=1)
log('done')
