"""Exact solve phase of the weight-four certification (after certify_w4b.py / certify_w4_part.py):
 tests = u(0)s (from tests_w4/ and tests_w4_part/) and e(-1)w, x(-1)w (computed here); exact membership in the span of
 the 6010 pivot words of N'_4 (piv_w4/); phase (d): every word of the five families of N'_4 lies in that span."""
import os, sys, pickle, json, time, gc
import lattice_voa as V, exact_engine as E, Kfield as K
t0 = time.time()
def log(*a): print(*a, '(%.0fs)' % (time.time() - t0), flush=True)
S = E.setup(0); T = S['T']
d = pickle.load(open('exact_S3_basis.pkl', 'rb'))
S3list = [v for vs in d['S3b'].values() for v in vs]; U2list = [v for vs in d['U2b'].values() for v in vs]
B4 = E.Blocks(T, 4); bidx = {mu: i for i, mu in enumerate(sorted(B4.blk))}
def load(path):
    out = []
    if not os.path.exists(path): return out
    with open(path, 'rb') as f:
        while True:
            try: out.append(pickle.load(f))
            except (EOFError, pickle.UnpicklingError): break
    return out
counts = json.load(open('piv_w4/DONE.json')); npiv = sum(counts.values())
# phase (d)
buf = {}; nd = [0, 0]
def flush():
    for mu, X in buf.items():
        ok, bad = E.solve_block(load('piv_w4/blk_%03d.pkl' % bidx[mu]), X, B4.pos[mu], len(B4.blk[mu]))
        nd[0] += len(X); nd[1] += len(bad)
    buf.clear(); gc.collect()
def collect(v):
    if v: buf.setdefault(B4.wt(v), []).append(v)
    if sum(len(x) for x in buf.values()) > 1500: flush()
O = E.build_op([(K.ONE, V.OMEGA)], 0, 3)
for s in S3list: collect(E.apply_op(O, s))
del O
for y in S['U1']:
    O = E.build_op([(K.ONE, y)], -1, 3)
    for s in S3list: collect(E.apply_op(O, s))
    O = E.build_op([(K.ONE, y)], -2, 2)
    for w in U2list: collect(E.apply_op(O, w))
    for y2 in S['U1k']: collect(E.kprod([(K.ONE, y)], -3, y2))
    del O; gc.collect()
log('  phase (d): affine and L(-1) families queued')
for name, idx, terms in S['mods']:
    O = E.build_op(terms, -1, 2)
    for w in U2list: collect(E.apply_op(O, w))
    del O
flush()
log('EXACT phase (d): words of the five families in the span of the %d pivots:' % npiv, nd[0], 'failures', nd[1])
json.dump({'phase_d': nd, 'N4_pivots': npiv}, open('phase_d_w4.json', 'w'), indent=1)
log('done')
