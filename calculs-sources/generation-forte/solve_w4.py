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
extra = {}
for name in ('e', 'x'):
    O = E.build_op(S['hw'][name], -1, 2)
    for j, w in enumerate(U2list):
        v = E.apply_op(O, w)
        if v: extra.setdefault(B4.wt(v), []).append((('u(-1)w', name, j), v))
seen = set(); ntests = 0; fails = []
for mu in sorted(B4.blk, key=lambda m: len(B4.blk[m])):
    i = bidx[mu]
    tests = load('tests_w4/blk_%03d.pkl' % i) + load('tests_w4_part/blk_%03d.pkl' % i) + extra.get(mu, [])
    tests = [(lab, v) for lab, v in tests if (lab not in seen and not seen.add(lab))]
    if not tests: continue
    piv = load('piv_w4/blk_%03d.pkl' % i)
    ok, bad = E.solve_block(piv, [v for _, v in tests], B4.pos[mu], len(B4.blk[mu]))
    assert ok, mu
    ntests += len(tests); fails += [str(tests[k][0]) for k in bad]
names = sorted({lab[1] for lab in seen if lab[0] == 'u(0)s'})
log('EXACT: u(0)S_3 for u in', names, 'and e(-1)U_2, x(-1)U_2 in N\'_4 :', ntests, 'tests, failures', len(fails))
res = {'N4_pivots': npiv, 'N4_by_family': counts, 'tests': ntests, 'n_failures': len(fails), 'failures': fails[:50], 'hw_tested': names}
json.dump(res, open('certify_w4.json', 'w'), indent=1)
log('done')
