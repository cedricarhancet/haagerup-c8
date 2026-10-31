"""Parallel helper for certify_w4b.py: exact test vectors u(0)s for a subset of highest weight vectors u,
stored per block in the directory given (same format as tests_w4/)."""
import sys, os, pickle, time, gc
import exact_engine as E
t0 = time.time()
names = sys.argv[1].split(';'); outdir = sys.argv[2]
S = E.setup(0); T = S['T']
d = pickle.load(open('exact_S3_basis.pkl', 'rb'))
S3list = [v for vs in d['S3b'].values() for v in vs]
B4 = E.Blocks(T, 4); bidx = {mu: i for i, mu in enumerate(sorted(B4.blk))}
os.makedirs(outdir, exist_ok=True); fh = {}
for name in names:
    terms = S['hw'][name]
    single = [(c, {k: x}) for c, v in terms for k, x in v.items()]
    acc = [dict() for _ in S3list]
    for c0 in range(0, len(single), 12):
        O = E.build_op(single[c0:c0 + 12], 0, 3)
        for j, s in enumerate(S3list): E.kv_add(acc[j], E.apply_op(O, s))
        del O; gc.collect()
    for j in range(len(S3list)):
        v = E.kv_clean(acc[j])
        if not v: continue
        i = bidx[B4.wt(v)]
        if i not in fh: fh[i] = open(outdir + '/blk_%03d.pkl' % i, 'ab')
        pickle.dump((('u(0)s', name, j), v), fh[i])
    for f in fh.values(): f.flush()
    del acc; gc.collect(); print('done', name, '%.0fs' % (time.time() - t0), flush=True)
for f in fh.values(): f.close()
