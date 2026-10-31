"""Exact products for the Zhu relations (weights <= 4), stored exactly on disk and modulo p in memory-light form.
 (1) a(n)b, n = 3, 2, 1, 0, -1, for a <= b in the lowering basis of U_2;
 (2) h(-1)z, T z (z in the lowering basis), weight 3;
 (3) h(-1)s, T s for the exact basis s of S_3 (956 vectors), weight 4;
 (4) T = L(-1) on V'_1, V'_2, V'_3 as sparse matrices modulo p.
Resumable: per-a files in prod/."""
import pickle, time, gc, os
import scipy.sparse as sp
from zhu_common import *
t0 = time.time()
def log(*a): print(*a, '(%.0fs)' % (time.time() - t0), flush=True)
OUT = ZD + '/prod'; os.makedirs(OUT, exist_ok=True)
S = E.setup(0); T = S['T']
d = pickle.load(open(ZD + '/U2_lowering_basis.pkl', 'rb')); Z = [b[3] for b in d['basis']]
POS = {w: {k: i for i, k in enumerate(V.basis(w))} for w in range(5)}
def mp(v, w): return modp_sparse(v, POS[w]) if v else (np.zeros(0, np.int64), np.zeros(0, np.int64))
# (4) translation matrices
fn = OUT + '/Tmat.pkl'
if not os.path.exists(fn):
    Tm = {}
    for w in (0, 1, 2, 3):
        O = E.build_op([(K.ONE, V.OMEGA)], 0, w); rows, cols, vals = [], [], []
        for k, img in O.items():
            for kk, a in img.items(): rows.append(POS[w + 1][kk]); cols.append(POS[w][k]); vals.append(E.kred(a))
        Tm[w] = sp.csr_matrix((np.array(vals, np.int64), (rows, cols)), shape=(len(POS[w + 1]), len(POS[w])))
        del O; gc.collect()
    pickle.dump(Tm, open(fn, 'wb')); log('T matrices done')
# (1) products in U_2
for a in range(120):
    fn = OUT + '/a%03d.pkl' % a
    if os.path.exists(fn): continue
    res = {}; exact = {}
    for n, wt in ((3, 0), (2, 1), (1, 2), (0, 3), (-1, 4)):
        O = E.build_op(E.kv_terms(Z[a]), n, 2)
        for b in range(a, 120):
            v = E.apply_op(O, Z[b]); res[b, n] = mp(v, wt)
            if n == -1: exact[b] = v
        del O; V._MODE.clear(); V._EMODE.clear(); gc.collect()
    pickle.dump(res, open(fn, 'wb')); pickle.dump(exact, open(OUT + '/exact_a%03d.pkl' % a, 'wb'))
    del res, exact; gc.collect()
    if a % 10 == 0: log('U2 products, a =', a)
# (2) currents on U_2
fn = OUT + '/hU2.pkl'
if not os.path.exists(fn):
    res = {}
    for h in range(6):
        O = E.build_op([(K.ONE, S['U1'][h])], -1, 2)
        for b in range(120): res['h', h, b] = mp(E.apply_op(O, Z[b]), 3)
    O = E.build_op([(K.ONE, V.OMEGA)], 0, 2)
    for b in range(120): res['T', b] = mp(E.apply_op(O, Z[b]), 3)
    pickle.dump(res, open(fn, 'wb')); log('currents on U_2 done')
# (3) S_3
fn = OUT + '/S3.pkl'
if not os.path.exists(fn):
    dd = pickle.load(open(os.path.join(GF, 'exact_S3_basis.pkl'), 'rb'))
    S3 = [v for vs in dd['S3b'].values() for v in vs]
    res = {'S3': [mp(s, 3) for s in S3]}
    for h in range(6):
        O = E.build_op([(K.ONE, S['U1'][h])], -1, 3)
        ex = []
        for i, s in enumerate(S3):
            v = E.apply_op(O, s); res['h', h, i] = mp(v, 4); ex.append(v)
        pickle.dump(ex, open(OUT + '/exact_hS3_%d.pkl' % h, 'wb')); del O, ex; gc.collect(); log('h(-1)S_3, h =', h)
    O = E.build_op([(K.ONE, V.OMEGA)], 0, 3)
    for i, s in enumerate(S3): res['T', i] = mp(E.apply_op(O, s), 4)
    pickle.dump(res, open(fn, 'wb')); log('S_3 done')
log('done')
