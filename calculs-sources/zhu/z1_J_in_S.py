"""Z1 (exact over K): the simple current J = L(3,3) x L(13,13) lies in S.
 The Cartan block (3,13) of V'_4 is one-dimensional (a single basis state); S_4 meets it (pivot of N'_4, an exact
 element of S); the vector is annihilated by e'(0), e''(0) and by all modes y(n), n >= 1, y in U_1 (exact)."""
import sys, os, pickle, json, time
sys.path.insert(0, os.path.abspath('.'))
import lattice_voa as V, exact_engine as E, Kfield as K
t0 = time.time()
S = E.setup(0); T = S['T']
B4 = E.Blocks(T, 4)
bidx = {mu: i for i, mu in enumerate(sorted(B4.blk))}
mu = (3, 13)
keys = B4.blk[mu]
print('Cartan block (3,13) of V\'_4: %d basis state(s):' % len(keys), keys)
fn = 'piv_w4/blk_%03d.pkl' % bidx[mu]
vecs = []
if os.path.exists(fn):
    with open(fn, 'rb') as f:
        while True:
            try: vecs.append(pickle.load(f))
            except EOFError: break
print('pivots of N\'_4 (exact elements of S_4) in this block:', len(vecs))
v = vecs[0]; print('vector:', {k: str(a) for k, a in v.items()})
res = {'block_dim': len(keys), 'n_pivots': len(vecs), 'key': str(keys[0])}
# singularity
tests = []
for name in ('e1', 'e2'):
    tests.append((name + '(0)', E.kprod([(K.ONE, T[name])], 0, v)))
for i, y in enumerate(S['U1']):
    for n in range(1, 5):
        tests.append(('y%d(%d)' % (i, n), E.kprod([(K.ONE, y)], n, v)))
bad = [lab for lab, w in tests if w]
print('singular-vector tests:', len(tests), 'nonzero:', bad)
# conformal weight check via Sugawara: L(0) eigenvalue 4 is the weight
res['singular_tests'] = len(tests); res['nonzero'] = bad
json.dump(res, open('../zhu/z1_J_in_S.json', 'w'), indent=1)
print('done (%.0fs)' % (time.time() - t0))
