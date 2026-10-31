"""Axis A, weight three, EXACT over K (choice 0 of the Griess subalgebra 9).
 (a) exact basis of S_3 = span{ x(-1)U_2, L(-1)U_2, x(-2)y, x(-3)1, u(0)U_2 (u = e, x, eight highest weight vectors) }
     closed under the zero modes of g' + g'' (pivots chosen modulo p, hence independent over K);
 (T1) every generating word lies in the span (exact); (T2) e_i(0), f_i(0) preserve the span (exact);
 (b) L(1) S_3 c U_2, L(2) S_3 c U_1, y(1) S_3 c U_2, y(2) S_3 c U_1 (y in U_1) (exact)."""
import time, pickle, json, sys
from gmpy2 import mpq
import lattice_voa as V
import exact_engine as E, Kfield as K
t0 = time.time()
def log(*a): print(*a, '(%.0fs)' % (time.time() - t0), flush=True)
S = E.setup(0); T = S['T']
B1, B2, B3 = E.Blocks(T, 1), E.Blocks(T, 2), E.Blocks(T, 3)
def byblock(B, vecs):
    out = {}
    for i, v in enumerate(vecs):
        if v: out.setdefault(B.wt(v), []).append((i, v))
    return out
def pivots_by_block(B, vecs, start=None):
    start = start or {}
    res = {mu: list(start.get(mu, [])) for mu in B.blk}
    for mu, lst in byblock(B, vecs).items():
        idx = E.select_pivots([v for _, v in lst], B.pos[mu], len(B.blk[mu]), start=res[mu])
        res[mu] += [lst[i][1] for i in idx]
    return {mu: v for mu, v in res.items() if v}
def check_in(B, basis, tests, label):
    bad = 0; ntest = 0
    for mu, lst in byblock(B, tests).items():
        ok, bads = E.solve_block(basis.get(mu, []), [v for _, v in lst], B.pos[mu], len(B.blk[mu]))
        assert ok, ('basis not independent', mu)
        bad += len(bads); ntest += len(lst)
    log('  %-40s tests %6d  failures %d' % (label, ntest, bad))
    return ntest, bad
res = {}
U2b = pivots_by_block(B2, S['U2span']); n2 = sum(len(v) for v in U2b.values())
U2list = [v for vs in U2b.values() for v in vs]
U1b = pivots_by_block(B1, S['U1k'])
log('U_2 exact basis:', n2, ' U_1:', sum(len(v) for v in U1b.values()))
res['U2_rank'] = n2
# weight-three words
ops = {}
words = []
for i, x in enumerate(S['U1']):
    O = E.build_op([(K.ONE, x)], -1, 2)
    words += [('x(-1)w', i, j, E.apply_op(O, w)) for j, w in enumerate(U2list)]
O = E.build_op([(K.ONE, V.OMEGA)], 0, 2); words += [('L(-1)w', 0, j, E.apply_op(O, w)) for j, w in enumerate(U2list)]
for i, x in enumerate(S['U1']):
    for j, y in enumerate(S['U1']):
        words.append(('x(-2)y', i, j, E.kprod([(K.ONE, x)], -2, S['U1k'][j])))
    words.append(('x(-3)1', i, 0, E.kprod([(K.ONE, x)], -3, {V.VAC: K.ONE})))
for name, terms in S['hw'].items():
    O = E.build_op(terms, 0, 2)
    words += [('u(0)w:' + name, 0, j, E.apply_op(O, w)) for j, w in enumerate(U2list)]
log('weight-three words:', len(words))
S3b = pivots_by_block(B3, [w[3] for w in words])
log('pivots before g-closure:', sum(len(v) for v in S3b.values()))
gops = [E.build_op([(K.ONE, T[k])], 0, 3) for k in ('e1', 'f1', 'e2', 'f2')]
for it in range(10):
    imgs = [E.apply_op(O, v) for O in gops for vs in S3b.values() for v in vs]
    before = sum(len(v) for v in S3b.values())
    S3b = pivots_by_block(B3, imgs, start=S3b)
    after = sum(len(v) for v in S3b.values())
    log('  g-closure iteration', it, ':', after)
    if after == before: break
res['S3_rank_lower_bound'] = after
pickle.dump({'S3b': S3b, 'U2b': U2b, 'U1b': U1b}, open('exact_S3_basis.pkl', 'wb'))
log('S_3 exact basis: %d independent vectors over K (minor nonzero mod p)' % after)
# T1, T2
res['T1'] = check_in(B3, S3b, [w[3] for w in words], 'T1: generating words in S_3')
res['T2'] = check_in(B3, S3b, [E.apply_op(O, v) for O in gops for vs in S3b.values() for v in vs], 'T2: g-stability of S_3')
# (b)
S3list = [v for vs in S3b.values() for v in vs]
for lab, terms, n, tgt, Bt, basis in [('L(1)S3 c U2', [(K.ONE, V.OMEGA)], 2, 2, B2, U2b), ('L(2)S3 c U1', [(K.ONE, V.OMEGA)], 3, 1, B1, U1b)] + \
        [('y%d(1)S3 c U2' % i, [(K.ONE, y)], 1, 2, B2, U2b) for i, y in enumerate(S['U1'])] + \
        [('y%d(2)S3 c U1' % i, [(K.ONE, y)], 2, 1, B1, U1b) for i, y in enumerate(S['U1'])]:
    O = E.build_op(terms, n, 3)
    res[lab] = check_in(Bt, basis, [E.apply_op(O, s) for s in S3list], lab)
json.dump(res, open('certify_w3.json', 'w'), indent=1)
log('done')
