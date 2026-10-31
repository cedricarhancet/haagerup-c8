"""Ranks modulo p of N_3 = L(-1)U_2 + U_1(-1)U_2 (+ x(-2)y, x(-3)1) and of the affine part A_2 of U_2, for the count
of strong generators (weight two: 120 - dim A_2 ; weight three: 956 - dim N_3)."""
import pickle, json
import lattice_voa as V, exact_engine as E, Kfield as K
S = E.setup(0); T = S['T']
d = pickle.load(open('exact_S3_basis.pkl', 'rb'))
U2list = [v for vs in d['U2b'].values() for v in vs]
B2, B3 = E.Blocks(T, 2), E.Blocks(T, 3)
def rank(B, vecs):
    ech = {mu: E.Ech(len(ks)) for mu, ks in B.blk.items()}; r = 0
    for v in vecs:
        if v: mu = B.wt(v); r += ech[mu].add(E.modp_vec(v, B.pos[mu], len(B.blk[mu])))
    return r
A2 = [E.kv_from_terms([(K.ONE, a)]) for a in S['A2'] if a]
words = []
for x in S['U1']:
    O = E.build_op([(K.ONE, x)], -1, 2); words += [E.apply_op(O, w) for w in U2list]
O = E.build_op([(K.ONE, V.OMEGA)], 0, 2); words += [E.apply_op(O, w) for w in U2list]
res = {'dim_A2_modp': rank(B2, A2), 'dim_N3_modp': rank(B3, words)}
print(res); json.dump(res, open('count_generators.json', 'w'))
