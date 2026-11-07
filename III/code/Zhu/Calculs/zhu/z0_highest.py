"""Exact highest weight vectors of U_2 and the lowering basis; checks: highest (exact), in U_2 and rank 120 (mod p)."""
import pickle, json, time
from zhu_common import *
t0 = time.time()
S = E.setup(0); T = S['T']
hv, cur = highest_vectors(S)
res = {'highest': []}
for lab, w, z in hv:
    bad = [n for n in ('e1', 'e2') if prod(cur[n], 0, z)]
    res['highest'].append({'label': lab, 'type': list(w), 'not_killed_by': bad, 'nterms': len(z)})
basis = lowering_basis(hv, cur)
B2 = E.Blocks(T, 2)
d = pickle.load(open(os.path.join(GF, 'exact_S3_basis.pkl'), 'rb'))
U2list = [v for vs in d['U2b'].values() for v in vs]
pos = {k: i for i, k in enumerate(B2.keys)}
def mat(vs):
    M = np.zeros((len(vs), len(B2.keys)), dtype=np.int64)
    for r, v in enumerate(vs):
        for k, a in v.items(): M[r, pos[k]] = E.kred(a)
    return M
from flint import nmod_mat
def rk(M): return nmod_mat(M.shape[0], M.shape[1], [int(x) for x in M.flatten()], P).rank()
r1 = rk(mat([b[3] for b in basis])); r2 = rk(mat(U2list + [b[3] for b in basis]))
res['lowering_basis_size'] = len(basis); res['rank_lowering_mod_p'] = r1; res['rank_with_U2_mod_p'] = r2
print(json.dumps(res, indent=1))
pickle.dump({'hv': [(l, w) for l, w, _ in hv], 'basis': basis}, open(ZD + '/U2_lowering_basis.pkl', 'wb'))
json.dump(res, open(ZD + '/z0_highest.json', 'w'), indent=1)
print('(%.0fs)' % (time.time() - t0))
