"""Axis D1 (modulo p3, coherent embedding): dimensions of C_n = Com_S(A)_n for n = 3, 4, as the common kernel of the
modes y_(k), k >= 0, y in U_1, on the g-invariant weight block (0,0) of S_n.  S_3: exact basis (certify_w3) reduced
mod p; S_4: checkpoint of closure_w4_modp (same prime and embedding)."""
import pickle, json, numpy as np
from gmpy2 import mpq
from flint import nmod_mat
import lattice_voa as V, exact_engine as E, cases as C
S = E.setup(0); T = S['T']; P = E.P; F = E.F
res = {}
def kernel_dim(basis_rows, keys_src, n):
    """dimension of the common kernel of y_(k) (y in U_1, 0 <= k <= n) on span(basis_rows), modulo p"""
    import scipy.sparse as sps
    outkeys = {}
    R_, C_, D_ = [], [], []
    for j, kk in enumerate(keys_src):
        for a, y in enumerate(S['U1']):
            for k in range(0, n + 1):
                r = V.mode_vec(y, k, {kk: mpq(1)})
                for k2, c in r.items():
                    if c:
                        key = (a, k, k2)
                        if key not in outkeys: outkeys[key] = len(outkeys)
                        R_.append(outkeys[key]); C_.append(j); D_.append(F.q(c))
    V._MODE.clear(); V._EMODE.clear()
    M = sps.csr_matrix((np.array(D_, dtype=np.int64), (np.array(R_), np.array(C_))), shape=(len(outkeys), len(keys_src)))
    Bm = np.array(basis_rows, dtype=np.int64).T
    Z = np.asarray(M @ Bm) % P
    Z = Z[np.any(Z, axis=1)]
    if Z.shape[0] == 0: return Bm.shape[1]
    return Bm.shape[1] - nmod_mat(Z.shape[0], Z.shape[1], [int(x) for x in Z.flatten()], P).rank()
d = pickle.load(open('exact_S3_basis.pkl', 'rb'))
B3 = E.Blocks(T, 3)
mu0 = (0, 0)
S3_00 = d['S3b'].get(mu0, [])
rows3 = [E.modp_vec(v, B3.pos[mu0], len(B3.blk[mu0])) for v in S3_00]
res['dim_S3_00'] = len(rows3); res['dim_C3'] = kernel_dim(rows3, B3.blk[mu0], 3)
print(res, flush=True)
z = np.load('../seconde-realisation/ckpt_S4_0_p8388571.npz')
B4 = E.Blocks(T, 4)
# the checkpoint blocks are indexed like closure_w5_modp (V.basis(4) order, blocks by weight) -> same order as B4.blk[mu]
rows4 = list(z['0_0'])
res['dim_S4_00'] = len(rows4); res['dim_C4'] = kernel_dim(rows4, B4.blk[mu0], 4)
# expected Vir(c1) x Vir(c2) vacuum descendants (c1, c2 generic): weight 3 -> 2, weight 4 -> 5
res['new_primaries_w3'] = res['dim_C3'] - 2; res['new_quasi_w4_upper'] = res['dim_C4'] - 5
print(res); json.dump(res, open('commutant_dims.json', 'w'), indent=1)
