"""Control of the predicted dimension of B: dimension of the g'+g''-invariants of V'_2 and of the A-primary
vectors of type (0,0) (= B), and of the weight-one A-primary vectors (R_b) by type."""
import sys
from gmpy2 import mpq
import lattice_voa as V, cases as C, commutant_griess as CG
for case in ('I', 'II'):
    T = C.build(case); T['case'] = case
    inv, _ = CG.annihilated(T, 2, (0, 0), ns=(0,))
    B, _ = CG.annihilated(T, 2, (0, 0), ns=(0, 1, 2))
    # weight-one A-primaries = highest weight vectors of V'_1 annihilated by e1(0), e2(0), x(1)
    types = {}
    for mu in sorted(set(C.weight_of_key(T, k) for k in V.basis(1))):
        if mu[0] < 0 or mu[1] < 0: continue
        W = [k for k in V.basis(1) if C.weight_of_key(T, k) == mu]
        # highest weight: e1_0 v = e2_0 v = 0 ; primary: x_1 v = 0
        import linalg as LA
        cols = []
        blocks = [(T['e1'], 0), (T['e2'], 0)] + [(x, 1) for x in CG.U1_basis(T)]
        keys = set(); imgs = []
        for x, n in blocks:
            imgs.append([V.clean(V.mode_vec(x, n, {k: mpq(1)})) for k in W])
        from flint import fmpq_mat
        allk = sorted(set((b, k) for b, blk in enumerate(imgs) for c in blk for k in c))
        idx = {kk: i for i, kk in enumerate(allk)}
        M = fmpq_mat(max(1, len(allk)), len(W))
        for b, blk in enumerate(imgs):
            for j, c in enumerate(blk):
                for k, x in c.items(): M[idx[(b, k)], j] = LA.to_fq(x)
        d = len(LA.nullspace(M))
        if d: types[mu] = d
    print('case', case, ': invariants of V\'_2:', len(inv), '; dim B = Com(A)_2:', len(B),
          '; weight-one highest weight vectors (incl. U_1) by (H\',H\'\') weight:', types)
