"""Axis B1-B2: A-primary highest weight vectors of weight two (spaces Pi_tau) and the linear conditions
imposed by C_2 = span(e, x) on the line p_tau = U_2 cap Pi_tau:
   (K) f_2 u = 0 for f in C_2   (u in U_2 primary => L_1 u, x_2 u are weight-one A-primaries of U, hence 0),
   (E) u is an eigenvector of x_1 on Pi_tau.
Coefficients of x may lie in Q(sqrt 13); computations with sympy over QQ<sqrt(13)>."""
import sys, json, time, pickle
import sympy as sp
from gmpy2 import mpq
from flint import fmpq_mat
import lattice_voa as V, cases as C, commutant_griess as CG, griess_subalgebras as GS, linalg as LA
PRIMARY_TYPES = {'I': {(0, 4): 1, (0, 8): 2, (0, 12): 1, (1, 1): 1, (1, 5): 1, (1, 7): 1, (1, 13): 1},
                 'II': {(0, 4): 1, (1, 1): 1, (1, 3): 1, (1, 5): 1, (1, 7): 1, (2, 2): 1, (2, 6): 1, (3, 3): 1}}
def primaries(T, wt, tau):
    """basis of A-primary highest weight vectors of weight wt and joint weight tau"""
    W = [k for k in V.basis(wt) if C.weight_of_key(T, k) == tau]
    blocks = [(T['e1'], 0), (T['e2'], 0)] + [(x, n) for x in CG.U1_basis(T) for n in range(1, wt + 1)]
    imgs = [[V.clean(V.mode_vec(x, n, {k: mpq(1)})) for k in W] for x, n in blocks]
    allk = sorted(set((b, k) for b, blk in enumerate(imgs) for c in blk for k in c))
    idx = {kk: i for i, kk in enumerate(allk)}
    M = fmpq_mat(max(1, len(allk)), len(W))
    for b, blk in enumerate(imgs):
        for j, c in enumerate(blk):
            for k, x in c.items(): M[idx[(b, k)], j] = LA.to_fq(x)
    ns = LA.nullspace(M)
    return [V.clean({W[j]: LA.to_mq(v[j]) for j in range(len(W)) if v[j] != 0}) for v in ns]
def to_sp(x): return sp.Rational(int(x.numerator), int(x.denominator))
def op_matrix(b, n, src, tgt):
    """matrix of b_n from span(src) to span(tgt) (exact rational); tgt must contain the images"""
    cols = []
    for u in src:
        im = V.clean(V.mode_vec(b, n, u))
        if not im:
            cols.append([0]*len(tgt)); continue
        cc = LA.coords(tgt, im)
        assert cc is not None, 'image not in target span'
        cols.append([sp.Rational(int(c.p), int(c.q)) for c in cc])
    return sp.Matrix(len(tgt), len(src), lambda i, j: cols[j][i]) if tgt else sp.zeros(0, len(src))
def setup(case):
    T, basis, prod, gram = GS.griess_data(case)      # basis = [e] + perp (rational)
    data = {'T': T, 'basis': basis, 'types': {}}
    for tau in PRIMARY_TYPES[case]:
        P = primaries(T, 2, tau)
        P1 = primaries(T, 1, tau)
        # operators: b_1 on P (into P), b_2 on P (into P1) for b in basis
        X1 = [op_matrix(b, 1, P, P) for b in basis]
        X2 = [op_matrix(b, 2, P, P1) for b in basis]
        data['types'][tau] = {'P': P, 'P1': P1, 'X1': X1, 'X2': X2}
    return data
if __name__ == '__main__':
    for case in sys.argv[1:] or ['I', 'II']:
        t0 = time.time()
        D = setup(case)
        print('case', case, ': dims of Pi_tau (weight 2) and of the weight-one primaries of type tau:')
        for tau, d in D['types'].items():
            print('   ', tau, 'dim Pi =', len(d['P']), '; weight-one primaries:', len(d['P1']))
        pickle.dump({'types': {str(t): {'X1': [sp.srepr(m) for m in d['X1']], 'X2': [sp.srepr(m) for m in d['X2']],
                                        'P': [{repr(k): str(v) for k, v in u.items()} for u in d['P']]}
                               for t, d in D['types'].items()}}, open(f'primaries_{case}.pkl', 'wb'))
        print('   time %.1fs' % (time.time() - t0))
