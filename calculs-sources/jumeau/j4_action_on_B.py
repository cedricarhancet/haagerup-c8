"""Axis B1: action of the 8 elements of Z on B = Com_{V'}(A)_2 (case II, basis of griess_B_II_basis.pkl).
A weight-two vector b of V' = L_e6(1) (x) L_sl3(1) is determined by the map x -> x_(1) b on V'_1 (no singular vector
in weight two), and z(b) is characterized by x_(1) z(b) = z((z^-1 x)_(1) b).  Since z commutes with U_1, z(B) = B, so
z(b_i) = sum_j c_ij b_j is found by a small linear system."""
import sys, os, pickle, itertools, json
HERE = os.path.abspath(os.path.dirname(__file__)); sys.path.insert(0, os.path.abspath(HERE + '/../generation-forte'))
os.chdir(HERE + '/../generation-forte')
from gmpy2 import mpq
import sympy as sp
import lattice_voa as V
Bb = pickle.load(open('griess_B_II_basis.pkl', 'rb'))
B1 = V.basis(1); n = len(B1); idx = {k: i for i, k in enumerate(B1)}
def tofrac(x): return sp.Rational(int(x.numerator), int(x.denominator))
def phi(b):
    M = sp.zeros(n, n)
    for j, k in enumerate(B1):
        r = V.clean(V.mode_vec({k: mpq(1)}, 1, b))
        for kk, x in r.items(): M[idx[kk], j] = tofrac(x)
    return M
Phi = [phi(b) for b in Bb]
print('phi computed; ranks', [p.rank() for p in Phi])
res = {}
for bits in itertools.product((1, 0), repeat=3):
    Z = sp.Matrix([[sp.sympify(x) for x in r] for r in pickle.load(open(HERE + '/j3_z_%d%d%d.pkl' % bits, 'rb'))])
    Zi = Z.inv()
    C = sp.zeros(5, 5)
    for i in range(5):
        target = Z*Phi[i]*Zi
        # solve sum_j c_j Phi_j = target
        A = sp.Matrix.hstack(*[p.reshape(n*n, 1) for p in Phi]); bvec = target.reshape(n*n, 1)
        sol, params = A.gauss_jordan_solve(bvec)
        assert not params.free_symbols if hasattr(params, 'free_symbols') else True
        for j in range(5): C[j, i] = sol[j]
        assert A*sol == bvec
    eps = tuple(2*b - 1 for b in bits)
    res[str(eps)] = [[str(x) for x in C.row(r)] for r in range(5)]
    print(eps, 'action on B (columns = images):', C.tolist())
json.dump(res, open(HERE + '/j4_action_on_B.json', 'w'), indent=1)
