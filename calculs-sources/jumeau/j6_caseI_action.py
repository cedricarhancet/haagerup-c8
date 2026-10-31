"""Axis B (case I): the 8 elements of the centralizer Z of U_1 (levels 1, 39) and their action on B_I = Com_{V'}(A)_2
(dim 4, basis of griess_subalgebras.griess_data('I')); which elements exchange the Griess subalgebras 5 and 6?
Groebner solution (j2_solve.py I): z_(0,4) = diag(e0, e1), z_(0,8) = e1, z_(1,5) = e2, z_(1,9) = z_(1,3) = e1 e2,
all other blocks identity."""
import sys, os, json, pickle, itertools
HERE = os.path.abspath(os.path.dirname(__file__)); sys.path.insert(0, os.path.abspath(HERE + '/../generation-forte'))
os.chdir(HERE + '/../generation-forte')
from gmpy2 import mpq
import sympy as sp
import lattice_voa as V, griess_subalgebras as GS
d = pickle.load(open(HERE + '/j1_equations_I.pkl', 'rb'))
label = [((tuple(l[0]), l[1], l[2], l[3]) if not isinstance(l[0], tuple) else tuple(l)) for l in d['label']]
P = sp.Matrix.hstack(*[sp.Matrix(c) for c in d['P']]); n = P.shape[0]; Pinv = P.inv()
types = pickle.load(open(HERE + '/j1_types_I.pkl', 'rb'))['types']
mults = {tuple(int(t) for t in k.strip('()').split(',')): len(v) for k, v in types.items()}
def zmat(e0, e1, e2):
    Zb = {mu: sp.eye(m) for mu, m in mults.items()}
    Zb[(0, 4)] = sp.diag(e0, e1); Zb[(0, 8)] = sp.Matrix([[e1]]); Zb[(1, 5)] = sp.Matrix([[e2]])
    Zb[(1, 9)] = sp.Matrix([[e1*e2]]); Zb[(1, 3)] = sp.Matrix([[e1*e2]])
    pos = {l: i for i, l in enumerate(label)}
    Zf = sp.zeros(n, n)
    for (mu, k, i, j), p in pos.items():
        for l in range(Zb[mu].shape[0]):
            Zf[pos[(mu, l, i, j)], p] = Zb[mu][l, k]
    return P*Zf*Pinv
T, basis, prod, gram = GS.griess_data('I')
B1 = V.basis(1); idx = {k: i for i, k in enumerate(B1)}
def tofrac(x): return sp.Rational(int(x.numerator), int(x.denominator))
def phi(b):
    M = sp.zeros(n, n)
    for j, k in enumerate(B1):
        for kk, x in V.clean(V.mode_vec({k: mpq(1)}, 1, b)).items(): M[idx[kk], j] = tofrac(x)
    return M
Phi = [phi(b) for b in basis]; m = len(basis)
A = sp.Matrix.hstack(*[p.reshape(n*n, 1) for p in Phi])
G = json.load(open('griess_subalgebras.json'))['I']['subalgebras']
X = [sp.Matrix([0] + [sp.sympify(c) for c in s['x']]) for s in G]
def same_line(u, v):
    M = sp.Matrix.hstack(u, v)
    return all(sp.simplify(M[i, 0]*M[j, 1] - M[j, 0]*M[i, 1]) == 0 for i in range(m) for j in range(m))
res = {}
for e0, e1, e2 in itertools.product((1, -1), repeat=3):
    Z = zmat(e0, e1, e2); Zi = Z.inv()
    # monomial / diagonal ?
    mono = all(sum(1 for i in range(n) if Z[i, j] != 0) == 1 for j in range(n))
    C = sp.zeros(m, m)
    for i in range(m):
        sol, params = A.gauss_jordan_solve((Z*Phi[i]*Zi).reshape(n*n, 1))
        for j in range(m): C[j, i] = sol[j]
    perm = {i: [j for j in range(len(X)) if same_line(C*X[i], X[j])] for i in range(len(X))}
    res[str((e0, e1, e2))] = {'monomial': mono, 'action_on_B': [[str(x) for x in C.row(r)] for r in range(m)], 'perm': perm}
    print((e0, e1, e2), 'monomial', mono, 'B-action', C.tolist(), 'perm', perm, flush=True)
    pickle.dump([[str(x) for x in Z.row(i)] for i in range(n)], open(HERE + '/j6_zI_%d%d%d.pkl' % ((e0 + 1)//2, (e1 + 1)//2, (e2 + 1)//2), 'wb'))
json.dump(res, open(HERE + '/j6_caseI_action.json', 'w'), indent=1)
