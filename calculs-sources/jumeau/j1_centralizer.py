"""Axis B1: centralizer Z of U_1 = g' + g'' in Aut(V') = Aut(e6 + sl3) (case II), and its action on B = Com_{V'}(A)_2.
z commutes with U_1, hence acts on V'_1 = sum_tau V(tau) (x) M_tau as id (x) Z_tau (multiplicity spaces = highest
weight vectors).  z is a Lie algebra automorphism iff  Z A_k = sum_l Z_{lk} A_l Z  for the ad-matrices A_k of a basis
m_k of the highest weight vectors (equivariance then gives all brackets).  Every automorphism of e6 + sl3 preserving
the normalized forms lifts uniquely to V' = L_e6(1,0) (x) L_sl3(1,0)."""
import sys, os, json, pickle, itertools
HERE = os.path.abspath(os.path.dirname(__file__)); sys.path.insert(0, os.path.abspath(HERE + '/../generation-forte'))
os.chdir(HERE + '/../generation-forte')
from gmpy2 import mpq
import sympy as sp
import lattice_voa as V, cases as C
CASE = sys.argv[1] if len(sys.argv) > 1 else 'II'
T = C.build(CASE)
B1 = V.basis(1); n = len(B1); idx = {k: i for i, k in enumerate(B1)}
def vec(d):
    v = [sp.Integer(0)]*n
    for k, x in d.items(): v[idx[k]] = sp.Rational(int(x.numerator), int(x.denominator))
    return sp.Matrix(v)
def admat(x):
    return sp.Matrix.hstack(*[vec(V.clean(V.mode_vec(x, 0, {k: mpq(1)}))) for k in B1])
ad = {s: admat(T[s]) for s in ('e1', 'f1', 'H1', 'e2', 'f2', 'H2')}
# invariant form on V'_1 : (u, v) = coefficient of vacuum in u_(1) v
def formmat():
    M = sp.zeros(n, n)
    for i, k in enumerate(B1):
        for j, l in enumerate(B1):
            c = V.mode_vec({k: mpq(1)}, 1, {l: mpq(1)}).get(V.VAC, 0)
            if c: M[i, j] = sp.Rational(int(c.numerator), int(c.denominator))
    return M
Fm = formmat()
# joint weights
wt = {k: C.weight_of_key(T, k) for k in B1}
types = {}
for mu in sorted(set(wt.values())):
    if mu[0] < 0 or mu[1] < 0: continue
    cols = [i for i, k in enumerate(B1) if wt[k] == mu]
    A = sp.Matrix.vstack(ad['e1'][:, cols], ad['e2'][:, cols])
    ns = A.nullspace()
    if ns:
        types[mu] = []
        for v in ns:
            full = sp.zeros(n, 1)
            for j, c in enumerate(cols): full[c] = v[j]
            types[mu].append(full)
print('multiplicities:', {k: len(v) for k, v in types.items()}, 'total dim', sum((a+1)*(b+1)*len(v) for (a, b), v in types.items()))
# put U_1 highest weight vectors first in their blocks, and orthogonalize the rest (w.r.t. a pairing with lowest vectors)
e1v, e2v = vec(T['e1']), vec(T['e2'])
def lowest(v, a, b):
    w = v
    for _ in range(a): w = ad['f1']*w
    for _ in range(b): w = ad['f2']*w
    return w
def hwpair(u, v, a, b):        # invariant pairing of highest weight vectors through lowest vectors
    return (u.T*Fm*lowest(v, a, b))[0]
for mu, vs in types.items():
    a, b = mu
    fixed = [x for x in (e1v, e2v) if any(sp.Matrix.hstack(*vs, x).rank() == len(vs) for _ in [0])]
    fixed = [x for x in fixed if sp.Matrix.hstack(*vs, x).rank() == len(vs)]
    rest = []
    # complement orthogonal to the fixed vectors
    if fixed:
        x = fixed[0]
        cand = [v - (hwpair(v, x, a, b)/hwpair(x, x, a, b))*x for v in vs]
        R = sp.Matrix.hstack(*cand); basis_rest = R.columnspace()
        types[mu] = [x] + basis_rest
        types[mu] = types[mu][:len(vs)]
    print(mu, 'fixed U_1 vector' if fixed else '', 'Gram:', sp.Matrix([[hwpair(u, v, a, b) for v in types[mu]] for u in types[mu]]))
pickle.dump({'types': {str(k): [list(v) for v in vs] for k, vs in types.items()}}, open(HERE + '/j1_types_%s.pkl' % CASE, 'wb'))
# full basis: f1^i f2^j m
basis = []; label = []
for mu, vs in types.items():
    a, b = mu
    for k, m in enumerate(vs):
        wi = m
        for i in range(a + 1):
            wj = wi
            for j in range(b + 1):
                basis.append(wj); label.append((mu, k, i, j)); wj = ad['f2']*wj
            wi = ad['f1']*wi
P = sp.Matrix.hstack(*basis); assert P.rank() == n, P.rank()
Pinv = P.inv()
# unknown block matrices
syms = []; Zb = {}
for mu, vs in types.items():
    m = len(vs)
    Mz = sp.zeros(m, m)
    for r in range(m):
        for c in range(m):
            s = sp.Symbol('z_%d%d_%d%d' % (mu[0], mu[1], r, c)); syms.append(s); Mz[r, c] = s
    if mu in ((2, 0), (0, 2)):     # the U_1 vector (first) is fixed
        for r in range(m): Mz[r, 0] = 1 if r == 0 else 0
    Zb[mu] = Mz
# full z in the basis P (block diagonal by (mu, i, j))
def zfull():
    Zf = sp.zeros(n, n)
    pos = {l: i for i, l in enumerate(label)}
    for (mu, k, i, j), p in pos.items():
        for l in range(len(types[mu])):
            Zf[pos[(mu, l, i, j)], p] = Zb[mu][l, k]
    return Zf
Zf = zfull()
# ad matrices of the highest weight vectors in the basis P
hw = [(mu, k) for mu, vs in types.items() for k in range(len(vs))]
Ahw = {}
for mu, k in hw:
    m = types[mu][k]
    # ad_m as a matrix on V'_1 : build from bracket with basis states
    Mx = sp.Matrix.hstack(*[vec(V.clean(V.mode_vec(dict((B1[i], mpq(int(sp.fraction(m[i])[0]), int(sp.fraction(m[i])[1]))) for i in range(n) if m[i] != 0), 0, {kk: mpq(1)}))) for kk in B1])
    Ahw[(mu, k)] = Pinv*Mx*P
eqs = set()
for mu, k in hw:
    if (mu, k) in (((2, 0), 0), ((0, 2), 0)): continue
    lhs = Zf*Ahw[(mu, k)]
    rhs = sum((Zb[mu][l, k]*Ahw[(mu, l)] for l in range(len(types[mu]))), sp.zeros(n, n))*Zf
    D = (lhs - rhs).applyfunc(sp.expand)
    for e in D:
        if e != 0: eqs.add(e)
print('equations:', len(eqs), 'unknowns:', len(Zf.free_symbols))
pickle.dump({'eqs': [str(e) for e in eqs], 'syms': [str(s) for s in Zf.free_symbols], 'label': label,
             'P': [list(P.col(i)) for i in range(n)]}, open(HERE + '/j1_equations_%s.pkl' % CASE, 'wb'))
