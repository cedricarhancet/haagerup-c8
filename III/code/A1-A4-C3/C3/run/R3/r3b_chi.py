"""R3b (exact): PCT operators of W compatible with the compact form of sl2^4, and their action on the four candidate
structures S x S (omega_{C x 1} = u_0..u_3, Groebner-certified list).
phi_chi(c mono e^beta) = (-1)^{#mono} chi(beta) conj(c) mono e^{-beta}, chi a positive character of L = E8+E8, is a PCT
operator of W (conjugate of the standard phi_0 by the positive torus element chi^{1/2}).  We solve chi(alpha) =
-d_alpha / c_alpha for alpha in the supports of E_i = sum c_alpha e^alpha, F_i = sum d_alpha e^{-alpha}, so that
phi_chi(E_i) = -F_i, phi_chi(H_i) = -H_i.  phi_chi preserves A x A and G = Com_W(A x A)_2, hence maps the Virasoro vector
omega_{C x 1} of a structure S x S to a solution of the same polynomial system: phi_chi permutes {u_0, ..., u_3}.
If phi_chi(u_k) = u_k, phi_chi preserves S x 1 = Com_W(1 x A, omega_CC - u_k) and 1 x S, hence S x S (generated in
weights 1, 2) and Prop. 7.12 (3) gives (H_U).  Output: r3b_chi.json"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pickle, json, time
from fractions import Fraction as Fr
from collections import defaultdict
import sympy as sp
from gmpy2 import mpq
import e8e8_voa as V
import common as CM
import qs6 as Q
t0 = time.time()
DG = pickle.load(open(os.path.join(CM.HERE, 'G.pkl'), 'rb'))
OC = pickle.load(open(os.path.join(CM.HERE, 'omegaC1.pkl'), 'rb'))
G = DG['G']
# lattice coordinates in a Z-basis: E8 simple roots in each copy
def lcoords(ib):
    y = V.VECS[ib]
    return V.coords8(y[:8]) + V.coords8(y[8:])
eqs = []
for nm in CM.NAMES:
    E, H, F = CM.SL[nm]
    for (b, m), c in E.items():
        assert not m
        d = F.get((V.NEGI[b], ()))
        assert d is not None, 'support of F is not -support of E'
        r = -Fr(int(d.numerator), int(d.denominator))/Fr(int(c.numerator), int(c.denominator))
        eqs.append((lcoords(b), r))
print('equations chi(alpha) = r_alpha:', len(eqs), '; all r > 0:', all(r > 0 for _, r in eqs))
# multiplicative system: solve prime by prime for the exponents
primes = sorted(set(p for _, r in eqs for p in list(sp.factorint(r.numerator)) + list(sp.factorint(r.denominator))))
A = sp.Matrix([list(c) for c, _ in eqs])
print('rank of the support lattice:', A.rank(), 'of 16')
expo = {}
for p in primes:
    b = sp.Matrix([sp.multiplicity(p, r.numerator) - sp.multiplicity(p, r.denominator) for _, r in eqs])
    sol, params = A.gauss_jordan_solve(b)
    assert len(params) == 0 or True
    sol = sol.subs({t: 0 for t in params})
    assert A*sol == b, 'inconsistent multiplicative system'
    expo[p] = sol
def chi(ib):
    c = sp.Matrix(lcoords(ib))
    val = sp.Integer(1)
    for p, e in expo.items():
        val *= sp.Integer(p)**((e.T*c)[0])
    return sp.nsimplify(val)
# check phi_chi(E) = -F
def phi_chi_vec(v):
    out = {}
    for (b, m), c in v.items():
        x = chi(b)
        out[(V.NEGI[b], m)] = (-1)**len(m)*x*sp.Rational(int(c.numerator), int(c.denominator))
    return out
ok = True
for nm in CM.NAMES:
    E, H, F = CM.SL[nm]
    pe = phi_chi_vec(E)
    ok &= all(pe.get(k, 0) == -sp.Rational(int(c.numerator), int(c.denominator)) for k, c in F.items()) and len(pe) == len(F)
print('phi_chi(E_i) = -F_i for the four sl2:', ok)
# action on G (rational basis): phi_chi^lin(g) for g in G, coordinates in G
cvals = {}
def phi_lin_G(g):
    out = {}
    for (b, m), c in g.items():
        if b not in cvals: cvals[b] = chi(b)
        x = cvals[b]
        assert x.is_rational, ('irrational chi on G support', x)
        out[(V.NEGI[b], m)] = mpq(int(x.p), int(x.q))*(-1)**len(m)*c
    return V.clean(out)
Mphi = []
for g in G:
    c = CM.span_coords(G, phi_lin_G(g)); assert c is not None, 'phi_chi does not preserve G'
    Mphi.append([CM.to_mq(x) for x in c])
n = len(G)
def apply(u):
    """phi_chi(sum u_i g_i) = sum conj(u_i) phi_lin(g_i)"""
    out = [Q.Z]*n
    for i in range(n):
        ci = Q.conj(u[i])
        if Q.iszero(ci): continue
        for k in range(n):
            if Mphi[i][k]: out[k] = Q.add(out[k], Q.scal(Mphi[i][k], ci))
    return out
cands = [c[0] for c in OC['cands']]
perm = []
for k, uk in enumerate(cands):
    img = apply(uk)
    j = [jj for jj, uj in enumerate(cands) if uj == img]
    perm.append(j[0] if j else None)
    print('phi_chi(u_%d) = u_%s' % (k, j[0] if j else 'NONE (?)'))
omCC = CM.span_coords(G, DG['omCC'])
print('phi_chi(omega_CC) = omega_CC:', apply([(CM.to_mq(x), mpq(0)) for x in omCC]) == [(CM.to_mq(x), mpq(0)) for x in omCC])
json.dump({'chi_rank': int(A.rank()), 'phiE_eq_minusF': bool(ok), 'perm': perm,
           'chi_on_simple_roots': {str(p): [str(x) for x in e] for p, e in expo.items()}}, open(os.path.join(CM.HERE, 'r3b_chi.json'), 'w'), indent=1)
print('time %.0fs' % (time.time() - t0))
