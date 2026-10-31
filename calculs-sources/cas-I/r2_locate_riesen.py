"""Axis A0/A1: locate Riesen's generators e, f, v (transported, r1) in our case (I) data:
 (1) e, f in B_I = Com_{V'}(A)_2 : coordinates in the basis of griess_data('I'); Griess subalgebra 5 or 6?
 (2) v in the type (1,13) primaries Pi; in the admissible space E (subalgebra 5 or 6); which isotropic line?
 (3) action of the centralizer Z (monomial elements, j6) on the two lines."""
import sys, os, json, pickle, itertools
HERE = os.path.abspath(os.path.dirname(__file__)); GF = os.path.abspath(HERE + '/../cas-II')
sys.path.insert(0, GF); os.chdir(GF)
import sympy as sp
from gmpy2 import mpq
import lattice_voa as V, griess_subalgebras as GS, isotropic_lines as IL
FULL = {n: {k: sp.sympify(c) for k, c in d.items()} for n, d in pickle.load(open(HERE + '/r1_riesen_transported.pkl', 'rb')).items()}
T, basis, prod, gram = GS.griess_data('I')
def rat(x): return sp.Rational(int(x.numerator), int(x.denominator))
from flint import fmpq_mat, fmpq
R13, R3, II = sp.symbols('R13 R3 II')
def split(x):
    """x in Q(sqrt13, sqrt3, i) -> dict monomial -> Rational (monomials in R13, R3, II reduced)"""
    e = sp.expand(sp.sympify(x)).xreplace({sp.sqrt(39): R13*R3, sp.sqrt(13): R13, sp.sqrt(3): R3, sp.I: II})
    e = sp.expand(e)
    P = sp.Poly(e, R13, R3, II)
    out = {}
    for (a, b, c), co in P.terms():
        f = sp.Integer(13)**(a//2)*sp.Integer(3)**(b//2)*(-1)**(c//2)
        key = (a % 2, b % 2, c % 2); out[key] = out.get(key, 0) + co*f
    return {k: v for k, v in out.items() if v != 0}
def coords(vecs, target, consts=None):
    """exact coordinates of target (coefficients in Q(sqrt13, sqrt3, i)) in the span of rational vectors"""
    keys = sorted(set().union(*[set(v) for v in vecs]) | set(target))
    def r_(v, k):
        x = v.get(k, 0)
        return sp.Rational(x) if isinstance(x, sp.Basic) else (rat(x) if k in v else sp.Integer(0))
    M = fmpq_mat(len(keys), len(vecs), [fmpq(int(r_(v, k).p), int(r_(v, k).q)) for k in keys for v in vecs])
    parts = {}
    for i, k in enumerate(keys):
        for mono, c in split(target.get(k, 0)).items():
            parts.setdefault(mono, [sp.Integer(0)]*len(keys))[i] = c
    sol = [sp.Integer(0)]*len(vecs)
    for mono, col in parts.items():
        b = fmpq_mat(len(keys), 1, [fmpq(int(sp.Rational(c).p), int(sp.Rational(c).q)) for c in col])
        A = fmpq_mat(len(keys), len(vecs) + 1, [M[i, j] if j < len(vecs) else b[i, 0] for i in range(len(keys)) for j in range(len(vecs) + 1)])
        Rr, rank = A.rref()
        x = [fmpq(0)]*len(vecs)
        for i in range(rank):
            j = next(j for j in range(len(vecs) + 1) if Rr[i, j] != 0)
            assert j < len(vecs), 'not in the span'
            x[j] = Rr[i, len(vecs)]
        # exact check
        for i in range(len(keys)):
            assert sum((M[i, j]*x[j] for j in range(len(vecs))), fmpq(0)) == b[i, 0], 'not in the span'
        m = R13**mono[0]*R3**mono[1]*II**mono[2]
        for j in range(len(vecs)): sol[j] += sp.Rational(int(x[j].p), int(x[j].q))*m
    return [sp.expand(c.subs({R13: sp.sqrt(13), R3: sp.sqrt(3), II: sp.I})) for c in sol]
res = {}
ce = coords(basis, FULL['e']); cf = coords(basis, FULL['f']); cw = coords(basis, FULL['omega_C'])
print('omega_C coords:', cw); print('e coords:', ce); print('f coords:', cf)
G = json.load(open('griess_subalgebras.json'))['I']['subalgebras']
X = [sp.Matrix([0] + [sp.sympify(c) for c in s['x']]) for s in G]
def in_sub(c, i):
    M = sp.Matrix.hstack(sp.Matrix([1] + [0]*(len(basis) - 1)), X[i], sp.Matrix(c))
    return sp.simplify(M.det()) == 0 if M.shape[0] == 3 else M.applyfunc(sp.radsimp).rank(simplify=True) == 2
res['e,f in subalgebra'] = [i for i in range(len(G)) if in_sub(ce, i) and in_sub(cf, i)]
print('e, f lie in Griess subalgebras:', res['e,f in subalgebra'])
ci = res['e,f in subalgebra'][0]
D, out = IL.admissible_spaces('I', ci)
tau = (1, 13)
P = out[tau]['P']; spaces = out[tau]['spaces']
print('type (1,13): dim Pi =', len(P), 'admissible spaces:', [(s['E'].shape[1], str(s['xi'])) for s in spaces])
cv = coords(P, FULL['v'], consts=(sp.sqrt(13), sp.sqrt(3), sp.sqrt(39)))
print('v in Pi, coordinates:', cv)
Em = spaces[0]['E']
cvec = sp.Matrix(cv)
rows = next(r for r in itertools.combinations(range(Em.rows), 2) if sp.simplify(Em[list(r), :].det()) != 0)
cvE = list((Em[list(rows), :].inv()*cvec[list(rows), :]).applyfunc(lambda z: sp.expand(sp.radsimp(z))))
chk = (Em*sp.Matrix(cvE) - cvec).applyfunc(lambda z: sp.simplify(sp.expand(z)))
print('v in E, coordinates:', cvE, '; exact residual zero:', chk == sp.zeros(Em.rows, 1))
res['v in E (exact)'] = chk == sp.zeros(Em.rows, 1)
L = IL.lines_in(spaces[0])
print('isotropic lines in E:', [list(l) for l in L])
which = [j for j, l in enumerate(L) if sp.simplify(cvE[0]*l[1] - cvE[1]*l[0]) == 0]
print('v is the line', which)
res.update({'subalgebra': ci, 'v_coords_E': [str(x) for x in cvE], 'lines': [[str(c) for c in l] for l in L], 'v_line': which})
json.dump(res, open(HERE + '/r2_locate_riesen.json', 'w'), indent=1)
pickle.dump({'E': [[str(x) for x in Em.row(i)] for i in range(Em.rows)], 'P': P, 'lines': [[str(c) for c in l] for l in L]}, open(HERE + '/r2_type113.pkl', 'wb'))
