"""Axis A0, control of r1: the transported conformal vector of V' and Sugawara vector of A equal ours (exact)."""
import sys, os, json, pickle
HERE = os.path.abspath(os.path.dirname(__file__)); sys.path.insert(0, HERE + '/../cas-II'); os.chdir(HERE + '/../cas-II')
import sympy as sp
import lattice_voa as V, cases as C, commutant_griess as CG
FULL = {n: {k: sp.sympify(c) for k, c in d.items()} for n, d in pickle.load(open(HERE + '/r1_riesen_transported.pkl', 'rb')).items()}
def rat(x): return sp.Rational(int(x.numerator), int(x.denominator))
def same(a, b): return all(sp.simplify(a.get(k, 0) - b.get(k, 0)) == 0 for k in set(a) | set(b))
om = {k: rat(x) for k, x in V.OMEGA.items()}
T = C.build('I'); T['case'] = 'I'; wA = {}; V.vadd(wA, CG.sugawara(T, 1)); V.vadd(wA, CG.sugawara(T, 2))
wA = {k: rat(x) for k, x in V.clean(wA).items()}
ef = {k: FULL['e'].get(k, 0) + FULL['f'].get(k, 0) for k in set(FULL['e']) | set(FULL['f'])}
res = {'psi_R(omega_W) = omega': same(FULL['omega_W'], om), 'psi_R(omega_A) = Sugawara(A)': same(FULL['omega_A'], wA),
       'psi_R(e) + psi_R(f) = psi_R(omega_C)': same(ef, FULL['omega_C'])}
print(res); json.dump(res, open(HERE + '/r1b_omega_check.json', 'w'), indent=1)
