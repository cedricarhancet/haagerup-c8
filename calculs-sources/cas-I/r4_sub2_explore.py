"""Axis A1 (case I, Griess subalgebra 2, charges 4/5 and 686/205): exploration.  Types with finitely many admissible
lines: (1,13) (h1 = 0, dim E = 2) and (1,1) (h1 = 2/5, dim E = 2).  Test (a): for u highest weight in the (1,13) line
and w' in the g'+g''-module of the (1,1) line, u_(2) w' must lie in U_1; since (1,13) x (1,1) contains no type (2,0)
or (0,2), this means u_(2) w' = 0.  Exact computation (rational spaces E, lines over quadratic fields)."""
import sys, os, json, pickle, itertools, time
HERE = os.path.abspath(os.path.dirname(__file__)); GF = os.path.abspath(HERE + '/../cas-II')
sys.path.insert(0, GF); os.chdir(GF)
import sympy as sp
from gmpy2 import mpq
import lattice_voa as V, isotropic_lines as IL, closure_explore as CE
t0 = time.time()
D, out = IL.admissible_spaces('I', 2)
T = D['T']
info = {}
for tau, v in out.items():
    print(tau, [(s['E'].shape[1], str(s['h1']), str(s['h2'])) for s in v['spaces']], flush=True)
def rat(x): return sp.Rational(int(x.numerator), int(x.denominator))
def vec_from(E, P, coef):
    """sum_k coef_k E[:,k] in the basis P -> dict key -> sympy"""
    out_ = {}
    for i in range(E.rows):
        c = sp.expand(sum(coef[k]*E[i, k] for k in range(E.cols)))
        if c != 0:
            for key, x in P[i].items(): out_[key] = out_.get(key, 0) + c*rat(x)
    return {k: sp.expand(x) for k, x in out_.items() if sp.expand(x) != 0}
def lines_of(tau, h1):
    sp_ = [s for s in out[tau]['spaces'] if str(s['h1']) == h1][0]
    L = IL.lines_in(sp_)
    return sp_, [list(l) for l in L]
def lower(v, f):
    return {k: sp.expand(x) for k, x in V.mode_vec_sym(f, 0, v).items()} if hasattr(V, 'mode_vec_sym') else None
# symbolic products: lattice_voa works with mpq coefficients; we expand u = sum_i c_i P_i and use bilinearity
def prod_sym(uP, ucoef, wP, wcoef, n):
    res = {}
    for a, ca in zip(uP, ucoef):
        if ca == 0: continue
        for b, cb in zip(wP, wcoef):
            if cb == 0: continue
            for k, x in V.mode_vec(a, n, b).items():
                res[k] = res.get(k, 0) + ca*cb*rat(x)
    return {k: sp.expand(x) for k, x in res.items() if sp.expand(x) != 0}
s113, L113 = lines_of((1, 13), '0'); s11, L11 = lines_of((1, 1), '2/5')
print('(1,13) lines:', L113); print('(1,1) lines:', L11, flush=True)
P113 = out[(1, 13)]['P']; P11 = out[(1, 1)]['P']
# basis vectors of E as dicts, and their f-lowerings for the (1,1) module
def Ebasis(spc, P):
    return [vec_from(spc['E'], P, [1 if j == k else 0 for j in range(spc['E'].cols)]) for k in range(spc['E'].cols)]
E113 = Ebasis(s113, P113); E11 = Ebasis(s11, P11)
def tompq(d): return {k: mpq(int(sp.Rational(x).p), int(sp.Rational(x).q)) for k, x in d.items()}
E113q = [tompq(d) for d in E113]; E11q = [tompq(d) for d in E11]
mods11 = [CE.lowerings(T, (1, 1), e) for e in E11q]        # mods11[k][w]
res = {}
for i, l1 in enumerate(L113):
    for j, l2 in enumerate(L11):
        worst = []
        nz = False
        for w in range(len(mods11[0])):
            z = prod_sym(E113q, l1, [m[w] for m in mods11], l2, 2)
            if z: nz = True; worst.append((w, len(z)))
        res['(1,13) line %d x (1,1) line %d' % (i, j)] = 'u_2 w != 0 (fails)' if nz else 'u_2 w = 0 (passes)'
        print('(1,13) line', i, 'x (1,1) line', j, ':', res['(1,13) line %d x (1,1) line %d' % (i, j)], worst[:3], flush=True)
json.dump({'lines_113': [[str(c) for c in l] for l in L113], 'lines_11': [[str(c) for c in l] for l in L11], 'test_a': res,
           'time': time.time() - t0}, open(HERE + '/r4_sub2_explore.json', 'w'), indent=1)
