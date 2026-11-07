import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pickle, itertools, time
import sympy as sp
import common as CM
D = pickle.load(open(os.path.join(CM.HERE, 'G_griess.pkl'), 'rb'))
prod, form, Mk, comp = D['prod'], D['form'], D['Mk'], D['comp']
n = len(form)
Qr = lambda x: sp.Rational(int(x.numerator), int(x.denominator))
rows = []; rhs = []
for kappa, secs in sorted(comp.items()):
    rows.append([sum(Qr(Mk[kappa][g][i][i]) for i in range(len(secs))) for g in range(n)])
    rhs.append(sp.Rational(str(sum(CM.omega_C1_eig(kappa, s) for s in secs))))
sol, params = sp.Matrix(rows).gauss_jordan_solve(sp.Matrix(rhs))
ts = sp.symbols('t0:%d' % len(params))
u = list(sol.subs(dict(zip(params, ts))))
P = {k: [Qr(x) for x in v] for k, v in prod.items()}
uu = [0]*n
for i in range(n):
    for j in range(n):
        if u[i] == 0 or u[j] == 0: continue
        c = sp.expand(u[i]*u[j])
        for k in range(n):
            if P[(i, j)][k]: uu[k] += c*P[(i, j)][k]
eqs = [sp.expand(uu[k] - 2*u[k]) for k in range(n)]
eqs.append(sp.expand(sum(Qr(form[i][j])*u[i]*u[j] for i in range(n) for j in range(n) if form[i][j]) - sp.Rational(9, 5)))
for kappa, secs in sorted(comp.items()):
    nonS = [z for z in secs if z not in ('S1x1', '1xS1')]
    if len(nonS) != 2: continue
    k = len(secs)
    M = sp.Matrix(k, k, lambda i, j: sp.expand(sum(Qr(Mk[kappa][g][i][j])*u[g] for g in range(n) if Mk[kappa][g][i][j])))
    lam = [CM.omega_C1_eig(kappa, z) for z in nonS]
    e2 = sum(M.extract(list(c), list(c)).det() for c in itertools.combinations(range(k), 2))
    eqs.append(sp.expand(e2 - sp.Rational(str(lam[0]*lam[1]))))
eqs = [e for e in eqs if e != 0]
t0 = time.time()
GB = sp.groebner(eqs, *ts, order='grevlex', domain='QQ')
print('sympy groebner (%.0fs):' % (time.time()-t0), len(GB.exprs))
for g in GB.exprs: print('  ', g)
print('zero-dimensional:', GB.is_zero_dimensional)
