"""Referee check: plug the 4 points of the Groebner variety into ALL equations of r3a_groebner.py (exactly, sympy,
over Q(sqrt(-6))) and compare with the 4 candidates of omegaC1.pkl (G-coordinates)."""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pickle, itertools
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
params = list(params); m = len(params)
ts = sp.symbols('t0:%d' % m)
u = list(sol.subs(dict(zip(params, ts))))
s6 = sp.sqrt(-6)
pts = []
for e7 in (1, -1):
    for e2 in (1, -1):
        t7 = sp.Rational(2, 39)*e7; t2 = e2*s6/78; t4 = 39*t2*t7
        pts.append({ts[0]: sp.Rational(-1, 39), ts[1]: sp.Rational(-8, 195), ts[2]: t2, ts[3]: 0, ts[4]: t4,
                    ts[5]: 0, ts[6]: t7/2, ts[7]: t7})
P = {k: [Qr(x) for x in v] for k, v in prod.items()}
F = [[Qr(x) for x in r] for r in form]
cands = pickle.load(open(os.path.join(CM.HERE, 'omegaC1.pkl'), 'rb'))['cands']
def cvec(c):
    return [Qr(a) + Qr(b)*s6 for a, b in c[0]]
CV = [cvec(c) for c in cands]
for ip, pt in enumerate(pts):
    uu = [sp.nsimplify(sp.expand(x.subs(pt))) for x in u]
    # idempotency
    prodv = [0]*n
    for i in range(n):
        for j in range(n):
            if uu[i] == 0 or uu[j] == 0: continue
            for k in range(n):
                if P[(i, j)][k]: prodv[k] += uu[i]*uu[j]*P[(i, j)][k]
    idem = all(sp.simplify(prodv[k] - 2*uu[k]) == 0 for k in range(n))
    nrm = sp.simplify(sum(F[i][j]*uu[i]*uu[j] for i in range(n) for j in range(n)))
    tr_ok = all(sp.simplify(sum(r[g]*uu[g] for g in range(n)) - b) == 0 for r, b in zip(rows, rhs))
    # full characteristic polynomial on every component (stronger than e_1, e_2)
    char_ok = True
    for kappa, secs in sorted(comp.items()):
        k = len(secs)
        M = sp.Matrix(k, k, lambda i, j: sum(Qr(Mk[kappa][g][i][j])*uu[g] for g in range(n) if Mk[kappa][g][i][j]))
        x = sp.Symbol('x')
        cp = sp.expand(M.charpoly(x).as_expr())
        target = sp.expand(sp.prod([x - sp.Rational(str(CM.omega_C1_eig(kappa, z))) for z in secs]))
        if sp.simplify(cp - target) != 0: char_ok = False
    match = [ic for ic, cv in enumerate(CV) if all(sp.simplify(a - b) == 0 for a, b in zip(uu, cv))]
    print('point', ip, ': idempotent', idem, '; (u,u) =', nrm, '; traces', tr_ok, '; full charpoly on all W_1 comps', char_ok,
          '; equals candidate', match)
