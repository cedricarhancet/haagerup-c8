"""Insert the explicit data of the chain U < V(E6+A2) < V(E8) into the weight-two relations.

For h in a Cartan subalgebra of s = U_1 < e6 + sl3, normalized by <h,h> = 1 (form of e8 at level one),
write h = (h1, h2), u = (h1,h1) in [0,1], (h2,h2) = 1 - u, and T = Tr_s (ad h)^2.
Lowest spaces (Axis A of the paper): R_b = (e6+sl3) - s, R_a = (27,3), R_c0 = (27,3bar), R_c1 = (1,3), R_c2 = (27,1).
Traces: Tr_{e6} ad^2 = 24 (h1,h1), Tr_{sl3} ad^2 = 6 (h2,h2), Tr_27 = 6 (h1,h1), Tr_3 = (h2,h2).
Unknown traces tau_j = Tr h(0)^2 on the lowest space of M_{dj} (tau_j >= 0 for real h).
Usage: python3 weight2_solve.py
"""
import json
import sympy as sp

u, T = sp.symbols('u T', real=True)
tau = sp.symbols('tau1:7', real=True)

def data(gamma):
    xd = [7 - gamma, 42 + 16*gamma, 1 if gamma == 1 else 119, 5 - 3*gamma, 13 - 5*gamma, 14 + 7*gamma]
    w = {}
    w['w0'] = sp.Rational(-1, 12)
    w['w0p'] = T - sp.Rational(1, 12)*((6 + 13*gamma) - 24)
    t_b = 24*u + 6*(1 - u) - T
    t_a = 18*u + 27*(1 - u)
    w['w_b'] = t_b - sp.Rational(80 - 13*gamma, 12)
    w['w_a'] = t_a - sp.Rational(81, 12)
    w['w_c0'] = t_a - sp.Rational(81, 12)
    w['w_c1'] = (1 - u) - sp.Rational(3, 12)
    w['w_c2'] = 6*u - sp.Rational(27, 12)
    for j in range(6):
        w['w_d%d' % (j + 1)] = tau[j] - sp.Rational(xd[j], 12)
    return w

results = {}
for gamma in (0, 1):
    R = json.load(open('weight2_relations_gamma%d.json' % gamma))
    w = data(gamma)
    eqs = [sp.nsimplify(sum(sp.Rational(c)*w[lab] for c, lab in zip(row, R['labels']))) for row in R['relations']]
    sol = sp.solve(eqs, list(tau) + [T], dict=True)
    print('gamma =', gamma)
    if not sol:
        # try solving for tau only, report the residual constraints on (u, T)
        sol_t = sp.solve(eqs, list(tau), dict=True)
        print('  no solution with T free; solving for tau only:', sol_t)
        lin = sp.linear_eq_to_matrix(eqs, list(tau))
        print('  system matrix rank:', lin[0].rank(), ' augmented rank:', lin[0].row_join(lin[1]).rank())
    for s in sol:
        for k in list(tau) + [T]:
            if k in s:
                print('   ', k, '=', sp.factor(s[k]))
    results[gamma] = sol
    # consistency check: remaining equations
    if sol:
        res = [sp.simplify(e.subs(sol[0])) for e in eqs]
        print('  residuals:', res)
