"""Theorem 6.1 (gamma = 1 is not realized), fully exact, from the EXACT weight-two relations
(weight2_relations_exact_gamma1.json, obstruction principle).  All steps are exact rational arithmetic.
(1) tau_j as affine functions of (u, T, tau5, tau6).
(2) Simple ideal g (indices i1, i2, level k = i1 + i2, u = i1/k, T = 2h/k): tau3 = 0 (h in [s,s], R_{d3} one-dim.)
    gives tau5 = T - 9u + 3; nonnegativity of tau1, tau4, tau5 gives 9u - 3 <= T <= 9u - 1, hence
    6 i1 <= 2h + 3 i2 and k <= (2h + 9 i2)/6; the finitely many cases are checked exactly.
(3) Central element (T = 0): exact elimination gives u >= 1/9.
(4) Enumeration of reductive s of dimension 19 with rank p1(s) <= 4: none.
"""
import json, itertools
from fractions import Fraction as Fr
import sympy as sp
u, T = sp.symbols('u T', real=True); tau = sp.symbols('tau1:7', real=True)
R = json.load(open('weight2_relations_exact_gamma1.json'))
xd = [6, 58, 1, 2, 8, 21]
w = {'w0': sp.Rational(-1, 12), 'w0p': T - sp.Rational(19 - 24, 12), 'w_b': 24*u + 6*(1 - u) - T - sp.Rational(67, 12),
     'w_a': 27 - 9*u - sp.Rational(81, 12), 'w_c0': 27 - 9*u - sp.Rational(81, 12), 'w_c1': 1 - u - sp.Rational(3, 12),
     'w_c2': 6*u - sp.Rational(27, 12)}
for j in range(6):
    w['w_d%d' % (j + 1)] = tau[j] - sp.Rational(xd[j], 12)
eqs = [sum(sp.Rational(c)*w[l] for c, l in zip(row, R['labels'])) for row in R['relations']]
sol = sp.solve(eqs, [tau[0], tau[1], tau[2], tau[3]], dict=True)[0]
assert all(sp.simplify(e.subs(sol)) == 0 for e in eqs)
full = [sp.expand(sol.get(t, t)) for t in tau]
print('(1) tau_j:'); [print('    tau%d =' % (j + 1), full[j]) for j in range(6)]
# (2) simple ideals
t5 = sp.solve(sp.Eq(full[2], 0), tau[4])[0]
print('(2) tau3 = 0  <=>  tau5 =', t5)
c1 = sp.expand(full[0].subs(tau[4], t5)); c4 = sp.expand(full[3].subs(tau[4], t5))
# tau1 >= 0 and tau4 >= 0 are affine in tau6 with opposite signs: eliminate tau6
a1 = c1.coeff(tau[5]); a4 = c4.coeff(tau[5])
assert a1 < 0 and a4 > 0
elim = sp.expand(c1*a4 - c4*a1)            # >= 0
print('    tau1 >= 0, tau4 >= 0  =>  %s >= 0 ;  tau5 >= 0  =>  %s >= 0' % (sp.factor(elim), t5))
SIMPLE = {'A1': (2, 3, 2), 'A2': (3, 8, 3), 'B2': (3, 10, 3), 'G2': (4, 14, 4), 'A3': (4, 15, 4)}   # h^vee, dim, rank? (rank below)
RANK = {'A1': 1, 'A2': 2, 'B2': 2, 'G2': 2, 'A3': 3}
allowed = []
for name, (hv, dim, _) in SIMPLE.items():
    for i2 in ((0, 1, 4) if name == 'A1' else (0,)):
        kmax = (2*hv + 9*i2)//6
        for i1 in range(1, kmax + 1):
            k = i1 + i2
            sub = {u: Fr(i1, k), T: Fr(2*hv, k)}
            if not (sp.Rational(elim.subs(sub)) >= 0 and sp.Rational(t5.subs(sub)) >= 0):
                continue
            # interval for tau6 with all tau_j >= 0
            lo, hi = Fr(0), None; ok = True
            for e in [full[j].subs(tau[4], t5) for j in range(6)] + [t5]:
                e = sp.expand(e.subs(sub)); a = sp.Rational(e.coeff(tau[5])); b = sp.Rational(e.subs(tau[5], 0))
                if a == 0:
                    ok = ok and b >= 0
                elif a > 0:
                    lo = max(lo, Fr(int((-b/a).p), int((-b/a).q)))
                else:
                    v = Fr(int((-b/a).p), int((-b/a).q)); hi = v if hi is None else min(hi, v)
            if ok and (hi is None or lo <= hi):
                allowed.append((name, i1, i2, k))
print('    allowed simple ideals (type, i1, i2, level):', allowed)
assert all(a[1] >= 1 for a in allowed)
# (3) central elements: T = 0, tau3, tau1, tau4 >= 0
c = [sp.expand(full[j].subs(T, 0)) for j in (2, 0, 3)]
# tau3 >= 0 is affine in tau5 only; tau1, tau4 involve tau5, tau6: eliminate tau6 then tau5
a1 = c[1].coeff(tau[5]); a4 = c[2].coeff(tau[5])
e14 = sp.expand(c[1]*a4 - c[2]*a1)       # >= 0, involves u, tau5
b3 = c[0].coeff(tau[4]); b14 = e14.coeff(tau[4])
assert b3 > 0 and b14 < 0
final = sp.expand(c[0]*(-b14) + e14*b3)  # >= 0, involves u only
print('(3) centre: eliminating tau6 and tau5 gives', sp.factor(final), '>= 0, i.e. u >=', sp.solve(sp.Eq(final, 0), u))
# (4) enumeration: rank p1(s) = sum rank g_i + dim z <= 4, at most one sl2 with i2 != 0, dim s = 19
types = sorted(set((a[0], SIMPLE[a[0]][1], RANK[a[0]], a[2] != 0) for a in allowed))
found = []
for counts in itertools.product(range(5), repeat=len(types)):
    if sum(cn for cn, t in zip(counts, types) if t[3]) > 1:
        continue
    rk = sum(cn*t[2] for cn, t in zip(counts, types)); dm = sum(cn*t[1] for cn, t in zip(counts, types))
    for z in range(0, 5):
        if rk + z <= 4 and dm + z == 19:
            found.append((counts, z))
print('(4) reductive s with dim 19 and rank p1(s) <= 4:', found if found else 'NONE')
assert not found
print('=> gamma = 1 is not realized (exact).')
