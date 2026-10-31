"""Test of the weight-four relations against the lowest spaces ANNOUNCED by Riesen (gamma = 0) (Axis B3).
Blocks V(a,b): a = highest weight for the level-1 sl2, b = for the level-39 sl2 (Riesen, working manuscript
2026-10-09, Theorem 11.4 and Section 11.1; same data as Calculs/gamma1/gamma0_riesen_check.py).
For h = H/sqrt(2k) in the level-k factor: tau = Tr H^2/(2k), Q = Tr H^4/(2k)^2, T = 4/k, T4 = 8/k^2,
u = i1/k with (k, i1) = (1, 1) and (39, 35).  The data are used ONLY as a test of the relations.
"""
from fractions import Fraction as Fr
import sympy as sp
import weight4_solve as W

blocks = {'d1': [(0, 6)], 'd2': [(0, 0), (0, 4), (0, 6), (0, 8), (1, 7), (1, 1)],
          'd3': [(0, 0), (0, 2), (0, 12), (1, 9), (1, 7), (1, 5), (1, 3), (1, 1), (0, 4), (0, 4), (0, 6), (0, 6), (0, 8), (0, 8)],
          'd4': [(0, 4)], 'd5': [(0, 4), (1, 3)], 'd6': [(0, 2), (0, 6), (1, 1)]}
def mom(n, p):
    return sum((n - 2*i)**p for i in range(n + 1))
eqs, tsol, R4 = W.weight4_system(0)
ok = True
for label, (k, i1) in {'level 1': (1, 1), 'level 39': (39, 35)}.items():
    tau_obs, Q_obs = [], []
    for j in range(1, 7):
        bl = blocks['d%d' % j]
        if k == 1:
            t2 = sum((b + 1)*mom(a, 2) for a, b in bl); t4 = sum((b + 1)*mom(a, 4) for a, b in bl)
        else:
            t2 = sum((a + 1)*mom(b, 2) for a, b in bl); t4 = sum((a + 1)*mom(b, 4) for a, b in bl)
        tau_obs.append(Fr(t2, 2*k)); Q_obs.append(Fr(t4, (2*k)**2))
    subs = {W.u: sp.Rational(i1, k), W.T: sp.Rational(4, k), W.T4: sp.Rational(8, k*k), W.tau[5]: sp.Rational(tau_obs[5].numerator, tau_obs[5].denominator)}
    for j in range(6):
        subs[W.Q[j]] = sp.Rational(Q_obs[j].numerator, Q_obs[j].denominator)
    res = [sp.simplify(e.subs(subs)) for e in eqs]
    print(label, ': announced Q =', [str(x) for x in Q_obs], ' residuals of the 9 relations:', res)
    ok = ok and all(r == 0 for r in res)
print('weight-four relations', 'AGREE' if ok else 'DISAGREE', 'with the announced lowest spaces')
