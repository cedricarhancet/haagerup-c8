"""Consistency of the weight-two relations (gamma = 0) with the lowest spaces ANNOUNCED by Riesen
(working manuscript of 2026-10-09, Theorem 11.4 and the profile table of Section 11.1; blocks V(a,b),
a = highest weight for the level-1 sl2, b = highest weight for the level-39 sl2, dim (a+1)(b+1)).
For h in the level-k factor normalized by <h,h> = 1, Tr_{V(a,b)} h(0)^2 = (a+1)(b+1)-weighted second
moment: for the level-1 factor sum over V(a) of m^2/2 times (b+1); for the level-39 factor (a+1) times
sum over V(b) of m^2/78.  Predicted values: tau_j from weight2 relations with u, T of each factor.
Usage: python3 gamma0_riesen_check.py
"""
from fractions import Fraction as Fr
import sympy as sp
import gamma1_exclusion as G     # reuses system(); importing runs its report, so silence it
blocks = {'d1': [(0, 6)], 'd2': [(0, 0), (0, 4), (0, 6), (0, 8), (1, 7), (1, 1)],
          'd3': [(0, 0), (0, 2), (0, 12), (1, 9), (1, 7), (1, 5), (1, 3), (1, 1), (0, 4), (0, 4), (0, 6), (0, 6), (0, 8), (0, 8)],
          'd4': [(0, 4)], 'd5': [(0, 4), (1, 3)], 'd6': [(0, 2), (0, 6), (1, 1)]}
dims = {'d1': 7, 'd2': 42, 'd3': 119, 'd4': 5, 'd5': 13, 'd6': 14}
def m2(n):   # sum of m^2 over weights n, n-2, ..., -n of the (n+1)-dim sl2 irrep
    return Fr(n*(n + 1)*(n + 2), 3)
for z, bl in blocks.items():
    assert sum((a + 1)*(b + 1) for a, b in bl) == dims[z], z
full, free = G.system(0)
for label, (k, i1, i2) in {'level 1': (1, 1, 0), 'level 39': (39, 35, 4)}.items():
    hv = 2
    subs = {G.T: Fr(2*hv, k), G.u: Fr(i1, k)}
    data = {}
    for z, bl in blocks.items():
        if k == 1:
            data[z] = sum((b + 1)*m2(a)/2 for a, b in bl)
        else:
            data[z] = sum((a + 1)*m2(b)/78 for a, b in bl)
    t6 = data['d6']
    pred = [sp.nsimplify(e.subs(subs).subs(G.tau[5], t6)) for e in full]
    obs = [data['d%d' % j] for j in range(1, 7)]
    print(label, ': observed (announced) tau =', [str(x) for x in obs])
    print('        predicted from relations =', [str(x) for x in pred])
    assert [sp.Rational(x.numerator, x.denominator) for x in obs] == pred
print('gamma = 0: the weight-two relations agree exactly with the announced lowest spaces')
