"""Is there an s-module U_3 with dim 956, the computed moments, L_{-1} U_2 <= U_3 <= V'_3 (multiplicities)?
Integer feasibility solved with scipy.optimize.milp (HiGHS); any solution found is verified exactly.
If the LP relaxation is infeasible, an exact Farkas certificate is produced and checked in rational arithmetic."""
import json
from collections import Counter
from fractions import Fraction as Fr
import numpy as np
from scipy.optimize import milp, LinearConstraint, Bounds, linprog
import sympy as sp
def M(n, p):
    return sum((n - 2*i)**p for i in range(n + 1))
PQ = [(2, 0), (0, 2), (4, 0), (0, 4), (2, 2), (6, 0), (0, 6), (4, 2), (2, 4)]
br = json.load(open('lattice_branching.json'))
u3 = json.load(open('u3_moments.json'))
u2e = json.load(open('u2_embedding_test.json'))
def parse(d):
    return {tuple(int(x) for x in k.strip('()').split(',')): v for k, v in d.items()}
out = {}
for name, full in (('A11+A139', 'A11+A139 (announced)'), ('A13+A113', 'A13+A113')):
    V3 = parse(br[name]['3'])
    U2 = parse(u2e[full]['solutions'][0])
    mom = {tuple(int(x) for x in k.strip('()').split(',')): int(sp.Rational(v)) for k, v in u3[full]['U3'].items()}
    types = sorted(V3)
    lower = np.array([U2.get(t, 0) for t in types]); upper = np.array([V3[t] for t in types])
    ok_bounds = bool(np.all(lower <= upper))
    Aeq = [[(t[0] + 1)*(t[1] + 1) for t in types]] + [[M(t[0], p)*M(t[1], q) for t in types] for p, q in PQ]
    beq = [956] + [mom[pq] for pq in PQ]
    Aeq = np.array(Aeq, dtype=float); beq = np.array(beq, dtype=float)
    res = milp(c=np.zeros(len(types)), constraints=LinearConstraint(Aeq, beq, beq), integrality=np.ones(len(types)),
               bounds=Bounds(lower, upper))
    line = '%s: bounds consistent %s; MILP status %s' % (name, ok_bounds, res.message)
    rec = {'status': res.message}
    if res.x is not None:
        x = np.round(res.x).astype(int)
        exact = all(sum(int(a)*int(b) for a, b in zip(row, x)) == int(bb) for row, bb in zip(Aeq.astype(int).tolist(), beq.astype(int).tolist()))
        sol = {str(t): int(v) for t, v in zip(types, x) if v}
        line += '; solution verified exactly: %s' % exact
        rec['solution'] = sol; rec['verified'] = exact
    else:
        lp = linprog(np.zeros(len(types)), A_eq=Aeq, b_eq=beq, bounds=list(zip(lower, upper)), method='highs')
        line += '; LP relaxation: %s' % lp.message
        rec['lp'] = lp.message
    print(line)
    if 'solution' in rec:
        print('    U_3 =', rec['solution'])
    out[name] = rec
json.dump(out, open('u3_embedding_test.json', 'w'), indent=1)
