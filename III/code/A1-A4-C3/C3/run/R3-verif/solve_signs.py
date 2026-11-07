import sys, json
sys.path.insert(0, '.')
from lat import *
from fractions import Fraction as Fr
import sympy as sp
from lie_bracket import br

def br2(x, y):
    out = {}
    for kx, cx in x.items():
        for ky, cy in y.items():
            r = br({kx: Fr(1)}, {ky: Fr(1)})
            for kk, vv in r.items():
                out[kk] = out.get(kk, 0) + sp.Rational(vv.numerator, vv.denominator)*cx*cy
    return {k: sp.together(v) for k, v in out.items() if sp.simplify(v) != 0}

result = {}
for cp in (0, 1):
    syms = {}
    EF = []
    for X in range(4):
        E = {}; F = {}
        for j, (a, k) in enumerate(support(X)):
            c = 0 if any(a[:8]) else 1
            if c != cp: continue
            s = sp.Symbol('c%d_%d' % (X, j)); syms[(X, j)] = s
            E[('e', a)] = s; F[('e', vneg(a))] = -k/s
        EF.append((E, F))
    eqs = []
    for X in range(4):
        for Y in range(4):
            if X == Y: continue
            E1, F1 = EF[X]; E2, F2 = EF[Y]
            for r in (br2(E1, F2), br2(E1, E2), br2(F1, F2)):
                for v in r.values():
                    eqs.append(sp.numer(sp.together(v)))
    eqs = list(set(sp.expand(e) for e in eqs))
    print('copy', cp+1, 'n eqs', len(eqs))
    # gauge: torus. find maximal set of support roots linearly independent -> set those c = 1
    keys = list(syms)
    vecs = []
    chosen = []
    M = sp.Matrix([[0]*16])
    for key in keys:
        X, j = key
        a = support(X)[j][0]
        Mn = M.col_join(sp.Matrix([list(a)]))
        if Mn.rank() > M.rank():
            M = Mn; chosen.append(key)
    print('gauge-fixed (set to 1):', chosen)
    subs = {syms[k]: 1 for k in chosen}
    eqs2 = [sp.expand(e.subs(subs)) for e in eqs]
    eqs2 = [e for e in eqs2 if e != 0]
    free = [syms[k] for k in keys if k not in chosen]
    print('free', free)
    for e in eqs2: print('   ', e)
    sol = sp.solve(eqs2, free, dict=True)
    print('solutions', sol)
    result[cp] = (chosen, sol)
