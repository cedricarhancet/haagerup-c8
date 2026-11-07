"""Proposition 8.4 (second part): among the partitions of the four lattice modules W_{e,a}, e*a != 0, into the classes
of their restrictions to S (two classes, by the first part), only {(1,1),(-1,-1)}, {(1,-1),(-1,1)} gives a matrix
Z = b b^T commuting with s^{V'}_{g,h} = (1/3) exp(-2 pi i (lambda_g, lambda_h)); the complex conjugate matrix gives the
same conclusion.  Exact arithmetic in Q(zeta_3) (sympy)."""
import sympy as sp, itertools, json
w = sp.exp(2 * sp.pi * sp.I / 3)
G = [(e, a) for e in range(3) for a in range(3)]
def B3(g, h): return (g[0] * h[0] * 4 + g[1] * h[1] * 2) % 3          # 3 (lambda_g, lambda_h) mod 3
idx = {g: i for i, g in enumerate(G)}
res = {}
for conj in (False, True):
    S = sp.Matrix(9, 9, lambda i, j: sp.Rational(1, 3) * w ** ((B3(G[i], G[j]) if conj else -B3(G[i], G[j])) % 3))
    four = [(1, 1), (2, 2), (1, 2), (2, 1)]
    parts = [[(1, 1), (2, 2)], [(1, 1), (1, 2)], [(1, 1), (2, 1)]]
    for p in parts:
        other = [g for g in four if g not in p]
        cols = [[(0, 0)], [(0, 0)], [(0, 1), (0, 2)], [(1, 0), (2, 0)], p, other]
        b = sp.zeros(9, len(cols))
        for j, gs in enumerate(cols):
            for g in gs: b[idx[g], j] += 1
        Z = b * b.T
        comm = sp.simplify(Z * S - S * Z)
        res[('conj' if conj else 'plain') + str(p)] = (comm == sp.zeros(9, 9))
print(json.dumps({k: bool(v) for k, v in res.items()}, indent=1))
json.dump({k: bool(v) for k, v in res.items()}, open('check_patterns.json', 'w'), indent=1)
