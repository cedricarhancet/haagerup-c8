"""Independent check of the certificates produced by e8_rank2_table.py.

For each line (det, 2a, b, 2c, N | u | v) of e8_rank2_certificates_DMAX36.txt, with u, v in
doubled coordinates (a vector x of E8 is stored as 2x), we check in exact integer arithmetic:
 (1) u and v belong to E8 = D8 u (D8 + (1/2,...,1/2)): in doubled coordinates, all entries even
     or all entries odd, and sum of entries = 0 mod 4;
 (2) the Gram matrix of (u, v) is [[2a, b], [b, 2c]] and has determinant det;
 (3) Zu + Zv is primitive in E8: the gcd of the 2x2 minors of the coordinate matrix of (u, v)
     in a Z-basis of E8 equals 1. We use the basis of Conway--Sloane type
       2e_1, -e_1+e_2, ..., -e_6+e_7, (1/2,...,1/2),
     which is different from the Bourbaki basis used in e8_rank2_table.py; we check that its
     Gram matrix has determinant one;
 (4) the number of roots of E8 orthogonal to u and v equals N, the 240 roots being listed
     directly as +-e_i +- e_j and (1/2)(+-1,...,+-1) with an even number of minus signs.
We also check that the list contains every Gauss-reduced form [[2a,b],[b,2c]], 0 <= b <= a <= c,
with determinant at most the largest one listed. Only N >= 14 is used in the paper.
Usage: python3 check_e8_certificates.py e8_rank2_certificates_DMAX36.txt
"""
import sys
from fractions import Fraction
from itertools import combinations, product
from math import gcd


def ip(x, y):
    """Inner product of two vectors given in doubled coordinates."""
    s = sum(a*b for a, b in zip(x, y))
    assert s % 4 == 0
    return s//4


def in_e8(x):
    par = {e % 2 for e in x}
    return len(par) == 1 and sum(x) % 4 == 0


# roots in doubled coordinates
roots = []
for i, j in combinations(range(8), 2):
    for si, sj in product((2, -2), repeat=2):
        r = [0]*8; r[i] = si; r[j] = sj; roots.append(tuple(r))
for signs in product((1, -1), repeat=8):
    if signs.count(-1) % 2 == 0:
        roots.append(signs)
assert len(roots) == 240 and all(ip(r, r) == 2 and in_e8(r) for r in roots)

# Z-basis of E8 in doubled coordinates
B = [[4, 0, 0, 0, 0, 0, 0, 0]]
for k in range(6):
    r = [0]*8; r[k] = -2; r[k+1] = 2; B.append(r)
B.append([1]*8)
assert all(in_e8(b) for b in B)


def det_frac(M):
    M = [[Fraction(x) for x in row] for row in M]
    n = len(M); d = Fraction(1)
    for i in range(n):
        p = next((k for k in range(i, n) if M[k][i] != 0), None)
        if p is None:
            return Fraction(0)
        if p != i:
            M[i], M[p] = M[p], M[i]; d = -d
        d *= M[i][i]
        for k in range(i+1, n):
            f = M[k][i]/M[i][i]
            M[k] = [M[k][j] - f*M[i][j] for j in range(n)]
    return d


assert det_frac([[ip(x, y) for y in B] for x in B]) == 1     # unimodular: B is a Z-basis of E8
Bt = [[Fraction(B[j][i]) for j in range(8)] for i in range(8)]  # columns = basis vectors


def coords(x):
    """Solve Bt c = x exactly; c must be integral since x is in E8."""
    M = [row[:] + [Fraction(x[i])] for i, row in enumerate(Bt)]
    n = 8
    for i in range(n):
        p = next(k for k in range(i, n) if M[k][i] != 0)
        M[i], M[p] = M[p], M[i]
        for k in range(n):
            if k != i and M[k][i] != 0:
                f = M[k][i]/M[i][i]
                M[k] = [M[k][j] - f*M[i][j] for j in range(n+1)]
    c = [M[i][n]/M[i][i] for i in range(n)]
    assert all(t.denominator == 1 for t in c)
    return [int(t) for t in c]


fname = sys.argv[1] if len(sys.argv) > 1 else 'e8_rank2_certificates_DMAX36.txt'
ok = True; count = 0; minN = None; seen = set()
for line in open(fname):
    if line.startswith('#') or not line.strip():
        continue
    head, us, vs = line.split('|')
    det, A, b, C, N = map(int, head.split())
    u = tuple(map(int, us.split())); v = tuple(map(int, vs.split()))
    c1 = in_e8(u) and in_e8(v)
    c2 = (ip(u, u), ip(u, v), ip(v, v)) == (A, b, C) and A*C - b*b == det
    cu, cv = coords(u), coords(v)
    g = 0
    for i, j in combinations(range(8), 2):
        g = gcd(g, cu[i]*cv[j] - cu[j]*cv[i])
    c3 = (g == 1)
    n = sum(1 for r in roots if ip(r, u) == 0 and ip(r, v) == 0)
    c4 = (n == N)
    good = c1 and c2 and c3 and c4 and N >= 14
    seen.add((A//2, b, C//2))
    ok &= good; count += 1; minN = n if minN is None else min(minN, n)
    print(f'det {det:3d}  K=[[{A},{b}],[{b},{C}]]'.ljust(30),
          f'E8:{c1} Gram:{c2} primitive:{c3} roots={n} (claimed {N}) ->', 'OK' if good else 'FAIL')
# completeness: every even positive-definite binary lattice of determinant <= DMAX is isometric
# to [[2a,b],[b,2c]] with 0 <= b <= a <= c (Gauss reduction); all such triples must occur
DMAX = max(4*a*c - b*b for a, b, c in seen)
need = {(a, b, c) for a in range(1, DMAX+1) for c in range(a, DMAX+1) for b in range(0, a+1)
        if 0 < 4*a*c - b*b <= DMAX}
complete = (need == seen)
ok &= complete
print('reduced forms with det <=', DMAX, ':', len(need), ' all present:', complete)
print(f'{count} certificates checked; minimum number of orthogonal roots = {minN};',
      'ALL OK' if ok else 'SOME FAILURES')
