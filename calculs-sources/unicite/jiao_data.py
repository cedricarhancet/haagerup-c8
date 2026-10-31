"""Jiao's explicit data ([Jia26], Section 3), in his root coordinates (edges 12,23,34,45,36,78), coefficients in
K = Q(t, z), t^2 = 13, z^2 = -3, stored as (a, b, c, d) = a + b t + c z + d tz (same basis as Kfield: 1, sqrt13,
sqrt-3, sqrt13 sqrt-3).  Lattice vectors: 'e' entries {momentum: K coefficient}; Heisenberg: 'h' (vector)."""
from fractions import Fraction as F
def k(a=0, b=0, c=0, d=0): return (F(a), F(b), F(c), F(d))
r1, r2, r3 = (1,1,1,0,0,0,0,0), (0,1,1,1,0,0,0,0), (0,0,1,1,1,0,0,0)
q = (1,2,3,2,1,2,0,0)
def A(i): v = [0]*8; v[i-1] = 1; return tuple(v)
def neg(v): return tuple(-x for x in v)
PE = {r1: k(1), r2: k(1), r3: k(1)}
PF = {neg(r1): k(-1), neg(r2): k(-1), neg(r3): k(-1)}
PH = (1,2,3,2,1,0,0,0)
BE = {A(1): k(2), A(2): k(1), A(4): k(2), A(5): k(1), q: k(1), A(7): k(2), A(8): k(1)}
BF = {neg(A(1)): k(-1), neg(A(2)): k(-2), neg(A(4)): k(-1), neg(A(5)): k(-2), neg(q): k(-1), neg(A(7)): k(-1), neg(A(8)): k(-2)}
BH = (3,4,3,4,3,2,2,2)
U26 = {(0,1,1,1,1,0,1,1): k(0, F(1,6), F(1,6)), (1,1,1,1,0,0,1,1): k(0, F(1,3), F(-1,3)),
       (1,1,1,1,1,0,0,1): k(F(1,6), 0, F(-1,6)), (1,1,1,1,1,0,1,0): k(F(1,3), 0, F(1,3)),
       (1,1,1,2,1,0,0,0): k(F(1,6), 0, 0, F(-1,6)), (1,2,1,1,1,0,0,0): k(F(1,12), 0, 0, F(1,12)),
       (2,3,4,3,2,2,0,0): k(1)}
U17 = {(0,0,-1,0,0,-1,1,1): k(0, 0, 1), (1,1,2,2,1,1,1,1): k(0, F(1,2), F(1,6)),
       (1,2,2,1,1,1,1,1): k(0, F(1,4), F(-1,12)), (1,2,2,2,1,1,0,1): k(F(1,2), 0, F(1,6)),
       (1,2,2,2,1,1,1,0): k(1, 0, F(-1,3)), (1,2,2,2,2,1,0,0): k(F(1,2), 0, 0, F(-1,6)),
       (2,2,2,2,1,1,0,0): k(1, 0, 0, F(1,3))}
GJ = [[2 if i == j else 0 for j in range(8)] for i in range(8)]
for a, b in [(0,1),(1,2),(2,3),(3,4),(2,5),(6,7)]: GJ[a][b] = GJ[b][a] = -1
def ipJ(a, b): return sum(a[i]*GJ[i][j]*b[j] for i in range(8) for j in range(8))
def epsJ(p, qq):   # (-1)^(sum p_i q_i + sum_{i>j} p_i G_ij q_j)
    e = sum(p[i]*qq[i] for i in range(8)) + sum(p[i]*GJ[i][j]*qq[j] for i in range(8) for j in range(i))
    return -1 if e % 2 else 1
