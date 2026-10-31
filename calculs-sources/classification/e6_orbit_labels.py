"""Independent check of the Bala--Carter labels used in Section 8 (added during the referee pass).

For a weighted Dynkin diagram (Bourbaki order, node 2 attached to node 4) with characteristic H,
dim of the centralizer of e = dim g_0 + dim g_1, orbit dim = 78 - that, and the Dynkin index of the
sl2-triple = (sum over roots of alpha(H)^2) / 48  (since Tr_{e6}(ad H)^2 = 2 h^vee (H,H)/... normalized
so that a long-root sl2 has index 1).  Orbit dimensions are compared with the table of
Collingwood--McGovern (E6): A1 22, 3A1 40, A2+2A1 50, 2A2+A1 54, A5 64.
"""
import e6_centralizers as E
from fractions import Fraction
TABLE = {'A1': 22, '3A1': 40, 'A2+2A1': 50, '2A2+A1': 54, 'A5': 64}
cases = [((0, 1, 0, 0, 0, 0), 'A1'), ((0, 0, 0, 1, 0, 0), '3A1'), ((1, 0, 0, 1, 0, 1), '2A2+A1'),
         ((2, 1, 1, 0, 1, 2), 'A5')]
for lab, name in cases:
    g = [sum(r[i]*lab[i] for i in range(6)) for r in E.roots]
    cent = 6 + g.count(0) + g.count(1)
    idx = Fraction(sum(x*x for x in g), 48)
    print(lab, 'centralizer of e:', cent, ' orbit dim:', 78 - cent, ' Dynkin index:', idx,
          ' label:', name, ' table dim:', TABLE[name], ' OK' if TABLE[name] == 78 - cent else ' MISMATCH')
print('index of the principal sl2 of a regular subalgebra = sum of components (A1: 1, A2: 4, A5: 35):')
print('  A2+2A1 -> 6, 2A2+A1 -> 9; hence the index-9 triple belongs to 2A2+A1, not A2+2A1.')
