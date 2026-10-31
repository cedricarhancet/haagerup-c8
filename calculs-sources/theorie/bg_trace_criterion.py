"""Bantay--Gannon criterion (arXiv:0705.2467, eq. (2.22)): the columns of Xi span the whole
module M(rho) of weakly holomorphic vector-valued modular functions iff
   Tr(Lambda) = 5d/12 + Tr(S)/4 + (2/(3 sqrt 3)) Re( exp(-pi i/6) Tr(U) ),  U = rho([[0,-1],[1,-1]]) = S T^{-1}.
"""
import mpmath as mp
from zh3_modular_data import zh3_S
mp.mp.dps = 40
from fractions import Fraction as Fr
Lam = [Fr(2,3)]*4 + [Fr(0), Fr(1,3), Fr(5,39), Fr(20,39), Fr(32,39), Fr(2,39), Fr(8,39), Fr(11,39)]
S = zh3_S()
d = 12
Tinv = mp.matrix(d, d)
for i in range(d):
    Tinv[i, i] = mp.exp(-2j*mp.pi*mp.mpf(Lam[i].numerator)/Lam[i].denominator)
U = S*Tinv
trS = sum(S[i, i] for i in range(d))
trU = sum(U[i, i] for i in range(d))
rhs = mp.mpf(5)*d/12 + trS/4 + 2/(3*mp.sqrt(3))*mp.re(mp.exp(-1j*mp.pi/6)*trU)
print('Tr S =', mp.nstr(trS, 20), ' Tr U =', mp.nstr(trU, 20))
print('rhs of (2.22) =', mp.nstr(rhs, 25), ' ; Tr Lambda =', sum(Lam))
