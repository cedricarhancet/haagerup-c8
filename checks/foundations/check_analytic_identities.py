#!/usr/bin/env python3
"""Independent exact checks of constants used by the analytical arguments.
Requires sympy and the supplied classification source directory. The proofs
in the manuscript do not depend on these checks.
"""
import argparse
import json
import os
import sys
from pathlib import Path
from fractions import Fraction as F

parser = argparse.ArgumentParser()
parser.add_argument('--sources', type=Path, required=True)
args = parser.parse_args()
src = args.sources.resolve()
sys.path.insert(0, str(src))
os.chdir(src)
import sympy as sp
import e6_tools as E
import weight4_solve as W

c=list(map(F,['4/3','1','5/3','2','4/3','2/3']))
assert [sum(E.C[i][j]*c[j] for j in range(6)) for i in range(6)] == [1,0,0,0,0,0]
assert E.form(c) == F(4,3)
mu=c.copy();mu[0]-=1
assert E.form(mu,c) == F(1,3)
assert any((2*x).denominator != 1 for x in c)
assert E.kappas() == (F(3), F(18))
eqs, _, _ = W.weight4_system(0)
solutions = sp.solve(eqs, list(W.Q[:5]), dict=True)
assert len(solutions) == 1
sol = solutions[0]
assert all(sp.expand(e.subs(sol)) == 0 for e in eqs)
assert sp.expand(sol[W.Q[4]]-(12*W.T+2*W.T4-27*W.u**2-90*W.u+105)/26) == 0
x,y,k,l=sp.symbols('x y k l')
q=2*k*x**2+2*l*y**2
rhs=12*q*8*(x*x+y*y)+2*32*(x**4+y**4)-12*q*q
mixed=sp.expand(rhs).coeff(x,2).coeff(y,2)/(26*6)
assert sp.expand(mixed-8*(2*(k+l)-k*l)/13) == 0
values={str((a,b)):str(mixed.subs({k:a,l:b})) for a,b in [(1,1),(1,3),(3,3)]}
assert values == {'(1, 1)':'24/13','(1, 3)':'40/13','(3, 3)':'24/13'}
print(json.dumps({'status':'PASS','omega1_norm':'4/3','mu_omega1_pairing':'1/3','quartic_constants':['3','18'],'quartic_relation_count':len(eqs),'all_quartic_residuals_zero':True,'mixed_d5_exclusions':values},indent=2))
