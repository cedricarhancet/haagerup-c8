"""Independent exhaustive norm <= 8 enumeration using rational LDL bounds.
No breadth-first/root-neighbour enumeration is used on the checking side.
"""
from pathlib import Path
import json,math,sys
import sympy as sp
work=Path(__file__).resolve().parent/'work';sys.path.insert(0,str(work))
import lattice_voa as V
G=sp.Matrix(V.G);L,D=G.LDLdecomposition(hermitian=False);U=L.T
assert G==U.T*D*U
assert all(D[i,i]>0 for i in range(8))
vecs=set();x=[0]*8

def rec(i,rem):
 if i<0:
  vecs.add(tuple(x));return
 c=sum((U[i,j]*x[j] for j in range(i+1,8)),sp.Integer(0));p,q=int(c.p),int(c.q)
 bound=sp.floor(rem*q*q/D[i,i]);h=math.isqrt(int(bound))
 # low = ceil((-h-p)/q), high = floor((h-p)/q)
 low=-((h+p)//q);high=(h-p)//q
 for t in range(low,high+1):
  cost=D[i,i]*(t+c)**2
  assert cost<=rem
  x[i]=t;rec(i-1,rem-cost)
rec(7,sp.Integer(8))
assert vecs==set(V.VECS),(len(vecs),len(V.VECS),vecs-set(V.VECS),set(V.VECS)-vecs)
counts={n:sum(V.ip(v,v)==n for v in vecs) for n in range(0,9,2)}
result={'status':'PASS','method':'independent rational LDL ellipsoid enumeration','norm_limit':8,'vector_count':len(vecs),'counts_by_norm':counts,'ambient_dimensions_weights_0_to_4':[len(V.basis(n)) for n in range(5)],'identical_to_root_BFS':True,'rational_diagonal':[str(D[i,i]) for i in range(8)]}
Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
