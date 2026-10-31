"""Standard-library verification of the stored 7x7 exact C4 matrices.
This checks the matrix identities, spectral multiplicities and independence of the
seven lattice vectors. It does not replace the VOA-product checks of check_commutant_c4.py.
Usage: python recheck_c4_matrices.py PATH_TO_C4_EXACT_WITNESSES_JSON
"""
from fractions import Fraction as Q
from pathlib import Path
import sys,json
path=Path(sys.argv[1]);d=json.loads(path.read_text());Z=(Q(0),)*4;ONE=(Q(1),Q(0),Q(0),Q(0));n=7
add=lambda a,b:tuple(x+y for x,y in zip(a,b))
neg=lambda a:tuple(-x for x in a)
def mul(a,b):
 r=[Q(0)]*4
 for i,x in enumerate(a):
  for j,y in enumerate(b):
   r[i^j]+=x*y*(13 if i&j&1 else 1)*(-3 if i&j&2 else 1)
 return tuple(r)
def mm(A,B):
 C=[[Z for j in range(n)] for i in range(n)]
 for i in range(n):
  for j in range(n):
   for k in range(n):C[i][j]=add(C[i][j],mul(A[i][k],B[k][j]))
 return C
I=[[ONE if i==j else Z for j in range(n)] for i in range(n)]
O=[[Z for j in range(n)] for i in range(n)]
def subscalar(A,x):
 return [[add(A[i][j],(-x,Q(0),Q(0),Q(0))) if i==j else A[i][j] for j in range(n)] for i in range(n)]
A=[[[*map(Q,x)] for x in row] for row in d['L1']];A=[[tuple(x) for x in row] for row in A]
B=[[[*map(Q,x)] for x in row] for row in d['L2']];B=[[tuple(x) for x in row] for row in B]
assert [[add(A[i][j],B[i][j]) for j in range(n)] for i in range(n)]==[[tuple(4*x for x in I[i][j]) for j in range(n)] for i in range(n)]
assert mm(A,B)==mm(B,A)
for M,values in [(A,[Q(0),Q(0),Q(2),Q(4),Q(4),Q(25,19),Q(50,19)]),(B,[Q(0),Q(0),Q(2),Q(4),Q(4),Q(26,19),Q(51,19)])]:
 product=I
 for x in sorted(set(values)):product=mm(product,subscalar(M,x))
 assert product==O
 power=I
 for k in range(5):
  tr=Z
  for i in range(n):tr=add(tr,power[i][i])
  assert tr==(sum(x**k for x in values),Q(0),Q(0),Q(0))
  power=mm(power,M)
p=d['prime'];r=d['sqrt13_modp'];q=d['sqrt_minus3_modp'];assert r*r%p==13 and q*q%p==p-3
def red(x):
 x=Q(x);return x.numerator*pow(x.denominator,-1,p)%p
rows=[]
for vec in d['basis']:
 row={}
 for lam,osc,a in vec:
  key=(tuple(lam),tuple(tuple(x) for x in osc));v=sum(red(x)*y for x,y in zip(a,[1,r,q,r*q]))%p
  if v:row[key]=v
 rows.append(row)
piv=[]
for row in rows:
 for key,b in piv:
  c=row.get(key,0)
  if c:
   for k,v in b.items():
    t=(row.get(k,0)-c*v)%p
    if t:row[k]=t
    else:row.pop(k,None)
 assert row
 k=min(row);inv=pow(row[k],-1,p);piv.append((k,{j:v*inv%p for j,v in row.items()}))
assert len(piv)==7
out={'result':'PASS','basis_rank_modp':7,'commuting_matrices':True,'sum_4I':True,'annihilating_polynomials':True,'spectral_multiplicities':True,'VOA_products_recomputed':False}
path.with_name('c4-matrix-recheck.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
