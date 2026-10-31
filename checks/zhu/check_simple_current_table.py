"""Exact elementary arithmetic for the seven-term simple-current identity.

The default run uses only Python's standard library.  --source-dir additionally
compares the two displayed vectors with exact_engine.setup(0); this comparison
is a normalization check, not a prerequisite for the hand calculation.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse,json,os,sys
parser=argparse.ArgumentParser()
parser.add_argument('--source-dir',type=Path)
parser.add_argument('--output',type=Path)
args=parser.parse_args()
out=args.output.resolve() if args.output else None
G=[[2,0,-1,0,0,0,0,0],[0,2,0,-1,0,0,0,0],[-1,0,2,-1,0,0,0,0],[0,-1,-1,2,-1,0,0,0],[0,0,0,-1,2,-1,0,0],[0,0,0,0,-1,2,0,0],[0,0,0,0,0,0,2,-1],[0,0,0,0,0,0,-1,2]]
alpha=[(1,1,2,3,2,1,1,0),(1,2,3,4,3,2,0,0),(1,1,2,3,2,1,0,1),(2,2,3,4,3,1,0,0),(0,1,1,2,2,1,1,1),(1,1,2,2,1,0,1,1),(1,0,1,1,1,1,1,1)]
w=(2,2,4,5,4,2,1,1)
beta=[tuple(x-y for x,y in zip(w,a)) for a in alpha]
a=[(F(39,5),F(1,5),F(-13,5),F(1,5)),(F(-104,5),0,0,F(-56,5)),(F(39,5),F(-1,5),F(13,5),F(1,5)),(-52,0,0,12),(0,F(-19,15),F(-13,15),0),(0,1,0,0),(0,F(3,5),F(-39,5),0)]
b=[(F(-13,4),F(3,4),F(13,4),F(1,4)),(0,1,0,0),(F(-13,4),F(-3,4),F(-13,4),F(1,4)),(0,F(-1,2),F(13,10),0),(F(39,2),0,0,F(9,2)),(13,0,0,7),(F(39,2),0,0,F(-3,2))]
signs=[1,-1,-1,-1,-1,1,-1]
labels1=(-1,-2,2,-1,2,-1,0,0)
labels2=(1,0,0,1,0,1,2,2)
def dot(x,y):return sum(xx*yy for xx,yy in zip(x,y))
def ip(x,y):return sum(x[i]*G[i][j]*y[j] for i in range(8) for j in range(8))
def eps(x,y):return -1 if (dot(x,y)+sum(x[i]*y[j]*G[i][j] for i in range(8) for j in range(i+1,8)))%2 else 1
def mul(x,y):
 z=[F(0)]*4
 for i,xx in enumerate(x):
  for j,yy in enumerate(y):z[i^j]+=xx*yy*(13 if i&j&1 else 1)*(-3 if i&j&2 else 1)
 return z
def solve(M,v):
 M=[[F(t) for t in row]+[F(v[i])] for i,row in enumerate(M)]
 for j in range(len(v)):
  k=next(i for i in range(j,len(v)) if M[i][j]);M[j],M[k]=M[k],M[j]
  r=M[j][j];M[j]=[t/r for t in M[j]]
  for i in range(len(v)):
   if i!=j:
    r=M[i][j];M[i]=[t-r*u for t,u in zip(M[i],M[j])]
 return tuple(row[-1] for row in M)
h1,h2=solve(G,labels1),solve(G,labels2)
assert ip(h1,h1)==6 and ip(h2,h2)==26 and ip(h1,h2)==0
assert tuple((x+y)/2 for x,y in zip(h1,h2))==w
assert ip(w,w)==8 and (dot(labels1,w),dot(labels2,w))==(3,13)
for i in range(7):
 assert ip(alpha[i],alpha[i])==ip(beta[i],beta[i])==4
 assert (dot(labels1,alpha[i]),dot(labels2,alpha[i]))==(1,7)
 assert (dot(labels1,beta[i]),dot(labels2,beta[i]))==(2,6)
 assert ip(alpha[i],beta[i])==0 and eps(alpha[i],beta[i])==signs[i]
 for j in range(7):
  if j!=i:assert ip(alpha[i],beta[j])>0
terms=[[sign*t for t in mul(x,y)] for x,y,sign in zip(a,b,signs)]
k=[26,104,26,104,65,65,117]
assert terms==[[F(0),F(t,5),F(7*t,5),F(0)] for t in k]
total=[sum(t[i] for t in terms) for i in range(4)]
assert total==[0,F(507,5),F(3549,5),0]
energy=lambda x,y:F(x*x,12)+F(y*y,52)
assert energy(5,13)==F(16,3) and energy(3,15)==F(66,13)
assert min(energy(x,y) for x,y in [(3,13),(1,13),(5,13),(3,11),(3,15)])==F(40,13)>3
compared=False
if args.source_dir:
 source=args.source_dir.resolve();sys.path.insert(0,str(source));os.chdir(source)
 import exact_engine as E
 data=E.setup(0)
 for name,states,coeffs in [('(1, 7)',alpha,a),('(2, 6)',beta,b)]:
  expected={(x,()):tuple(F(t) for t in c) for x,c in zip(states,coeffs)}
  found={key:tuple(F(str(t)) for t in c) for key,c in E.kv_from_terms(data['hw'][name]).items()}
  assert expected==found
 compared=True
result={'status':'PASS','basis':['1','sqrt(13)','sqrt(-3)','sqrt(-39)'],'coefficient':[str(t) for t in total],'seven_pair_coefficients':[str(F(t,5)) for t in k],'other_pairs_positive_inner_product':42,'cartan_norms':[6,26],'minimal_positive_mode_energy':'40/13','compared_with_exact_engine':compared}
if out:out.write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result))
