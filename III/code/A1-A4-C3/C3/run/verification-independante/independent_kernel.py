"""Independent lattice-current mode computation of the commutant at weight 2.
This audit imports only our independent root-system code, standard library, and FLINT for rank.
No original VOA module is imported. Formula:
Y(e^a,z)e^b=eps(a,b) z^(a,b) exp(sum_{m>0} a(-m) z^m/m)e^(a+b).
Commuting e^a(n) past h_i(-m) gives h_i(-m)e^a(n)-a_i e^a(n-m).
"""
import sys,runpy,time
from pathlib import Path
from functools import lru_cache
from collections import defaultdict
from flint import nmod_mat
GLOB=runpy.run_path(str(Path(__file__).with_name('independent_triples.py')))
globals().update({k:GLOB[k] for k in ('Q','Z','dot','eps','triples','norm4G','weight','d','zero')})
t0=time.monotonic()
def add(out,v,c=Q(1)):
 for k,a in v.items():
  z=out.get(k,0)+a*c
  if z:out[k]=z
  elif k in out:del out[k]
 return out
def oscillator(v,p):return {(a,tuple(sorted(mono+(p,)))):c for (a,mono),c in v.items()}
@lru_cache(None)
def exponential(a,degree):
 if degree==0:return {():Q(1)}
 out={}
 # Recurrence n S_n = sum_{m=1}^n a(-m) S_{n-m}.
 for m in range(1,degree+1):
  for mono,c in exponential(a,degree-m).items():
   for i,ai in enumerate(a):
    if ai:
     mm=tuple(sorted(mono+((m,i),)));out[mm]=out.get(mm,0)+Q(ai,2)*c/degree
 return {k:c for k,c in out.items() if c}
@lru_cache(None)
def root_mode(a,n,state):
 b,mono=state
 if mono:
  m,i=mono[0];rest=(b,mono[1:])
  ans=oscillator(root_mode(a,n,rest),(m,i))
  return add(ans,root_mode(a,n-m,rest),-Q(a[i],2))
 degree=-n-1-Q(dot(a,b),4)
 if degree<0:return {}
 assert degree.denominator==1
 gamma=tuple(x+y for x,y in zip(a,b))
 return {(gamma,m):eps(a,b)*c for m,c in exponential(a,int(degree)).items()}
def h_mode(h,n,state):
 a,mono=state
 if n<0:
  return {(a,tuple(sorted(mono+((-n,i),)))):c for i,c in enumerate(h) if c}
 if n==0:
  c=Q(dot(h,a),2);return {state:c} if c else {}
 ans={}
 for p,(m,i) in enumerate(mono):
  if m==n:
   mm=mono[:p]+mono[p+1:];k=(a,mm);ans[k]=ans.get(k,0)+n*h[i]
 return {k:c for k,c in ans.items() if c}
def current_mode(X,n,state):
 h,rs=X;out=h_mode(h,n,state)
 for a,c in rs.items():add(out,root_mode(a,n,state),c)
 return out
basis=[(Z,((2,i),)) for i in range(16)]+[(Z,((1,i),(1,j))) for i in range(16) for j in range(i,16)]+[(a,()) for a in sorted(norm4G) if weight(a)==zero]
assert len(basis)==644
currents=[x for t in triples for x in t]
columns=[]
for state in basis:
 col={}
 for i,x in enumerate(currents):
  for n in range(3):
   for target,c in current_mode(x,n,state).items():col[(i,n,target)]=c
 columns.append(col)
rows=sorted(set(k for c in columns for k in c));ri={k:i for i,k in enumerate(rows)}
print('Independent 36-mode system: columns',len(columns),'rows',len(rows),'elapsed',round(time.monotonic()-t0,2),flush=True)
bi={s:i for i,s in enumerate(basis)}
for j,v in enumerate(d['basis_G']):
 accum={}
 for a,mon,c in v:
  state=(tuple(a),tuple(tuple(x) for x in mon))
  assert state in bi
  add(accum,columns[bi[state]],Q(c))
 assert not accum,j
print('PASS all 16 certificate basis vectors annihilated exactly by all 36 modes',flush=True)
# A single modular minor proves rational linear independence of the 16 vectors.
p=1000000007
B=nmod_mat(len(basis),16,p)
for j,v in enumerate(d['basis_G']):
 for a,mon,c in v:
  state=(tuple(a),tuple(tuple(x) for x in mon));x=Q(c)
  B[bi[state],j]=int(x.numerator)*pow(int(x.denominator),-1,p)%p
assert B.rank()==16
print('PASS 16 rational null vectors are linearly independent: modular rank 16',flush=True)
for p in (1000000007,1000000009):
 M=nmod_mat(len(rows),len(columns),p)
 for j,c in enumerate(columns):
  for k,x in c.items():M[ri[k],j]=int(x.numerator)*pow(int(x.denominator),-1,p)%p
 rank=M.rank();assert rank==628
 print('PASS modular rank',rank,'mod',p,'=> rational rank >=628',flush=True)
print('PASS exact dim G = 644 - 628 = 16, using 16 independently certified rational null vectors.')
print('Total kernel time',round(time.monotonic()-t0,2),'seconds')
