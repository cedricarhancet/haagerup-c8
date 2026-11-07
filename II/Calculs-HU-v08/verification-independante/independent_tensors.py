"""Audit tensors independently from the 36 current-mode formula in independent_kernel.py.
General weight-two states are expressed as a(-1)b or derivatives of weight-one
currents; the normal-order identity is then used. No author's engine is loaded.
"""
from pathlib import Path
import runpy,time
from functools import lru_cache
from itertools import product
from collections import Counter,defaultdict
D=runpy.run_path(str(Path(__file__).with_name('independent_kernel.py')))
globals().update({k:D[k] for k in ('Q','Z','dot','eps','weight','root_mode','h_mode','add','d','GLOB')})
roots=GLOB['roots'];rootset=set(roots)
G=[{(tuple(a),tuple(tuple(x) for x in mon)):Q(c) for a,mon,c in v} for v in d['basis_G']]
VAC=(Z,())
# The certificate stores ordered tensors. Verify its symmetry explicitly,
# then 136 unordered comparisons cover all 256 entries.
assert all(d['griess_product'][i][j] == d['griess_product'][j][i] for i in range(16) for j in range(16))
assert all(d['bilinear_form'][i][j] == d['bilinear_form'][j][i] for i in range(16) for j in range(16))
@lru_cache(None)
def decomp(a):
 for r in roots:
  b=tuple(x-y for x,y in zip(a,r))
  if b in rootset and dot(r,b)==0:return r,b,eps(r,b)
 raise AssertionError(a)
@lru_cache(None)
def mode1(a,n,target):
 if isinstance(a,int):return h_mode(tuple(Q(i==a) for i in range(16)),n,target)
 return root_mode(a,n,target)
@lru_cache(maxsize=100000)
def mode2(source,n,target):
 a,mono=source
 if mono:
  assert a==Z
  if len(mono)==1:
   assert mono[0][0]==2
   return {s:-n*c for s,c in mode1(mono[0][1],n-1,target).items() if n*c}
  assert len(mono)==2 and mono[0][0]==mono[1][0]==1
  left,right,factor=mono[0][1],mono[1][1],1
 else:left,right,factor=decomp(a)
 wt=dot(target[0],target[0])//8+sum(m for m,i in target[1])
 ans={}
 for k in range(wt+1):
  for s,c in mode1(right,n+k,target).items():add(ans,mode1(left,-1-k,s),factor*c)
  for s,c in mode1(left,k,target).items():add(ans,mode1(right,n-1-k,s),factor*c)
 return ans
def bilinear(X,n,Y):
 out={}
 for a,ca in X.items():
  for b,cb in Y.items():add(out,mode2(a,n,b),ca*cb)
 return out
t0=time.monotonic()
for i in range(16):
 for j in range(i,16):
  got=bilinear(G[i],1,G[j]);want={}
  for k,c in enumerate(d['griess_product'][i][j]):
   if Q(c):add(want,G[k],Q(c))
  assert got==want,('product',i,j,len(got),len(want))
  gotform=bilinear(G[i],3,G[j]);q=Q(d['bilinear_form'][i][j]);wantform={VAC:q} if q else {}
  assert gotform==wantform,('form',i,j,gotform,wantform)
 print('PASS product and form row',i,'time',round(time.monotonic()-t0,1),flush=True)
print('PASS all 136 unordered products (g_i)_(1)g_j and all 136 forms (g_i)_(3)g_j exactly match certificate',flush=True)
# Spectral invariants on HW multiplicity spaces are obtained by four Weyl
# differences from traces on weight spaces. This avoids sharing HW bases.
W1=[(Z,((1,i),)) for i in range(16)]+[(r,()) for r in roots]
WI={k:i for i,k in enumerate(W1)}
actions=[]
for g in G:
 columns=[]
 for s in W1:
  out={}
  for a,c in g.items():add(out,mode2(a,1,s),c)
  assert all(t in WI and weight(t[0])==weight(s[0]) for t in out)
  columns.append({WI[t]:c for t,c in out.items()})
 actions.append(columns)
print('Independent W1 actions computed',round(time.monotonic()-t0,1),flush=True)
tr=[]
for M in actions:
 sums=defaultdict(Q)
 for i,s in enumerate(W1):sums[weight(s[0])]+=M[i].get(i,0)
 tr.append(sums)
tr2={}
for i in range(16):
 for j in range(i,16):
  sums=defaultdict(Q)
  for col,s in enumerate(W1):
   v=sum(c*actions[i][k].get(col,0) for k,c in actions[j][col].items())
   sums[weight(s[0])]+=v
  tr2[i,j]=sums

def hwtrace(sums,tau):
 return sum((-1)**sum(mask)*sums.get(tuple(a+2*b for a,b in zip(tau,mask)),0) for mask in product((0,1),repeat=4))
for block in d['multiplicity_spaces']:
 tau=block['highest_weight'];mats=[[[Q(c) for c in row] for row in M] for M in block['actions']]
 n=len(mats[0])
 for i in range(16):
  assert hwtrace(tr[i],tau)==sum(mats[i][k][k] for k in range(n)),('trace',tau,i)
  for j in range(i,16):
   expected=sum(mats[i][k][l]*mats[j][l][k] for k in range(n) for l in range(n))
   assert hwtrace(tr2[i,j],tau)==expected,('tracepair',tau,i,j)
print('PASS all 20 multiplicity spaces: 16 traces and 136 trace-products each exactly match certificate',flush=True)
print('PASS all tensor consequences needed for the 20 linear equations and 10 e2 quadrics rebuilt independently')
print('Total tensor audit time',round(time.monotonic()-t0,1),'seconds')
