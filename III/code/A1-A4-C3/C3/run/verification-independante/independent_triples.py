"""Independent exact audit of the manuscript's explicit E,H,F, cocycle and characters.
Only fractions/standard library are used. No author module is imported.
"""
from fractions import Fraction as Q
from collections import Counter
from functools import lru_cache
from itertools import product,combinations
from pathlib import Path
import json,sys,hashlib
if len(sys.argv)!=2: raise SystemExit('Usage: python3 -I '+Path(__file__).name+' <Calculs-HU directory>')
ROOT=Path(sys.argv[1]).resolve()
for name in ('omega_certificate.json','data/t12_W1.json'):
 print('Input SHA256',name,hashlib.sha256((ROOT/name).read_bytes()).hexdigest(),flush=True)
Z=(0,)*16
# Roots in doubled orthonormal coordinates, independently enumerated.
roots8=[]
for i,j in combinations(range(8),2):
 for s,t in product((-2,2),repeat=2):
  x=[0]*8;x[i]=s;x[j]=t;roots8.append(tuple(x))
roots8 += [x for x in product((-1,1),repeat=8) if sum(x)%4==0]
assert len(roots8)==240
simple=[(1,-1,-1,-1,-1,-1,-1,1),(2,2,0,0,0,0,0,0),(-2,2,0,0,0,0,0,0),(0,-2,2,0,0,0,0,0),(0,0,-2,2,0,0,0,0),(0,0,0,-2,2,0,0,0),(0,0,0,0,-2,2,0,0),(0,0,0,0,0,-2,2,0)]
# Rational elimination for coordinates, deliberately no floating point inverse.
aug=[[Q(simple[j][i]) for j in range(8)]+[Q(i==j) for j in range(8)] for i in range(8)]
for k in range(8):
 r=next(i for i in range(k,8) if aug[i][k]);aug[k],aug[r]=aug[r],aug[k]
 a=aug[k][k];aug[k]=[v/a for v in aug[k]]
 for i in range(8):
  if i!=k:
   a=aug[i][k];aug[i]=[x-a*y for x,y in zip(aug[i],aug[k])]
inv=[r[8:] for r in aug]
edges=[(0,2),(1,3),(2,3),(3,4),(4,5),(5,6),(6,7)]
def dot(a,b):return sum(x*y for x,y in zip(a,b))
@lru_cache(None)
def coord(a):
 c=tuple(dot(r,a) for r in inv);assert all(x.denominator==1 for x in c);return tuple(int(x) for x in c)
@lru_cache(None)
def eps8(a,b):
 x,y=coord(a),coord(b);return (-1)**((dot(x,y)+sum(x[i]*y[j] for i,j in edges))%2)
def eps(a,b):return eps8(a[:8],b[:8])*eps8(a[8:],b[8:])
for a,b in product(roots8,repeat=2):
 assert eps8(a,b)*eps8(b,a)==(-1)**(dot(a,b)//4%2)
for a in roots8:assert eps8(a,a)==-1
print('PASS exact cocycle diagonal and skew identities on all 57,600 root pairs of E8')
# Lie vectors are (Cartan tuple in actual coordinates, root -> coefficient).
def lv(terms=(),cart=Z):return (tuple(map(Q,cart)),{a:Q(c) for c,a in terms})
def scaled(c,X):return (tuple(c*a for a in X[0]),{a:c*v for a,v in X[1].items() if c*v})
def bracket(X,Y):
 hx,rx=X;hy,ry=Y;h=[Q(0)]*16;r={}
 def add(a,c):r[a]=r.get(a,0)+c
 for a,c in ry.items():add(a,c*dot(hx,a)/2)
 for a,c in rx.items():add(a,-c*dot(hy,a)/2)
 for a,ca in rx.items():
  for b,cb in ry.items():
   s=tuple(x+y for x,y in zip(a,b));p=dot(a,b)//4
   if s==Z:
    c=ca*cb*eps(a,b);h=[x+c*y/2 for x,y in zip(h,a)]
   elif p==-1:add(s,ca*cb*eps(a,b))
 return (tuple(h),{a:c for a,c in r.items() if c})
def rt(copy,terms):
 a=[0]*16
 for j,c in terms:a[8*copy+j-1]=2*c
 return tuple(a)
def swap(X):return (X[0][8:]+X[0][:8],{a[8:]+a[:8]:c for a,c in X[1].items()})
def phi(X):
 return (tuple(-x for x in X[0]),{tuple(-x for x in a):c*Q(2)**(-sum(a[i] for i in (6,7,14,15))//2) for a,c in X[1].items()})
E1=lv([(-1,rt(0,[(1,1),(2,-1)])),(-1,rt(0,[(1,1),(2,1)])),(-1,rt(1,[(3,1),(4,1)]))])
H1=lv(cart=(2,0,0,0,0,0,0,0,0,0,1,1,0,0,0,0))
E2=lv([(1,rt(0,[(5,1),(2,-1)])),(1,rt(0,[(2,1),(5,1)])),(2,rt(0,[(7,1),(6,-1)])),(-2,rt(0,[(6,1),(8,1)])),(1,rt(1,[(3,1),(4,-1)])),(-1,rt(1,[(6,1),(2,-1)])),(-1,rt(1,[(2,1),(6,1)])),(-2,rt(1,[(7,1),(5,-1)])),(1,rt(1,[(5,1),(8,-1)]))])
H2=lv(cart=(0,0,0,0,2,0,2,2,0,0,1,-1,0,2,2,-2))
triples=[(E,H,scaled(-1,phi(E))) for E,H in [(E1,H1),(E2,H2),(swap(E1),swap(H1)),(swap(E2),swap(H2))]]
for i,(E,H,F) in enumerate(triples):
 assert bracket(E,F)==H
 assert bracket(H,E)==scaled(2,E)
 assert bracket(H,F)==scaled(-2,F)
 assert phi(F)==scaled(-1,E)
 assert phi(H)==scaled(-1,H)
 indices=tuple(dot(x,x)/2 for x in (H[0][:8],H[0][8:]))
 assert indices==[(2,1),(6,7),(1,2),(7,6)][i]
 print('PASS triple',i+1,'indices',indices)
for i,j in combinations(range(4),2):
 for x,y in product(triples[i],triples[j]):assert bracket(x,y)==lv()
print('PASS all 54 cross brackets vanish')
d=json.loads((ROOT/'omega_certificate.json').read_text())
def fromcert(v):
 h=[Q(0)]*16;r={}
 for a,mono,c in v:
  if not mono:r[tuple(a)]=Q(c)
  else:
   assert a==list(Z) and len(mono)==1 and mono[0][0]==1
   h[mono[0][1]]+=Q(c)
 return (tuple(h),r)
for t,name in zip(triples,['g1p','g1pp','g2p','g2pp']):assert t==tuple(fromcert(v) for v in d['sl2_triples'][name])
print('PASS manuscript triples exactly equal portable certificate triples')
H=[t[1][0] for t in triples]
roots=[a+(0,)*8 for a in roots8]+[(0,)*8+a for a in roots8]
def weight(a):return tuple(dot(h,a)/2 for h in H)
actual=Counter(weight(a) for a in roots);actual[(0,0,0,0)]+=16
expected=Counter()
for key,m in json.loads((ROOT/'data/t12_W1.json').read_text())['W1'].items():
 tau=tuple(map(int,key.split(',')))
 for mu in product(*(range(-a,a+1,2) for a in tau)):expected[mu]+=m
assert expected==actual
print('PASS exact W1 character under manuscript H equals all 20 T12 highest-weight types; dimension',sum(actual.values()))
norm4=set(tuple(x+y for x,y in zip(a,b)) for a,b in combinations(roots8,2) if dot(a,b)==0)
assert len(norm4)==2160
# Exhaustive direct coordinate enumeration: squared doubled norm = 16.
# Even coordinates: one +-4 or four +-2; odd: one +-3 and seven +-1.
norm4_direct=set()
for i in range(8):
 for sign in (-4,4):
  a=[0]*8;a[i]=sign;norm4_direct.add(tuple(a))
for positions in combinations(range(8),4):
 for signs in product((-2,2),repeat=4):
  a=[0]*8
  for i,sign in zip(positions,signs):a[i]=sign
  norm4_direct.add(tuple(a))
for i in range(8):
 for signs in product((-1,1),repeat=8):
  a=list(signs);a[i]*=3
  if sum(a)%4==0:norm4_direct.add(tuple(a))
assert len(norm4_direct)==2160 and norm4==norm4_direct
print('PASS exhaustive direct E8 norm-four enumeration (16+1120+1024) equals orthogonal root sums')
norm4G=[a+(0,)*8 for a in norm4]+[(0,)*8+a for a in norm4]+[a+b for a,b in product(roots8,repeat=2)]
assert len(norm4G)==61920
zero=(0,)*4
n0roots=sum(weight(a)==zero for a in roots)
n0norm4=sum(weight(a)==zero for a in norm4G)
dim0=16+136+16*n0roots+n0norm4
assert dim0==644
print('PASS direct weight count dim W2^0 =',dim0,'(oscillators 152, roots',n0roots,'norm-four',n0norm4,')')
