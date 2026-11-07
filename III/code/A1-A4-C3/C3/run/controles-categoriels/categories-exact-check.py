from fractions import Fraction as Q
from pathlib import Path

# Exact arithmetic in Q[z]/(1+z+...+z^12), independently of supplied scripts.
class C:
    def __init__(self,a=0):
        if isinstance(a,C):self.a=a.a;return
        if isinstance(a,(int,Q)):a=[a]
        a=list(map(Q,a));a += [Q(0)]*max(0,12-len(a))
        for k in range(len(a)-1,11,-1):
            t=a[k]
            for j in range(13):a[k-12+j]-=t
        self.a=tuple(a[:12])
    def __add__(self,o):o=C(o);return C([a+b for a,b in zip(self.a,o.a)])
    __radd__=__add__
    def __neg__(self):return C([-a for a in self.a])
    def __sub__(self,o):return self+-C(o)
    def __rsub__(self,o):return C(o)+-self
    def __mul__(self,o):
        o=C(o);a=[Q(0)]*23
        for i,x in enumerate(self.a):
            for j,y in enumerate(o.a):a[i+j]+=x*y
        return C(a)
    __rmul__=__mul__
    def __truediv__(self,o):return C([x/Q(o) for x in self.a])
    def __pow__(self,k):
        r=C(1)
        for _ in range(k):r=r*self
        return r
    def __eq__(self,o):return self.a==C(o).a
z=C([0,1]); powers=[z**i for i in range(13)]
s=sum((1 if pow(i,6,13)==1 else -1)*powers[i] for i in range(1,13))
assert s*s==13
x=(13-3*s)/26;y=3*s/13
M=[[C(0) for _ in range(12)] for _ in range(12)]
M[0]=[x,1-x]+[C(1)]*4+[y]*6
M[1]=[1-x,x]+[C(1)]*4+[-y]*6
p=[0,1,2,3,5,4]
for i in range(2,6):
    M[i][0]=M[i][1]=C(1)
    for j in range(2,6):M[i][j]=C(2 if p[i]==j else -1)
for i in range(6,12):
    M[i][0]=y;M[i][1]=-y
    for j in range(6,12):
        v=((i-5)*(j-5))%13
        M[i][j]=-y*(powers[v]+powers[(-v)%13])
for i in range(12):
    for j in range(12):assert sum(M[i][k]*M[k][j] for k in range(12))==9*(i==j)
lam=[(11+3*s)/2,(11-3*s)/2]+[C(1)]*4+[C(-1)]*6
Nb=[[sum(M[i][k]*lam[k]*M[k][j] for k in range(12))/9 for j in range(12)] for i in range(12)]
for i in range(12):
    if i==0: want=[int(j==1) for j in range(12)]
    elif i==1: want=[1]*12
    elif i<6:want=[int(j!=0)+int(j==i) for j in range(12)]
    else:want=[int(j not in (0,i)) for j in range(12)]
    assert Nb[i]==[C(t) for t in want],i
assert sum(powers[6*l*l%13] for l in range(1,7))==(-1-s)/2
assert sum(powers[12*l*l%13] for l in range(1,7))==(-1+s)/2
out='PASS: exact S^2=I and all 144 coefficients of N_b match the three fusion identities added in v07.\nPASS: the two quadratic Gauss sums at m=1 equal (-1-sqrt(13))/2 and (-1+sqrt(13))/2.\n'
print(out,end='')
