import numpy as np
r13=np.sqrt(13); x=(13-3*r13)/26; y=3/r13
c=lambda j: -2*y*np.cos(2*np.pi*j/13)
S=np.zeros((12,12))
S[0,:6]=[x,1-x,1,1,1,1]; S[0,6:]=y
S[1,:6]=[1-x,x,1,1,1,1]; S[1,6:]=-y
S[2,:6]=[1,1,2,-1,-1,-1]; S[3,:6]=[1,1,-1,2,-1,-1]
S[4,:6]=[1,1,-1,-1,-1,2]; S[5,:6]=[1,1,-1,-1,2,-1]
for l in range(1,7):
    S[5+l,0]=y; S[5+l,1]=-y
    for m in range(1,7): S[5+l,5+m]=c(l*m)
S/=3
print("symmetric:",np.allclose(S,S.T)," unitary:",np.allclose(S@S.T,np.eye(12)))
names=['0','b','a','c0','c1','c2']+['d%d'%l for l in range(1,7)]
N=lambda i,j,k: sum(S[i,m]*S[j,m]*S[k,m]/S[0,m] for m in range(12))
print("b x b =",{names[k]:round(N(1,1,k),6) for k in range(12)})
print("min/max over all N_ijk:",min(N(i,j,k) for i in range(12) for j in range(12) for k in range(12)))
for i in (2,3,4,5,6):
    print(names[i],"x",names[i],"=",{names[k]:int(round(N(i,i,k))) for k in range(12) if abs(N(i,i,k))>1e-9})
print("dims:",[round(S[0,k]/S[0,0],4) for k in range(12)])
