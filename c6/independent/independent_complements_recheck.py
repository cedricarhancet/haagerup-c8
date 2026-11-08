import sympy as S, json
from pathlib import Path
G=2*S.eye(8)
for i,j in [(0,2),(1,3),(2,3),(3,4),(4,5),(6,7)]:G[i,j]=G[j,i]=-1
g=S.Matrix([4,3,5,6,4,2,2,1])/3
K=S.Matrix([[2,-1,0,0,0,-1],[-1,2,0,0,0,0],[0,0,2,-1,0,-1],[0,0,-1,2,0,0],[0,0,0,0,2,-1],[-1,0,-1,0,-1,4]])
BI=S.Matrix([[1,1,-2,1,1,0],[0,0,0,-3,0,0],[2,-1,-1,-1,-1,3],[0,0,0,-3,0,3],[-2,1,1,-2,1,3],[-1,2,2,-1,-1,0],[-1,-1,-1,2,-1,0],[1,-2,1,1,1,-3]])/3
BII=S.Matrix([[1,1,-1,2,-2,1],[0,0,0,3,0,-3],[-1,2,-2,4,-1,-1],[0,3,0,3,0,-6],[1,1,2,2,1,-5],[-1,2,1,1,2,-1],[-1,-1,1,-2,-1,2],[1,-2,-1,-1,1,1]])/3
hs=[([0,0,0,1,0,0,0,0],[8,10,14,19,14,8,2,2]),([0,-1,1,0,1,0,0,0],[4,5,7,10,7,4,2,2])]
out={}
for nm,B,(a,b) in zip(['I','II'],[BI,BII],hs):
 h1,h2=S.Matrix(a),S.Matrix(b);L=h1.row_join((h1+h2)/2)
 assert B.T*G*B==K
 assert B.T*G*h1==S.zeros(6,1) and B.T*G*h2==S.zeros(6,1)
 assert all(any(all(t.q==1 for t in B[:,j]-k*g) for k in range(3)) for j in range(6))
 assert (L.T*G*L).det()==39
 out[nm]={'gram_L':str(L.T*G*L),'gram_check':True,'orthogonality':True,'all_basis_in_E8plus':True}
assert G.det()==9 and (g.T*G*g)[0]%2==0
out['E8plus_integrality']=all(t.q==1 for t in G*g)
out['E8plus_norm_glue']=str((g.T*G*g)[0])
print(json.dumps(out,indent=2));Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2))
