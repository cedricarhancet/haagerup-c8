import sympy as S, itertools as I, math, json
from pathlib import Path
K=S.Matrix([[2,-1,0,0,0,-1],[-1,2,0,0,0,0],[0,0,2,-1,0,-1],[0,0,-1,2,0,0],[0,0,0,0,2,-1],[-1,0,-1,0,-1,4]])
F=S.Matrix([[2,-1,-1,1,-1,-1],[-1,2,1,0,0,1],[-1,1,2,0,0,0],[1,0,0,2,-1,-1],[-1,0,0,-1,2,0],[-1,1,0,-1,0,8]])
N=S.Matrix([[4,-1,-2,-2,-1,-1],[-1,4,-1,-1,2,2],[-2,-1,4,1,-1,-1],[-2,-1,1,4,1,1],[-1,2,-1,1,4,1],[-1,2,-1,1,1,4]])
def enumerate_to(G,top):
 inv=G.inv(); bounds=[math.isqrt(int(S.floor(top*inv[i,i]))) for i in range(6)]
 g=[[int(G[i,j]) for j in range(6)] for i in range(6)]
 shells={i:[] for i in range(top+1)}
 for x in I.product(*[range(-a,a+1) for a in bounds]):
  q=sum(g[i][i]*x[i]*x[i] for i in range(6))+2*sum(g[i][j]*x[i]*x[j] for i in range(6) for j in range(i))
  if q<=top:shells[q].append(x)
 return shells
def aut_count(G,shells):
 g=[[int(G[i,j]) for j in range(6)] for i in range(6)]
 candidates=[shells[g[i][i]] for i in range(6)]
 gc=[{x:tuple(sum(g[i][j]*x[j] for j in range(6)) for i in range(6)) for x in c} for c in candidates]
 count=0
 def rec(chosen):
  nonlocal count
  k=len(chosen)
  if k==6:
   count+=1;return
  for x in candidates[k]:
   gx=gc[k][x]
   if all(sum(y[i]*gx[i] for i in range(6))==g[j][k] for j,y in enumerate(chosen)):
    rec(chosen+[x])
 rec([]);return count
out={}
for name,G in [('K',K),('Kflat',F),('N',N)]:
 sh=enumerate_to(G,8 if name=='Kflat' else 6)
 out[name]={'det':int(G.det()),'sylvester':[int(G[:i,:i].det()) for i in range(1,7)],'shells':{q:len(v) for q,v in sh.items() if len(v)},'root_rank':S.Matrix(sh[2]).rank()}
 if name!='N':out[name]['O']=aut_count(G,sh)
 if name!='N':out[name]['qgen']=str(G.inv()[0,0]/2 if name=='K' else G.inv()[1,1]/2)
 print(name,out[name],flush=True)
ker=[]
for v in I.product(range(3),repeat=6):
 if not any(v):continue
 if all(sum(int(N[i,j])*v[j] for j in range(6))%3==0 for i in range(6)) and sum(int(N[i,j])*v[i]*v[j] for i in range(6) for j in range(6))%18==0:
  if v<tuple((-x)%3 for x in v):ker.append(v)
counts={str(v):0 for v in ker}
g=[[int(N[i,j]) for j in range(6)] for i in range(6)]
for z in I.product(range(-3,4),repeat=6):
 q=sum(g[i][j]*z[i]*z[j] for i in range(6) for j in range(6))
 if q!=18:continue
 zz=tuple(x%3 for x in z)
 for v in ker:
  if zz==(0,)*6 or zz==v or zz==tuple((-x)%3 for x in v):counts[str(v)]+=1
out['index3_lines']=ker;out['overlattice_roots']=counts
B3=lambda x:x**3-S.Rational(3,2)*x*x+x/2
bern=39**2*sum(S.jacobi_symbol(a,39)*B3(S.Rational(a,39)) for a in range(1,40))
out['B3chi']=str(bern);out['mass_from_formula']=str(bern/S.Integer(138240));out['mass_from_orders']=str(1/S.Integer(out['K']['O'])+1/S.Integer(out['Kflat']['O']))
mins={}
for m,n in I.product(range(-12,13),range(-25,26)):
 if (m-n)%2:continue
 e=S.Rational(m*m,12)+S.Rational(n*n,52);p=(26*m+21*n)%39
 if p not in mins or e<mins[p]['energy']:mins[p]={'energy':e,'momenta':[(m,n)]}
 elif e==mins[p]['energy']:mins[p]['momenta'].append((m,n))
out['energy_minima']={str(p):{'energy':str(d['energy']),'momenta':d['momenta']} for p,d in sorted(mins.items())}
out['over_one']={str(p):{'energy':str(d['energy']),'momenta':d['momenta']} for p,d in sorted(mins.items()) if d['energy']>=1}
out['pointed_isometries']=[u for u in range(39) if math.gcd(u,39)==1 and (u*u-1)%39==0]
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2))
print(json.dumps({k:v for k,v in out.items() if k!='energy_minima'},indent=2),flush=True)
