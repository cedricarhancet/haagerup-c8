"""Exact projective-chart exclusions in case I, using rational polynomial equations.
Usage: python check_caseI_exclusions.py PATH_TO_GENERATION_FORTE
This does not factor eigenvalue polynomials or use numerical eigenvectors.
"""
import sys,os,json,hashlib,time
from pathlib import Path
source=Path(sys.argv[1]).resolve();os.chdir(source);sys.path.insert(0,str(source))
import sympy as sp
import griess_subalgebras as GS, primaries as PR,isotropy_test as IT
T,basis,product,gram=GS.griess_data('I');data=PR.setup('I');griess=json.loads(Path('griess_subalgebras.json').read_text())['I']['subalgebras'];n=len(basis);gram=sp.Matrix(gram)
result={'case':'I','coefficient_field':'Q','candidates':{}}
for ci in (0,1,3,4):
 xs=[sp.sympify(c) for c in griess[ci]['x']];assert all(x.is_Rational for x in xs)
 xv=sp.Matrix([0]+xs);ev=sp.Matrix([1]+[0]*(n-1));ys=sp.Matrix.vstack(ev.T*gram,xv.T*gram).nullspace()
 found=False
 for tau in [(1,13),(1,1),(1,7),(0,12)]:
  primary=data['types'][tau];d=len(primary['P']);u=sp.symbols('u0:%d'%d);v=sp.Matrix(u);z=sp.Symbol('z');P=IT.pairing(data['T'],tau,primary['P'])
  X=sum((xv[i]*primary['X1'][i] for i in range(n)),sp.zeros(d,d));Y=[sum((y[i]*primary['X1'][i] for i in range(n)),sp.zeros(d,d)) for y in ys]
  eqs=[sp.expand((v.T*P*y*v)[0]) for y in Y];Xv=X*v
  eqs +=[sp.expand(Xv[i]*u[j]-Xv[j]*u[i]) for i in range(d) for j in range(i+1,d)]
  r2=primary['X2'][0].rows
  if r2:
   for A in [primary['X2'][0],sum((xv[i]*primary['X2'][i] for i in range(n)),sp.zeros(r2,d))]:eqs +=list(A*v)
  norm=sp.expand((v.T*P*v)[0]);charts=[];allunit=True
  for k in range(d):
   sub={u[j]:0 for j in range(k)};sub[u[k]]=1;polys=[sp.expand(e.subs(sub)) for e in eqs]+[sp.expand(z*norm.subs(sub)-1)];vars=list(u[k+1:])+[z]
   gb=sp.groebner(polys,*vars,order='grevlex',domain=sp.QQ);unit=list(gb)==[1];allunit=allunit and unit
   charts.append({'first_nonzero_coordinate':k,'variables':list(map(str,vars)),'polynomials':list(map(str,polys)),'unit':unit,'groebner_basis':list(map(str,gb))})
  print('candidate',ci,'type',tau,'dim',d,'all_unit',allunit,flush=True)
  if allunit:
   result['candidates'][str(ci)]={'primary_type':list(tau),'dimension':d,'charts':charts};found=True;break
 assert found,ci
result['source_sha256']={p:hashlib.sha256(Path(p).read_bytes()).hexdigest() for p in ['primaries_I.pkl','griess_subalgebras.json','isotropy_test.py','primaries.py','griess_subalgebras.py']}
result['verified_projective_charts']=sum(len(v['charts']) for v in result['candidates'].values());Path(__file__).with_suffix('.json').write_text(json.dumps(result,indent=2));print('PASS',result['verified_projective_charts'],flush=True)
