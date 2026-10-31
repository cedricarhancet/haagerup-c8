from pathlib import Path
import os,sys,json,itertools,math
gf=Path(sys.argv[1]).resolve()
sys.path.insert(0,str(gf));os.chdir(gf)
import cases as C,lattice_voa as V
import sympy as sp
from sympy.matrices.normalforms import smith_normal_form
out={}
for case in ('I','II'):
 T=C.build(case)
 support=sorted({k[0] for nm in ('e1','e2','f1','f2') for k in T[nm] if k[0]!=V.ZERO})
 positive=sorted({k[0] for nm in ('e1','e2') for k in T[nm] if k[0]!=V.ZERO})
 zero=[a for a in V.ROOTS if C.hval(T['c1'],a)==C.hval(T['c2'],a)==0]
 M=sp.Matrix(positive);D=smith_normal_form(M,domain=sp.ZZ)
 minors=[]
 for ids in itertools.combinations(range(len(positive)),8):
  det=int(M[list(ids),:].det())
  if abs(det)==2:minors.append({'rows':list(ids),'det':det});break
 print(case,'positive',positive,'zero roots',zero,'SNF',[int(D[i,i]) for i in range(8)],'minor',minors,flush=True)
 out[case]={'H1':[str(x) for x in T['c1']],'H2':[str(x) for x in T['c2']],'raising_supports':positive,'zero_roots':zero,'smith_diagonal':[int(D[i,i]) for i in range(8)],'det2_minor':minors,'kernel_parity_verified':all(sum(a[i] for i in ((1,2,4) if case=='I' else (0,3,5)))%2==0 for a in support)}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2))
