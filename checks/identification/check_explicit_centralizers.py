"""Direct exact checks for the analytic centralizer proof, without solving bracket equations.
Usage: python check_explicit_centralizers.py --sources PATH [--membership]
PATH is the computational source directory containing generation-forte and cas-I.
The explicit lattice lifts are constructed from the formulas in the paper, without j3/j6 matrices.
"""
from pathlib import Path
import argparse,sys,os,json,pickle,itertools
p=argparse.ArgumentParser();p.add_argument('--sources',required=True);p.add_argument('--membership',action='store_true');args=p.parse_args()
sources=Path(args.sources).resolve();gf=sources/'generation-forte';sys.path.insert(0,str(gf));os.chdir(gf)
from gmpy2 import mpq
import sympy as sp
import lattice_voa as V,cases as C
out={}
def lift(perm,cs):
 cs=list(map(mpq,cs)); simples=[tuple(int(i==j) for i in range(8)) for j in range(8)]
 def gv(a):
  r=[0]*8
  for i in range(8):r[perm[i]]=a[i]
  return tuple(r)
 delta=[[V.eps(a,b)*V.eps(gv(a),gv(b)) for b in simples] for a in simples]
 assert all(delta[i][i]==1 for i in range(8))
 assert all(V.ip(gv(a),gv(b))==V.ip(a,b) for a in simples for b in simples)
 def coeff(a):
  c=mpq(1)
  for i in range(8):
   c*=cs[i]**a[i]
   for j in range(i+1,8):
    if delta[i][j]==-1 and (a[i]*a[j])%2:c=-c
  return c
 def apply(v,ring='Q'):
  result={}
  for (lam,osc),a in v.items():
   key=(gv(lam),tuple(sorted((n,perm[i]) for n,i in osc)))
   val=coeff(lam)
   if ring=='K':
    import Kfield as K
    result[key]=K.add(result.get(key,K.Z),K.scal(val,a))
   else:
    if ring=='S':val=sp.Rational(int(val.numerator),int(val.denominator))
    result[key]=result.get(key,0)+val*a
  return {k:a for k,a in result.items() if a} if ring!='K' else {k:a for k,a in result.items() if not K.iszero(a)}
 return apply
identity=list(range(8));gE=[5,1,4,3,2,0,6,7];gA=[0,1,2,3,4,5,7,6];gEA=[5,1,4,3,2,0,7,6]
B1=[{k:mpq(1)} for k in V.basis(1)]
for case in ('I','II'):
 T=C.build(case);cs=['4/5','1','-5/3','-1','-3/5','5/4','1','1'] if case=='II' else ['4/5','-1','1/3','1','3','5/4','1','1']
 rho_cs=['-1','1','1','-1','1','-1','1','1'] if case=='II' else ['1','-1','-1','1','-1','1','1','1']
 ze,za,rho=lift(gE,cs),lift(gA,['1']*8),lift(identity,rho_cs)
 maps={'zE':ze,'zA':za,'rho':rho}
 checks={}
 checks['six_currents_fixed']={name:all(z(T[n])==T[n] for n in ('e1','f1','H1','e2','f2','H2')) for name,z in maps.items()}
 checks['involutions']={name:all(z(z(b))==b for b in B1) for name,z in maps.items()}
 checks['commuting']={a+' '+b:all(maps[a](maps[b](x))==maps[b](maps[a](x)) for x in B1) for a,b in itertools.combinations(maps,2)}
 group=[]
 for bits in itertools.product((0,1),repeat=3):
  images=[]
  for x in B1:
   for bit,z in zip(bits,maps.values()):
    if bit:x=z(x)
   images.append(tuple(sorted((k,str(v)) for k,v in x.items())))
  group.append(tuple(images))
 checks['group_order']=len(set(group));assert checks['group_order']==8
 assert all(all(x.values()) for x in [checks['six_currents_fixed'],checks['involutions'],checks['commuting']])
 out[case]=checks
 print(case,checks,flush=True)
 if case=='I':
  full={n:{k:sp.sympify(c) for k,c in d.items()} for n,d in pickle.load(open(sources/'cas-I/r1_riesen_transported.pkl','rb')).items()}
  assert all(not a.has(sp.Float) for vec in full.values() for a in vec.values())
  def eq(a,b):return all(sp.simplify(sp.expand(a.get(k,0)-b.get(k,0)))==0 for k in set(a)|set(b))
  def galois(v,s13,si):
   r,s,i=sp.symbols('r s i');result={}
   for k,a in v.items():
    a=sp.expand(a).xreplace({sp.sqrt(39):r*s,sp.sqrt(13):r,sp.sqrt(3):s,sp.I:i})
    result[k]=sp.expand(a.xreplace({r:s13*sp.sqrt(13),s:sp.sqrt(3),i:si*sp.I}))
   return result
  vprime={k:sp.expand(sp.I*a) for k,a in full['v'].items()}
  import Kfield as K
  for vec in [full['e'],full['f'],vprime]:
   for a in vec.values():K.from_sympy(a)
  checks['normalized_generators_in_K']=True
  z1=lambda v:ze(za(v,'S'),'S');tau=lambda v:z1(galois(v,1,-1));sig=lambda v:ze(galois(v,-1,1),'S')
  checks['descent']={}
  for name,psi in [('sigma',sig),('tau',tau)]:
   for n in ('e','f'):checks['descent'][name+' fixes '+n]=eq(psi(full[n]),full[n])
   checks['descent'][name+' vprime = -vprime']=eq(psi(vprime),{k:-a for k,a in vprime.items()})
  print('I descent',checks['descent'],flush=True)
  assert all(checks['descent'].values())
 if case=='II':
  import exact_engine as E,Kfield as K
  bb=pickle.load(open('griess_B_II_basis.pkl','rb'))
  sb=[dict(b) for b in bb];V.vadd(sb[4],bb[1],mpq(-1,12));V.vadd(sb[4],bb[2],mpq(5,24));V.vadd(sb[4],bb[3],mpq(1,16));V.vadd(sb[4],bb[4],mpq(-2));sb[4]=V.clean(sb[4])
  checks['B_action']={'zE=s':all(ze(b)==t for b,t in zip(bb,sb)),'zA=s':all(za(b)==t for b,t in zip(bb,sb)),'rho=id':all(rho(b)==b for b in bb)}
  assert all(checks['B_action'].values());print('II B',checks['B_action'],flush=True)
  if args.membership:
   S=E.setup(0);B2=E.Blocks(S['T'],2);vs={}
   for v in S['U2span']:vs.setdefault(B2.wt(v),[]).append(v)
   basis={mu:[v[j] for j in E.select_pivots(v,B2.pos[mu],len(B2.blk[mu]))] for mu,v in vs.items()}
   sig=lambda v:{k:(a[0],-a[1],a[2],-a[3]) for k,a in v.items()}
   tau=lambda v:{k:(a[0],a[1],-a[2],-a[3]) for k,a in v.items()}
   psis={'sigma':lambda v:ze(sig(v),'K'),'tau':lambda v:ze(za(tau(v),'K'),'K')}
   checks['descent_membership']={}
   for name,psi in psis.items():
    imgs={}
    for v in S['U2span']:
     w=psi(v);imgs.setdefault(B2.wt(w),[]).append(w)
    ok=True
    for mu,ws in imgs.items():
     good,bad=E.solve_block(basis[mu],ws,B2.pos[mu],len(B2.blk[mu]));ok=ok and good and not bad
    checks['descent_membership'][name]=bool(ok);assert ok;print('II descent',name,ok,flush=True)
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2));print('PASS',flush=True)
