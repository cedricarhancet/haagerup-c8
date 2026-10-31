"""Independently re-evaluate recorded F4 zero products without operator matrices.
Run against the completed manifests of a fresh regeneration. Exact Cartan
energy excludes impossible weights; all remaining products use kprod directly.
The recorded zero labels are checked, not inferred from missing block files.
"""
from pathlib import Path
import argparse, hashlib, json, os, pickle, resource, sys, time
from fractions import Fraction
if not __debug__:
    raise RuntimeError("Python optimized mode disables required exact-check assertions")
ap=argparse.ArgumentParser();ap.add_argument('--work',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);ap.add_argument('--seconds',type=float,default=240);ap.add_argument('--existing',type=Path)
a=ap.parse_args();work=a.work.resolve();out=a.output.resolve();existing=a.existing.resolve() if a.existing else None;os.chdir(work);sys.path.insert(0,str(work))
t0=time.monotonic()
import exact_engine as E, Kfield as K, lattice_voa as V
S=E.setup(0)
with open('exact_S3_basis.pkl','rb') as f:d=pickle.load(f)
S3=[v for vs in d['S3b'].values() for v in vs];U2=[v for vs in d['U2b'].values() for v in vs]
assert len(S3)==956 and len(U2)==120
mods={(n,i):terms for n,i,terms in S['mods']}
def word(lab):
 kind=lab[0]
 if kind=='u(0)s':return S['hw'][lab[1]],0,S3[lab[2]]
 if kind=='u(-1)w':return S['hw'][lab[1]],-1,U2[lab[2]]
 if kind=='L(-1)S3':return [(K.ONE,V.OMEGA)],0,S3[lab[1]]
 if kind=='y(-1)S3':return [(K.ONE,S['U1'][lab[1]])],-1,S3[lab[2]]
 if kind=='y(-2)U2':return [(K.ONE,S['U1'][lab[1]])],-2,U2[lab[2]]
 if kind=='y(-3)U1':return [(K.ONE,S['U1'][lab[1]])],-3,S['U1k'][lab[2]]
 if kind=='m(-1)U2':return mods[lab[1],lab[2]],-1,U2[lab[3]]
 raise ValueError(lab)
def wt(v):
 weights={E.C.weight_of_key(S['T'],key) for key,c in v.items() if not K.iszero(c)}
 assert len(weights)==1
 return next(iter(weights))
paths=sorted(Path('.').glob('f4_tests_*/COMPLETE.json'))
if Path('f4_families/COMPLETE.json').is_file():paths.append(Path('f4_families/COMPLETE.json'))
manifests=[(p,json.loads(p.read_text())) for p in paths]
hashes={str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p,_ in manifests}
all_expected={(str(p),tuple(lab)) for p,m in manifests for lab in m['zero_labels']}
assert len(all_expected)==sum(m['zero'] for _,m in manifests)
total=len(all_expected);records=[];done=set();basehash=hashlib.sha256(Path('exact_S3_basis.pkl').read_bytes()).hexdigest()
if existing:
 old=json.loads(existing.read_text());assert old['S3_basis_sha256']==basehash
 assert all(hashes.get(p)==h for p,h in old['manifest_hashes'].items())
 for r in old['records']:
  key=(r['manifest'],tuple(r['label']));assert key in all_expected and key not in done
  done.add(key);records.append(r)
for p,m in manifests:
 labels=[tuple(x) for x in m['zero_labels']];assert len(labels)==len(set(labels))==m['zero']
 for lab in labels:
  if (str(p),lab) in done:continue
  if time.monotonic()-t0>a.seconds:break
  terms,n,v=word(lab);uw=wt(E.kv_from_terms(terms));vw=wt(v);mu=tuple(x+y for x,y in zip(uw,vw))
  impossible=Fraction(mu[0]**2,12)+Fraction(mu[1]**2,52)>4
  if not impossible:
   assert not E.kprod(terms,n,v),('nonzero alleged zero',str(p),lab)
   V._MODE.clear();V._EMODE.clear()
  records.append({'manifest':str(p),'label':list(lab),'method':'Cartan energy > 4' if impossible else 'direct exact kprod'})
  done.add((str(p),lab))
  if len(records)%20==0:print('verified',len(records),'elapsed',round(time.monotonic()-t0,2),'maxrss_kb',resource.getrusage(resource.RUSAGE_SELF).ru_maxrss,flush=True)
  if resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>900000:break
 if time.monotonic()-t0>a.seconds or resource.getrusage(resource.RUSAGE_SELF).ru_maxrss>900000:break
assert done <= all_expected
result={'status':'PASS' if len(records)==total else 'PARTIAL','scope':'completed manifests present at start; zero records only','total_zero_labels_in_manifests_at_start':total,'verified':len(records),'records':records,'manifest_hashes':hashes,'S3_basis_sha256':hashlib.sha256(Path('exact_S3_basis.pkl').read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-t0,'maxrss_kb':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss}
out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k not in ['records','manifest_hashes']}),flush=True)
