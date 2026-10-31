if not __debug__:
 raise RuntimeError("Run without -O or PYTHONOPTIMIZE; assertions are part of the verifier")
"""Verify full F4 with strict coverage and exact all-row membership checks.
Usage: python3 verify_f4_strict.py [--part INDEX/TOTAL].
Inputs are only COMPLETE.json manifests made by regenerate_f4.py. Missing,
truncated, duplicate or conflicting records are fatal. Each recorded pivot
must be exactly equal to its generating family word. All manifest labels,
including explicitly computed zero products, are checked against the complete
mathematical domains (9800 products and 18368 normal-product family words).
"""
import argparse, hashlib, json, os, pickle, time
from pathlib import Path
import exact_engine as E

ap=argparse.ArgumentParser();ap.add_argument('--part',default='0/1');a=ap.parse_args()
part,nparts=map(int,a.part.split('/'));assert 0<=part<nparts
root=Path('.'); t0=time.monotonic()
source_manifest=json.loads(Path('SOURCE-MANIFEST.json').read_text())
for source in source_manifest['archived_files']+source_manifest['new_verifiers']:
 path=Path(source['name'])
 assert path.stat().st_size==source['bytes'] and hashlib.sha256(path.read_bytes()).hexdigest()==source['sha256'],('source mismatch',path)
def log(*v):print(*v, f'({time.monotonic()-t0:.1f}s)',flush=True)
S=E.setup(0);B4=E.Blocks(S['T'],4);weights=sorted(B4.blk)
basehash=hashlib.sha256(Path('exact_S3_basis.pkl').read_bytes()).hexdigest()
def tuplify(a):return tuple(a)
def expected_families():
 x={('L(-1)S3',j) for j in range(956)}
 x|={('y(-1)S3',i,j) for i in range(6) for j in range(956)}
 x|={('y(-2)U2',i,j) for i in range(6) for j in range(120)}
 x|={('y(-3)U1',i,j) for i in range(6) for j in range(6)}
 x|={('m(-1)U2',name,idx,j) for name,idx,_ in S['mods'] for j in range(120)}
 assert len(x)==18368
 return x

def read_manifest(dirname):
 path=root/dirname;f=path/'COMPLETE.json'
 m=json.loads(f.read_text()); assert m['schema']==1
 assert m['S3_basis_sha256']==basehash,(dirname,'source hash')
 nz=list(map(tuplify,m['nonzero_labels']));zz=list(map(tuplify,m['zero_labels']))
 assert len(nz)==len(set(nz)) and len(zz)==len(set(zz))
 assert not set(nz)&set(zz)
 assert (len(nz),len(zz),len(nz)+len(zz))==(m['nonzero'],m['zero'],m['total'])
 blocks={b['block']:b for b in m['blocks']};assert len(blocks)==len(m['blocks'])
 assert sum(b['count'] for b in blocks.values())==len(nz)
 for i,b in blocks.items():
  p=path/f'blk_{i:03}.pkl'
  assert b['weight']==list(weights[i])
  assert p.stat().st_size==b['bytes'] and hashlib.sha256(p.read_bytes()).hexdigest()==b['sha256'],str(p)
 actual={p.name for p in path.glob('blk_*.pkl')}
 assert actual=={f'blk_{i:03}.pkl' for i in blocks},(dirname,'unexpected/missing files')
 return {'path':path,'meta':m,'blocks':blocks,'nonzero':set(nz),'zero':set(zz),'all':set(nz)|set(zz),'seen':set()}

piv=read_manifest('f4_pivots');families=read_manifest('f4_families')
assert piv['meta']['total']==6010 and not piv['zero']
assert families['all']==expected_families()
assert piv['all'] <= families['nonzero']
tests=[];alltests=set()
for name in S['hw']:
 tag=name.replace(' ','').replace('(','').replace(')','').replace(',','_')
 d=read_manifest('f4_tests_'+tag)
 expected={('u(0)s',name,j) for j in range(956)}
 if name in ('e','x'):expected|={('u(-1)w',name,j) for j in range(120)}
 assert d['all']==expected,(name,'domain')
 assert not alltests&d['all'];alltests|=d['all'];tests.append(d)
assert len(alltests)==9800
log('Domains checked',len(alltests),len(families['all']),'; fresh pivot count',len(piv['all']))

def read_block(data,i):
 if i not in data['blocks']:return []
 b=data['blocks'][i];p=data['path']/f'blk_{i:03}.pkl';out=[];seen=set()
 with p.open('rb') as f:
  for j in range(b['count']):
   try:record=pickle.load(f)
   except Exception as exc:raise ValueError(('truncated/invalid record',str(p),j)) from exc
   assert isinstance(record,tuple) and len(record)==2
   lab,v=record;assert isinstance(lab,tuple) and isinstance(v,dict) and v
   assert lab not in seen and lab in data['nonzero'];seen.add(lab)
   for k in v:assert k in B4.pos[weights[i]],('wrong block',p,lab,k)
   assert lab not in data['seen'],('duplicate across blocks',data['path'],lab)
   data['seen'].add(lab)
   out.append((lab,v))
  assert not f.read(1),('extra bytes',p)
 return out

outdir=Path('f4_verified');outdir.mkdir(exist_ok=True)
results=[]
for i,mu in sorted(enumerate(weights),key=lambda q:len(B4.blk[q[1]])):
 if i%nparts!=part:continue
 P=read_block(piv,i);F=read_block(families,i);TT=[x for d in tests for x in read_block(d,i)]
 if not P and not F and not TT:continue
 fmap=dict(F);assert len(fmap)==len(F)
 for lab,v in P:assert fmap.get(lab)==v,('pivot provenance',i,lab)
 start=time.monotonic()
 rank,bad=E.solve_block([v for _,v in P],[v for _,v in F+TT],B4.pos[mu],len(B4.blk[mu]))
 assert rank,('dependent pivot block',mu)
 assert not bad,('membership failure',mu,[(F+TT)[j][0] for j in bad[:10]])
 rec={'block':i,'weight':list(mu),'ambient':len(B4.blk[mu]),'pivots':len(P),'family_nonzero':len(F),'tests_nonzero':len(TT),'all_rows_verified':True,'bad':0,'elapsed_seconds':time.monotonic()-start}
 q=outdir/f'block_{i:03}.json';q.with_suffix('.tmp').write_text(json.dumps(rec,indent=2));q.with_suffix('.tmp').replace(q)
 results.append(rec);log('BLOCK VERIFIED',i,mu,'rank',len(P),'family',len(F),'tests',len(TT))
if nparts==1:
 for d in [piv,families]+tests:
  assert d['seen']==d['nonzero'],('missing nonzero records',d['path'])
meta={'source_manifest_sha256':hashlib.sha256(Path('SOURCE-MANIFEST.json').read_bytes()).hexdigest(),'seen_labels':{str(d['path']):[list(x) for x in sorted(d['seen'],key=str)] for d in [piv,families]+tests},'part':part,'parts':nparts,'domain_products':9800,'domain_families':18368,'pivots':6010,'input_manifest_hashes':{str(d['path']):hashlib.sha256((d['path']/'COMPLETE.json').read_bytes()).hexdigest() for d in [piv,families]+tests},'blocks':results,'elapsed_seconds':time.monotonic()-t0,'status':'PASS'}
q=outdir/f'PART_{part:03}_OF_{nparts:03}.json';q.with_suffix('.tmp').write_text(json.dumps(meta,indent=2));q.with_suffix('.tmp').replace(q)
log('EXACT VERIFICATION PART COMPLETE',part,nparts)
