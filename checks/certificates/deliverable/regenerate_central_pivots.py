"""Regenerate 280 exact N4' pivot words of Cartan weight (0,0).
Only source coordinates actually appearing in the required weight blocks are
used to construct each operator. This is exact linearity, with no omitted
input coordinate. This auxiliary job supplies C4 while full F4 continues.
"""
import gc,hashlib,json,pickle,time
from pathlib import Path
from gmpy2 import mpq
import exact_engine as E,Kfield as K,lattice_voa as V
start=time.monotonic()
def log(*x):print(*x,round(time.monotonic()-start,2),flush=True)
S=E.setup(0);T=S['T'];B=E.Blocks(T,4);mu=(0,0);pos=B.pos[mu];n=len(B.blk[mu]);ech=E.Ech(n)
d=pickle.loads(Path('exact_S3_basis.pkl').read_bytes());S3=[v for vs in d['S3b'].values() for v in vs];U2=[v for vs in d['U2b'].values() for v in vs]
weights={k:m for m,keys in B.blk.items() for k in keys}
def wt(v):return E.C.weight_of_key(T,next(iter(v)))
def sourcevs(vs,m):return [(j,v) for j,v in enumerate(vs) if wt(v)==m]
def op(terms,mode,vs):
 keys=set().union(*(set(v) for _,v in vs));out={}
 for key in keys:
  acc={}
  for c,u in terms:
   r=V.mode_vec(u,mode,{key:mpq(1)})
   for k,x in r.items():
    if x:acc[k]=K.add(acc.get(k,K.Z),K.scal(x,c))
  acc=E.kv_clean(acc)
  if acc:out[key]=acc
 V._MODE.clear();V._EMODE.clear();return out
piv=[];labels=[];counts={};words=0
class Finished(Exception):pass
def offer(label,v):
 global words
 words+=1
 if not v:return
 assert all(k in pos for k in v)
 if ech.add(E.modp_vec(v,pos,n)):
  piv.append(v);labels.append(label);counts[label[0]]=counts.get(label[0],0)+1
  if len(piv)==280:raise Finished
try:
 vs=sourcevs(S3,mu);O=op([(K.ONE,V.OMEGA)],0,vs)
 for j,v in vs:offer(('L(-1)S3',j),E.apply_op(O,v))
 del O;gc.collect();log('L family',len(piv))
 for i,y in enumerate(S['U1']):
  w=wt(y);opp=tuple(-a for a in w);vs=sourcevs(S3,opp);O=op([(K.ONE,y)],-1,vs)
  for j,v in vs:offer(('y(-1)S3',i,j),E.apply_op(O,v))
  del O;gc.collect()
 log('y(-1)',len(piv))
 for i,y in enumerate(S['U1']):
  opp=tuple(-a for a in wt(y));vs=sourcevs(U2,opp);O=op([(K.ONE,y)],-2,vs)
  for j,v in vs:offer(('y(-2)U2',i,j),E.apply_op(O,v))
  for j,v in sourcevs(S['U1k'],opp):offer(('y(-3)U1',i,j),E.kprod([(K.ONE,y)],-3,v))
  del O;gc.collect()
 for name,idx,terms in S['mods']:
  field=E.kv_from_terms(terms);opp=tuple(-a for a in wt(field));vs=sourcevs(U2,opp)
  if not vs:continue
  O=op(terms,-1,vs)
  for j,v in vs:offer(('m(-1)U2',name,idx,j),E.apply_op(O,v))
  del O;gc.collect()
  log('primary',name,idx,len(piv))
except Finished:pass
assert len(piv)==280
p=Path('pivot-block-00-legacy.pkl')
with p.open('wb') as f:
 for v in piv:pickle.dump(v,f,protocol=5)
metadata={'status':'PASS','scope':'280 independent N4 normal-product words in the central Cartan block; not a full F4 verification','weight':[0,0],'ambient_dimension':n,'rank':280,'word_labels':labels,'by_family':counts,'words_evaluated':words,'prime':E.P,'sqrt13':E.F.s13,'sqrt_minus3':E.F.s3,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'S3_basis_sha256':hashlib.sha256(Path('exact_S3_basis.pkl').read_bytes()).hexdigest(),'elapsed_seconds':time.monotonic()-start}
Path('pivot-block-00.json').write_text(json.dumps(metadata,indent=2));log('CENTRAL PIVOTS COMPLETE',len(piv),p.resolve())
