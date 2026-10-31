"""Strict exact regeneration of F4 data; never consumes incomplete legacy streams.
Run from a fresh work directory copied from generation-forte. Each phase has
its own directory and atomic completion manifest; zero products have explicit
labels. Coordinates use the historical exact arithmetic engine, unchanged.
"""
import argparse, gc, hashlib, json, os, pickle, sys, time
from pathlib import Path
import exact_engine as E, Kfield as K, lattice_voa as V

ap=argparse.ArgumentParser();ap.add_argument('phase',choices=['pivots','tests']);ap.add_argument('--names');args=ap.parse_args()
t0=time.monotonic()
def log(*a): print(*a, f'({time.monotonic()-t0:.1f}s)',flush=True)
S=E.setup(0);T=S['T']
with open('exact_S3_basis.pkl','rb') as f:d=pickle.load(f)
S3=[v for vs in d['S3b'].values() for v in vs];U2=[v for vs in d['U2b'].values() for v in vs]
assert len(S3)==956 and len(U2)==120
B4=E.Blocks(T,4);bidx={mu:i for i,mu in enumerate(sorted(B4.blk))}
key_weights={k:mu for mu,keys in B4.blk.items() for k in keys}
assert len(key_weights)==len(B4.keys)
def weight(v):
    ws={key_weights[k] for k,a in v.items() if not K.iszero(a)}
    if len(ws)!=1:raise ValueError(('nonhomogeneous or zero',ws))
    return next(iter(ws))
class Store:
 def __init__(self,path):
  self.path=Path(path);self.path.mkdir(exist_ok=False);self.handles={};self.records={};self.labels=[];self.zero=[]
 def add(self,label,v):
  label=tuple(label)
  if label in self.records:raise ValueError(('duplicate label',label))
  if not v:self.records[label]=None;self.zero.append(label);return
  mu=weight(v);i=bidx[mu]
  if i not in self.handles:self.handles[i]=open(self.path/f'blk_{i:03}.pkl','ab')
  pickle.dump((label,v),self.handles[i],protocol=5)
  self.records[label]=i;self.labels.append(label)
 def finish(self,extra):
  for f in self.handles.values():f.close()
  blocks=[]
  for i in sorted(self.handles):
   p=self.path/f'blk_{i:03}.pkl';blocks.append({'block':i,'weight':list(sorted(B4.blk)[i]),'count':sum(x==i for x in self.records.values()),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
  data={'schema':1,'phase':args.phase,'elapsed_seconds':time.monotonic()-t0,'blocks':blocks,'nonzero_labels':self.labels,'zero_labels':self.zero,'total':len(self.records),'nonzero':len(self.labels),'zero':len(self.zero),'S3_basis_sha256':hashlib.sha256(Path('exact_S3_basis.pkl').read_bytes()).hexdigest(),**extra}
  tmp=self.path/'manifest.tmp';tmp.write_text(json.dumps(data,indent=1));tmp.replace(self.path/'COMPLETE.json');return data

def families():
 O=E.build_op([(K.ONE,V.OMEGA)],0,3)
 for j,s in enumerate(S3):yield ('L(-1)S3',j),E.apply_op(O,s)
 del O;gc.collect();log('L(-1) family')
 for i,y in enumerate(S['U1']):
  O=E.build_op([(K.ONE,y)],-1,3)
  for j,s in enumerate(S3):yield ('y(-1)S3',i,j),E.apply_op(O,s)
  del O;gc.collect()
 log('y(-1) family')
 for i,y in enumerate(S['U1']):
  O=E.build_op([(K.ONE,y)],-2,2)
  for j,w in enumerate(U2):yield ('y(-2)U2',i,j),E.apply_op(O,w)
  for j,y2 in enumerate(S['U1k']):yield ('y(-3)U1',i,j),E.kprod([(K.ONE,y)],-3,y2)
  del O;gc.collect()
 log('y(-2), y(-3) families')
 for name,idx,terms in S['mods']:
  O=E.build_op(terms,-1,2)
  for j,w in enumerate(U2):yield ('m(-1)U2',name,idx,j),E.apply_op(O,w)
  del O;gc.collect()
 log('primary families')

if args.phase=='pivots':
 store=Store('f4_pivots');allstore=Store('f4_families')
 ech={mu:E.Ech(len(ks)) for mu,ks in B4.blk.items()};counts={}
 for label,v in families():
  allstore.add(label,v)
  if not v:continue
  mu=weight(v)
  if ech[mu].add(E.modp_vec(v,B4.pos[mu],len(B4.blk[mu]))):store.add(label,v);counts[label[0]]=counts.get(label[0],0)+1
 result=store.finish({'independent_pivots':len(store.records),'by_family':counts,'prime':E.P,'sqrt13':E.F.s13,'sqrt_minus3':E.F.s3})
 assert result['total']==6010,result['total']
 assert len(allstore.records)==18368,len(allstore.records)
 allstore.finish({'families_exhaustive':True})
 log('PIVOT AND FAMILY GENERATION COMPLETE',result['total'],counts)
else:
 names=args.names.split(';') if args.names else list(S['hw'])
 if len(set(names))!=len(names) or any(n not in S['hw'] for n in names):raise ValueError(names)
 for name in names:
  tag=name.replace(' ','').replace('(','').replace(')','').replace(',','_')
  store=Store('f4_tests_'+tag)
  terms=S['hw'][name];single=[(c,{k:x}) for c,v in terms for k,x in v.items()]
  acc=[{} for _ in S3]
  for start in range(0,len(single),6):
   O=E.build_op(single[start:start+6],0,3)
   for j,s in enumerate(S3):E.kv_add(acc[j],E.apply_op(O,s))
   del O;gc.collect()
  for j,v in enumerate(acc):store.add(('u(0)s',name,j),E.kv_clean(v))
  del acc;gc.collect()
  if name in ('e','x'):
   O=E.build_op(S['hw'][name],-1,2)
   for j,w in enumerate(U2):store.add(('u(-1)w',name,j),E.apply_op(O,w))
   del O;gc.collect()
  result=store.finish({'name':name,'names_available':list(S['hw'])})
  assert result['total']==956+(120 if name in ('e','x') else 0)
  log('TEST GENERATOR COMPLETE',name,result['nonzero'],result['zero'])
