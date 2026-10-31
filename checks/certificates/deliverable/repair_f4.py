"""Recover explicitly missing labelled words and rebuild every pivot by provenance.
Source streams are immutable snapshots. Their hashes and every existing record
are checked; no absent record is silently treated as zero. Missing nonzero words
are recomputed exactly with operators restricted to all source coordinates used
by those words. Final SafeStore manifests require complete rereading.
"""
if not __debug__:raise RuntimeError('Run without -O or PYTHONOPTIMIZE')
import argparse,gc,hashlib,json,pickle,time
from pathlib import Path
from gmpy2 import mpq
import exact_engine as E,Kfield as K,lattice_voa as V
from safe_store import SafeStore
p=argparse.ArgumentParser();p.add_argument('phase',choices=['families','tests']);p.add_argument('--names');a=p.parse_args()
t0=time.monotonic()
def log(*x):print(*x,round(time.monotonic()-t0,1),flush=True)
S=E.setup(0);B=E.Blocks(S['T'],4);weights=sorted(B.blk)
basehash=hashlib.sha256(Path('exact_S3_basis.pkl').read_bytes()).hexdigest();d=pickle.loads(Path('exact_S3_basis.pkl').read_bytes())
S3=[v for vs in d['S3b'].values() for v in vs];U2=[v for vs in d['U2b'].values() for v in vs]
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
def group(lab):return lab[:-1]
def restricted_op(terms,n,vs):
 keys=set().union(*(set(v) for v in vs));O={}
 for key in keys:
  out={}
  for c,u in terms:
   for kk,x in V.mode_vec(u,n,{key:mpq(1)}).items():
    if x:out[kk]=K.add(out.get(kk,K.Z),K.scal(x,c))
  out=E.kv_clean(out)
  if out:O[key]=out
 V._MODE.clear();V._EMODE.clear();return O

def recover(name):
 source=Path('input-streams')/name;m=json.loads((source/'COMPLETE.json').read_text())
 assert m['S3_basis_sha256']==basehash
 nz=list(map(tuple,m['nonzero_labels']));zero=list(map(tuple,m['zero_labels']))
 assert len(nz)==len(set(nz))==m['nonzero'] and len(zero)==len(set(zero))==m['zero']
 assert not set(nz)&set(zero) and len(nz)+len(zero)==m['total']
 store=SafeStore(name,B,basehash);expected=set(nz);present=set();source_blocks=[]
 for rec in m['blocks']:
  i=rec['block'];q=source/f'blk_{i:03}.pkl';assert rec['weight']==list(weights[i])
  assert q.stat().st_size==rec['bytes'] and hashlib.sha256(q.read_bytes()).hexdigest()==rec['sha256'],('source snapshot changed',q)
  count=0;truncation=None;last_complete=0
  with q.open('rb') as f:
   while f.tell()<q.stat().st_size:
    at=f.tell()
    try:lab,v=pickle.load(f)
    except (EOFError,pickle.UnpicklingError) as exc:
     truncation={'exception':type(exc).__name__,'offset':at,'tail_bytes':q.stat().st_size-at}
     break
    lab=tuple(lab);assert lab in expected and lab not in present and v
    assert all(k in B.pos[weights[i]] for k in v)
    present.add(lab);count+=1;last_complete=f.tell();store.add(lab,v)
  assert count<=rec['count'],('surplus records',q)
  source_blocks.append({'block':i,'expected':rec['count'],'present':count,'sha256':rec['sha256'],'truncation':truncation,'last_complete_offset':last_complete})
 missing=[lab for lab in nz if lab not in present];groups={}
 for lab in missing:groups.setdefault(group(lab),[]).append(lab)
 log(name,'present',len(present),'missing',len(missing),'groups',len(groups))
 for key,labs in groups.items():
  terms,n,_=word(labs[0]);vs=[word(lab)[2] for lab in labs];O=restricted_op(terms,n,vs)
  for lab,v in zip(labs,vs):
   z=E.apply_op(O,v);assert z,('previous nonzero label now zero',lab);store.add(lab,z)
  del O;gc.collect();log('recovered',key,len(labs))
 for lab in zero:store.add(lab,{})
 assert set(store.records)==expected|set(zero)
 meta=store.finish({'producer':'repair_f4.py','input_snapshot_manifest_sha256':hashlib.sha256((source/'COMPLETE.json').read_bytes()).hexdigest(),'preserved_nonzero_records':len(present),'recomputed_nonzero_records':len(missing),'recomputed_labels':missing,'source_blocks':source_blocks})
 log('RECOVERY COMPLETE',name,meta['total'])
 return store,meta

if a.phase=='families':
 store,meta=recover('f4_families');assert meta['total']==18368
 old=Path('input-streams/f4_pivots/COMPLETE.json');pm=json.loads(old.read_text());plabs=list(map(tuple,pm['nonzero_labels']))
 assert len(plabs)==len(set(plabs))==6010 and not pm['zero_labels']
 byblock={}
 for lab in plabs:
  assert lab in store.records and store.records[lab] is not None
  byblock.setdefault(store.records[lab],[]).append(lab)
 piv=SafeStore('f4_pivots',B,basehash);ranks={}
 for i,labs in sorted(byblock.items()):
  wanted=set(labs);data={};q=Path('f4_families')/f'blk_{i:03}.pkl'
  with q.open('rb') as f:
   while f.tell()<q.stat().st_size:
    lab,v=pickle.load(f)
    if lab in wanted:data[lab]=v
  assert set(data)==wanted
  rows=[E.modp_vec(data[lab],B.pos[weights[i]],len(B.blk[weights[i]])) for lab in labs]
  rank=E.rank_mod(rows,len(B.blk[weights[i]]));assert rank==len(labs)
  ranks[str(weights[i])]=rank
  for lab in labs:piv.add(lab,data[lab])
  del rows,data;gc.collect()
 assert sum(ranks.values())==6010
 piv.finish({'producer':'repair_f4.py','independent_pivots':6010,'prime':E.P,'sqrt13':E.F.s13,'sqrt_minus3':E.F.s3,'rank_by_block':ranks,'provenance':'Exact labelled family words from the recovered family streams; independently reranked modulo the stated prime','input_pivot_labels_sha256':hashlib.sha256(old.read_bytes()).hexdigest()})
 log('PIVOT RECONSTRUCTION COMPLETE',sum(ranks.values()))
else:
 for name in a.names.split(';'):
  tag=name.replace(' ','').replace('(','').replace(')','').replace(',','_')
  _,meta=recover('f4_tests_'+tag);assert meta['total']==956+(120 if name in ('e','x') else 0)
