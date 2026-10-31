"""Closed-record streams with mandatory full reread before a COMPLETE manifest.
A missing record, duplicate, changed block, truncation or unexpected file aborts.
The whole production should run in a private temporary directory; only completed
files are copied to the deliverable. No assumption about the cause of earlier
stream losses is needed by this protocol.
"""
if not __debug__:raise RuntimeError('Run without -O or PYTHONOPTIMIZE')
from pathlib import Path
import hashlib,json,pickle,time
class SafeStore:
 def __init__(self,path,blocks,basehash):
  self.path=Path(path);self.path.mkdir(exist_ok=False)
  self.blocks=blocks;self.weights=sorted(blocks.blk);self.idx={m:i for i,m in enumerate(self.weights)}
  self.key_weights={k:m for m,ks in blocks.blk.items() for k in ks};self.basehash=basehash
  self.records={};self.zero=[];self.nonzero=[];self.started=time.monotonic()
 def add(self,label,v):
  label=tuple(label)
  assert label not in self.records,('duplicate label',label)
  if not v:self.records[label]=None;self.zero.append(label);return
  ws={self.key_weights[k] for k in v};assert len(ws)==1,('nonhomogeneous',label)
  i=self.idx[next(iter(ws))]
  with (self.path/f'blk_{i:03}.pkl').open('ab') as f:pickle.dump((label,v),f,protocol=5)
  self.records[label]=i;self.nonzero.append(label)
 def finish(self,extra=None):
  byblock={}
  for lab,i in self.records.items():
   if i is not None:byblock.setdefault(i,set()).add(lab)
  actual_files={int(p.stem.split('_')[1]) for p in self.path.glob('blk_*.pkl')}
  assert actual_files==set(byblock)
  blocks=[];allseen=set()
  for i,expected in sorted(byblock.items()):
   p=self.path/f'blk_{i:03}.pkl';seen=set()
   with p.open('rb') as f:
    for _ in range(len(expected)):
     lab,v=pickle.load(f)
     assert isinstance(lab,tuple) and lab in expected and lab not in seen
     assert isinstance(v,dict) and v and all(k in self.blocks.pos[self.weights[i]] for k in v)
     seen.add(lab)
    assert not f.read(1),('extra bytes',p)
   assert seen==expected and not allseen&seen
   allseen|=seen
   blocks.append({'block':i,'weight':list(self.weights[i]),'count':len(seen),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
  assert allseen==set(self.nonzero)
  data={'schema':1,'S3_basis_sha256':self.basehash,'total':len(self.records),'nonzero':len(self.nonzero),'zero':len(self.zero),'nonzero_labels':self.nonzero,'zero_labels':self.zero,'blocks':blocks,'post_write_full_reread':True,'elapsed_seconds':time.monotonic()-self.started,**(extra or {})}
  p=self.path/'COMPLETE.json';q=self.path/'COMPLETE.tmp';q.write_text(json.dumps(data,indent=1));q.replace(p)
  return data
