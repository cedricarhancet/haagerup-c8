"""Fault injection for the actual F4 reader functions, on tiny isolated fixtures.
The functions are extracted unchanged from the verifier's AST, avoiding its
large top-level initialization. No production input or output is modified.
"""
from pathlib import Path
from types import SimpleNamespace
import argparse, ast, hashlib, json, pickle, tempfile
ap=argparse.ArgumentParser();ap.add_argument('--verifier',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args()
source=a.verifier.read_text();tree=ast.parse(source)
functions=[n for n in tree.body if isinstance(n,ast.FunctionDef) and n.name in ('tuplify','read_manifest','read_block')]
assert len(functions)==3
key=((0,)*8,());value={key:(1,0,0,0)};lab=('u(0)s','e',0);other=('u(0)s','e',1)
results=[]
with tempfile.TemporaryDirectory(prefix='f4-reader-negative-') as tmp:
 root=Path(tmp);basehash=hashlib.sha256(b'tiny source fixture').hexdigest()
 env={'root':root,'basehash':basehash,'weights':[(0,0)],'B4':SimpleNamespace(pos={(0,0):{key:0}}),'hashlib':hashlib,'json':json,'pickle':pickle}
 exec(compile(ast.Module(body=functions,type_ignores=[]),str(a.verifier),'exec'),env)
 def make(name,records=None,labels=None):
  path=root/name;path.mkdir();records=records if records is not None else [(lab,value)];labels=labels if labels is not None else [lab]
  data=b''.join(pickle.dumps(x,protocol=5) for x in records);p=path/'blk_000.pkl';p.write_bytes(data)
  m={'schema':1,'S3_basis_sha256':basehash,'nonzero_labels':[list(x) for x in labels],'zero_labels':[],'nonzero':len(labels),'zero':0,'total':len(labels),'blocks':[{'block':0,'weight':[0,0],'count':len(records),'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()}]}
  save=lambda:(path/'COMPLETE.json').write_text(json.dumps(m))
  save();return path,p,m,save
 def refresh(p,m):m['blocks'][0].update(bytes=p.stat().st_size,sha256=hashlib.sha256(p.read_bytes()).hexdigest())
 def exercise(name,mutate=None,records=None,labels=None,reject=True):
  path,p,m,save=make(name,records,labels)
  if mutate:mutate(path,p,m);save()
  try:
   d=env['read_manifest'](name);env['read_block'](d,0)
  except Exception as exc:
   assert reject,(name,'valid fixture rejected',repr(exc))
   results.append({'case':name,'status':'rejected','exception':type(exc).__name__})
  else:
   assert not reject,(name,'bad fixture accepted')
   results.append({'case':name,'status':'accepted'})
 exercise('valid',reject=False)
 exercise('missing_block',lambda path,p,m:p.unlink())
 exercise('changed_hash',lambda path,p,m:p.write_bytes(p.read_bytes()[:-1]+bytes([p.read_bytes()[-1]^1])))
 def truncate(path,p,m):p.write_bytes(p.read_bytes()[:-7]);refresh(p,m)
 exercise('truncated_record_with_updated_hash',truncate)
 exercise('duplicate_record_with_valid_hash',records=[(lab,value),(lab,value)],labels=[lab,other])
 exercise('unannounced_record_label',records=[(other,value)])
 def trailing(path,p,m):p.write_bytes(p.read_bytes()+b'junk');refresh(p,m)
 exercise('extra_trailing_bytes_with_updated_hash',trailing)
 exercise('missing_second_record',lambda path,p,m:m['blocks'][0].update(count=2),records=[(lab,value)],labels=[lab,other])
 exercise('unexpected_block',lambda path,p,m:(path/'blk_001.pkl').write_bytes(b'x'))
 exercise('extra_noncanonical_block_name',lambda path,p,m:(path/'blk_00.pkl').write_bytes(p.read_bytes()))
 exercise('wrong_weight_coordinate',records=[(lab,{((1,)+(0,)*7,()):(1,0,0,0)})])
 exercise('duplicate_manifest_label',lambda path,p,m:m.update(nonzero_labels=[list(lab),list(lab)],nonzero=2,total=2))
 exercise('overlap_zero_nonzero',lambda path,p,m:m.update(zero_labels=[list(lab)],zero=1,total=2))
 exercise('wrong_source_hash',lambda path,p,m:m.update(S3_basis_sha256='0'*64))
result={'status':'PASS','verifier_sha256':hashlib.sha256(source.encode()).hexdigest(),'method':'unmodified reader function AST; temporary fixtures only','valid_fixtures':1,'rejected_fixtures':sum(r['status']=='rejected' for r in results),'results':results}
a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
