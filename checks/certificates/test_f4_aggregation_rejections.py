"""Synthetic complete receipts test the F4 aggregator's failure paths.
No F4 mathematics is claimed by these fixtures. All writes are temporary.
"""
from pathlib import Path
import argparse,hashlib,json,subprocess,sys,tempfile
ap=argparse.ArgumentParser();ap.add_argument('--aggregator',required=True,type=Path);ap.add_argument('--output',required=True,type=Path);a=ap.parse_args()
source=a.aggregator.read_bytes();results=[]
def write(p,x):p.write_text(json.dumps(x,separators=(',',':')))
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def build(root):
 (root/'f4_verified').mkdir();(root/'exact_S3_basis.pkl').write_bytes(b'synthetic S3 basis, never decoded')
 (root/'source_stub.py').write_bytes(b'synthetic source bytes')
 sm={'archived_files':[{'name':'source_stub.py','bytes':len((root/'source_stub.py').read_bytes()),'sha256':sha(root/'source_stub.py')}],'new_verifiers':[]}
 write(root/'SOURCE-MANIFEST.json',sm)
 specs=[('f4_pivots',6010,0),('f4_families',18368,174)]+[(f'f4_tests_{i}',980,712 if i==0 else 0) for i in range(10)]
 manifests={};hashes={};zero_records=[]
 for name,total,zero in specs:
  nz=[[name,'nonzero',j] for j in range(total-zero)];zz=[[name,'zero',j] for j in range(zero)]
  m={'total':total,'nonzero':len(nz),'zero':len(zz),'nonzero_labels':nz,'zero_labels':zz,'blocks':[{'block':0},{'block':1}]}
  directory=root/name;directory.mkdir();write(directory/'COMPLETE.json',m);manifests[name]=m;hashes[name]=sha(directory/'COMPLETE.json')
  for lab in zz:zero_records.append({'manifest':name+'/COMPLETE.json','label':lab,'method':'direct exact kprod'})
 parts=[]
 for i in range(2):
  r={'parts':2,'part':i,'status':'PASS','source_manifest_sha256':sha(root/'SOURCE-MANIFEST.json'),'input_manifest_hashes':dict(hashes),'seen_labels':{name:m['nonzero_labels'][i::2] for name,m in manifests.items()},'blocks':[{'block':i,'all_rows_verified':True,'bad':0,'pivots':3005,'family_nonzero':9097,'tests_nonzero':4544}]}
  parts.append(r);write(root/'f4_verified'/f'PART_{i:03}_OF_002.json',r)
 z={'status':'PASS','verified':886,'records':zero_records,'manifest_hashes':{name+'/COMPLETE.json':hashes[name] for name in manifests if name!='f4_pivots'},'S3_basis_sha256':sha(root/'exact_S3_basis.pkl')}
 write(root/'zero-check.json',z);return parts,z
with tempfile.TemporaryDirectory(prefix='f4-aggregation-negative-') as tmp:
 temp=Path(tmp);script=temp/'aggregate-under-test.py';script.write_bytes(source)
 cases=['valid','missing_part','duplicate_part_number','inconsistent_part_manifest_hash','missing_seen_label','duplicate_seen_label','missing_block_with_preserved_counts','failed_all_rows','duplicate_zero_record_with_preserved_count','missing_zero_manifest','changed_source_bytes','optimized_python']
 for name in cases:
  root=temp/name;root.mkdir();parts,z=build(root)
  if name=='missing_part':(root/'f4_verified/PART_001_OF_002.json').unlink()
  elif name=='duplicate_part_number':parts[1]['part']=0
  elif name=='inconsistent_part_manifest_hash':parts[1]['input_manifest_hashes']['f4_pivots']='0'*64
  elif name=='missing_seen_label':parts[0]['seen_labels']['f4_families'].pop()
  elif name=='duplicate_seen_label':parts[1]['seen_labels']['f4_families'].append(parts[0]['seen_labels']['f4_families'][0])
  elif name=='missing_block_with_preserved_counts':
   for k in ('pivots','family_nonzero','tests_nonzero'):parts[0]['blocks'][0][k]*=2
   parts[1]['blocks']=[]
  elif name=='failed_all_rows':parts[0]['blocks'][0]['all_rows_verified']=False
  elif name=='duplicate_zero_record_with_preserved_count':z['records'][-1]=z['records'][0]
  elif name=='missing_zero_manifest':z['manifest_hashes'].pop('f4_tests_9/COMPLETE.json')
  elif name=='changed_source_bytes':(root/'source_stub.py').write_bytes(b'changed source bytes!!')
  if name!='missing_part':
   for i,r in enumerate(parts):write(root/'f4_verified'/f'PART_{i:03}_OF_002.json',r)
  write(root/'zero-check.json',z)
  cmd=[sys.executable]+(['-O'] if name=='optimized_python' else [])+[str(script),'--zeros',str(root/'zero-check.json')]
  run=subprocess.run(cmd,cwd=root,capture_output=True,text=True)
  accepted=run.returncode==0
  results.append({'case':name,'accepted':accepted,'expected_acceptance':name=='valid','returncode':run.returncode})
result={'status':'PASS' if all(r['accepted']==r['expected_acceptance'] for r in results) else 'FAIL','aggregator_sha256':hashlib.sha256(source).hexdigest(),'scope':'synthetic receipt tests only, no mathematical F4 verification','results':results}
a.output.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
if result['status']!='PASS':raise SystemExit(1)
