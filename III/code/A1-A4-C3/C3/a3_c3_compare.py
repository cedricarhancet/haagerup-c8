"""Compare C3 to the fixed public Article II release, without unpickling data."""
from pathlib import Path
import argparse,gzip,hashlib,json,sys
HERE=Path(__file__).resolve().parent
EXPECTED_TAG='article-II-v14-supplements-v03'
EXPECTED_COMMIT='134ce1476c7d4ce9b336ea1d167d707870650e33'
def sha(data):return hashlib.sha256(data).hexdigest()
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--reference',type=Path,required=True);p.add_argument('--report',type=Path,default=Path('c3-II-comparison.json'));a=p.parse_args()
 ref=a.reference.resolve()
 if ref.name=='Calculs-HU-v08':ref=ref.parent
 if (ref/'II/MANIFEST-II.json').is_file():ref=ref/'II'
 if load(ref/'VERSION-CONCORDANCE.json')['supplement_tag']!=EXPECTED_TAG:raise RuntimeError('Use the fixed public II release '+EXPECTED_TAG)
 for n,r in load(ref/'MANIFEST-II.json').items():
  rel=Path(n)
  if rel.is_absolute() or '..' in rel.parts:raise ValueError('Unsafe manifest path')
  b=(ref/rel).read_bytes()
  if len(b)!=r['bytes'] or sha(b)!=r['sha256']:raise RuntimeError('II manifest mismatch: '+n)
 rows=[];only=[];errors=[];run=HERE/'run';calc=ref/'Calculs-HU-v08'
 def compare(n,x,y,kind):
  equal=x==y;canonical=None
  if not equal and n.endswith('.json'):canonical=json.loads(x.decode('utf-8'))==json.loads(y.decode('utf-8'))
  ok=equal or canonical is True
  rows.append({'file':n,'kind':kind,'C3_sha256':sha(x),'II_sha256':sha(y),'identical_bytes':equal,'canonical_json_equal':canonical,'pass':ok})
  if not ok:errors.append(n)
 for x in sorted(run.rglob('*')):
  if not x.is_file() or '__pycache__' in x.parts or x.suffix=='.pyc':continue
  n=x.relative_to(run).as_posix();y=calc/n
  if y.is_file():compare(n,x.read_bytes(),y.read_bytes(),'direct shared file')
  elif n=='R3-verif/s2_space_0_all.pkl':
   g=ref/'checks/release-II-v14-v02/generated/R3-verif/s2_space_0_all.pkl.gz'
   compare(n,x.read_bytes(),gzip.decompress(g.read_bytes()),'secondary space delivered compressed by II')
  elif n=='omega_certificate-reexport-C3.json':compare(n,x.read_bytes(),(calc/'omega_certificate.json').read_bytes(),'re-export vs published portable certificate')
  else:only.append(n)
 if not rows:raise RuntimeError('No common files compared')
 result={'status':'FAIL' if errors else 'PASS','II_tag':EXPECTED_TAG,'II_commit':EXPECTED_COMMIT,'shared_files':rows,'III_only_files':only,'errors':errors,'scope':'Shared delivered source/data/log files are compared byte for byte, or as JSON objects when line endings differ. The secondary space is decompressed from the pinned II artifact. III-only diagnostics are listed explicitly and are not described as comparisons to II. This does not rerun the full 69 archived C3 production commands.'}
 a.report.parent.mkdir(parents=True,exist_ok=True);a.report.write_text(json.dumps(result,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
 print('PASS' if not errors else 'FAIL',len(rows),'shared files;',len(only),'III-only files;',len(errors),'differences')
 if errors:print('\n'.join(errors))
 return 1 if errors else 0
if __name__=='__main__':raise SystemExit(main())
