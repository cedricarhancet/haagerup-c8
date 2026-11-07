"""Check the delivered v02 receipts, counters and generated artifact hashes."""
from pathlib import Path
import gzip,hashlib,json,re,os
ROOT=Path(__file__).resolve().parent
CHECKS=ROOT/'checks/release-II-v14-v02'
def sha(data):return hashlib.sha256(data).hexdigest()
def main():
 if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):raise RuntimeError('Run without -O/-OO/PYTHONOPTIMIZE')
 summary=json.loads((CHECKS/'SUMMARY.json').read_text(encoding='utf-8'))
 assert summary['status']=='PASS' and summary['new_commands']==15
 root=json.loads((CHECKS/'extended/RECEIPT.json').read_text(encoding='utf-8'))
 assert root['status']=='PASS' and root['source_unchanged'] and not root['errors'] and len(root['results'])==5
 count=0
 for label in ['CI0','CI1','CI2','CI3','second-engine']:
  p=CHECKS/'extended'/label;r=json.loads((p/'RECEIPT.json').read_text(encoding='utf-8'))
  assert r['status']=='PASS' and len(r['checks'])==3
  for c in r['checks']:
   assert c['exit_code']==0
   for ext,key in [('out','stdout_sha256'),('err','stderr_sha256')]:assert sha((p/'logs'/(c['name']+'.'+ext)).read_bytes())==c[key]
  v=r['validated_results']
  if label.startswith('CI'):
   ci=int(label[-1]);n=1200 if ci==0 else 180
   assert v['dimension']==120 and v['closure_pairs']==n and v['modes_per_pair']==3 and v['PCT_pass']
   assert not any(v['closure_failures'].values()) and v['commutant_failures']==0
   d=json.loads((ROOT/f'Calculs-HU-v08/R3/S2full_{ci}.json').read_text())
   assert d['dim']==120 and d['U2_match'] and d['commutant_failures']==0 and d['closure_pairs']==n and d['full']==(ci==0) and not any(d['closure_failures'].values())
   pct=json.loads((ROOT/f'Calculs-HU-v08/R3/r3b_verify_{ci}.json').read_text())
   assert pct['currents'] and pct['primaries'] and pct['omegaC1_fixed'] and pct['S2_stable_failures']==0
   assert pct['auto_samples'][1]==0 and pct['positivity_samples_fail']==0
   dims={'d4':5,'d1':7,'d5':13,'d6':14,'c2':27}
   assert set(pct['slices'])==set(dims)
   assert all(x['proportional'] and x['inertia']==[dims[k],0,0] for k,x in pct['slices'].items())
   count+=n
  else:
   assert v['dimension']==120 and v['weight_multiplicities_equal'] and v['PCT_stable_all_blocks']
   assert v['hermitian'] and v['positive_definite'] and v['principal_minors']==['50/507']*5
   assert v['twist_rational_choices']==16 and v['twist_other_choices']==48 and v['twist_failures']==0 and v['twist_coefficients_tested']>0
   log=(p/'logs/S2-check.out').read_text(encoding='utf-8')
   assert 'positive definite: True' in log and 'Hermitian: True' in log and 'phi_chi maps (S x 1) cap W_2 into itself (all blocks): True' in log
   twist=(p/'logs/twist.out').read_text(encoding='utf-8');assert re.search(r'coefficients\): failures 0',twist)
 for row in json.loads((CHECKS/'ARTIFACTS.json').read_text(encoding='utf-8')):
  p=ROOT/row['path'];data=p.read_bytes();assert len(data)==row['bytes'] and sha(data)==row['sha256']
  if p.suffix=='.gz':
   raw=gzip.decompress(data);assert len(raw)==row['uncompressed_bytes'] and sha(raw)==row['uncompressed_sha256']
 for row in json.loads((CHECKS/'MATH-CODE-IDENTITY.json').read_text(encoding='utf-8')):
  assert sha((ROOT/'Calculs-HU-v08'/row['path']).read_bytes())==row['sha256']
 allow=json.loads((ROOT/'ALLOWLIST-COMPUTATIONS.json').read_text(encoding='utf-8'))
 for row in allow['files']:
  data=(ROOT/row['path']).read_bytes();assert len(data)==row['bytes'] and sha(data)==row['sha256']
 assert count==1740
 print('PASS: 15 completed commands,',count,'pairs across the stated scopes,',count*3,'mode checks; secondary PCT/positivity and generated hashes verified.')
if __name__=='__main__':main()
