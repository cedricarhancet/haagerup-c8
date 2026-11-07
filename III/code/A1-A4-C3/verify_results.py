"""Check the delivered receipts, counters and public log hashes (no long replay)."""
from pathlib import Path
import hashlib,json,os,re
ROOT=Path(__file__).resolve().parent;BASE=ROOT/'checks/release-III-v10-v01'
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def check_receipt(path,count):
 r=load(path)
 if r['status']!='PASS' or len(r['checks'])!=count:raise RuntimeError('Incomplete receipt: '+str(path))
 for row in r['checks']:
  if row['exit_code']!=0 or row.get('reported_failed_check'):raise RuntimeError('Failed command: '+row['name'])
  for kind in ['stdout','stderr']:
   p=path.parent/row[kind+'_file']
   if hashlib.sha256(p.read_bytes()).hexdigest()!=row[kind+'_sha256']:raise RuntimeError('Changed public log: '+str(p))
 return r
def main():
 if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):raise RuntimeError('Run without Python optimization')
 s=load(BASE/'SUMMARY.json');assert s['status']=='PASS' and s['source_I_foundation_files_byte_identical']==19
 if (ROOT/'C3').is_dir():
  check_receipt(BASE/'A1-A4/RECEIPT.json',10);check_receipt(BASE/'C3/RECEIPT.json',4);check_receipt(BASE/'reference-v08/RECEIPT.json',4)
  for p in (BASE/'A1-A4').glob('a4_item*.out'):
   m=re.search(r'(\d+) checks, (\d+) FAUX',p.read_text(encoding='utf-8'));assert m and int(m[2])==0
  assert sum(x['checks'] for x in s['A4_detailed_checks'])==320
  d=load(BASE/'C3/comparison-II.json');assert d['status']=='PASS' and not d['errors'] and len(d['shared_files'])==131 and all(x['pass'] for x in d['shared_files'])
  assert not load(BASE/'C3/tables.json')['discrepancies']
  assert load(BASE/'reference-v08/comparison-v08.json')['summary']=={'regenerated_with_content_difference':0,'logs_differing':0,'missing':0}
  print('PASS: ten A1/A4 commands, 320 detailed A4 checks, four current C3 checks and the original-reference comparison; public log hashes verified.')
 else:
  check_receipt(BASE/'RECEIPT.json',4)
  v=load(BASE/'verify/verif/verif_exact.json');assert v['C_overlaps']==131698 and v['C_unresolved']==0 and v['D_inclusions']==151 and v['D_unresolved']==0
  assert v['B_generators_not_reducing_to_zero']==0 and v['F_dimension']==762 and v['G_non_integral']==0 and v['H_rules_with_different_reduction']==0
  f=load(BASE/'verify/Calculs/zhu-exact/specfib/z4_upper_algebra_modp.json');assert (f['dimension'],f['trace_form_rank'],f['center_dim'])==(762,762,12)
  assert load(BASE/'verify/Calculs/zhu-exact/exactZ.json')==load(ROOT/'reference/exactZ.json')
  assert load(BASE/'verify/Calculs/zhu-exact/check_patterns.json')==load(ROOT/'reference/check_patterns.json')
  assert load(BASE/'projection/result.json')=={'mode':'projection','status':'passed'}
  assert load(BASE/'projection/Calculs/zhu-exact/z3x.json')==load(ROOT/'reference/z3x.json')
  assert load(BASE/'modes10/result.json')=={'mode':'sample','status':'passed'}
  m=load(BASE/'MODE-SAMPLE.json');old=load(ROOT/'reexecution/modes/mode_sample.json')
  assert m['rows']==10 and m['pairs']==1155 and len(m['files'])==10 and m['all_ten_rows_byte_identical_to_historical_sample']
  assert all(x['sha256']==old['sample_digests'][x['file']] for x in m['files'])
  log=(BASE/'verify/membership_b_targeted.log').read_text(encoding='utf-8');assert re.search(r'554 of 554 failures 0\b',log)
  profile=(BASE/'verify.out').read_text(encoding='utf-8')
  for message in ['Vertex 12 in the original ideal','Vertex 13 in the original ideal','46 live eliminations in the original ideal','5105 initial rules in the original ideal','554 completion rules in the original ideal']:assert message in profile
  print('PASS: full confluence counters, both ideal-containment profile, special fibre, exact weights/branching, projection and ten-row mode sample; public log hashes verified.')
if __name__=='__main__':main()
