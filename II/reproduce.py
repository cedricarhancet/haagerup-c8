"""Reproduce the published checks on a fresh working copy; keep input files intact."""
from pathlib import Path
import argparse,sys,os,json,shutil,subprocess,time,datetime,hashlib,importlib.metadata

ROOT=Path(__file__).resolve().parent
CALC=ROOT/'Calculs-HU-v08'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):raise RuntimeError('Run without -O, -OO or PYTHONOPTIMIZE')
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--group',choices=['smoke','certificate','modes','originals','extended','auxiliary','all'],default='smoke')
 p.add_argument('--output',type=Path,required=True)
 p.add_argument('--seconds-per-command',type=int,default=21600)
 a=p.parse_args();out=a.output.resolve()
 if a.group=='extended':
  from extended_replay import run_extended
  raise SystemExit(run_extended(ROOT,out,workers=1,seconds=a.seconds_per_command))
 delivery_root=ROOT.parent if (ROOT.parent/'MANIFEST-PUBLIC.json').is_file() else ROOT
 if out==delivery_root or out.is_relative_to(delivery_root):raise ValueError('Choose an output directory outside the delivery repository')
 if out.exists() and any(out.iterdir()):raise ValueError('Output directory must be new or empty')
 work=out/'work';logs=out/'logs';logs.mkdir(parents=True)
 shutil.copytree(CALC,work)
 env=os.environ.copy();env['PYTHONIOENCODING']='utf-8';env['PYTHONDONTWRITEBYTECODE']='1';env.pop('PYTHONOPTIMIZE',None)
 input_hashes={p.relative_to(CALC).as_posix():sha(p) for p in CALC.rglob('*') if p.is_file()}
 jobs=[]
 def add(name,cwd,script,*args):jobs.append((name,cwd,[script,*args]))
 if a.group in ['smoke','certificate','all']:
  add('categories','controles-categoriels','categories-exact-check.py')
  if a.group!='smoke':add('certificate','.','verify_omega_certificate.py','omega_certificate.json')
 if a.group in ['modes','all']:add('independent-modes','verification-independante','independent_tensors.py','..')
 if a.group in ['originals','all']:
  add('table-I','data','t12_table.py',str(ROOT/'inputs-I/Haagerup-c8-article-I-v14-arxiv.tex'))
  add('character','data','t12_W1.py','t12_table.json')
  add('split','data','t12_e8split.py','t12_W1.json')
  for s in ['test_e8e8','sl2x4','weights_W','t12b_optionB','r3a_G','r3a_Griess','r3a_groebner','r3a_omega_exact','r3b_chi']:
   add('R3-'+s,'R3',s+'.py')
  for i in range(4):add('R3-PCT-'+str(i),'R3','r3b_verify.py',str(i))
  add('reexport','.','export_omega_certificate.py','R3','omega_certificate-reexport.json')
  for s in ['test_voa','check_W1_char','solve_signs','embed','g_space','griess','setup_system','solve_groebner','verify_solutions','phi_chi','compare_R3']:
   add('R3-verif-'+s,'R3-verif',s+'.py')
 if a.group in ['auxiliary','all']:
  for s in ['check_moments','check_galois','check_m','verify_simplification']:add('section6-'+s,'section6',s+'.py')
  for s in ['e8facts','enum_sl2','sectors']:add('F1-'+s,'scripts-rapporteurs/F1',s+'.py')
  for s in ['chk','verlinde']:add('K1-'+s,'scripts-rapporteurs/K1',s+'.py')
  for s in ['gb_sympy','check_gb']:
   shutil.copy2(work/'scripts-rapporteurs/R3P'/(s+'.py'),work/'R3'/(s+'.py'))
   add('R3P-'+s,'R3',s+'.py')
 versions={}
 for name in ['sympy','python-flint','gmpy2','numpy']:
  try:versions[name]=importlib.metadata.version(name)
  except importlib.metadata.PackageNotFoundError:versions[name]=None
 receipt={'group':a.group,'python':sys.version,'dependencies':versions,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':[],'status':'RUNNING'}
 def save(): (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2),encoding='utf-8')
 save();failed=False
 for name,cwd,args in jobs:
  cmd=[sys.executable,'-B',*args];start=time.monotonic();print('RUN',name,flush=True)
  try:
   with (logs/(name+'.out')).open('wb') as stdout,(logs/(name+'.err')).open('wb') as stderr:
    cp=subprocess.run(cmd,cwd=work/cwd,env=env,stdout=stdout,stderr=stderr,timeout=a.seconds_per_command)
   code=cp.returncode
  except subprocess.TimeoutExpired:code='TIMEOUT'
  receipt['checks'].append({'name':name,'cwd':cwd,'command':['python','-B',*args],'exit_code':code,'seconds':round(time.monotonic()-start,3),'stdout_sha256':sha(logs/(name+'.out')),'stderr_sha256':sha(logs/(name+'.err'))})
  save();print('DONE',name,code,receipt['checks'][-1]['seconds'],flush=True)
  if code!=0:failed=True;break
 # Regenerations are compared to the selected original artifacts. Exact equality
 # of canonical JSON is stronger than matching printed success messages.
 comparisons=[]
 if a.group in ['originals','all'] and not failed:
  for f in ['data/t12_table.json','data/t12_W1.json','data/t12_e8split.json','R3/omegaC1.json','R3/r3b_chi.json',*[f'R3/r3b_verify_{i}.json' for i in range(4)]]:
   eq=json.loads((work/f).read_text())==json.loads((CALC/f).read_text())
   comparisons.append({'file':f,'canonical_json_equal':eq});failed|=not eq
  reexport=json.loads((work/'omega_certificate-reexport.json').read_text())
  cert=json.loads((CALC/'omega_certificate.json').read_text())
  eq=reexport==cert;comparisons.append({'file':'omega_certificate-reexport.json','canonical_json_equal':eq});failed|=not eq
 unchanged=all(sha(CALC/f)==h for f,h in input_hashes.items())
 receipt.update(status='FAIL' if failed or not unchanged else 'PASS',finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),original_delivery_unchanged=unchanged,regenerated_json_comparisons=comparisons)
 save();print(receipt['status'],len(receipt['checks']),'commands',flush=True)
 raise SystemExit(1 if receipt['status']=='FAIL' else 0)
if __name__=='__main__':main()
