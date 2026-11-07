"""Replay A1/A4 and the fixed C3 comparisons in a fresh external working copy."""
from pathlib import Path
import argparse,datetime,hashlib,json,os,re,shutil,subprocess,sys,time
ROOT=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def main():
 if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):raise RuntimeError('Use Python without optimization')
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--group',choices=['A','C3','quick'],default='quick');p.add_argument('--output',type=Path,required=True);p.add_argument('--zhu',type=Path);p.add_argument('--II',type=Path);p.add_argument('--artI',type=Path);p.add_argument('--artIII',type=Path);p.add_argument('--seconds-per-command',type=int,default=21600);a=p.parse_args();out=a.output.resolve()
 if out==ROOT or out.is_relative_to(ROOT):raise ValueError('Use an output directory outside the archive')
 if out.exists():raise ValueError('Use a new output directory')
 if a.group in ['C3','quick'] and any(x is None for x in [a.zhu,a.II,a.artI,a.artIII]):p.error('C3 comparisons require --zhu, --II, --artI and --artIII')
 out.mkdir(parents=True);work=out/'work';shutil.copytree(ROOT,work);logs=out/'logs';logs.mkdir()
 source={p.relative_to(ROOT).as_posix():sha(p) for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc'}
 receipt={'status':'RUNNING','group':a.group,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'python':sys.version,'checks':[]}
 def save():(out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
 jobs=[]
 if a.group in ['A','quick']:
  jobs += [(n,['A1/'+n+'.py']) for n in ['a1_lattices','a1b_extra']]
  jobs += [(n,['A4/'+n+'.py']) for n in ['a4_kfield','a4_smoke','a4_item1_lattice','a4_item2_section5','a4_item3_section6','a4_item4_section7','a4_item5_prop81','a4_item6_riesen']]
 if a.group in ['C3','quick']:
  jobs += [('certificate',['C3/run/verify_omega_certificate.py','C3/run/omega_certificate.json']),('categories',['C3/run/controles-categoriels/categories-exact-check.py']),('tables',['C3/a3_c3_table_check.py','--artI',str(a.artI.resolve()),'--artIII',str(a.artIII.resolve()),'--zhu',str(a.zhu.resolve()),'--out',str(out/'tables.json')]),('comparison-II',['C3/a3_c3_compare.py','--reference',str(a.II.resolve()),'--report',str(out/'comparison-II.json')])]
 env=os.environ.copy();env['PYTHONIOENCODING']='utf-8';env['PYTHONDONTWRITEBYTECODE']='1';env.pop('PYTHONOPTIMIZE',None);save()
 for label,args in jobs:
  start=time.monotonic();print('RUN',label,flush=True);log=logs/(label+'.out');err=logs/(label+'.err')
  try:
   with log.open('wb') as so,err.open('wb') as se:r=subprocess.run([sys.executable,'-I','-X','utf8','-B','-u',*args],cwd=work,env=env,stdout=so,stderr=se,timeout=a.seconds_per_command)
   code=r.returncode
  except subprocess.TimeoutExpired:code='TIMEOUT'
  text=log.read_text(encoding='utf-8');bad=re.search(r'^FAUX\b',text,re.M) is not None
  row={'name':label,'command':['python','-I','-X','utf8','-B','-u',*args],'exit_code':code,'seconds':round(time.monotonic()-start,3),'stdout_sha256':sha(log),'stderr_sha256':sha(err),'reported_failed_check':bad};receipt['checks'].append(row);save();print('DONE',label,code,flush=True)
  if code!=0 or bad:receipt.update(status='FAIL');save();return 1
  if label in ['a1_lattices','a1b_extra'] and 'ALL CHECKS PASSED' not in text:raise RuntimeError('Missing final A1 validation')
 if a.group in ['C3','quick']:
  if load(out/'tables.json')['discrepancies'] or load(out/'comparison-II.json')['status']!='PASS':raise RuntimeError('C3 result validation failed')
 receipt.update(status='PASS',finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),source_unchanged=all(sha(ROOT/n)==h for n,h in source.items()))
 if not receipt['source_unchanged']:raise RuntimeError('The original archive changed')
 save();print('PASS:',len(jobs),'commands; source package unchanged');return 0
if __name__=='__main__':raise SystemExit(main())
