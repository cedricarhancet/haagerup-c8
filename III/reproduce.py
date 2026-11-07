"""Download the fixed packages and run selected checks in an external workspace."""
from pathlib import Path
import argparse,datetime,hashlib,json,os,subprocess,sys,time
ROOT=Path(__file__).resolve().parent
sys.path.insert(0,str(ROOT))
from fetch_supplements import fetch
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):raise RuntimeError('Run without -O, -OO or PYTHONOPTIMIZE')
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--group',choices=['smoke','A','C3','zhu-verify','zhu-confluence','zhu-membership','zhu-fibre','zhu-weights','zhu-projection','zhu-modes10','zhu-regenerate','all'],default='smoke');p.add_argument('--jobs',type=int,default=1);p.add_argument('--rebuild-rules',action='store_true');p.add_argument('--local-assets',type=Path);a=p.parse_args()
 out=a.output.resolve();repo=ROOT.parent if (ROOT.parent/'MANIFEST-PUBLIC.json').is_file() else ROOT
 if out==repo or out.is_relative_to(repo):raise ValueError('Use a directory outside the repository')
 if out.exists():raise FileExistsError('Use a fresh output directory')
 if a.jobs<1:raise ValueError('--jobs must be positive')
 out.mkdir(parents=True)
 keys=['calculations'] if a.group in ['smoke','A'] else (['zhu'] if a.group.startswith('zhu-') else ['calculations','zhu','II','I_source','III_source'])
 packages={k:Path(v) for k,v in fetch(keys,out/'packages',a.local_assets).items()};jobs=[]
 if a.group=='smoke':jobs.append(('A1-smoke',[str(packages['calculations']/'A1/a1_lattices.py')]))
 if a.group in ['A','C3','all']:
  args=[str(packages['calculations']/'replay.py'),'--group',('quick' if a.group=='all' else a.group),'--output',str(out/'calculation-replay')]
  if a.group!='A':args+=['--zhu',str(packages['zhu']),'--II',str(packages['II']),'--artI',str(packages['I_source']),'--artIII',str(packages['III_source'])]
  jobs.append(('calculations',args))
 modes=['verify','projection','sample'] if a.group=='all' else ([a.group[4:]] if a.group.startswith('zhu-') else [])
 for mode in modes:
  m='sample' if mode=='modes10' else mode
  args=[str(packages['zhu']/'run_hj.py'),m,str(out/('zhu-'+mode))]
  if m in ['sample','regenerate']:args+=['--jobs',str(a.jobs)]
  if m=='sample':args+=['--rows','10']
  if m=='regenerate' and a.rebuild_rules:args+=['--rebuild-rules']
  jobs.append(('zhu-'+mode,args))
 env=os.environ.copy();env['PYTHONIOENCODING']='utf-8';env['PYTHONDONTWRITEBYTECODE']='1';env.pop('PYTHONOPTIMIZE',None)
 receipt={'status':'RUNNING','group':a.group,'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'checks':[]};logs=out/'logs';logs.mkdir()
 def save():(out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
 save()
 for name,args in jobs:
  print('RUN',name,flush=True);start=time.monotonic();log=logs/(name+'.out');err=logs/(name+'.err')
  with log.open('wb') as so,err.open('wb') as se:r=subprocess.run([sys.executable,'-I','-X','utf8','-B','-u',*args],env=env,stdout=so,stderr=se)
  row={'name':name,'command':['python','-I','-X','utf8','-B','-u',*args],'exit_code':r.returncode,'seconds':round(time.monotonic()-start,3),'stdout_sha256':sha(log),'stderr_sha256':sha(err)};receipt['checks'].append(row);save();print('DONE',name,r.returncode,flush=True)
  if r.returncode:receipt['status']='FAIL';save();return r.returncode
 receipt.update(status='PASS',finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());save();print('PASS:',len(jobs),'completed commands');return 0
if __name__=='__main__':raise SystemExit(main())
