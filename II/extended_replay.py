"""Complete the bounded S2 diagnostics in isolated work copies, with strict results."""
from pathlib import Path
import argparse,concurrent.futures,datetime,hashlib,importlib.metadata,json,os,re,shutil,subprocess,sys,threading,time

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def utc():return datetime.datetime.now(datetime.timezone.utc).isoformat()
def normalized(text):
 return '\n'.join(s.strip() for s in text.replace('\r\n','\n').splitlines() if not re.match(r'^time\s+[\d.]+s$',s.strip()))

def run_extended(source,output,workers=2,seconds=21600):
 if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):raise RuntimeError('Run without -O/-OO/PYTHONOPTIMIZE')
 source=Path(source).resolve();output=Path(output).resolve();calc=source/'Calculs-HU-v08'
 delivery_root=source.parent if (source.parent/'MANIFEST-PUBLIC.json').is_file() else source
 if output==delivery_root or output.is_relative_to(delivery_root):raise ValueError('Output must be outside the delivery repository')
 if output.exists() and any(output.iterdir()):raise ValueError('Use a fresh output directory')
 output.mkdir(parents=True,exist_ok=True)
 source_hashes={p.relative_to(calc).as_posix():sha(p) for p in sorted(calc.rglob('*')) if p.is_file() and '__pycache__' not in p.parts}
 versions={n:importlib.metadata.version(n) for n in ['sympy','python-flint','gmpy2','numpy']}
 state={'status':'RUNNING','started_utc':utc(),'python':sys.version,'dependencies':versions,'workers':workers,'source_hashes':source_hashes,'active':{},'results':[],'scope':'CI0 full closure (1200 pairs); CI1-3 archived scope (180 pairs each); secondary S2 reconstruction, PCT and rational-form diagnostic.'}
 lock=threading.Lock()
 def save():
  p=output/'RECEIPT.json';tmp=output/'RECEIPT.json.tmp';tmp.write_text(json.dumps(state,indent=2)+'\n',encoding='utf-8');tmp.replace(p)
 with lock:save()
 env=os.environ.copy();env['PYTHONIOENCODING']='utf-8';env['PYTHONDONTWRITEBYTECODE']='1';env.pop('PYTHONOPTIMIZE',None)

 def job(label,ci=None):
  dest=output/label;work=dest/'work';logs=dest/'logs';logs.mkdir(parents=True)
  shutil.copytree(calc,work)
  receipt={'name':label,'status':'RUNNING','started_utc':utc(),'checks':[]}
  def save_job(): (dest/'RECEIPT.json').write_text(json.dumps(receipt,indent=2)+'\n',encoding='utf-8')
  if ci is not None:
   jobs=[('S2','R3',['r3a_S2.py',str(ci)]),('S2full','R3',['r3a_S2_full.py',str(ci),*(['full'] if ci==0 else [])]),('PCT','R3',['r3b_verify.py',str(ci)])]
  else:jobs=[('S2-space','R3-verif',['s2_space.py','0','all']),('S2-check','R3-verif',['s2_check.py','0']),('twist','R3-verif',['twist_check.py'])]
  save_job()
  for name,cwd,args in jobs:
   start=time.monotonic()
   with lock:state['active'][label]={'command':args,'started_utc':utc()};save()
   print('RUN',label,name,flush=True)
   try:
    with (logs/(name+'.out')).open('wb') as so,(logs/(name+'.err')).open('wb') as se:
     cp=subprocess.run([sys.executable,'-B','-u',*args],cwd=work/cwd,env=env,stdout=so,stderr=se,timeout=seconds)
    code=cp.returncode
   except subprocess.TimeoutExpired:code='TIMEOUT'
   receipt['checks'].append({'name':name,'command':['python','-B','-u',*args],'cwd':cwd,'exit_code':code,'seconds':round(time.monotonic()-start,3),'stdout_sha256':sha(logs/(name+'.out')),'stderr_sha256':sha(logs/(name+'.err'))});save_job()
   print('DONE',label,name,code,receipt['checks'][-1]['seconds'],flush=True)
   if code!=0:
    receipt.update(status='FAIL',finished_utc=utc());save_job();raise RuntimeError(label+' '+name+' failed: '+str(code))
  try:
   comparisons=[]
   if ci is not None:
    f=f'R3/S2full_{ci}.json';d=json.loads((work/f).read_text())
    assert d['dim']==120 and d['U2_match'] is True and d['commutant_failures']==0,d
    assert not any(d['closure_failures'].values()),d
    assert d['full']==(ci==0) and d['closure_pairs']==(1200 if ci==0 else 180),d
    pct=json.loads((work/f'R3/r3b_verify_{ci}.json').read_text())
    assert pct['currents'] and pct['primaries'] and pct['omegaC1_fixed']
    assert pct['S2_stable_failures']==0 and pct['auto_samples'][1]==0 and pct['positivity_samples_fail']==0,pct
    dims={'d4':5,'d1':7,'d5':13,'d6':14,'c2':27}
    assert set(pct['slices'])==set(dims)
    assert all(v['proportional'] and v['inertia']==[dims[k],0,0] for k,v in pct['slices'].items())
    names=[f'R3/S2_{ci}.pkl',f'R3/S2full_{ci}.pkl',f'R3/S2full_{ci}.json',f'R3/r3b_verify_{ci}.json']
    for n in names:
     p=work/n;q=calc/n;row={'path':n,'bytes':p.stat().st_size,'sha256':sha(p),'baseline_sha256':sha(q),'bytes_equal':p.read_bytes()==q.read_bytes()}
     if p.suffix=='.json':row['canonical_json_equal']=json.loads(p.read_text())==json.loads(q.read_text());assert row['canonical_json_equal'],n
     else:assert row['bytes_equal'],n
     comparisons.append(row)
    receipt['validated_results']={'dimension':120,'commutant_failures':0,'closure_pairs':d['closure_pairs'],'modes_per_pair':3,'closure_failures':d['closure_failures'],'full':d['full'],'PCT_pass':True}
   else:
    txt=(logs/'S2-check.out').read_text(encoding='utf-8')
    required=['u0 : dim_K (S x 1) cap W_2 = 120','phi_chi maps (S x 1) cap W_2 into itself (all blocks): True','Hermitian: True','positive definite: True','leading principal minors: [50/507, 50/507, 50/507, 50/507, 50/507]']
    assert all(s in txt for s in required),txt
    assert re.search(r'weight multiplicities equal those of U2 = .* : True',txt),txt
    assert normalized(txt)==normalized((calc/'R3-verif/s2_check_0.out').read_text(encoding='utf-8'))
    twist=(logs/'twist.out').read_text(encoding='utf-8')
    assert re.search(r'coordinates in Q: 16 ; in Q\(sqrt\(-6\)\) but not Q: 48',twist),twist
    assert 'sl2^4 generators untouched by the twist: True' in twist
    assert all('u%d field generators used: [\'1\']'%i in twist for i in range(4)),twist
    match=re.search(r'twisted Q-form closed under products .*?, (\d+) coefficients\): failures (\d+)',twist)
    assert match and int(match[1])>0 and int(match[2])==0,twist
    for name in ['R3-verif/s2_space_0_all.pkl']:
     p=work/name;comparisons.append({'path':name,'bytes':p.stat().st_size,'sha256':sha(p),'baseline':'not shipped as a pickle in v01; generated in this run'})
    receipt['validated_results']={'dimension':120,'weight_multiplicities_equal':True,'PCT_stable_all_blocks':True,'hermitian':True,'positive_definite':True,'principal_minors':['50/507']*5,'twist_rational_choices':16,'twist_other_choices':48,'twist_coefficients_tested':int(match[1]),'twist_failures':0,'s2_check_matches_historical_output_except_time':True,'twist_matches_historical_output':normalized(twist)==normalized((calc/'R3-verif/twist_check.out').read_text(encoding='utf-8'))}
   receipt.update(status='PASS',finished_utc=utc(),comparisons=comparisons);save_job()
  except Exception as e:
   receipt.update(status='FAIL_VALIDATION',finished_utc=utc(),error=str(e));save_job();raise
  with lock:state['active'].pop(label,None);state['results'].append({'name':label,'status':'PASS','receipt':label+'/RECEIPT.json','validated_results':receipt['validated_results']});save()
  return receipt

 futures=[];errors=[]
 with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as pool:
  for label,ci in [('CI0',0),('second-engine',None),('CI1',1),('CI2',2),('CI3',3)]:futures.append(pool.submit(job,label,ci))
  for fut in concurrent.futures.as_completed(futures):
   try:fut.result()
   except Exception as e:errors.append(str(e));print('JOB ERROR',str(e),flush=True)
 unchanged={n:h for n,h in source_hashes.items() if not (calc/n).is_file() or sha(calc/n)!=h}
 state.update(status='FAIL' if errors or unchanged else 'PASS',finished_utc=utc(),errors=errors,source_unchanged=not unchanged)
 with lock:save()
 print('FINAL',state['status'],len(state['results']),'groups',flush=True)
 return 1 if state['status']!='PASS' else 0

if __name__=='__main__':
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--source',type=Path,default=Path(__file__).resolve().parent);p.add_argument('--output',type=Path,required=True);p.add_argument('--workers',type=int,choices=[1,2],default=2);p.add_argument('--seconds-per-command',type=int,default=21600)
 a=p.parse_args();raise SystemExit(run_extended(a.source,a.output,a.workers,a.seconds_per_command))
