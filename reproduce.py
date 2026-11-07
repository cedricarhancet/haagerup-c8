#!/usr/bin/env python3
"""Relancer les controles publics de I v14 sans ecraser les sorties archivees."""
import argparse, datetime, hashlib, importlib.metadata, json, os, platform, shutil, subprocess, sys, time
from pathlib import Path

ROOT=Path(__file__).resolve().parent
def main():
    if not __debug__ or os.environ.get('PYTHONOPTIMIZE'):
        raise SystemExit('Executer sans -O et sans PYTHONOPTIMIZE.')
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--group',choices=['smoke','foundations','modes','lattice','quick'],default='smoke')
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    if a.group in ['smoke','quick']:
        shutil.copy2(ROOT/'checks/identification/c4/c4-exact-witnesses.json',out/'c4-exact-witnesses.json')
    groups={
      'smoke':[
        ['checks/identification/recheck_c4_matrices.py',str(out/'c4-exact-witnesses.json')],
        ['verification-certificates/check_closure_coverage.py']],
      'foundations':[
        ['checks/foundations/check_analytic_identities.py','--sources','calculs-sources/classification'],
        ['verification-foundations/verify_modular_certificates.py','--output',str(out/'modular-certificates.json')],
        ['verification-foundations/verify_u2_moments.py']],
      'modes':[
        ['checks/independent-modes/independent_modes.py','--engine','checks/independent-modes/reference/lattice_voa.py','--output',str(out/'basis-modes.json')],
        ['checks/independent-modes/replay_generators.py','--engine','checks/independent-modes/reference/lattice_voa.py','--inputs','checks/independent-modes/generator-inputs.json','--record','checks/independent-modes/generator-modes-results.json','--output',str(out/'generator-modes.json')]],
      'lattice':[
        ['checks/certificates/deliverable/check_lattice_ball.py','--sources','calculs-sources/generation-forte','--output',str(out/'lattice-ball.json')],
        ['verification-certificates/check_simple_current_word.py','--source-dir','calculs-sources/generation-forte','--output',str(out/'simple-current-word.json')]]}
    selected=['smoke','foundations','modes','lattice'] if a.group=='quick' else [a.group]
    receipt={'date_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'python':sys.version,'platform':platform.platform(),'group':a.group,'dependencies':{},'checks':[]}
    for name in ['sympy','numpy','gmpy2','python-flint','mpmath']:
        try:receipt['dependencies'][name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:pass
    status=0
    for group in selected:
        for argv in groups[group]:
            name=Path(argv[0]).stem;log=out/(name+'.log');t=time.monotonic()
            print('RUN',name,flush=True)
            with log.open('w',encoding='utf-8') as f:
                r=subprocess.run([sys.executable,'-u']+argv,cwd=ROOT,stdout=f,stderr=subprocess.STDOUT)
            public_argv=[s.replace(str(out),'<OUTPUT>') for s in argv]
            receipt['checks'].append({'id':name,'command':['python']+public_argv,'exit_code':r.returncode,'seconds':round(time.monotonic()-t,3),'log':log.name,'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest()})
            print('PASS' if r.returncode==0 else 'FAIL',name,flush=True)
            if r.returncode:status=1;break
        if status:break
    receipt['status']='PASS' if status==0 else 'FAIL'
    (out/'RECEIPT.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
    return status
if __name__=='__main__':sys.exit(main())
