#!/usr/bin/env python3
"""Case II, candidate 9: exhaustive eigenspaces, exact weights and binary lines.
The first stage examines nonlinear factors and refuses undecided exclusions.
The second is called only after the first has certified that the admissible
space is unique for each type. Original VOA/source formulas are reused.
"""
import argparse,json,os,sys,time
from pathlib import Path
import sympy as sp
ap=argparse.ArgumentParser();ap.add_argument('--source-dir',required=True);ap.add_argument('--output',required=True);a=ap.parse_args()
source=Path(a.source_dir).resolve();output=Path(a.output).resolve();sys.path.insert(0,str(source));os.chdir(source)
import exact_isotropy_s13 as ES, isotropic_lines as IL, exact_isotropy as EI
start=time.time();r=ES.run('II',9)
for typ,row in r['types'].items():
    kept=[v for v in row['eigen'] if v['verdict']=='exists']
    assert len(kept)==1 and kept[0].get('dimE')==2 and 'xi' in kept[0],(typ,kept)
    assert not any('undecided' in v.get('why','') or v['verdict']=='undecided' for v in row['eigen'])
D,spaces=IL.admissible_spaces('II',9)
cand=json.loads((source/'griess_subalgebras.json').read_text())['II']['subalgebras'][9]
b=1/sp.sqrt(sp.sympify(cand['D']));a0=(1-b*sp.sympify(cand['lambda']))/2
expected={(0,4):('112/95','8/19'),(1,1):('111/95','12/19'),(1,3):('147/95','1/19'),(1,5):('7/95','68/57'),(1,7):('56/95','4/19'),(2,2):('56/95','50/57'),(2,6):('56/95','4/19'),(3,3):('7/19','12/19')}
checks=[]
for tau,row in spaces.items():
    assert len(row['spaces'])==1;s=row['spaces'][0];assert s['E'].shape[1]==2
    kept=next(v for v in r['types'][str(tau)]['eigen'] if v['verdict']=='exists')
    assert sp.simplify(s['xi']-sp.sympify(kept['xi']))==0
    h1=sp.simplify(a0*(2-EI.h_tau('II',tau))+2*b*s['xi']);h2=sp.simplify(2-EI.h_tau('II',tau)-h1)
    assert sp.simplify(h1-sp.sympify(expected[tau][0]))==0
    assert sp.simplify(h2-sp.sympify(expected[tau][1]))==0
    lines=IL.lines_in(s);assert isinstance(lines,list) and len(lines)==2
    norms=[]
    for u in lines:
        assert all(sp.simplify((u.T*Q*u)[0])==0 for Q in s['Qs'])
        norm=sp.simplify((u.T*s['Q0']*u)[0]);assert norm.is_zero is False,(tau,norm)
        norms.append(str(norm))
    determinant=sp.simplify(sp.Matrix.hstack(*lines).det());assert determinant.is_zero is False
    checks.append({'type':str(tau),'dimension':2,'h1':str(h1),'h2':str(h2),'norms':norms,'line_determinant':str(determinant)})
    print(tau,'PASS',h1,h2,flush=True)
output.write_text(json.dumps({'status':'PASS','elapsed_seconds':time.time()-start,'all_eigenspaces':r,'surviving_spaces':checks},indent=2)+'\n')
