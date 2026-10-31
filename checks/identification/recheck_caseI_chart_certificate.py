"""Recompute the 13 stored unit ideals in lex order, without the VOA engine."""
from pathlib import Path
import json
import sympy as sp
source=Path(__file__).with_name('check_caseI_exclusions.json');data=json.loads(source.read_text())
assert data['case']=='I' and data['coefficient_field']=='Q'
assert set(data['candidates'])=={'0','1','3','4'}
count=0
for c in data['candidates'].values():
 assert [x['first_nonzero_coordinate'] for x in c['charts']]==list(range(c['dimension']))
 for x in c['charts']:
  variables=list(map(sp.Symbol,x['variables']));polys=list(map(sp.sympify,x['polynomials']))
  assert list(sp.groebner(polys,*variables,domain=sp.QQ,order='lex'))==[1]
  count+=1
assert count==13
out={'result':'PASS','charts':count,'VOA_engine_used':False,'term_order':'lex'}
Path(__file__).with_suffix('.json').write_text(json.dumps(out,indent=2));print(out)
