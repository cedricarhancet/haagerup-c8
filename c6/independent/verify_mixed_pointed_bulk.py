"""Exact monodromy and pointed-relevant labels for the mixed extension."""
from pathlib import Path
from fractions import Fraction as F
import json

BASE=Path(__file__).resolve().parent
def q(p):return F(17*p*p,39)%1
def monodromy(p,s):return (q(p+s)-q(p)-q(s))%1
local=[p for p in range(39) if monodromy(p,13)==monodromy(p,26)==0]
assert local==list(range(0,39,3))
assert all((p+q_)%39 in local and (-p)%39 in local for p in local for q_ in local)
assert all(q(3*r)==(-F(r*r,13))%1 for r in range(13))
relevant={3:F(12,13),36:F(12,13),13:F(2,3),26:F(2,3),
          4:F(38,39),35:F(38,39),17:F(38,39),22:F(38,39)}
pointed_relevant={p:h for p,h in relevant.items() if p in local}
assert pointed_relevant=={3:F(12,13),36:F(12,13)}
result={'status':'PASS','arithmetic':'standard-library Fraction',
        'mixed_support_pointed_degrees':[0,13,26],
        'local_pointed_labels':local,'order':len(local),
        'restricted_metric':'q(3r)=-r^2/13 modulo1',
        'pointed_relevant_labels':{str(p):str(h) for p,h in pointed_relevant.items()},
        'scalar_scaling_dimension':'24/13','complex_dimension':2,'Hermitian_real_dimension':2,
        'premises':'The complete pointed lowest-weight list and one-dimensional tops are imported from the established Omega calculation.',
        'scope':'Checks finite monodromy arithmetic. Free-module simplicity and the full-center type-I formula are separately proved in pointed-bulk-mixed.tex. No mixed Haagerup action is assumed.'}
(BASE/'mixed-pointed-bulk-arithmetic.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result))
