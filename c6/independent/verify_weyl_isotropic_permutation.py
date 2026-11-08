"""Exact F3 check of the four isotropic lines and the Weyl permutation."""
from itertools import product
from pathlib import Path
from math import gcd
import json

BASE = Path(__file__).resolve().parent

def q(v):
    x,y,z = v
    return (x*y+2*z*z)%3

def normalize(v):
    first = next(x for x in v if x)
    inverse = 1 if first == 1 else 2
    return tuple((inverse*x)%3 for x in v)

def w(v):
    return (v[0],v[1],(-v[2])%3)

vectors = list(product(range(3),repeat=3))
assert all(q(w(v)) == q(v) and w(w(v)) == v for v in vectors)
lines = sorted({normalize(v) for v in vectors if any(v) and q(v)==0})
assert lines == [(0,1,0),(1,0,0),(1,1,1),(1,1,2)]
permutation = {str(v):normalize(w(v)) for v in lines}
assert normalize(w((0,1,0))) == (0,1,0)
assert normalize(w((1,0,0))) == (1,0,0)
assert normalize(w((1,1,1))) == (1,1,2)
assert normalize(w((1,1,2))) == (1,1,1)
# A_Q has q(e,t)=(2e^2+t^2)/3. The coordinate change is invertible.
change = {((t+e)%3,(t-e)%3) for e,t in product(range(3),repeat=2)}
assert len(change) == 9
assert all(((t+e)*(t-e)-2*e*e-t*t)%3==0
           for e,t in product(range(3),repeat=2))
pointed_isometries=[u for u in range(39) if gcd(u,39)==1 and (u*u-1)%39==0]
assert pointed_isometries==[1,14,25,38]
assert all((14*(26*m+21*n)-(-26*m+21*n))%39==0
           and (25*(26*m+21*n)-(26*m-21*n))%39==0
           for m,n in product(range(39),repeat=2))
assert (14*14)%39==(25*25)%39==1 and (14*25)%39==38
result = {'status':'PASS','field':'F3','vectors_checked':len(vectors),
          'quadratic_form_numerator':'xy+2z^2',
          'isotropic_lines':lines,'weyl_action':'(x,y,z) -> (x,y,-z)',
          'permutation':permutation,'pure_lines_fixed':2,'mixed_lines_exchanged':2,
          'pointed_isometries':pointed_isometries,
          'individual_weyl_multipliers':[14,25],'product_multiplier':38,
          'momentum_pairs_checked':39**2,
          'scope':'Finite metric-group arithmetic only; the actual VOA/graph action is proved and audited separately.'}
(BASE/'weyl-isotropic-permutation.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
print(json.dumps(result))
