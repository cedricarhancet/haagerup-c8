"""a4_smoke.py -- smoke test of the engine: sl2 relations, Virasoro relations and central charges."""
import os as _os, sys as _sys
_sys.path.insert(0, _os.path.dirname(_os.path.abspath(__file__)))
from a4_common import *
from a4_voa import sl2_checks, virasoro_checks, sugawara, commute_checks

L, voa = make_voa()
t = Timer()
g = huang_generators(voa)
j = jiao_generators(voa)
print(t('built'))
for r in sl2_checks(voa, g['EL'], g['FL'], g['HL'], 1, 'L: '):
    print(r)
for r in sl2_checks(voa, g['ER'], g['FR'], g['HR'], 39, 'R: '):
    print(r)
print(commute_checks(voa, [g['EL'], g['FL'], g['HL']], [g['ER'], g['FR'], g['HR']]))
for r in sl2_checks(voa, j['PE'], j['PF'], j['PH'], 3, 'P: '):
    print(r)
for r in sl2_checks(voa, j['BE'], j['BF'], j['BH'], 13, 'B: '):
    print(r)
print(commute_checks(voa, [j['PE'], j['PF'], j['PH']], [j['BE'], j['BF'], j['BH']]))
print(t('sl2'))
om = voa.omega()
res, c = virasoro_checks(voa, om, 'omega: ')
print(res, c)
print(t('omega'))
for name, (E, F, H, k) in {"wL": (g['EL'], g['FL'], g['HL'], 1), "wR": (g['ER'], g['FR'], g['HR'], 39),
                           "wP": (j['PE'], j['PF'], j['PH'], 3), "wB": (j['BE'], j['BF'], j['BH'], 13)}.items():
    w = sugawara(voa, E, F, H, k)
    res, c = virasoro_checks(voa, w, name + ': ')
    print(name, res, c, len(w))
    # omega - w is Virasoro too (commutant)
    res2, c2 = virasoro_checks(voa, om - w, 'omega-' + name + ': ')
    print(name, res2, c2)
    print(t(name))
