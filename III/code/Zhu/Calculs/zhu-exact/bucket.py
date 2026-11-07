"""Split exact products prodx/a%03d.pkl by Cartan weight of the pair: prodw/<kind>_<w0>_<w1>/a%03d.pkl, kind in B, N.
Small products (p1, p2, Tp2) are collected in prodw/small.pkl."""
import pickle, os, glob, gc
from collections import defaultdict
d = pickle.load(open('/home/claude/hj/work/Calculs/zhu/U2_lowering_basis.pkl', 'rb'))
HV = d['hv']; LB = [(a, i, j) for a, i, j, _ in d['basis']]; TYPE = [w for _, w in HV]
def wt(b): a, i, j = LB[b]; r, l = TYPE[a]; return (r - 2 * i, l - 2 * j)
OUT = 'prodw'; os.makedirs(OUT, exist_ok=True)
small = {}
for a in range(120):
    res = pickle.load(open('prodx/a%03d.pkl' % a, 'rb'))
    groups = defaultdict(dict)
    for b, r in res.items():
        w = (wt(a)[0] + wt(b)[0], wt(a)[1] + wt(b)[1])
        small[a, b] = {'p2': r['p2'], 'p1': r['p1'], 'Tp2': r['Tp2']}
        groups['B', w][a, b] = r['B']
        if r['N'] is not None: groups['N', w][a, b] = r['N']
    for (kind, w), dct in groups.items():
        dd = '%s/%s_%d_%d' % (OUT, kind, w[0], w[1]); os.makedirs(dd, exist_ok=True)
        pickle.dump(dct, open(dd + '/a%03d.pkl' % a, 'wb'))
    del res, groups; gc.collect()
pickle.dump(small, open(OUT + '/small.pkl', 'wb'))
print('done', len(small))
