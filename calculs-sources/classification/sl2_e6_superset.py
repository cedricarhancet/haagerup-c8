"""sl2-subalgebras of e6: a computed SUPERSET of the possible neutral elements (Axis C3).

By Dynkin (see Collingwood--McGovern, Thm 3.5.4), the neutral element H of an sl2-triple of e6 is conjugate
to a dominant element with alpha_i(H) in {0,1,2}.  Necessary conditions tested here, for each of the
3^6 labellings: the eigenvalues of ad H on e6 and of H on the 27 form the weight multiset of a genuine
sl2-module (multiplicities m(j) are non-increasing in |j| within each parity class).  The surviving
labellings form a superset of the 20 nonzero nilpotent orbits; we report the Dynkin index (H,H)/2,
the decomposition of 27 and of e6 into sl2-irreducibles, and the dimension of the centralizer of the sl2.
Using a superset is safe for exclusion arguments.
"""
import itertools, json
from fractions import Fraction as Fr
import e6_tools as E

Cinv = None
def H_from_labels(lab):
    # H = sum c_i alpha_i^vee with alpha_j(H) = sum_i c_i C[i][j] = lab_j  ->  c = lab C^{-1}
    import sympy as sp
    M = sp.Matrix(E.C)
    c = sp.Matrix([lab]) * M.inv()
    return [Fr(int(sp.fraction(x)[0]), int(sp.fraction(x)[1])) for x in c]

def decompose(eigs):
    """multiset of integer eigenvalues -> dict {highest weight n: multiplicity} or None if not an sl2 module"""
    from collections import Counter
    cnt = Counter(eigs)
    if any(Fr(x).denominator != 1 for x in cnt):
        return None
    cnt = Counter({int(k): v for k, v in cnt.items()})
    if any(cnt[j] != cnt[-j] for j in cnt):
        return None
    top = max(abs(j) for j in cnt)
    dec = {}
    for n in range(top, -1, -1):
        m = cnt.get(n, 0) - cnt.get(n + 2, 0)
        if m < 0:
            return None
        if m:
            dec[n] = m
    return dec

out = []
for lab in itertools.product((0, 1, 2), repeat=6):
    if not any(lab):
        continue
    c = H_from_labels(lab)
    ad = [E.ev(c, r) for r in E.ROOTS] + [0]*6
    w = [E.ev(c, x) for x in E.W27]
    dad = decompose(ad); d27 = decompose(w)
    if dad is None or d27 is None:
        continue
    index = E.form(c)/2
    assert index.denominator == 1
    out.append({'labels': lab, 'index': int(index), 'dec27': d27, 'dec_ad': dad,
                'centralizer_dim': dad.get(0, 0)})
idx = sorted(set(o['index'] for o in out))
print('labellings passing the necessary test:', len(out))
print('Dynkin indices occurring:', idx)
for o in sorted(out, key=lambda o: o['index']):
    print('  index %3d  labels %s  27 = %s  centralizer dim %d' % (o['index'], o['labels'], o['dec27'], o['centralizer_dim']))
json.dump([{**o, 'labels': list(o['labels']), 'dec27': {str(k): v for k, v in o['dec27'].items()},
            'dec_ad': {str(k): v for k, v in o['dec_ad'].items()}} for o in out], open('sl2_e6_superset.json', 'w'), indent=0)
