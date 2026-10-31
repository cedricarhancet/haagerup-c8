"""Axis G (v12).  (1) Classes and genera of the even binary lattices of determinant 39 (discriminant -39), exact.
(2) Character test of the genus picture U = (+)_{alpha in L'/L} C^alpha (x) V_{alpha+L}: from the dimensions of the
blocks of fixed h-weight mu = (m', m'') (eigenvalues of the coroots h'(0), h''(0)) of U_n, n <= 3, each case gives the
graded dimensions f_alpha(x) of the C-modules C^alpha (x = conformal weight), by
    sum_x f_alpha(x) q^(x + |mu|^2/2) = prod_m (1-q^m)^2 * sum_n dim U_n^mu q^n,
with |mu|^2 = m'^2/(2k') + m''^2/(2k''), alpha = mu mod L, L = <(2k',0), (0,2k''), (k',k'')>, L' = {m' = m'' mod 2}.
Checks: (a) consistency inside each case (different mu in the same coset); (b) comparison of the two cases through
each isometry psi of the discriminant forms (Z/39, q).  Case II: exact bases (U_2, U_3, all blocks); case I: U^R_2,
U^R_3 modulo 7999993 (all blocks, ranks = 120, 956).  Usage: python3 s20_genus_check.py"""
import sys, os, json, pickle, itertools
from fractions import Fraction as F
import numpy as np
HERE = os.path.abspath(os.path.dirname(__file__)); GF = os.path.abspath(HERE + '/../generation-forte')
res = {}
# (1) binary forms of discriminant -39: reduced forms (a,b,c), b^2-4ac = -39, |b| <= a <= c, b >= 0 if |b| = a or a = c
forms = []
for a in range(1, 40):
    for b in range(-a + 1, a + 1):
        if (b*b + 39) % (4*a): continue
        c = (b*b + 39)//(4*a)
        if c < a or (c == a and b < 0): continue
        forms.append((a, b, c))
def dform_values(a, b, c):
    G = np.array([[2*a, b], [b, 2*c]]); D = 39
    Gi = [[F(2*c, D), F(-b, D)], [F(-b, D), F(2*a, D)]]
    vals = {}
    for x in range(D):
        for y in range(D):
            v = (Gi[0][0]*x + Gi[0][1]*y, Gi[1][0]*x + Gi[1][1]*y)
            key = (v[0] % 1, v[1] % 1)
            if key not in vals:
                q = (v[0]*v[0]*2*a + 2*v[0]*v[1]*b + v[1]*v[1]*2*c)/2
                vals[key] = q % 1
    return sorted(str(t) for t in vals.values())
genera = {}
for f in forms: genera.setdefault(tuple(dform_values(*f)), []).append(f)
res['reduced_forms_disc_-39'] = [list(f) for f in forms]
res['genera (as lists of reduced forms; (2,1,5) and (2,-1,5) are improperly equivalent)'] = [[list(f) for f in g] for g in genera.values()]
print('reduced forms of discriminant -39:', forms); print('genera:', list(genera.values()))
# (2) block dimensions
def blocks_caseII():
    d = pickle.load(open(GF + '/exact_S3_basis.pkl', 'rb'))
    out = {0: {(0, 0): 1}, 1: {mu: len(v) for mu, v in d['U1b'].items() if v}, 2: {mu: len(v) for mu, v in d['U2b'].items() if v},
           3: {mu: len(v) for mu, v in d['S3b'].items() if v}}
    return out
def blocks_caseI():
    SCR = '/home/claude/work/scratch-v10/cas-I'
    out = {0: {(0, 0): 1}, 1: {(0, 0): 2, (2, 0): 1, (-2, 0): 1, (0, 2): 1, (0, -2): 1}}
    for n in (2, 3):
        z = np.load(SCR + '/b1_S%d_allblocks_p7999993.npz' % n)
        out[n] = {}
        for k in z.files:
            m1, m2 = map(int, k.split('_')); r = z[k].shape[0]
            if r: out[n][(m1, m2)] = r
    return out
CASES = {'II': (3, 13, blocks_caseII()), 'I': (1, 39, blocks_caseI())}
for c_, (k1, k2, bl) in CASES.items():
    print('case', c_, 'total dims', [sum(bl[n].values()) for n in range(4)])
def norm(mu, k1, k2): return F(mu[0]**2, 2*k1) + F(mu[1]**2, 2*k2)
def coset(mu, k1, k2):
    r1, r2 = mu[0] % (2*k1), mu[1] % (2*k2)
    s1, s2 = (r1 + k1) % (2*k1), (r2 + k2) % (2*k2)
    return min((r1, r2), (s1, s2))
P2 = [1, -2, -1, 2, 1]          # prod_{m>=1} (1-q^m)^2 up to q^4
def fseries(c_):
    k1, k2, bl = CASES[c_]
    mus = set(mu for n in bl for mu in bl[n])
    out = {}
    for mu in mus:
        assert (mu[0] - mu[1]) % 2 == 0
        h = norm(mu, k1, k2)/2
        if h > 3: continue
        D = [bl[n].get(mu, 0) for n in range(4)]
        Fq = [sum(P2[j]*D[n - j] for j in range(n + 1)) for n in range(4)]
        al = coset(mu, k1, k2)
        for n in range(4):
            x = n - h
            if x < 0:
                assert Fq[n] == 0, (c_, mu, n)
                continue
            out.setdefault(al, {}).setdefault(x, set()).add(Fq[n])
    bad = {str(a): {str(x): sorted(v) for x, v in d.items() if len(v) > 1} for a, d in out.items()}
    bad = {a: d for a, d in bad.items() if d}
    return {a: {x: next(iter(v)) for x, v in d.items()} for a, d in out.items()}, bad
f = {}
for c_ in CASES:
    f[c_], bad = fseries(c_)
    res['case %s: internal inconsistencies' % c_] = bad
    print('case', c_, 'cosets seen:', len(f[c_]), 'internal inconsistencies:', bad)
# discriminant groups and isometries
def group(c_):
    k1, k2, _ = CASES[c_]
    els = {}
    for m1 in range(-2*k1, 2*k1 + 1):
        for m2 in range(-2*k2, 2*k2 + 1):
            if (m1 - m2) % 2: continue
            a = coset((m1, m2), k1, k2)
            if a not in els: els[a] = (m1, m2)
    return els
G = {c_: group(c_) for c_ in CASES}
assert len(G['I']) == len(G['II']) == 39
def add(c_, a, b):
    k1, k2, _ = CASES[c_]; x, y = G[c_][a], G[c_][b]
    return coset((x[0] + y[0], x[1] + y[1]), k1, k2)
def mult(c_, n, a):
    r = coset((0, 0), *CASES[c_][:2])
    for _ in range(n % 39): r = add(c_, r, a)
    return r
def q(c_, a): k1, k2, _ = CASES[c_]; return (norm(G[c_][a], k1, k2)/2) % 1
gII = next(a for a in G['II'] if len({mult('II', n, a) for n in range(39)}) == 39)
isos = []
for gI in G['I']:
    if len({mult('I', n, gI) for n in range(39)}) != 39: continue
    psi = {mult('II', n, gII): mult('I', n, gI) for n in range(39)}
    if all(q('II', a) == q('I', psi[a]) for a in psi): isos.append(psi)
print('isometries of discriminant forms II -> I:', len(isos))
res['number_of_isometries'] = len(isos)
comp = []
for i, psi in enumerate(isos):
    agree, disagree, compared = 0, [], 0
    for a, d in f['II'].items():
        b = psi[a]
        if b not in f['I']: continue
        for x, v in d.items():
            if x in f['I'][b]:
                compared += 1
                if f['I'][b][x] == v: agree += 1
                else: disagree.append((str(a), str(b), str(x), v, f['I'][b][x]))
    comp.append({'isometry': i, 'compared': compared, 'agree': agree, 'disagreements': disagree[:10]})
    print('isometry', i, 'compared', compared, 'agree', agree, 'disagree', len(disagree))
res['cross_case_comparison'] = comp
res['f_series'] = {c_: {str(a): {str(x): v for x, v in sorted(d.items())} for a, d in f[c_].items()} for c_ in f}
json.dump(res, open(HERE + '/s20_genus_check.json', 'w'), indent=1)
