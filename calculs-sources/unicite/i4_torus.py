"""Axis I: torus element exp(h(0)) (a character t of L) matching the transported currents E with ours:
t(p) c_Jiao(p) = c_ours(p) on the supports of E' and E''.  Ratios are positive rationals; solve in log-coordinates
over the primes involved (exactly, with sympy), then check F', F'' and the two seeds after the full map."""
import sys, os, json, pickle
HERE = os.path.abspath('.'); sys.path.insert(0, HERE); sys.path.insert(0, os.path.abspath('../generation-forte')); os.chdir('../generation-forte')
import sympy as sp
from gmpy2 import mpq
import lattice_voa as V, exact_engine as E, Kfield as K
TR = pickle.load(open(HERE + '/jiao_transported.pkl', 'rb'))
S = E.setup(0); T = S['T']
ours = {n: E.kv_from_terms([(K.ONE, T[n])]) for n in ('e1', 'f1', 'e2', 'f2')}
eqs = []
for a, b in (('PE', 'e1'), ('BE', 'e2')):
    for k, c in TR[a].items():
        d = ours[b][k]; assert not any(c[1:]) and not any(d[1:])
        eqs.append((k[0], sp.Rational(int(d[0].numerator), int(d[0].denominator)) / sp.Rational(int(c[0].numerator), int(c[0].denominator))))
primes = sorted({p for _, r in eqs for p in sp.factorint(abs(r.p) * abs(r.q))})
print('ratios', [(str(p), str(r)) for p, r in eqs], 'primes', primes)
assert all(r > 0 for _, r in eqs)
y = sp.symbols('y0:%d' % (8 * len(primes)))
sol_all = {}
system = []
for p, r in eqs:
    fr = sp.factorint(r.p); fq = sp.factorint(r.q)
    for j, pr in enumerate(primes):
        e = fr.get(pr, 0) - fq.get(pr, 0)
        system.append(sum(p[i] * y[8 * j + i] for i in range(8)) - e)
sol = sp.solve(system, y, dict=True)
print('solvable:', bool(sol))
sol = sol[0]
yv = [sol.get(s, 0).subs({s2: 0 for s2 in y}) for s in y]       # free parameters set to 0
def t(p):   # value of the character on p: prod_primes pr^(sum p_i y_{j,i}) (exponents may be rational -> radicals)
    val = sp.Integer(1)
    for j, pr in enumerate(primes):
        ex = sum(p[i] * yv[8 * j + i] for i in range(8)); val *= sp.Integer(pr) ** ex
    return sp.nsimplify(val)
def kfrom(c):   # K element -> sympy
    return sum(sp.Rational(int(x.numerator), int(x.denominator)) * b for x, b in zip(c, (1, sp.sqrt(13), sp.sqrt(-3), sp.sqrt(13) * sp.sqrt(-3))))
def mapped(name):
    return {k: sp.simplify(t(k[0]) * kfrom(c)) for k, c in TR[name].items()}
res = {'character_exponents': [str(v) for v in yv], 'primes': primes}
for a, b in (('PE', 'e1'), ('BE', 'e2'), ('PF', 'f1'), ('BF', 'f2')):
    m = mapped(a); o = {k: kfrom(c) for k, c in ours[b].items()}
    res[a + ' -> ' + b] = (set(m) == set(o)) and all(sp.simplify(m[k] - o[k]) == 0 for k in m)
for s in ('U26', 'U17'):
    m = mapped(s)
    res[s + ' coefficients'] = {str(k[0]): str(v) for k, v in m.items()}
print(json.dumps(res, indent=1))
json.dump(res, open(HERE + '/i4_torus.json', 'w'), indent=1)
pickle.dump({'yv': [str(v) for v in yv], 'primes': primes}, open(HERE + '/i4_torus.pkl', 'wb'))
