"""Axis A3 (exact over Q(sqrt13)): is there an automorphism of V' normalizing A = U_1 mapping the Griess subalgebra 9
to 10?  Such an automorphism preserves B = Com_{V'}(A)_2, e = omega_C/2, the product and the form, hence induces an
algebra automorphism of B mapping the idempotent f of charge 32/19 in span{e, x_9} to the one in span{e, x_10}.  The
characteristic polynomial of the multiplication L(f) on B (dim 5) is an invariant.  Idempotents: f = a e + b x with
x.x = lambda x + mu e, b^2 D = 1, a = (1 - b lambda)/2, D = lambda^2 + 4 mu, charge c_C a.  Data from griess_cache.py."""
import json
import sympy as sp
d = json.load(open('griess_B_II.json'))
n = len(d['gram'])
prod = {tuple(map(int, k.split(','))): [sp.sympify(c) for c in v] for k, v in d['prod'].items()}
gram = sp.Matrix([[sp.sympify(c) for c in r] for r in d['gram']])
def mult(a, b):
    out = [0]*n
    for i in range(n):
        if a[i] == 0: continue
        for j in range(n):
            if b[j] == 0: continue
            c = a[i]*b[j]
            for k in range(n): out[k] += c*prod[(i, j)][k]
    return [sp.radsimp(sp.expand(x)) for x in out]
G = json.load(open('griess_subalgebras.json'))['II']['subalgebras']
t = sp.symbols('t'); res = {}
cC = sp.Rational(18, 5)
for k in (9, 10):
    x = [0] + [sp.sympify(c) for c in G[k]['x']]
    e = [1] + [0]*(n - 1)
    xx = mult(x, x)
    lam = sp.sympify(G[k]['lambda']); mu = sp.sympify(G[k]['mu']); D = sp.sympify(G[k]['D'])
    assert all(sp.simplify(xx[i] - lam*x[i] - mu*e[i]) == 0 for i in range(n)), 'x.x relation'
    sD = sp.sqrt(sp.radsimp(D)); sD = sp.sqrtdenest(sD)
    out = []
    for s in (1, -1):
        b = s/sD; a = (1 - b*lam)/2
        f = [sp.radsimp(sp.expand(a*e[i] + b*x[i])) for i in range(n)]
        ff = mult(f, f)
        assert all(sp.simplify(ff[i] - f[i]) == 0 for i in range(n)), 'not idempotent'
        ch = sp.nsimplify(sp.radsimp(cC*a))
        Lf = sp.Matrix([[mult(f, [1 if i == j else 0 for i in range(n)])[r] for j in range(n)] for r in range(n)])
        cp = sp.Poly(sp.expand((Lf - t*sp.eye(n)).det()), t)
        cpc = [sp.radsimp(sp.nsimplify(c)) for c in cp.all_coeffs()]
        ev = [sp.radsimp(z) for z in sp.Poly(sum(c*t**(n - i) for i, c in enumerate(cpc)), t).all_roots()] if False else None
        out.append({'charge': str(ch), 'charpoly_coeffs': [str(c) for c in cpc]})
        print('subalgebra', k, 'charge', ch, ': charpoly of L(f) =', cpc, flush=True)
    res[k] = out
json.dump(res, open('a3_griess_invariants.json', 'w'), indent=1)
