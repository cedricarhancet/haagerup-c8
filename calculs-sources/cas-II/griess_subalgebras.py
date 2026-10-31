"""Axis A2: all two-dimensional subalgebras of the Griess algebra B = Com_{V'}(A)_2 containing e = omega_C/2.
Such a subalgebra is span(e, x) with x in e^perp (invariant form), and x.x in span(e, x), i.e. x is a
projective eigenvector of the quadratic map q(x) = (x.x)^perp on e^perp: q(x) = lambda x.  Then
x.x = lambda x + mu e with mu = <x,x>/<e,e>, D = lambda^2 + 4 mu, and if D != 0 the subalgebra is C+C,
spanned by orthogonal idempotents f, e - f, giving commuting Virasoro vectors 2f, 2(e-f) (Miyamoto) with
central charges c_C (1 -/+ lambda / sqrt(D)) / 2.  If D = 0 the subalgebra is C[eps]/eps^2."""
import sys, json, time
import sympy as sp
from gmpy2 import mpq
import lattice_voa as V, cases as C, commutant_griess as CG, linalg as LA
def griess_data(case):
    T = C.build(case); T['case'] = case
    Bv, _ = CG.annihilated(T, 2, (0, 0))
    wA = {}; V.vadd(wA, CG.sugawara(T, 1)); V.vadd(wA, CG.sugawara(T, 2))
    omC = V.clean(V.vadd(dict(V.OMEGA), wA, -1))
    e = {k: x/2 for k, x in omC.items()}
    # orthogonal complement of e in B
    gB = [V.form(e, b) for b in Bv]
    ee = V.form(e, e)
    perp = []
    for b, g in zip(Bv, gB):
        perp.append(V.clean(V.vadd(dict(b), e, -g/ee)))
    M, keys = LA.matrix_from_cols(perp); Rr, r = M.rref()
    piv = [next(j for j in range(M.ncols()) if Rr[i, j] != 0) for i in range(r)]
    perp = [perp[j] for j in piv]
    basis = [e] + perp
    n = len(basis)
    prod = {}
    for i in range(n):
        for j in range(i, n):
            p = V.clean(V.mode_vec(basis[i], 1, basis[j]))
            cc = LA.coords(basis, p); assert cc is not None
            prod[(i, j)] = prod[(j, i)] = [sp.Rational(int(c.p), int(c.q)) for c in cc]
    gram = [[sp.Rational(int(V.form(basis[i], basis[j]).numerator), int(V.form(basis[i], basis[j]).denominator)) for j in range(n)] for i in range(n)]
    return T, basis, prod, gram
def analyse(case):
    T, basis, prod, gram = griess_data(case)
    n = len(basis); m = n - 1
    cC = 8*gram[0][0]
    assert prod[(0, 0)] == [1] + [0]*m, 'e is not idempotent'
    for i in range(1, n):
        assert prod[(0, i)] == [0]*i + [1] + [0]*(n - 1 - i), 'e is not a unit'
    xs = sp.symbols('x1:%d' % n); lam = sp.Symbol('lam')
    # q(x) = perp-part of x.x
    q = [sp.expand(sum(xs[i-1]*xs[j-1]*prod[(i, j)][k] for i in range(1, n) for j in range(1, n))) for k in range(1, n)]
    mu_num = [sp.expand(sum(xs[i-1]*xs[j-1]*prod[(i, j)][0] for i in range(1, n) for j in range(1, n)))]
    sols = []
    # charts: x_t = 1, x_s = 0 for s < t
    for t in range(m):
        sub = {xs[s]: 0 for s in range(t)}; sub[xs[t]] = 1
        eqs = [sp.expand((q[k] - lam*xs[k]).subs(sub)) for k in range(m)]
        vars_ = [xs[s] for s in range(t + 1, m)] + [lam]
        Gb = sp.groebner([e_ for e_ in eqs if e_ != 0], *vars_, order='lex', domain='QQ')
        if list(Gb) == [1]:
            continue
        S = sp.solve(list(Gb), vars_, dict=True)
        for s in S:
            x = [sp.nsimplify(0)]*m
            for u in range(m):
                if u < t: x[u] = 0
                elif u == t: x[u] = 1
                else: x[u] = sp.simplify(s[xs[u]])
            L = sp.simplify(s[lam])
            sols.append((x, L))
    out = []
    ee = gram[0][0]
    for x, L in sols:
        xx = sp.simplify(sum(x[i-1]*x[j-1]*gram[i][j] for i in range(1, n) for j in range(1, n)))
        mu = sp.simplify(xx/ee)
        D = sp.simplify(L**2 + 4*mu)
        if D == 0:
            out.append({'x': [str(v) for v in x], 'lambda': str(L), 'mu': str(mu), 'D': '0', 'type': 'nilpotent'})
            continue
        r = sp.nsimplify(sp.simplify(L/sp.sqrt(D)))
        c1 = sp.nsimplify(sp.simplify(cC*(1 - r)/2)); c2 = sp.nsimplify(sp.simplify(cC*(1 + r)/2))
        out.append({'x': [str(v) for v in x], 'lambda': str(L), 'mu': str(mu), 'D': str(D),
                    'c': sorted([str(c1), str(c2)], key=lambda s: float(sp.sympify(s))),
                    'cnum': sorted([float(c1), float(c2)])})
    return cC, n, out
if __name__ == '__main__':
    res = {}
    for case in sys.argv[1:] or ['I', 'II']:
        t0 = time.time()
        cC, n, out = analyse(case)
        print('case', case, ': c_C =', cC, '; dim B =', n, '; two-dimensional subalgebras containing e:', len(out))
        for o in out:
            print('   ', o.get('type', 'C+C'), ' central charges', o.get('c'), '~', o.get('cnum'))
        print('   time %.1fs' % (time.time() - t0))
        res[case] = {'c_C': str(cC), 'dimB': n, 'subalgebras': out}
    json.dump(res, open('griess_subalgebras.json', 'w'), indent=1)
