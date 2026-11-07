"""R3a step 2e (certification of the candidate list): exact Groebner basis over Q (python-flint, Buchberger) of the
polynomial system satisfied by ANY omega_{C x 1} of an embedding S x S in W with (S x S)_1 = our sl2^4:
   u = u_0 + N t (affine space of the trace conditions), u.u = 2u in G, (u,u) = 9/5, and on each sl2^4-component of
   W_1 with two non-S sectors: e_2(o(u)) = lambda_z lambda_z' (2 x 2 minors sum; the S-line is killed by every o(g)).
The ideal is expected to be zero-dimensional with exactly the 4 points of omegaC1.pkl.
Output: groebner.json"""
import sys, os; sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import pickle, json, time, itertools
import sympy as sp
import flint
import common as CM
t0 = time.time()
D = pickle.load(open(os.path.join(CM.HERE, 'G_griess.pkl'), 'rb'))
prod, form, Mk, comp = D['prod'], D['form'], D['Mk'], D['comp']
n = len(form)
Qr = lambda x: sp.Rational(int(x.numerator), int(x.denominator))
rows = []; rhs = []
for kappa, secs in sorted(comp.items()):
    rows.append([sum(Qr(Mk[kappa][g][i][i]) for i in range(len(secs))) for g in range(n)])
    rhs.append(sp.Rational(str(sum(CM.omega_C1_eig(kappa, s) for s in secs))))
sol, params = sp.Matrix(rows).gauss_jordan_solve(sp.Matrix(rhs))
params = list(params); m = len(params)
ts = sp.symbols('t0:%d' % m)
u = list(sol.subs(dict(zip(params, ts))))
P = {k: [Qr(x) for x in v] for k, v in prod.items()}
uu = [0]*n
for i in range(n):
    if u[i] == 0: continue
    for j in range(n):
        if u[j] == 0: continue
        c = sp.expand(u[i]*u[j])
        for k in range(n):
            if P[(i, j)][k]: uu[k] += c*P[(i, j)][k]
eqs = [sp.expand(uu[k] - 2*u[k]) for k in range(n)]
eqs.append(sp.expand(sum(Qr(form[i][j])*u[i]*u[j] for i in range(n) for j in range(n) if form[i][j]) - sp.Rational(9, 5)))
for kappa, secs in sorted(comp.items()):
    nonS = [z for z in secs if z not in ('S1x1', '1xS1')]
    if len(nonS) != 2: continue
    k = len(secs)
    M = sp.Matrix(k, k, lambda i, j: sp.expand(sum(Qr(Mk[kappa][g][i][j])*u[g] for g in range(n) if Mk[kappa][g][i][j])))
    lam = [CM.omega_C1_eig(kappa, z) for z in nonS]
    target = sp.Rational(str(lam[0]*lam[1]))
    e2 = sum(M.extract(list(c), list(c)).det() for c in itertools.combinations(range(k), 2))
    eqs.append(sp.expand(e2 - target))
eqs = [e for e in eqs if e != 0]
print('variables', m, '; equations', len(eqs), '(%.0fs)' % (time.time() - t0), flush=True)
ctx = flint.fmpz_mpoly_ctx.get(['t%d' % i for i in range(m)], 'degrevlex')
polys = []
for e in eqs:
    num, den = sp.fraction(sp.together(e))
    pn = sp.Poly(sp.expand(num), *ts)
    terms = {tuple(mon): int(c) for mon, c in pn.terms()}
    polys.append(ctx.from_dict(terms))
V = flint.fmpz_mpoly_vec(polys, ctx)
GB = V.buchberger_naive()
print('Groebner basis computed (%.0fs): %d elements ; is_groebner %s' % (time.time() - t0, len(GB), GB.is_groebner()), flush=True)
GBr = GB.autoreduction()
print('reduced basis: %d elements' % len(GBr))
out = [str(GBr[i]) for i in range(len(GBr))]
for g in out[:30]: print('   ', g[:300])
json.dump({'vars': m, 'equations': len(eqs), 'groebner': out}, open(os.path.join(CM.HERE, 'groebner.json'), 'w'), indent=1)
print('time %.0fs' % (time.time() - t0))
