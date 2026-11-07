# Step 3c: solve the polynomial system exactly (linear elimination + Groebner basis)
import sys, pickle, time
sys.path.insert(0, '.')
import sympy as sp
t0 = time.time()
S = pickle.load(open('system.pkl', 'rb'))
lin, nl, xs = S['lin'], S['nl'], list(S['xs'])
sol = sp.solve(lin, xs, dict=True)
assert len(sol) == 1
sol = sol[0]
free = [x for x in xs if x not in sol]
print('eliminated by linear equations:', sorted(map(str, sol)), ' free:', free)
polys = []
for e in nl:
    p = sp.expand(e.subs(sol))
    if p != 0:
        polys.append(p)
polys = list(dict.fromkeys(polys))
print('nonlinear polys after substitution:', len(polys), ' degrees:', sorted(sp.Poly(p, *free).total_degree() for p in polys))
GB = sp.groebner(polys, *free, order='grevlex', domain='QQ')
print('grevlex GB computed: %d polys, %.1fs' % (len(GB.exprs), time.time()-t0), flush=True)
print('zero-dimensional:', GB.is_zero_dimensional)
# count standard monomials
lms = [sp.Poly(g, *free).monoms(order='grevlex')[0] for g in GB.exprs]
def divisible(m, l): return all(a >= b for a, b in zip(m, l))
# enumerate monomials not divisible by any leading monomial (BFS)
std = []
frontier = [tuple([0]*len(free))]
seen = set(frontier)
while frontier:
    m = frontier.pop()
    if any(divisible(m, l) for l in lms):
        continue
    std.append(m)
    for i in range(len(free)):
        mm = list(m); mm[i] += 1; mm = tuple(mm)
        if mm not in seen:
            seen.add(mm); frontier.append(mm)
    if len(std) > 10000: break
print('number of standard monomials (solutions with multiplicity):', len(std))
GBl = sp.groebner(GB.exprs, *free, order='lex', domain='QQ')
print('lex GB: %.1fs' % (time.time()-t0))
for g in GBl.exprs:
    print('   ', sp.factor(g))
pickle.dump({'sol_lin': sol, 'free': free, 'gb_grevlex': GB.exprs, 'gb_lex': GBl.exprs}, open('groebner.pkl', 'wb'))
print('time %.1fs' % (time.time()-t0))
