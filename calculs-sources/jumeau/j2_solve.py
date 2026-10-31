import sys, pickle, sympy as sp, time
CASE = sys.argv[1] if len(sys.argv) > 1 else 'II'
d = pickle.load(open('j1_equations_%s.pkl' % CASE, 'rb'))
syms = sorted([sp.Symbol(s) for s in d['syms']], key=str)
eqs = [sp.sympify(e) for e in d['eqs']]
# normalize: make primitive, dedupe
norm = set()
for e in eqs:
    p = sp.Poly(e, *syms)
    p = p.primitive()[1]
    lc = p.LC()
    p = sp.Poly(p.as_expr()/lc, *syms)
    norm.add(p.as_expr())
eqs = sorted(norm, key=lambda e: (sp.Poly(e, *syms).total_degree(), len(str(e))))
print('distinct', len(eqs), 'by degree', {dg: sum(1 for e in eqs if sp.Poly(e, *syms).total_degree() == dg) for dg in (1, 2, 3)})
for e in eqs[:40]: print('  ', e)
t = time.time()
G = sp.groebner(eqs, *syms, order='grevlex')
print('groebner grevlex', time.time() - t, 's, size', len(G.exprs))
for g in G.exprs[:60]: print('  ', g)
pickle.dump([str(g) for g in G.exprs], open('j2_groebner_%s.pkl' % CASE, 'wb'))
