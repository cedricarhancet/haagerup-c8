"""Exact characteristic polynomials of x_1 on Pi_tau for every candidate (factored over Q or Q(sqrt13))."""
import sys, json, pickle
import sympy as sp
import griess_subalgebras as GS, primaries as PR
def load(case):
    d = pickle.load(open(f'primaries_{case}.pkl', 'rb'))
    return {eval(t): {'X1': [eval(m, sp.__dict__) for m in v['X1']], 'X2': [eval(m, sp.__dict__) for m in v['X2']], 'n': len(v['P'])} for t, v in d['types'].items()}
t = sp.Symbol('t')
if __name__ == '__main__':
    case = sys.argv[1]; which = [int(a) for a in sys.argv[2:]]
    G = json.load(open('griess_subalgebras.json'))[case]['subalgebras']
    types = load(case)
    for ci in which:
        xs = [sp.sympify(c) for c in G[ci]['x']]
        print('candidate', ci, G[ci]['c'])
        for tau, d in types.items():
            n = d['n']
            X = sum((xs[i]*d['X1'][i + 1] for i in range(len(xs))), sp.zeros(n, n))
            cp = X.charpoly(t).as_expr()
            ext = [sp.sqrt(13)] if any(sp.sympify(c).has(sp.sqrt(13)) for c in G[ci]['x']) else []
            print('   ', tau, sp.factor(cp, extension=ext) if ext else sp.factor(cp))
