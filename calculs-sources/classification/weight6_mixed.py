"""Mixed weight-six relations for U_1 = g' + g'' (coefficients of x^4 y^2 and x^2 y^4 in W6(x h' + y h'')).

Leading coefficient of the x^4 y^2 form on a lowest space R:
   15 Tr h'^4 h''^2 - (5/4)(Tr h'^4 + 6 Tr h'^2 h''^2) + (5/16)(Tr h''^2 + 2 Tr h'^2) - 15 dim R/576,
and on the vacuum: w0 = -15/576, w0p = -(5/4) T4' + (5/16)(T'' + 2 T') + 990/576.
Known sectors from e8_joint.joint_sums (explicit e8); d-sectors: Tr h'^2, h'^4, h''^2, h''^4, h'^2 h''^2 known
(single-ideal and mixed quartic steps), Tr h'^4 h''^2 =: P_j unknown.  The exact relations
(weight6_relations_exact_gamma0.json) give conditions; we report them for both candidate pairs.
"""
import json, sympy as sp
import e8_joint as J, sl2_pair_feasibility as S, weight6_single as W6

DIMS = [7, 42, 119, 5, 13, 14]
R6 = [[sp.Rational(c) for c in r] for r in json.load(open('weight6_relations_exact_gamma0.json'))['relations']]
labels = json.load(open('weight6_relations_exact_gamma0.json'))['labels']
def M(n, p):
    return sum((n - 2*i)**p for i in range(n + 1))

CANDS = {'A11+A139 (announced)': ((1, 0), (35, 4), (2, 1, 1, 0, 1, 2)),
         'A13+A113': ((3, 0), (9, 4), (1, 0, 0, 1, 0, 1))}

def analyse(name, p_key, q_key, lab):
    d1, d2 = S.single[p_key], S.single[q_key]          # d1 = g' (i2 = 0), d2 = g'' (diagonal)
    k1, k2 = d1['k'], d2['k']
    js = J.joint_sums(lab, p_key[0])
    sol4 = S.mixed_solution(d1, d2)
    st6 = [s for s in S.joint_structures(14, (d1['t2'][5], d1['t4'][5]), (d2['t2'][5], d2['t4'][5]))
           if max(a for a, b in s[0]) <= k1 and max(b for a, b in s[0]) <= k2]
    assert len(st6) == 1
    X6 = sp.Rational(st6[0][1], 4*k1*k2)
    Xn = [sp.Rational(sol4[S.X[j]].subs(S.X[5], X6)) for j in range(5)] + [X6]
    def norm(p, q):
        return sp.Integer(2*k1)**sp.Rational(p, 2)*sp.Integer(2*k2)**sp.Rational(q, 2)
    sec = {}
    for (p, q), v in js.items():
        sec[(p, q)] = {'b': (v['e6'] + v['sl3'] - v['s'])/norm(p, q), 'a': v['27x3']/norm(p, q), 'c0': v['27x3']/norm(p, q),
                       'c1': v['1x3']/norm(p, q), 'c2': v['27']/norm(p, q)}
    # consistency with the quartic data already used
    for z in ('b', 'a', 'c1', 'c2'):
        pass
    tau1 = [sp.Rational(t, 2*k1) for t in d1['t2']]; tau2 = [sp.Rational(t, 2*k2) for t in d2['t2']]
    Q1 = [sp.Rational(t, 4*k1*k1) for t in d1['t4']]; Q2 = [sp.Rational(t, 4*k2*k2) for t in d2['t4']]
    T1, T2 = sp.Rational(4, k1), sp.Rational(4, k2); T41, T42 = sp.Rational(8, k1*k1), sp.Rational(8, k2*k2)
    dims = {'b': 80, 'a': 81, 'c0': 81, 'c1': 3, 'c2': 27}
    P = sp.symbols('P1:7')
    results = {}
    for form, (pp, qq) in (('x4y2', (4, 2)), ('x2y4', (2, 4))):
        if form == 'x4y2':
            Qa, ta, tb, T4a, Ta, Tb = Q1, tau1, tau2, T41, T1, T2
        else:
            Qa, ta, tb, T4a, Ta, Tb = Q2, tau2, tau1, T42, T2, T1
        w = {'w0': sp.Rational(-15, 576), 'w0p': -sp.Rational(5, 4)*T4a + sp.Rational(5, 16)*(Tb + 2*Ta) + sp.Rational(990, 576)}
        for z, dz in dims.items():
            s6 = sec[(pp, qq)][z]; s4 = sec[(pp - 2, qq - 2) if False else (pp, 0)][z] if False else None
            # Tr h_a^4 on known sectors, Tr h'^2 h''^2, Tr h_a^2, Tr h_b^2
            q4 = sec[(4, 0)][z] if form == 'x4y2' else sec[(0, 4)][z]
            x22 = sec[(2, 2)][z]
            t2a = sec[(2, 0)][z] if form == 'x4y2' else sec[(0, 2)][z]
            t2b = sec[(0, 2)][z] if form == 'x4y2' else sec[(2, 0)][z]
            w['w_' + z] = 15*s6 - sp.Rational(5, 4)*(q4 + 6*x22) + sp.Rational(5, 16)*(t2b + 2*t2a) - sp.Rational(15*dz, 576)
        for j in range(6):
            w['w_d%d' % (j + 1)] = 15*P[j] - sp.Rational(5, 4)*(Qa[j] + 6*Xn[j]) + sp.Rational(5, 16)*(tb[j] + 2*ta[j]) - sp.Rational(15*DIMS[j], 576)
        eqs = [sp.expand(sum(c*w[l] for c, l in zip(row, labels))) for row in R6]
        Mx, b = sp.linear_eq_to_matrix(eqs, list(P))
        conds = [sp.nsimplify(sp.expand((kv.T*b)[0])) for kv in Mx.T.nullspace()]
        sol = sp.solve(eqs, list(P), dict=True)
        results[form] = {'rank': Mx.rank(), 'conditions': conds, 'solution': sol[0] if sol else None}
        print('%s, form %s: rank %d; consistency conditions (must vanish): %s' % (name, form, Mx.rank(), conds))
    return results, d1, d2, Xn

if __name__ == '__main__':
    out = {}
    for name, (p, q, lab) in CANDS.items():
        res, d1, d2, Xn = analyse(name, p, q, lab)
        out[name] = {f: {'rank': r['rank'], 'conditions': [str(c) for c in r['conditions']],
                         'solution': {str(k): str(v) for k, v in r['solution'].items()} if r['solution'] else None} for f, r in res.items()}
    json.dump(out, open('weight6_mixed.json', 'w'), indent=1)
