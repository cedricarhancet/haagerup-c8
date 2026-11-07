#!/usr/bin/env python3
"""Verify the finite algebraic stage of Lemma omega-C and its PCT check.

Usage: python3 verify_omega_certificate.py [omega_certificate.json]
Dependency: SymPy (tested with 1.14.0); no pickle, gmpy2, FLINT or VOA engine.
The JSON supplies the exact Griess tensor and multiplicity-space actions.
This verifies their algebraic consequences, not their derivation from modes,
nor that the displayed 16 independent vectors exhaust the commutant kernel.
"""
import itertools
import json
import sys
from pathlib import Path

import sympy as sp


def q(value):
    return sp.Rational(value)


def vec(entries):
    return {(tuple(e[0]), tuple(tuple(x) for x in e[1])): q(e[2]) for e in entries}


def clean(v):
    return {k: sp.expand(x) for k, x in v.items() if sp.expand(x) != 0}


def add(v, w, c=1):
    ans = dict(v)
    for k, x in w.items():
        ans[k] = ans.get(k, 0) + c*x
    return clean(ans)


def phi(v):
    ans = {}
    for (y, mono), c in v.items():
        exponent = -sp.Rational(y[6] + y[7] + y[14] + y[15], 2)
        assert exponent.is_integer
        ans[(tuple(-a for a in y), mono)] = (
            (-1)**len(mono) * sp.Integer(2)**exponent * sp.conjugate(c))
    return clean(ans)


def eigenvalue(component, sector, heights):
    if sector in ('S1x1', '1xS1'):
        return sp.Integer(0)
    a, b, _, _ = component
    return heights[sector] - sp.Rational(a*(a+2), 20) - sp.Rational(b*(b+2), 60)


def algebraic_data(data):
    """Rebuild the trace system and quadratic equations from tensors."""
    n = len(data['basis_G'])
    ts = sp.symbols('t0:8')
    heights = {k: q(v) for k, v in data['lowest_weights'].items()}
    rows, rhs = [], []
    blocks = []
    for block in data['multiplicity_spaces']:
        component, sectors = block['highest_weight'], block['sectors']
        mats = [sp.Matrix([[q(x) for x in row] for row in M]) for M in block['actions']]
        assert len(mats) == n
        assert all(M.shape == (len(sectors), len(sectors)) for M in mats)
        lam = [eigenvalue(component, z, heights) for z in sectors]
        rows.append([M.trace() for M in mats])
        rhs.append(sum(lam))
        blocks.append((sectors, mats, lam))
    A, b = sp.Matrix(rows), sp.Matrix(rhs)
    sol, parameters = A.gauss_jordan_solve(b)
    assert len(parameters) == len(ts)
    u = sp.Matrix(sol.subs(dict(zip(parameters, ts))))
    prod = [[[q(c) for c in row] for row in plane] for plane in data['griess_product']]
    form = sp.Matrix([[q(c) for c in row] for row in data['bilinear_form']])
    eqs = [sp.expand(sum(prod[i][j][k]*u[i]*u[j] for i in range(n)
                         for j in range(n)) - 2*u[k]) for k in range(n)]
    eqs.append(sp.expand((u.T*form*u)[0] - sp.Rational(9, 5)))
    e2_count = 0
    for sectors, mats, lam in blocks:
        positions = [i for i, z in enumerate(sectors) if z not in ('S1x1', '1xS1')]
        if len(positions) != 2:
            continue
        M = sum((u[g]*mats[g] for g in range(n)), sp.zeros(len(sectors)))
        e2 = sum(M.extract(list(c), list(c)).det()
                 for c in itertools.combinations(range(len(sectors)), 2))
        eqs.append(sp.expand(e2 - lam[positions[0]]*lam[positions[1]]))
        e2_count += 1
    return ts, u, A, b, [e for e in eqs if e != 0], e2_count


def polynomial(terms, ts):
    return sum(q(c)*sp.prod(x**k for x, k in zip(ts, mon)) for mon, c in terms)


def main(path):
    data = json.loads(Path(path).read_text())
    assert data['format'] == 'haagerup-omega-certificate-v1'
    basis = [vec(x) for x in data['basis_G']]
    n = len(basis)
    assert n == 16
    support = sorted(set().union(*(set(v) for v in basis)))
    B = sp.Matrix([[v.get(k, 0) for v in basis] for k in support])
    assert B.rank() == 16
    for y, mono in support:
        assert len(y) == 16
        assert sp.Rational(sum(a*a for a in y), 8) + sum(k for k, i in mono) == 2
        assert all(k > 0 and 0 <= i < 16 for k, i in mono)
    print('16 linearly independent explicit weight-two vectors supplied.')

    ts, u, A, b, eqs, e2_count = algebraic_data(data)
    assert A.shape == (20, 16) and A.rank() == 8
    assert e2_count == 10 and len(eqs) == 25
    affine = data['affine_coordinates']
    c = sp.Matrix([q(x) for x in affine['constant']])
    N = sp.Matrix([[q(x) for x in row] for row in affine['linear']])
    assert N.shape == (16, 8)
    assert u == c + N*sp.Matrix(ts)
    assert A*c == b and A*N == sp.zeros(20, 8)
    assert [sp.expand(polynomial(p, ts)-e) for p, e in zip(data['quadratic_equations'], eqs)] == [0]*len(eqs)
    assert len(data['quadratic_equations']) == len(eqs)
    print('20 trace equations have rank 8; the affine parametrization is exact.')
    print('The 25 nonzero quadrics match the supplied multiplication, form and 10 e2 conditions.')

    t0, t1, t2, t3, t4, t5, t6, t7 = ts
    claimed = [39*t0+1, 195*t1+8, t3, t5,
               t7+26*t2*t4, 2*t6-t7, 1014*t2*t2+1, 507*t4*t4+2]
    GI = sp.groebner(eqs, *ts, order='grevlex', domain=sp.QQ)
    GJ = sp.groebner(claimed, *ts, order='grevlex', domain=sp.QQ)
    assert GI.exprs == GJ.exprs
    assert GI.is_zero_dimensional
    print('Exact equality of ideals verified by equal reduced Groebner bases over QQ.')
    for p in GI.exprs:
        print('  ', p)
    # Eliminate the six determined variables.  The quotient is
    # Q[t2,t4]/(1014*t2^2+1, 507*t4^2+2), basis 1,t2,t4,t2*t4.
    assert sp.gcd(1014*t2**2+1, 2028*t2) == 1
    assert sp.gcd(507*t4**2+2, 1014*t4) == 1
    print('Quotient dimension 4; the two univariate factors are separable, so the quotient is reduced.')

    s = sp.sqrt(6)*sp.I
    points = []
    for a, d in itertools.product((-1, 1), repeat=2):
        p = [-sp.Rational(1, 39), -sp.Rational(8, 195), a*s/78, 0,
             d*s/39, 0, sp.Rational(a*d, 39), sp.Rational(2*a*d, 39)]
        assert all(sp.expand(e.subs(dict(zip(ts, p)))) == 0 for e in eqs)
        points.append(p)
        up = u.subs(dict(zip(ts, p)))
        vector = {}
        for g, coeff in zip(basis, up):
            vector = add(vector, g, coeff)
        assert phi(vector) == vector
    assert len({tuple(p) for p in points}) == 4
    print('All four geometric points lie in Q(sqrt(-6)) and are fixed by the displayed PCT involution.')

    for name, triple in data['sl2_triples'].items():
        E, H, F = [vec(x) for x in triple]
        assert phi(E) == clean({k: -c for k, c in F.items()})
        assert phi(H) == clean({k: -c for k, c in H.items()})
    print('The displayed PCT sends E_k to -F_k and H_k to -H_k for all four supplied triples.')
    print('PASS (algebraic stage only; mode-derived input tensors and commutant dimension remain separate).')


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).with_name('omega_certificate.json'))
