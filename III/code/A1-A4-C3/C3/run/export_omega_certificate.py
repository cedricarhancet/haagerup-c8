#!/usr/bin/env python3
"""Export exact data from the trusted original R3 supplement to portable JSON.

Usage: python3 export_omega_certificate.py /path/to/Calculs-HU-v04/R3 [output.json]
This exporter imports e8e8_voa for the original state-index convention and
unpickles three original supplement files.  Use only trusted source pickles.
It needs SymPy and gmpy2; the subsequent verifier needs only SymPy.
"""
import hashlib
import importlib
import json
import pickle
import sys
from pathlib import Path

import sympy as sp

from verify_omega_certificate import algebraic_data


def rational(value):
    return str(sp.Rational(int(value.numerator), int(value.denominator)))


def main(source, destination):
    source = Path(source).resolve()
    sys.path.insert(0, str(source))
    V = importlib.import_module('e8e8_voa')
    with (source/'G.pkl').open('rb') as f:
        G = pickle.load(f)
    with (source/'G_griess.pkl').open('rb') as f:
        D = pickle.load(f)
    with (source/'sl2x4.pkl').open('rb') as f:
        SL = pickle.load(f)

    def vector(v):
        return [[list(V.VECS[b]), [list(x) for x in mono], rational(c)]
                for (b, mono), c in sorted(v.items())]

    n = len(G['G'])
    sources = ['G.pkl', 'G_griess.pkl', 'sl2x4.pkl', 'e8e8_voa.py',
               'common.py', 'r3a_G.py', 'r3a_Griess.py', 'r3a_groebner.py',
               'sl2x4.py', 'r3b_chi.py']
    data = {
        'format': 'haagerup-omega-certificate-v1',
        'scope': ('Finite algebraic consequences of supplied Griess structure constants and '
                  'multiplicity-space actions. Does not recompute lattice VOA modes or the '
                  'full commutant kernel. Exported from original first implementation R3.'),
        'source_sha256': {name: hashlib.sha256((source/name).read_bytes()).hexdigest() for name in sources},
        'state_convention': {
            'entry': '[doubled_lattice_vector, oscillator_monomial, rational_coefficient]',
            'doubled_lattice_vector': '16 integers y; actual vector is y/2 in E8 plus E8',
            'oscillator_monomial': 'list of [n,i] pairs representing eps_i(-n), indices 0 through 15',
            'basis_order': 'basis_G[j] is g_j, with j from 0 to 15; u=sum_j x_j*g_j',
            'cocycle': 'epsilon(a,b)=(-1)^(c(a)^T U c(b)) in each E8; simple roots and U are below',
            'simple_roots_doubled': [list(a) for a in V.E8_SIMPLE],
            'cocycle_U': V.U8,
        },
        'basis_G': [vector(g) for g in G['G']],
        'sl2_triples': {name: [vector(x) for x in triple] for name, triple in SL.items()},
        'griess_product': [[[rational(x) for x in D['prod'][(i, j)]]
                            for j in range(n)] for i in range(n)],
        'bilinear_form': [[rational(x) for x in row] for row in D['form']],
        'lowest_weights': {'c1':'1/3', 'c2':'2/3', 'd1':'6/13', 'd5':'7/13', 'd4':'5/13', 'd6':'8/13'},
        'multiplicity_spaces': [
            {'highest_weight': list(k), 'sectors': D['comp'][k],
             'actions': [[[rational(x) for x in row] for row in M] for M in D['Mk'][k]]}
            for k in sorted(D['comp'])],
        'tensor_convention': {
            'griess_product': '[i][j][k] = coefficient of g_k in (g_i)_(1) g_j',
            'bilinear_form': '[i][j] = vacuum coefficient of (g_i)_(3) g_j',
            'actions': 'actions[g][i][j] = row i, column j of (g_g)_(1) on the indicated highest-weight multiplicity space',
        },
    }
    ts, u, A, b, eqs, e2_count = algebraic_data(data)
    c = u.subs(dict.fromkeys(ts, 0))
    N = u.jacobian(ts)
    data['affine_coordinates'] = {
        'definition': 'x=c+N*t, in the order g_0,...,g_15 and t_0,...,t_7',
        'constant': [str(x) for x in c],
        'linear': [[str(N[i, j]) for j in range(8)] for i in range(n)],
        'free_coordinate_indices': [i for i in range(n) if c[i] == 0 and sum(N[i, j] != 0 for j in range(8)) == 1
                                    and any(N[i, j] == 1 for j in range(8))],
    }
    data['quadratic_equations'] = [[[list(m), str(c)] for m, c in sp.Poly(e, *ts).terms()] for e in eqs]
    Path(destination).write_text(json.dumps(data, separators=(',', ':'))+'\n')
    print('Exported', destination)
    print('Affine coordinates, x_j:')
    for j, entry in enumerate(u):
        print(j, entry)
    print('Input: 16 basis vectors, 20 multiplicity spaces, trace rank', A.rank(),
          ',', e2_count, 'e2 conditions,', len(eqs), 'nonzero quadrics.')


if __name__ == '__main__':
    main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else Path(__file__).with_name('omega_certificate.json'))
