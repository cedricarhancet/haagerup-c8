#!/usr/bin/env python3
"""Exact preprocessing checks for the identification proof.

Only the Python standard library is used. This checks the fourteen radical
coefficients, the chosen residue roots and the elementary nonvanishing in
case (I). It does not rerun the lattice membership or rank computations.
"""
from __future__ import annotations

import ast
from fractions import Fraction as F
import hashlib
import json
from math import isqrt
from pathlib import Path

RADICANDS = (13, 3, 5, -1)
K_BASIS = (0, 1, 10, 11)


def mono(mask, coefficient=1):
    return {mask: F(coefficient)} if coefficient else {}


def add(a, b):
    out = dict(a)
    for mask, coefficient in b.items():
        out[mask] = out.get(mask, F(0)) + coefficient
    return {mask: coefficient for mask, coefficient in out.items() if coefficient}


def mul(a, b):
    out = {}
    for mask_a, coefficient_a in a.items():
        for mask_b, coefficient_b in b.items():
            coefficient = coefficient_a * coefficient_b
            for bit, radicand in enumerate(RADICANDS):
                if (mask_a & mask_b) & (1 << bit):
                    coefficient *= radicand
            mask = mask_a ^ mask_b
            out[mask] = out.get(mask, F(0)) + coefficient
    return {mask: coefficient for mask, coefficient in out.items() if coefficient}


def inv_monomial(a):
    assert len(a) == 1, a
    mask, coefficient = next(iter(a.items()))
    square = coefficient * coefficient
    for bit, radicand in enumerate(RADICANDS):
        if mask & (1 << bit):
            square *= radicand
    return mono(mask, coefficient / square)


def radical(n):
    assert n > 0
    mask = 0
    for bit, radicand in enumerate(RADICANDS[:3]):
        if n % radicand == 0:
            mask |= 1 << bit
            n //= radicand
    assert n == 1, n
    return mono(mask)


def evaluate(node):
    if isinstance(node, ast.Constant) and isinstance(node.value, int):
        return mono(0, node.value)
    if isinstance(node, ast.Name) and node.id == 'I':
        return mono(8)
    if isinstance(node, ast.UnaryOp) and isinstance(node.op, ast.USub):
        return mul(mono(0, -1), evaluate(node.operand))
    if isinstance(node, ast.Call):
        assert isinstance(node.func, ast.Name) and node.func.id == 'sqrt'
        assert len(node.args) == 1 and not node.keywords
        n = ast.literal_eval(node.args[0])
        return radical(n)
    if isinstance(node, ast.BinOp):
        a, b = evaluate(node.left), evaluate(node.right)
        if isinstance(node.op, ast.Add):
            return add(a, b)
        if isinstance(node.op, ast.Sub):
            return add(a, mul(mono(0, -1), b))
        if isinstance(node.op, ast.Mult):
            return mul(a, b)
        if isinstance(node.op, ast.Div):
            return mul(a, inv_monomial(b))
    raise ValueError(ast.dump(node))


def parse(text):
    return evaluate(ast.parse(text, mode='eval').body)


def main():
    path = Path(__file__).with_name('jiao-i4-torus-input.json')
    raw = path.read_bytes()
    data = json.loads(raw)
    expected = {
        'U26 coefficients': [
            ('0', '-5/12', '-5/12', '0'),
            ('0', '1/3', '-1/3', '0'),
            ('5/3', '0', '-5/3', '0'),
            ('5/3', '0', '5/3', '0'),
            ('10/3', '0', '0', '-10/3'),
            ('-5/2', '0', '0', '-5/2'),
            ('-10', '0', '0', '0'),
        ],
        'U17 coefficients': [
            ('3', '0', '0', '0'),
            ('1/6', '0', '0', '-1/6'),
            ('1/8', '0', '0', '1/8'),
            ('-1', '0', '1', '0'),
            ('1', '0', '1', '0'),
            ('0', '-5', '-5', '0'),
            ('0', '4', '-4', '0'),
        ],
    }
    coefficients = {}
    for name, tuples in expected.items():
        scale = parse('sqrt(5)*I') if name.startswith('U17') else mono(0)
        assert len(data[name]) == len(tuples) == 7
        results = {}
        for (state, expression), numbers in zip(data[name].items(), tuples):
            exact = parse(expression)
            normalized = mul(exact, inv_monomial(scale))
            reconstructed = {mask: F(number) for mask, number in zip(K_BASIS, numbers) if F(number)}
            assert normalized == reconstructed, (state, expression)
            assert mul(reconstructed, scale) == exact
            results[state] = list(numbers)
        coefficients[name] = results
    p = 8388571
    roots = {'sqrt13': 1649360, 'sqrt-3': 2036966}
    assert roots['sqrt13'] ** 2 % p == 13
    assert roots['sqrt-3'] ** 2 % p == p - 3
    assert all(0 < root < p // 2 for root in roots.values())
    # Elementary primality check, sufficient at this size.
    assert all(p % d for d in range(2, isqrt(p) + 1))
    difference = 25 * 4403 - 19758
    assert difference == 90317 and difference > 0
    print(json.dumps({
        'input_sha256': hashlib.sha256(raw).hexdigest(),
        'verified_coefficients': 14,
        'basis': ['1', 'sqrt(13)', 'sqrt(-3)', 'sqrt(-39)'],
        'normalized_coefficients': coefficients,
        'prime': p,
        'residue_roots': roots,
        'case_I_radical_difference': difference,
        'membership_and_rank_rerun': False,
        'result': 'PASS',
    }, indent=2))


if __name__ == '__main__':
    main()
