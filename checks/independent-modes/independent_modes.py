"""Independent low-weight audit by simultaneous normal ordering.

The reference engine is used only for comparison and selection of sample states.
This evaluator does not call its pairing, cocycle, Schur, Heisenberg or mode code.
Run with --engine PATH; gmpy2 and sympy are reference-engine dependencies only.
"""
from fractions import Fraction as F
from functools import lru_cache
from math import comb, factorial
from collections import Counter
import argparse, hashlib, importlib.util, json, random, time
from pathlib import Path

ZERO = (0,) * 8
EDGES = {(0, 2), (1, 3), (2, 3), (3, 4), (4, 5), (6, 7)}

def gram(i, j):
    return 2 if i == j else -1 if tuple(sorted((i, j))) in EDGES else 0

def dot(a, b):
    return 2 * sum(x*y for x, y in zip(a, b)) - sum(
        a[i]*b[j] + a[j]*b[i] for i, j in EDGES)

def pairing(a, i):
    return sum(a[j]*gram(j, i) for j in range(8))

def cocycle(a, b):
    return (-1)**((sum(x*y for x, y in zip(a, b)) -
                   sum(a[i]*b[j] for i, j in EDGES)) % 2)

def weight(state):
    a, mono = state
    return dot(a, a)//2 + sum(n for n, i in mono)

def accumulate(out, key, c):
    if c:
        out[key] = out.get(key, F(0)) + c
        if not out[key]:
            del out[key]

@lru_cache(maxsize=8192)
def exponential(a, degree):
    """Direct multivariate exponential: product (a_i/n)^m / m!."""
    if degree < 0:
        return ()
    variables = [(n, i) for n in range(1, degree+1)
                 for i in range(8) if a[i]]
    out = []
    def rec(pos, remaining, mono, coef):
        if not remaining:
            out.append((mono, coef))
            return
        if pos == len(variables):
            return
        n, i = variables[pos]
        for mult in range(remaining//n + 1):
            rec(pos+1, remaining-mult*n, mono+((n, i),)*mult,
                coef * F(a[i], n)**mult / factorial(mult))
    rec(0, degree, (), F(1))
    return tuple(out)

def annihilate(i, mode, state):
    b, mono = state
    if mode == 0:
        q = pairing(b, i)
        return {state: F(q)} if q else {}
    out = {}
    for entry, multiplicity in Counter(mono).items():
        n, j = entry
        if n == mode and gram(i, j):
            rem = list(mono)
            rem.remove(entry)
            accumulate(out, (b, tuple(rem)), F(multiplicity*n*gram(i, j)))
    return out

def direct_mode(u, k, v):
    """Coefficient z^(-k-1) of :prod d^(n-1)h_i/(n-1)! Y(e^a,z):."""
    a, left = u
    b, right = v
    ab = tuple(x+y for x, y in zip(a, b))
    output_weight = weight(u) + weight(v) - k - 1
    bound = output_weight - dot(ab, ab)//2
    if bound < 0:
        return {}
    out = {}

    def finish(state, creators, power, coef):
        _, mono = state
        # E^+(a,z) acts by h_j(-s) -> h_j(-s) - (a,h_j) z^-s.
        terms = [((), 0, coef*cocycle(a, b))]
        for n, j in mono:
            expanded = []
            for rem, shift, cc in terms:
                expanded.append((rem+((n, j),), shift, cc))
                q = pairing(a, j)
                if q:
                    expanded.append((rem, shift-n, -q*cc))
            terms = expanded
        for rem, shift, cc in terms:
            degree = -k-1-power-dot(a, b)-shift
            for expmono, ec in exponential(a, degree):
                newmono = tuple(sorted(creators+rem+expmono))
                assert sum(n for n, i in newmono) == bound
                accumulate(out, (ab, newmono), cc*ec)

    def normal_order(pos, state, creators, power, coef, created_weight):
        if pos == len(left):
            finish(state, creators, power, coef)
            return
        n, i = left[pos]
        # Negative Heisenberg modes, placed to the left of Y(e^a,z).
        for r in range(n, bound-created_weight+1):
            normal_order(pos+1, state, creators+((r, i),), power+r-n,
                         coef*comb(r-1, n-1), created_weight+r)
        # Nonnegative modes, placed to its right.
        for m in sorted({0} | {s for s, j in state[1]}):
            for newstate, cc in annihilate(i, m, state).items():
                normal_order(pos+1, newstate, creators, power-m-n,
                             coef*cc*(-1)**(n-1)*comb(m+n-1, n-1),
                             created_weight)
    normal_order(0, v, (), 0, F(1), 0)
    return out

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--engine", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    started = time.time()
    spec = importlib.util.spec_from_file_location("reference_lattice", args.engine)
    ref = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ref)
    counts, nonzero = Counter(), Counter()
    def check(u, k, v, group):
        expected = {s: F(str(c)) for s, c in ref.mode(u, k, v).items() if c}
        actual = direct_mode(u, k, v)
        if actual != expected:
            raise AssertionError((group, u, k, v, expected, actual))
        counts[group] += 1
        nonzero[group] += bool(actual)
    basis1 = ref.basis(1)
    for u in basis1:
        for v in basis1:
            for target in range(3):
                check(u, 1-target, v, "all_weight_one_pairs")
    print("Weight-one products passed", dict(counts), flush=True)
    rng = random.Random(20261031)
    pools = {w: ref.basis(w) for w in range(1, 5)}
    # Explicitly oversample oscillator states and products retaining the charge.
    oscillators = {w: [x for x in states if x[1]] for w, states in pools.items()}
    for j in range(2500):
        wu, wv = rng.randrange(1, 5), rng.randrange(1, 5)
        u = rng.choice(oscillators[wu] if j % 2 else pools[wu])
        v = rng.choice(oscillators[wv] if j % 3 else pools[wv])
        if j % 4 == 0:
            v = (tuple(-t for t in u[0]), v[1])
            if weight(v) > 4:
                continue
        target = j % 5
        check(u, weight(u)+weight(v)-target-1, v, "stratified_weights_one_to_four")
    print("Stratified products passed", dict(counts), flush=True)
    # Out-of-window counterexample: alpha=2 alpha_0, beta=alpha_0.
    # (e^(2 alpha_0))_(-5)e^(alpha_0)=e^(3 alpha_0), of weight 9.
    u = ((2, 0, 0, 0, 0, 0, 0, 0), ())
    v = ((1, 0, 0, 0, 0, 0, 0, 0), ())
    actual = direct_mode(u, -5, v)
    expected_key = ((3, 0, 0, 0, 0, 0, 0, 0), ())
    assert actual == {expected_key: F(1)}
    assert ref.mode(u, -5, v) == {}
    result = {
        "status": "PASS", "seed": 20261031, "counts": dict(counts),
        "nonzero": dict(nonzero), "seconds": round(time.time()-started, 2),
        "engine_sha256": hashlib.sha256(args.engine.read_bytes()).hexdigest(),
        "audit_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "outside_window": {"input_weights": [4, 1], "mode": -5,
                           "output_weight": 9, "reference_silently_zero": True,
                           "independent_coefficient": "1"},
        "scope": "Exact sampled mode comparisons, not a regeneration of all F4 coordinates."
    }
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))

if __name__ == "__main__":
    main()
