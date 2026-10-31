"""Re-evaluate selected actual generators with independent_modes.direct_mode.

Generator coefficients are read from the existing setup. Their extraction is
not independent; field arithmetic and mode evaluation are independent.
"""
import argparse, hashlib, json, os, sys, time
from pathlib import Path
from fractions import Fraction as F
from independent_modes import direct_mode, ZERO

Z = (F(0),)*4

def mul(a, b):
    out = [F(0)]*4
    for i, ai in enumerate(a):
        for j, bj in enumerate(b):
            factor = (13 if i & j & 1 else 1) * (-3 if i & j & 2 else 1)
            out[i ^ j] += factor*ai*bj
    return tuple(out)

def rationalize(v):
    return {s: tuple(F(str(x)) for x in c) for s, c in v.items() if any(c)}

def product(u, n, v):
    out = {}
    for ku, cu in u.items():
        for kv, cv in v.items():
            coefficient = mul(cu, cv)
            for state, c in direct_mode(ku, n, kv).items():
                out[state] = tuple(a+c*b for a, b in zip(out.get(state, Z), coefficient))
    return {s: c for s, c in out.items() if any(c)}

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--source", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    output = args.output.resolve()
    started = time.time()
    src = args.source.resolve()
    sys.path.insert(0, str(src))
    os.chdir(src)
    import exact_engine as E
    import Kfield as K
    setup = E.setup(0)
    hw = {name: rationalize(E.kv_from_terms(t)) for name, t in setup["hw"].items()}
    currents = [rationalize(u) for u in setup["U1k"]]
    count = 0
    records = []
    for a in E.BAS:
        for b in E.BAS:
            assert mul(a, b) == K.mul(a, b)
    def check(u, n, v, name):
        nonlocal count
        reference = rationalize(E.kprod(E.kv_terms(u), n, v))
        actual = product(u, n, v)
        assert actual == reference, (name, n)
        count += 1
        records.append({"name": name, "mode": n, "support": len(actual)})
    for i, u in enumerate(currents):
        for name, v in hw.items():
            for n in range(4):
                check(u, n, v, "current_"+str(i)+" / "+name)
    print("240 current actions passed", flush=True)
    for name, u in hw.items():
        for n in range(-1, 4):
            check(u, n, u, name+" / "+name)
    names = list(hw)
    for j, name in enumerate(names):
        for n in (0, 1):
            other = names[(j+1) % len(names)]
            check(hw[name], n, hw[other], name+" / "+other)
    omega = rationalize(E.kv_from_terms([(K.ONE, E.V.OMEGA)]))
    vac = {(ZERO, ()): (F(1), F(0), F(0), F(0))}
    assert product(omega, 3, omega) == {k: tuple(4*x for x in v) for k, v in vac.items()}
    assert product(omega, 2, omega) == {}
    assert product(omega, 1, omega) == {k: tuple(2*x for x in v) for k, v in omega.items()}
    encoded = {}
    for name, v in list(hw.items()) + [("current_"+str(i), v) for i, v in enumerate(currents)]:
        encoded[name] = [[list(k[0]), [list(p) for p in k[1]], [str(x) for x in c]]
                         for k, c in sorted(v.items())]
    input_path = output.with_name("generator-inputs.json")
    input_path.write_text(json.dumps(encoded, indent=2)+"\n")
    result = {"status": "PASS", "actual_generator_products": count,
              "field_basis_products": 16, "virasoro_expected_identities": 3,
              "seconds": round(time.time()-started, 2),
              "input_sha256": hashlib.sha256(input_path.read_bytes()).hexdigest(),
              "records": records,
              "limitation": "The input generator coefficients come from E.setup(0); only their products are independently reevaluated."}
    output.write_text(json.dumps(result, indent=2)+"\n")
    print({k: v for k, v in result.items() if k != "records"})

if __name__ == "__main__":
    main()
