"""Portable replay of the 310 recorded products; no original setup needed."""
import argparse, importlib.util, json
from pathlib import Path
from fractions import Fraction as F
from generator_modes import product, mul, Z
from checked_modes import checked_mode

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--engine", type=Path, required=True)
    p.add_argument("--inputs", type=Path, required=True)
    p.add_argument("--record", type=Path, required=True)
    p.add_argument("--output", type=Path, required=True)
    args = p.parse_args()
    spec = importlib.util.spec_from_file_location("reference_lattice", args.engine)
    engine = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(engine)
    raw = json.loads(args.inputs.read_text())
    inputs = {name: {(tuple(a), tuple(map(tuple, mono))): tuple(F(x) for x in coef)
                     for a, mono, coef in rows} for name, rows in raw.items()}
    record = json.loads(args.record.read_text())
    count = 0
    for row in record["records"]:
        a, b = row["name"].split(" / ")
        u, v, n = inputs[a], inputs[b], row["mode"]
        expected = {}
        for ku, cu in u.items():
            for kv, cv in v.items():
                c = mul(cu, cv)
                for s, r in checked_mode(engine, ku, n, kv).items():
                    expected[s] = tuple(x+F(str(r))*y for x, y in zip(expected.get(s, Z), c))
        expected = {s: c for s, c in expected.items() if any(c)}
        actual = product(u, n, v)
        assert expected == actual
        assert len(actual) == row["support"]
        count += 1
    u = ((2, 0, 0, 0, 0, 0, 0, 0), ())
    v = ((1, 0, 0, 0, 0, 0, 0, 0), ())
    try:
        checked_mode(engine, u, -5, v)
    except ValueError:
        pass
    else:
        raise AssertionError("Guard failed to reject the output of weight 9")
    result = {"status": "PASS", "replayed_products": count,
              "all_reference_calls_within_weight_window": True,
              "outside_window_rejected": True}
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(result)

if __name__ == "__main__":
    main()
