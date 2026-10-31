#!/usr/bin/env python3
"""Check the primary-line reductions and an int64 accumulation bound."""
import argparse
import json
from pathlib import Path
import sys


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source", type=Path, nargs="?",
                        default=Path(__file__).resolve().parent.parent / "calculs-sources" / "generation-forte")
    args = parser.parse_args()
    source = args.source.resolve()
    sys.path.insert(0, str(source))
    import sympy as sp
    import Kfield as K
    import lattice_voa as V
    import modp

    field = modp.Fp(2**23)
    lines = json.loads((source / "isotropic_lines_II_9.json").read_text())
    choice = json.loads((source / "closure_explore_II_9.json").read_text())["combos"][0]
    keys = sorted(lines, key=lambda key: tuple(map(int, key.strip("()").split(","))))
    checked = 0
    for key, index in zip(keys, choice):
        line = sp.sympify(lines[key][0]["lines"][index])
        for coefficient in line:
            exact = K.from_sympy(coefficient)
            powers = [1, field.s13, field.s3, field.s13*field.s3]
            residue = sum(field.q(exact[i])*powers[i] for i in range(4)) % field.p
            if residue != field.alg(coefficient):
                raise ValueError(f"Incoherent reduction: {key}, {coefficient}")
            checked += 1
    if checked != 16:
        raise ValueError(f"Unexpected primary coefficient count: {checked}")
    dimensions = {str(w): len(V.basis(w)) for w in range(5)}
    bound = max(dimensions.values())*(field.p-1)**2
    if bound >= 2**63:
        raise ValueError("Potential int64 overflow in modular matrix products")
    if field.s13**2 % field.p != 13 or field.s3**2 % field.p != field.p-3:
        raise ValueError("Invalid residue-field roots")
    print(json.dumps({"status": "PASS", "p": int(field.p),
                      "sqrt13_image": int(field.s13), "sqrt_minus3_image": int(field.s3),
                      "chosen_primary_line_coefficients_checked": checked,
                      "ambient_dimensions": dimensions,
                      "global_dot_product_bound": int(bound), "int64_max": 2**63-1}, indent=2))


if __name__ == "__main__":
    main()
