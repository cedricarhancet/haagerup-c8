"""Cross-check the exported quiver, order and block tables of Article III."""
import json
import lzma
from fractions import Fraction
from pathlib import Path
import pickle

import gmpy2
from sympy import isprime

ROOT = Path(__file__).resolve().parent / "Calculs" / "zhu-exact"
PRIME = 8388571


class RationalUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if module == "gmpy2.gmpy2" and name == "from_binary":
            return gmpy2.from_binary
        raise pickle.UnpicklingError(f"Unexpected global {module}.{name}")


def read_pickle(name):
    path = ROOT / name
    open_file = lzma.open if path.suffix == ".xz" else Path.open
    with open_file(path, "rb") as stream:
        return RationalUnpickler(stream).load()


def check(test, what):
    if not test:
        raise AssertionError(what)


def shape_of(word, vertices, ep):
    if not word:
        check(ep[0] == ep[1], "a vertex term in a non-loop corner")
        return ep
    check(all(vertices[a]["source"] == vertices[b]["target"] for a, b in zip(word, word[1:])),
          "noncomposable word")
    return vertices[word[0]]["target"], vertices[word[-1]]["source"]


def integral(value):
    return all(q.denominator % PRIME != 0 for q in value)


def main():
    check(isprime(PRIME) and pow(1649360, 2, PRIME) == 13
          and pow(2036966, 2, PRIME) == PRIME - 3, "chosen prime and field embedding")
    raw = read_pickle("rows_exact.pkl.xz")
    sys = read_pickle("upper_exact.pkl")
    arrows = raw["ARROWS"]
    live = set(sys["verts"])
    active = set(sys["active"])
    elim = sys["alllin"]
    rules = sys["rules"]
    dead = {a for a, arrow in enumerate(arrows)
            if arrow["source"] not in live or arrow["target"] not in live}
    check(len(arrows) == 515 and len(raw["ROWS"]) == 19981, "counts of arrows and relations")
    check(active == set(range(len(arrows))) - dead - set(elim)
          and live == set(range(12)) and len(active) == 239
          and len(set(elim) - dead) == 46, "partition of vertices and arrows")
    for endpoints, _, row in raw["ROWS"]:
        check(all(shape_of(word, arrows, endpoints) == endpoints and integral(value)
                  for word, value in row.items()), "projected row endpoints and integrality")
    for a, rhs in elim.items():
        ep = arrows[a]["target"], arrows[a]["source"]
        check(all(len(word) <= 1 and (not word or word[0] > a)
                  and shape_of(word, arrows, ep) == ep and integral(value)
                  for word, value in rhs.items()), "triangular elimination")
    for lead, rhs in rules.items():
        check(set(lead) <= active and len(lead) in (2, 3), "active leading word")
        ep = shape_of(lead, arrows, None)
        check(all(set(word) <= active and (-len(lead), lead) < (-len(word), word)
                  and shape_of(word, arrows, ep) == ep and integral(value)
                  for word, value in rhs.items()), "rule order, endpoints, integrality and active support")
    quadratic = sum(len(w) == 2 for w in rules)
    cubic = sum(len(w) == 3 for w in rules)
    check((quadratic, cubic) == (5221, 438), "rewriting rule counts")

    special = json.loads((ROOT / "specfib/z4_upper_algebra_modp.json").read_text(encoding="utf-8"))
    constituent_sizes = [(r + 1) * (l + 1) for r, l in raw["CONS"][:12]]
    sizes = []
    top_dimensions = []
    for block in special["blocks"]:
        ranks = {int(k): v for k, v in block["vertex_ranks"].items()}
        check(set(ranks) <= live and sum(ranks.values()) == block["size"], "block vertex ranks")
        sizes.append(block["size"])
        top_dimensions.append(sum(constituent_sizes[v] * rank for v, rank in ranks.items()))
    check(sorted(sizes) == [1, 1, 2, 2, 2, 3, 5, 7, 11, 12, 12, 16]
          and sum(n * n for n in sizes) == 762, "corner block sizes")
    check(sorted(top_dimensions) == [1, 3, 5, 7, 13, 14, 27, 42, 80, 81, 81, 119]
          and sum(n * n for n in top_dimensions) == 36625, "full Zhu algebra dimensions")
    weights = [Fraction(w) for w in ("0", "1/3", "2/3", "1", "5/13", "6/13",
                                     "7/13", "8/13", "11/13", "15/13")]
    residues = {w.numerator * pow(w.denominator, -1, PRIME) % PRIME for w in weights}
    check(len(residues) == len(weights), "distinct reductions of ten weights")
    squared = {x * x % 13 for x in range(1, 13)}
    nonsquares = {2, 5, 6, 7, 8, 11}
    check({5 * x % 13 for x in squared} == nonsquares, "Galois orbit of nonresidue twists")
    print(json.dumps({
        "status": "passed", "projected_rows": len(raw["ROWS"]), "arrows": len(arrows),
        "quadratic_rules": quadratic, "cubic_rules": cubic,
        "corner_sizes": sorted(sizes), "top_dimensions": sorted(top_dimensions),
        "dimension_corner": sum(n*n for n in sizes),
        "dimension_Zhu": sum(n*n for n in top_dimensions),
        "distinct_weight_reductions": len(residues),
    }, indent=2))


if __name__ == "__main__":
    main()
