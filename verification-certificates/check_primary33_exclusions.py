"""Exact exclusion of thirteen case-II Griess candidates using Pi_(3,3).

Usage: python3 check_primary33_exclusions.py SOURCE_DIRECTORY OUTPUT_JSON
Every calculation is over Q.  Projective charts are disjoint: the first
nonzero coordinate is set to one.  The equation z*P(u,u)-1 excludes
isotropic/zero vectors.  Unit Groebner bases on all three charts certify
that no eigenline satisfying all necessary equations exists.
"""
import hashlib
import json
import os
from pathlib import Path
import sys
import time

source = Path(sys.argv[1]).resolve()
output = Path(sys.argv[2]).resolve()
sys.path.insert(0, str(source))
os.chdir(source)
import sympy as sp
import griess_subalgebras as GS
import primaries as PR
import isotropy_test as IT

T, basis, product, gram = GS.griess_data("II")
data = PR.setup("II")
primary = data["types"][(3, 3)]
assert len(primary["P"]) == 3
P = IT.pairing(data["T"], (3, 3), primary["P"])
gram_matrix = sp.Matrix(gram)
griess = json.loads(Path("griess_subalgebras.json").read_text())["II"]["subalgebras"]
n = len(basis)
u = sp.symbols("u0:3")
v = sp.Matrix(u)
z = sp.Symbol("z")
result = {"coefficient_field": "Q", "primary_type": [3, 3], "candidates": {}}
for ci in range(15):
    if ci in (9, 10):
        continue
    start = time.monotonic()
    xs = [sp.sympify(c) for c in griess[ci]["x"]]
    assert all(c.is_Rational for c in xs)
    xv = sp.Matrix([0] + xs)
    ev = sp.Matrix([1] + [0] * (n - 1))
    Ys = sp.Matrix.vstack(ev.T * gram_matrix, xv.T * gram_matrix).nullspace()
    X = sum((xv[i] * primary["X1"][i] for i in range(n)), sp.zeros(3, 3))
    Yops = [sum((y[i] * primary["X1"][i] for i in range(n)), sp.zeros(3, 3)) for y in Ys]
    equations = [sp.expand((v.T * P * Y * v)[0]) for Y in Yops]
    Xv = X * v
    equations += [sp.expand(Xv[i] * u[j] - Xv[j] * u[i]) for i in range(3) for j in range(i + 1, 3)]
    norm = sp.expand((v.T * P * v)[0])
    charts = []
    for k in range(3):
        substitution = {u[j]: 0 for j in range(k)}
        substitution[u[k]] = 1
        polys = [sp.expand(e.subs(substitution)) for e in equations]
        polys.append(sp.expand(z * norm.subs(substitution) - 1))
        variables = list(u[k + 1:]) + [z]
        gb = sp.groebner(polys, *variables, order="grevlex", domain=sp.QQ)
        is_unit = list(gb) == [1]
        assert is_unit, (ci, k, list(gb))
        charts.append({"first_nonzero_coordinate": k, "variables": list(map(str, variables)),
                       "polynomials": list(map(str, polys)), "groebner_basis": ["1"]})
    result["candidates"][str(ci)] = {
        "x_coordinates": list(map(str, xs)), "X": str(X), "P": str(P),
        "Y_operators": list(map(str, Yops)), "charts": charts,
    }
    print("candidate", ci, "all three charts excluded over Q", round(time.monotonic() - start, 3), flush=True)
assert len(result["candidates"]) == 13
result["verified_projective_charts"] = 39
result["source_sha256"] = {name: hashlib.sha256(Path(name).read_bytes()).hexdigest()
                           for name in ("griess_subalgebras.json", "primaries_II.pkl", "primaries.py", "griess_subalgebras.py", "isotropy_test.py")}
output.parent.mkdir(parents=True, exist_ok=True)
output.write_text(json.dumps(result, indent=2) + "\n")
print("PASS: 13 candidates, 39 unit Groebner bases; output:", output, flush=True)
