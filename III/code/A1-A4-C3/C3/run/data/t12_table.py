"""T12 -- parse Article I, Table tab:lowest-spaces (case II), read-only, and export data.

Writes t12_table.json : for each label z, the g'+g''-module R_z as {(a,b): mult},
the conformal weights h_z (gamma = 0) and dim R_z, with consistency checks.
Usage: python3 -I t12_table.py <LaTeX source of Article I [Arh26], version v14>
(file Haagerup-c8-article-I-v14-arxiv.tex, included in ../../inputs-I; see ../../REPRODUCIBILITY.md)
"""
import sys, re, json
from fractions import Fraction as F

src = open(sys.argv[1], encoding="utf-8").read()
i_lab = src.index(r"\label{tab:lowest-spaces}")
i_beg = src.rindex(r"\begin{tabular}", 0, i_lab)
i_end = src.index(r"\end{tabular}", i_beg)
body = src[i_beg:i_end].split("\n")

rows = {}          # label -> list of strings (case II column)
cur = None
for line in body:
    if "&" not in line:
        continue
    parts = line.split("&")
    if len(parts) != 3:
        continue
    head, c1, c2 = parts
    head = head.strip()
    m = re.match(r"\$(R_\{?[a-z0-9_]+\}?(?:,\s*R_\{?[a-z0-9_]+\}?)*|U_2)\$", head)
    if m:
        cur = m.group(1)
        rows[cur] = [c2]
    elif head == "" and cur is not None:
        rows[cur].append(c2)

def parse(cell_list):
    txt = " ".join(cell_list).replace("\\\\", " ").replace("$", " ").replace("\\ ", " ")
    txt = txt.replace("\\hline", " ")
    d = {}
    for mult, a, b in re.findall(r"(\d*)\((\d+),(\d+)\)", txt.replace(" ", "")):
        k = (int(a), int(b))
        d[k] = d.get(k, 0) + (int(mult) if mult else 1)
    return d

labels = {}
for key, cells in rows.items():
    mod = parse(cells)
    if key == "U_2":
        labels["U2"] = mod
        continue
    for nm in re.findall(r"R_\{?([a-z0-9_]+)\}?", key):
        nm = nm.replace("_", "")
        labels[nm] = mod

# conformal weights and dimensions from eq:lowest-spaces with gamma = 0
h = {"0": F(0), "b": F(1), "a": F(1), "c0": F(1), "c1": F(1, 3), "c2": F(2, 3),
     "d1": F(6, 13), "d2": F(11, 13), "d3": F(15, 13), "d4": F(5, 13), "d5": F(7, 13), "d6": F(8, 13)}
dimR = {"0": 1, "b": 80, "a": 81, "c0": 81, "c1": 3, "c2": 27,
        "d1": 7, "d2": 42, "d3": 119, "d4": 5, "d5": 13, "d6": 14}
# check the eq:lowest-spaces block in the source too (gamma = 0 substitution)
blk = src[src.index(r"\label{eq:lowest-spaces}"):]
blk = blk[:blk.index(r"\end{equation}")]
assert r"\frac{6}{13}" in blk and r"\frac{11}{13}" in blk and r"\frac{2}{13}+\epsi" in blk
labels["0"] = {(0, 0): 1}
out = {}
ok = True
for z in h:
    mod = labels[z]
    d = sum((a + 1) * (b + 1) * m for (a, b), m in mod.items())
    if d != dimR[z]:
        ok = False
    # level bound (3,13)
    for (a, b) in mod:
        assert a <= 3 and b <= 13
    out[z] = {"h": str(h[z]), "dimR": d, "R": {f"{a},{b}": m for (a, b), m in sorted(mod.items())}}
    print(f"{z:3s} h={str(h[z]):6s} dimR={d:4d} (expected {dimR[z]:4d})  R={dict(sorted(mod.items()))}")
dU2 = sum((a + 1) * (b + 1) * m for (a, b), m in labels["U2"].items())
print("U2 dim =", dU2, "(expected 120)")
out["U2"] = {f"{a},{b}": m for (a, b), m in sorted(labels["U2"].items())}
print("all dims ok:", ok and dU2 == 120)
json.dump(out, open("t12_table.json", "w"), indent=1)
