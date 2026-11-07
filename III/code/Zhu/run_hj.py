"""Reproduce the Zhu calculations for Article III in a fresh working directory."""
import argparse
import hashlib
import json
import lzma
import os
from pathlib import Path
import pickle
import re
import shutil
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

PACKAGE = Path(__file__).resolve().parent
REFERENCE = PACKAGE / "reference"
EXPECTED_ROWS = "d4a4390e2e80cd9ba7fd2b4b26f44a1bf7fbdf8806ddac31aef68b7946a91541"
PATH_LITERAL = re.compile(r"(['\"])/home/claude/(?:hj/)?work/([^'\"\n]*)\1")
ENTRY = (
    "from pathlib import Path; import runpy,sys;"
    "sys.argv=sys.argv[1:]; sys.path.insert(0,str(Path(sys.argv[0]).parent));"
    "runpy.run_path(sys.argv[0],run_name='__main__')"
)


def sha(path):
    with path.open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def json_read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def load(path):
    # Inputs are the locally archived data verified against SOURCE-SHA256.txt.
    with path.open("rb") as stream:
        return pickle.load(stream)


def require(test, what):
    if not test:
        raise RuntimeError(what)
    print(what, flush=True)


def prepare(work, from_modes):
    for line in (PACKAGE / "SOURCE-SHA256.txt").read_text(encoding="utf-8").splitlines():
        expected, relative = line.split("  ", 1)
        if sha(PACKAGE / relative) != expected:
            raise RuntimeError("Source checksum mismatch: " + relative)
    if work.exists():
        raise FileExistsError("Choose a fresh run directory: " + str(work))
    shutil.copytree(PACKAGE / "Calculs", work / "Calculs")
    shutil.copytree(PACKAGE / "verif", work / "verif")
    (work / "logs").mkdir()
    for path in work.rglob("*.py"):
        content = path.read_text(encoding="utf-8-sig")
        relocated = PATH_LITERAL.sub(lambda m: repr((work / m.group(2)).as_posix()), content)
        if content != relocated:
            path.write_text(relocated, encoding="utf-8")
    exact = work / "Calculs" / "zhu-exact"
    if not from_modes:
        with lzma.open(exact / "rows_exact.pkl.xz", "rb") as compressed:
            with (exact / "rows_exact.pkl").open("wb") as expanded:
                shutil.copyfileobj(compressed, expanded)
        require(sha(exact / "rows_exact.pkl") == EXPECTED_ROWS, "Expanded exact rows checksum")
        shutil.copyfile(exact / "rows_exact_reduced_modp.pkl", exact / "specfib" / "corner_rows_modp.pkl")
    return exact


def run(work, script, *args, cwd=None, name=None, local_pdir=None):
    cmd = [sys.executable, "-I", "-B", "-u", "-c", ENTRY, str(script), *map(str, args)]
    logfile = work / "logs" / ((name or script.stem) + ".log")
    env = os.environ.copy()
    for key in ("HAAG_PSTART", "HAAG_PDIR"):
        env.pop(key, None)
    if local_pdir is not None:
        env["HAAG_PDIR"] = str(local_pdir)
    print("Running", script.name, *args, flush=True)
    with logfile.open("w", encoding="utf-8") as out:
        result = subprocess.run(cmd, cwd=cwd or script.parent, env=env, stdout=out, stderr=subprocess.STDOUT)
    if result.returncode:
        raise RuntimeError(f"{script.name} returned {result.returncode}; see {logfile}")
    return logfile


def verify_confluence(work, exact, overlap_limit=None):
    script = work / "verif" / "verif_exact.py"
    expected = json_read(REFERENCE / "verif/verif_exact.json")
    if overlap_limit is not None:
        source = script.read_text(encoding="utf-8")
        before, middle = source.split("# ---------------- C. overlaps ----------------", 1)
        middle, after = middle.split("# ---------------- D. inclusions ----------------", 1)
        middle = middle.replace("for q1, r1 in RULES.items():",
                                f"for q1, r1 in RULES.items():\n    if nC >= {overlap_limit}: break")
        middle = middle.replace("    for k in range(1, len(q1)):",
                                f"    for k in range(1, len(q1)):\n        if nC >= {overlap_limit}: break")
        middle = middle.replace("        for q2 in LEADS_BY_FIRST.get(suf[0], ()):",
                                f"        for q2 in LEADS_BY_FIRST.get(suf[0], ()):\n            if nC >= {overlap_limit}: break")
        script = script.with_name("verif_exact_sample.py")
        script.write_text(before + middle + after, encoding="utf-8")
        expected["C_overlaps"] = overlap_limit
    run(work, script)
    got = json_read(work / "verif" / "verif_exact.json")
    require(got == expected,
            "Exact critical pair checks, forward containment, integrality and reduction match")


def verify_membership(work, exact):
    membership = work / "verif" / "membership"
    log = run(work, membership / "membership_0_vertices_E.py").read_text(encoding="utf-8")
    for v in (12, 13):
        require(re.search(rf"^vertex {v} .*IN SPAN, exact identity verified", log, re.M) is not None,
                f"Vertex {v} in the original ideal")
    run(work, membership / "membership_a.py")
    data = load(membership / "G0.pkl")
    rules = load(exact / "upper_exact.pkl")
    killed = set(rules["killed"])
    arrows = load(exact / "rows_exact.pkl")["ARROWS"]
    live_elim = {
        a for a in rules["alllin"]
        if arrows[a]["source"] not in killed and arrows[a]["target"] not in killed
    }
    require(set(data["cert_elim"]) == live_elim and len(live_elim) == 46,
            "46 live eliminations in the original ideal")
    require(len(set(data["G0"])) == 5105 and set(data["G0"]) <= rules["rules"].keys(),
            "5105 initial rules in the original ideal")
    log = run(work, membership / "membership_b_targeted.py").read_text(encoding="utf-8")
    require(re.search(r"554 of 554 failures 0\b", log) is not None,
            "554 completion rules in the original ideal")


def verify_fibre(work, exact):
    specfib = exact / "specfib"
    run(work, specfib / "z4_upper_algebra_modp-ArticleI.py", local_pdir=specfib)
    data = json_read(specfib / "z4_upper_algebra_modp.json")
    blocks = sorted((b["size"], b["dim"], sorted(b["vertex_ranks"].items())) for b in data["blocks"])
    saved = json_read(PACKAGE / "Calculs/zhu-exact/specfib/z4_upper_algebra_modp.json")
    expected = sorted((b["size"], b["dim"], sorted(b["vertex_ranks"].items())) for b in saved["blocks"])
    require((data["dimension"], data["trace_form_rank"], data["center_dim"], blocks)
            == (762, 762, 12, expected), "Special fibre dimension, trace form and split blocks")
    new = load(specfib / "z4_upper_algebra_modp.pkl")
    old = load(PACKAGE / "Calculs/zhu-exact/specfib/z4_upper_algebra_modp.pkl")
    require(all(new[k] == old[k] for k in ("rules", "active", "verts")),
            "Recomputed special fibre rewriting rules equal the archived reductions")
    run(work, exact / "specfib_Z_from_rules.py", name="specfib_Z")
    run(work, exact / "specfib_blocks.py")
    sorted_blocks = lambda entries: sorted(entries, key=lambda e: e["root"])
    require(sorted_blocks(json_read(specfib / "specfib_blocks.json")) ==
            sorted_blocks(json_read(PACKAGE / "Calculs/zhu-exact/specfib/specfib_blocks.json")),
            "Special fibre weights and corner dimensions match")


def verify_weights(work, exact):
    run(work, exact / "exactZ.py")
    require(json_read(exact / "exactZ.json") == json_read(REFERENCE / "exactZ.json"),
            "Exact central polynomial and occurrence of every listed root")
    run(work, exact / "check_patterns.py")
    require(json_read(exact / "check_patterns.json") == json_read(REFERENCE / "check_patterns.json"),
            "Both lattice branchings checked against all three possible partitions")


def project_relations(work, exact):
    run(work, exact / "z3x.py")
    archived_rows = PACKAGE / "Calculs/zhu-exact/rows_exact.pkl.xz"
    with lzma.open(archived_rows, "rb") as stream:
        old_rows = pickle.load(stream)
    require(load(exact / "rows_exact.pkl") == old_rows, "All 19981 projected exact rows agree")
    require(load(exact / "rows_exact_reduced_modp.pkl")
            == load(PACKAGE / "Calculs/zhu-exact/rows_exact_reduced_modp.pkl"),
            "Regenerated reduction agrees")


def regenerate(work, exact, jobs, rebuild_rules, sample_rows=None):
    zhu = exact.parent / "zhu"
    run(work, zhu / "z0_highest.py")
    basis = load(zhu / "U2_lowering_basis.pkl")
    require(basis == load(PACKAGE / "Calculs/zhu/U2_lowering_basis.pkl"),
            "Exact lowering basis regenerated")
    ranks = json_read(zhu / "z0_highest.json")
    require(ranks["rank_lowering_mod_p"] == ranks["rank_with_U2_mod_p"] == 120
            and not any(v["not_killed_by"] for v in ranks["highest"]),
            "Highest vectors and 120-dimensional lowering basis checked")
    product_script = exact / "prod_exact.py"
    if sample_rows is not None:
        source = product_script.read_text(encoding="utf-8")
        sample = product_script.with_name("prod_exact_sample.py")
        sample.write_text(source.replace("for a in range(120):", f"for a in range({sample_rows}):"),
                          encoding="utf-8")
        product_script = sample
    with ThreadPoolExecutor(max_workers=jobs) as pool:
        futures = [pool.submit(run, work, product_script, k, jobs, name=f"products-{k}")
                   for k in range(jobs)]
        for future in futures:
            future.result()
    count = 120 if sample_rows is None else sample_rows
    for a in range(count):
        data = load(exact / f"prodx/a{a:03}.pkl")
        if set(data) != set(range(a, 120)):
            raise RuntimeError(f"Incomplete product row {a}")
    print(f"{count} product rows regenerated", flush=True)
    if sample_rows is not None:
        return
    run(work, exact / "bucket.py")
    run(work, exact / "z2x.py")
    for name in ("rel_exact.pkl", "rel_exact_reduced_modp.pkl"):
        require(load(exact / name) == load(PACKAGE / "Calculs/zhu-exact" / name),
                "Regenerated exact data agree: " + name)
    project_relations(work, exact)
    if rebuild_rules:
        run(work, exact / "z4x.py")
        require(load(exact / "upper_exact.pkl") == load(PACKAGE / "Calculs/zhu-exact/upper_exact.pkl"),
                "Regenerated exact rewriting system agrees")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("confluence", "membership", "fibre", "weights", "verify", "regenerate", "projection", "sample", "confluence-sample"))
    parser.add_argument("run_dir", type=Path, help="fresh directory for this run")
    parser.add_argument("--jobs", type=int, default=1, help="product workers in regenerate mode")
    parser.add_argument("--rebuild-rules", action="store_true", help="also recompute the exact completion")
    parser.add_argument("--rows", type=int, default=10, help="initial product rows in sample mode")
    args = parser.parse_args()
    if args.jobs < 1:
        parser.error("--jobs must be positive")
    work = args.run_dir.resolve()
    if not 1 <= args.rows <= 120:
        parser.error("--rows must lie in 1..120")
    exact = prepare(work, args.mode in ("regenerate", "sample", "projection"))
    modes = {"confluence": verify_confluence, "membership": verify_membership,
             "fibre": verify_fibre, "weights": verify_weights}
    if args.mode in ("regenerate", "sample"):
        regenerate(work, exact, args.jobs, args.rebuild_rules,
                   args.rows if args.mode == "sample" else None)
    elif args.mode == "projection":
        project_relations(work, exact)
    elif args.mode == "confluence-sample":
        verify_confluence(work, exact, overlap_limit=1024)
    elif args.mode == "verify":
        for fn in modes.values():
            fn(work, exact)
    else:
        modes[args.mode](work, exact)
    (work / "result.json").write_text(json.dumps({"mode": args.mode, "status": "passed"}) + "\n",
                                     encoding="utf-8")
    print("Verified", args.mode, flush=True)


if __name__ == "__main__":
    main()
