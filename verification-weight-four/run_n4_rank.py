#!/usr/bin/env python3
"""Regenerate the N4 modular rank in a disposable source copy.

The complete run can take several minutes. The original source files and
recorded results are not overwritten. Install the requirements first.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import pickle
import shutil
import subprocess
import sys
import tempfile

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def export_npz(source, destination):
    record = pickle.loads((source / "N4_rank_blocks.pkl").read_bytes())
    arrays = {}
    names = {}
    for i, (weight, matrix) in enumerate(sorted(record["reduced_basis"].items())):
        name = f"block_{i:03d}"
        if np.any(matrix < 0) or np.any(matrix >= record["p"]):
            raise ValueError(f"Noncanonical residue in block {weight}")
        arrays[name] = matrix.astype(np.uint32)
        names[name] = str(weight)
    np.savez_compressed(destination / "N4_rank_blocks.npz", **arrays)
    (destination / "N4_rank_block_names.json").write_text(json.dumps(names, indent=2))
    shutil.copy2(source / "N4_rank_summary.json", destination)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE / "rerun")
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    source = ROOT / "calculs-sources" / "generation-forte"
    manifest = json.loads((ROOT / "calculs-sources" / "MANIFEST-sources.json").read_text())
    for entry in manifest:
        path = ROOT / entry["path"]
        if path.parent != source:
            continue
        if path.stat().st_size != entry["bytes"] or hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
            raise ValueError(f"Source changed: {path}")
    with (output / "check_n4_inputs.json").open("w") as log:
        subprocess.run([sys.executable, str(HERE / "check_n4_inputs.py"), str(source)],
                       stdout=log, check=True)
    with tempfile.TemporaryDirectory(prefix="haagerup-n4-rank-") as temporary:
        work = Path(temporary) / "generation-forte"
        shutil.copytree(source, work)
        shutil.copy2(HERE / "closure_w4_rank_export.py", work)
        environment = dict(os.environ, MODE="C2", PSTART="8388608")
        with (output / "N4_rank_run.out").open("w") as log:
            subprocess.run([sys.executable, "-u", "closure_w4_rank_export.py", "0"],
                           cwd=work, env=environment, stdout=log,
                           stderr=subprocess.STDOUT, check=True)
        export_npz(work, output)
        shutil.copy2(work / "closure_c2_modp_0_p8388571.json", output)
    with (output / "check_n4_rank_blocks.json").open("w") as log:
        subprocess.run([sys.executable, str(HERE / "check_n4_rank_blocks.py"), str(output)],
                       stdout=log, check=True)
    print(f"N4 rank generation and independent rank check completed: {output}")


if __name__ == "__main__":
    main()
