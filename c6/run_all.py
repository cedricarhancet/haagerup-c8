#!/usr/bin/env python3
"""Replay every computation of the c6 supplements in a fresh directory.

Usage: python -I -X utf8 -B run_all.py --output ../../c6-run

The scripts are copied to the output directory before running, because some of
them write their JSON certificates next to themselves. The new text outputs are
compared with the recorded ones, and the new JSON certificates are compared,
as parsed JSON, with the recorded certificates.
"""
import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent

TEXT_RUNS = [
    ('lattice/verify_c6_section6_v01.py', 'lattice/verify_c6_section6_v01-output.txt'),
    ('lattice/verify_c6_complements_v12.py', 'lattice/verify_c6_complements_v12-output.txt'),
]
JSON_RUNS = [
    ('independent/independent_lattice_recheck.py', 'independent/independent_lattice_recheck.json'),
    ('independent/independent_complements_recheck.py', 'independent/independent_complements_recheck.json'),
    ('independent/verify_weyl_isotropic_permutation.py', 'independent/weyl-isotropic-permutation.json'),
    ('independent/verify_mixed_pointed_bulk.py', 'independent/mixed-pointed-bulk-arithmetic.json'),
]


def run(script, cwd):
    proc = subprocess.run([sys.executable, '-I', '-X', 'utf8', '-B', str(script)],
                          cwd=cwd, capture_output=True, text=True, encoding='utf-8')
    if proc.returncode != 0:
        raise SystemExit('FAILED: %s\n%s' % (script.name, proc.stderr))
    return proc.stdout


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', required=True, help='new directory outside the repository')
    out = Path(parser.parse_args().output).resolve()
    if out.exists():
        raise SystemExit('The output directory must not exist: %s' % out)
    if HERE in out.parents or out == HERE:
        raise SystemExit('Choose an output directory outside the supplements folder.')
    for sub in ('lattice', 'independent'):
        shutil.copytree(HERE / sub, out / sub, ignore=shutil.ignore_patterns('*.json', '*-output.txt'))
    for script, recorded in TEXT_RUNS:
        new = run(out / script, (out / script).parent)
        (out / recorded).write_text(new, encoding='utf-8')
        old = (HERE / recorded).read_text(encoding='utf-8')
        if new.replace('\r\n', '\n') != old.replace('\r\n', '\n'):
            raise SystemExit('Output differs from the recorded one: %s' % recorded)
        print('PASS', script)
    for script, recorded in JSON_RUNS:
        run(out / script, (out / script).parent)
        new = json.loads((out / recorded).read_text(encoding='utf-8'))
        old = json.loads((HERE / recorded).read_text(encoding='utf-8'))
        if new != old:
            raise SystemExit('Certificate differs from the recorded one: %s' % recorded)
        print('PASS', script)
    print('ALL SUPPLEMENTS REPRODUCED in', out)


if __name__ == '__main__':
    main()
