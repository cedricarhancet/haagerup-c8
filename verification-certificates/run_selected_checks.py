#!/usr/bin/env python3
"""Rerun the corrected, selected checks in a disposable copy of the source data.

This is not the full proof pipeline. Install requirements-selected.txt first.
The original snapshot and stored outputs are never overwritten.
"""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
SOURCES = ROOT / 'calculs-sources'


def run(command, cwd, log):
    print('Running:', ' '.join(command), 'in', cwd.name, flush=True)
    with log.open('w') as stream:
        subprocess.run(command, cwd=cwd, stdout=stream, stderr=subprocess.STDOUT, check=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check', choices=['primary', 'closure', 'jiao', 'e6', 'all'], default='all')
    parser.add_argument('--output', type=Path, default=HERE / 'rerun')
    args = parser.parse_args()
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((SOURCES / 'MANIFEST-sources.json').read_text())
    for entry in manifest:
        source = ROOT / entry['path']
        assert source.stat().st_size == entry['bytes'], source
        assert hashlib.sha256(source.read_bytes()).hexdigest() == entry['sha256'], source
    print('Original snapshot verified:', len(manifest), 'files', flush=True)
    with tempfile.TemporaryDirectory(prefix='haagerup-selected-checks-') as temporary:
        work = Path(temporary)
        requested = ['primary', 'closure', 'jiao', 'e6'] if args.check == 'all' else [args.check]
        needed = {'primary': ['cas-II'], 'closure': ['cas-II'], 'jiao': ['unicite', 'zhu', 'generation-forte'], 'e6': ['gamma1']}
        for name in sorted({name for check in requested for name in needed[check]}):
            shutil.copytree(SOURCES / name, work / name)
        if 'primary' in requested:
            run([sys.executable, str(HERE / 'check_primary33_exclusions.py'), str(work / 'cas-II'),
                 str(output / 'primary33-exclusions.json')], work / 'cas-II',
                output / 'primary33-exclusions.out')
        if 'closure' in requested:
            target = work / 'cas-II'
            shutil.copy2(HERE / 'closure_certify_exact_intervals.py', target)
            run([sys.executable, 'closure_certify_exact_intervals.py', 'II', '9'], target,
                output / 'closure-certify-exact-intervals.out')
            shutil.copy2(target / 'closure_certify_II_9.json', output / 'closure-certify-exact-intervals.json')
            shutil.copy2(HERE / 'check_closure_coverage.py', output)
            run([sys.executable, str(output / 'check_closure_coverage.py')], output,
                output / 'check-closure-coverage.json')
        if 'jiao' in requested:
            target = work / 'unicite'
            shutil.copy2(ROOT / 'checks' / 'check-jiao-normalization.py', target / 'exact_radicals.py')
            shutil.copy2(HERE / 'i5_membership_exact.py', target / 'i5_membership.py')
            shutil.copy2(HERE / 'i6_generation_certificate.py', target)
            run([sys.executable, 'i5_membership.py'], target, output / 'i5-membership-rerun.out')
            run([sys.executable, 'i6_generation_certificate.py'], target, output / 'i6-certificate-rerun.out')
            for name in ['i5_membership.json', 'i6_generation.json', 'i6_rank_certificate.json']:
                shutil.copy2(target / name, output / name)
            shutil.copy2(HERE / 'check_i6_rank_certificate.py', output)
            run([sys.executable, str(output / 'check_i6_rank_certificate.py')], output,
                output / 'check-i6-rank-certificate.json')
        if 'e6' in requested:
            run([sys.executable, 'e6_rank_bound.py'], work / 'gamma1', output / 'e6-rank-bound-rerun.out')
    print('Selected checks completed. This does not certify the complete manuscript.', flush=True)


if __name__ == '__main__':
    main()
