"""Exercise the real recovery control flow with a tiny mocked mode operator.
This tests coverage and failure paths, not lattice-VOA identities.
"""
from pathlib import Path
import argparse
import ast
import gc
import hashlib
import importlib.util
import json
import os
import pickle
import tempfile
from types import SimpleNamespace

parser = argparse.ArgumentParser()
parser.add_argument('--repair', required=True, type=Path)
parser.add_argument('--store', required=True, type=Path)
parser.add_argument('--output', required=True, type=Path)
args = parser.parse_args()
repair = args.repair.resolve()
store_source = args.store.resolve()
output = args.output.resolve()
spec = importlib.util.spec_from_file_location('fixture_safe_store', store_source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
tree = ast.parse(repair.read_text())
selected = [node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name in ('group', 'recover')]
assert len(selected) == 2
code = compile(ast.Module(body=selected, type_ignores=[]), str(repair), 'exec')
blocks = SimpleNamespace(blk={(0, 0): ['a'], (1, 1): ['c']}, pos={(0, 0): {'a': 0}, (1, 1): {'c': 0}})
labels = [('u(0)s', 'fixture', i) for i in range(4)]
zero = ('u(0)s', 'fixture', 4)
results = []
old_cwd = Path.cwd()

def dump_records(path, records, suffix=b''):
    with path.open('wb') as handle:
        for record in records:
            pickle.dump(record, handle, protocol=5)
        handle.write(suffix)

with tempfile.TemporaryDirectory(prefix='f4-recovery-fixtures-') as tmp:
    root = Path(tmp)
    for case in ('missing_records', 'truncated_record', 'duplicate_label', 'unexpected_label', 'wrong_block', 'changed_source', 'duplicate_manifest_label'):
        work = root / case
        source = work / 'input-streams' / 'fixture'
        source.mkdir(parents=True)
        path = source / 'blk_000.pkl'
        records = [(labels[0], {'a': 2}), (labels[1], {'a': 4})]
        suffix = b''
        if case == 'truncated_record':
            suffix = pickle.dumps((labels[2], {'a': 6}), protocol=5)[:-3]
        elif case == 'duplicate_label':
            records[1] = records[0]
        elif case == 'unexpected_label':
            records[1] = (('u(0)s', 'fixture', 99), {'a': 4})
        elif case == 'wrong_block':
            records[1] = (labels[1], {'c': 4})
        dump_records(path, records, suffix)
        nonzero = [list(label) for label in labels]
        if case == 'duplicate_manifest_label':
            nonzero[1] = nonzero[0]
        manifest = {'schema': 1, 'S3_basis_sha256': 'base', 'nonzero_labels': nonzero,
                    'zero_labels': [list(zero)], 'nonzero': 4, 'zero': 1, 'total': 5,
                    'blocks': [{'block': 0, 'weight': [0, 0], 'count': 4,
                                'bytes': path.stat().st_size, 'sha256': hashlib.sha256(path.read_bytes()).hexdigest()}]}
        (source / 'COMPLETE.json').write_text(json.dumps(manifest))
        if case == 'changed_source':
            with path.open('ab') as handle:
                handle.write(b'changed')
        calls = []
        def word(label):
            return [('mocked', 1)], 0, {'a': label[-1] + 1}
        def apply_op(operator, vector):
            calls.append(vector['a'] - 1)
            return {'a': 2 * vector['a']}
        namespace = {'Path': Path, 'json': json, 'hashlib': hashlib, 'pickle': pickle,
                     'basehash': 'base', 'B': blocks, 'weights': sorted(blocks.blk),
                     'SafeStore': module.SafeStore, 'word': word, 'gc': gc,
                     'restricted_op': lambda terms, mode, vectors: object(),
                     'E': SimpleNamespace(apply_op=apply_op), 'log': lambda *values: None}
        exec(code, namespace)
        os.chdir(work)
        try:
            store, final = namespace['recover']('fixture')
        except (AssertionError, EOFError, pickle.UnpicklingError, FileNotFoundError, KeyError) as exc:
            assert case not in ('missing_records', 'truncated_record'), (case, exc)
            assert not (work / 'fixture' / 'COMPLETE.json').exists()
            results.append({'case': case, 'rejected': True, 'exception': type(exc).__name__})
        else:
            assert case in ('missing_records', 'truncated_record')
            assert calls == [2, 3], calls
            assert final['preserved_nonzero_records'] == 2
            assert final['recomputed_nonzero_records'] == 2
            assert final['recomputed_labels'] == labels[2:]
            assert set(store.records) == set(labels) | {zero}
            assert final['post_write_full_reread'] is True
            truncation = final['source_blocks'][0]['truncation']
            assert bool(truncation) == (case == 'truncated_record')
            if truncation:
                assert truncation['tail_bytes'] == len(suffix)
                assert truncation['offset'] == final['source_blocks'][0]['last_complete_offset']
            results.append({'case': case, 'recovered_exact_missing_labels': True,
                            'preserved': 2, 'recomputed': 2})
        finally:
            os.chdir(old_cwd)

result = {'status': 'PASS', 'repair_sha256': hashlib.sha256(repair.read_bytes()).hexdigest(),
          'store_sha256': hashlib.sha256(store_source.read_bytes()).hexdigest(),
          'recovered': 2, 'rejected': 5, 'results': results,
          'scope': 'actual recovery control flow with mocked mode operator; not a VOA calculation'}
output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: v for k, v in result.items() if k != 'results'}))
