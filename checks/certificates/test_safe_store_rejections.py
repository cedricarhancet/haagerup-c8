"""Exercise the actual SafeStore in small temporary fixtures, without VOA data."""
from pathlib import Path
import argparse
import hashlib
import importlib.util
import json
import pickle
import subprocess
import sys
import tempfile
from types import SimpleNamespace

parser = argparse.ArgumentParser()
parser.add_argument('--source', required=True, type=Path)
parser.add_argument('--output', required=True, type=Path)
args = parser.parse_args()
source = args.source.resolve()
spec = importlib.util.spec_from_file_location('actual_safe_store', source)
module = importlib.util.module_from_spec(spec)
spec.loader.exec_module(module)
blocks = SimpleNamespace(
    blk={(0, 0): ['a', 'b'], (1, 1): ['c']},
    pos={(0, 0): {'a': 0, 'b': 1}, (1, 1): {'c': 0}},
)
results = []

def fixture(root, name):
    store = module.SafeStore(root / name, blocks, 'fixed-source-hash')
    store.add(('first', 0), {'a': 2})
    store.add(('second', 0), {'b': 3})
    store.add(('third', 1), {'c': 4})
    store.add(('zero', 0), {})
    return store

def rewrite(path, records):
    with path.open('wb') as handle:
        for record in records:
            pickle.dump(record, handle, protocol=5)

with tempfile.TemporaryDirectory(prefix='f4-safestore-fixtures-') as tmp:
    root = Path(tmp)
    store = fixture(root, 'valid')
    manifest = store.finish()
    assert manifest['total'] == 4 and manifest['nonzero'] == 3 and manifest['zero'] == 1
    assert manifest['post_write_full_reread'] is True
    results.append({'case': 'valid', 'accepted': True})

    cases = [
        'missing_block', 'missing_record', 'duplicate_record',
        'unexpected_label', 'wrong_block', 'truncated_record',
        'trailing_bytes', 'unexpected_block', 'replaced_with_prefix',
        'duplicate_add', 'mixed_weight_add',
    ]
    for name in cases:
        store = fixture(root, name)
        path = store.path / 'blk_000.pkl'
        try:
            if name == 'missing_block':
                path.unlink()
            elif name in ('missing_record', 'replaced_with_prefix'):
                destination = path if name == 'missing_record' else path.with_suffix('.tmp')
                rewrite(destination, [(('first', 0), {'a': 2})])
                if destination != path:
                    destination.replace(path)
            elif name == 'duplicate_record':
                rewrite(path, [(('first', 0), {'a': 2}), (('first', 0), {'a': 2})])
            elif name == 'unexpected_label':
                rewrite(path, [(('first', 0), {'a': 2}), (('unannounced', 0), {'b': 3})])
            elif name == 'wrong_block':
                rewrite(path, [(('first', 0), {'a': 2}), (('second', 0), {'c': 3})])
            elif name == 'truncated_record':
                path.write_bytes(path.read_bytes()[:-3])
            elif name == 'trailing_bytes':
                with path.open('ab') as handle:
                    handle.write(b'junk')
            elif name == 'unexpected_block':
                rewrite(store.path / 'blk_002.pkl', [(('extra', 2), {'c': 1})])
            elif name == 'duplicate_add':
                store.add(('first', 0), {'a': 9})
            elif name == 'mixed_weight_add':
                store.add(('mixed', 0), {'a': 1, 'c': 1})
            store.finish()
        except (AssertionError, EOFError, pickle.UnpicklingError, FileNotFoundError, KeyError) as exc:
            assert not (store.path / 'COMPLETE.json').exists()
            results.append({'case': name, 'rejected': True, 'exception': type(exc).__name__})
        else:
            raise AssertionError(('invalid fixture accepted', name))

    optimized = subprocess.run([sys.executable, '-O', str(source)], capture_output=True, text=True)
    assert optimized.returncode != 0 and 'RuntimeError' in optimized.stderr
    results.append({'case': 'python_optimized', 'rejected': True, 'exception': 'RuntimeError'})

result = {'status': 'PASS', 'source_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
          'accepted': 1, 'rejected': len(results) - 1, 'results': results,
          'scope': 'temporary I/O protocol fixtures only; no F4 mathematical identities'}
args.output.write_text(json.dumps(result, indent=2) + '\n')
print(json.dumps({k: v for k, v in result.items() if k != 'results'}))
