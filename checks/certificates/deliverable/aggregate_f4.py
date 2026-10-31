if not __debug__:
 raise RuntimeError("Run without -O or PYTHONOPTIMIZE; assertions are part of the verifier")
"""Accept F4 only after complete exact verification and an independent zero check."""
from pathlib import Path
import argparse,hashlib,json
p=argparse.ArgumentParser();p.add_argument('--zeros',type=Path,required=True);a=p.parse_args()
root=Path('.');smhash=hashlib.sha256((root/'SOURCE-MANIFEST.json').read_bytes()).hexdigest()
source_manifest=json.loads((root/'SOURCE-MANIFEST.json').read_text())
for source in source_manifest['archived_files']+source_manifest['new_verifiers']:
 q=root/source['name']
 assert q.stat().st_size==source['bytes'] and hashlib.sha256(q.read_bytes()).hexdigest()==source['sha256'],('source changed after verification',q)
parts=[json.loads(q.read_text()) for q in sorted((root/'f4_verified').glob('PART_*_OF_*.json'))]
assert parts,'No completed exact verification parts'
n=parts[0]['parts'];assert len(parts)==n and {r['part'] for r in parts}==set(range(n))
assert all(r['parts']==n and r['status']=='PASS' and r['source_manifest_sha256']==smhash for r in parts)
manifests={str(p.parent):json.loads(p.read_text()) for p in sorted(root.glob('f4_*/COMPLETE.json'))}
assert all(set(manifests)==set(r['input_manifest_hashes'])==set(r['seen_labels']) for r in parts)
for name,m in manifests.items():
 h=hashlib.sha256((root/name/'COMPLETE.json').read_bytes()).hexdigest()
 assert all(r['input_manifest_hashes'].get(name)==h for r in parts)
 expected={tuple(x) for x in m['nonzero_labels']};seen=set()
 for r in parts:
  xs=[tuple(x) for x in r['seen_labels'][name]]
  assert len(xs)==len(set(xs)) and not seen&set(xs),('duplicate part labels',name)
  seen.update(xs)
 assert seen==expected,('incomplete part coverage',name)
blocks=[b for r in parts for b in r['blocks']]
assert len({b['block'] for b in blocks})==len(blocks)
expected_blocks={b['block'] for m in manifests.values() for b in m['blocks']}
assert {b['block'] for b in blocks}==expected_blocks,('missing or extra checked blocks')
assert all(b['all_rows_verified'] and b['bad']==0 for b in blocks)
assert sum(b['pivots'] for b in blocks)==6010
assert sum(b['family_nonzero'] for b in blocks)==manifests['f4_families']['nonzero']
tests=[(name,m) for name,m in manifests.items() if name.startswith('f4_tests_')]
assert len(tests)==10 and sum(m['total'] for _,m in tests)==9800
assert manifests['f4_families']['total']==18368
assert sum(b['tests_nonzero'] for b in blocks)==sum(m['nonzero'] for _,m in tests)
zeros=json.loads(a.zeros.read_text());assert zeros['status']=='PASS'
zero_names=['f4_families']+[name for name,m in tests]
expected_zero_manifest_hashes={name+'/COMPLETE.json':hashlib.sha256((root/name/'COMPLETE.json').read_bytes()).hexdigest() for name in zero_names}
assert zeros['manifest_hashes']==expected_zero_manifest_hashes,('zero check does not cover all current input manifests')
assert zeros['verified']==sum(manifests[name]['zero'] for name in zero_names)
expected_zeros={(name+'/COMPLETE.json',tuple(label)) for name in zero_names for label in manifests[name]['zero_labels']}
zero_records=zeros['records'];assert len(zero_records)==zeros['verified']
actual_zeros=[(r['manifest'],tuple(r['label'])) for r in zero_records]
assert len(actual_zeros)==len(set(actual_zeros)) and set(actual_zeros)==expected_zeros
assert all(r['method'] in ('direct exact kprod','Cartan energy > 4') for r in zero_records)
assert zeros['S3_basis_sha256']==hashlib.sha256(Path('exact_S3_basis.pkl').read_bytes()).hexdigest()
result={'status':'PASS','scope':'Full exact F4 reconstruction and all-row memberships; independent re-evaluation of every zero input','pivot_words':6010,'cartan_blocks_checked':len(blocks),'family_inputs':18368,'family_nonzero':manifests['f4_families']['nonzero'],'family_zero':manifests['f4_families']['zero'],'product_inputs':9800,'product_nonzero':sum(m['nonzero'] for _,m in tests),'product_zero':sum(m['zero'] for _,m in tests),'independently_verified_zero_inputs':zeros['verified'],'source_manifest_sha256':smhash,'input_manifest_hashes':parts[0]['input_manifest_hashes'],'zero_check_sha256':hashlib.sha256(a.zeros.read_bytes()).hexdigest(),'parts':n,'errors':0}
p=Path('COMPLETE-F4.json');p.with_suffix('.tmp').write_text(json.dumps(result,indent=2));p.with_suffix('.tmp').replace(p)
print(json.dumps(result,indent=2))
