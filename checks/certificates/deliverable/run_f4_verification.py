#!/usr/bin/env python3
"""Replay the exact F4 verification from archived coordinates, or regenerate.

Default: verify all input domains, source hashes, pivot provenance against
family records, and all-row exact memberships. With --regenerate, additionally
regenerate every word and every zero record from the lattice VOA formulas.
Use a fresh --output directory; no existing certificate is overwritten.
"""
if not __debug__:
 raise RuntimeError("Run without -O or PYTHONOPTIMIZE")
import argparse,hashlib,json,os,shutil,subprocess,sys,zipfile
from pathlib import Path
HERE=Path(__file__).resolve().parent

def source_default():
 for p in [HERE]+list(HERE.parents):
  candidate=p/'calculs-sources'/'generation-forte'
  if candidate.is_dir():return candidate
 raise FileNotFoundError('Pass --sources with the generation-forte directory')

def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument('--sources',type=Path)
 ap.add_argument('--archive',type=Path,nargs='+')
 ap.add_argument('--output',type=Path,required=True)
 ap.add_argument('--regenerate',action='store_true')
 args=ap.parse_args();sources=(args.sources or source_default()).resolve();out=args.output.resolve()
 out.mkdir(parents=True,exist_ok=False)
 manifest=json.loads((HERE/'SOURCE-MANIFEST.json').read_text())
 for row in manifest['archived_files']:
  p=sources/row['name'];data=p.read_bytes()
  assert len(data)==row['bytes'] and hashlib.sha256(data).hexdigest()==row['sha256'],('Source differs',p)
  shutil.copy2(p,out/p.name)
 for row in manifest['new_verifiers']:
  name=row['name'];q=HERE/name
  assert q.stat().st_size==row['bytes'] and hashlib.sha256(q.read_bytes()).hexdigest()==row['sha256']
  shutil.copy2(q,out/name)
 shutil.copy2(HERE/'SOURCE-MANIFEST.json',out/'SOURCE-MANIFEST.json')
 env=dict(os.environ)
 if args.regenerate:
  commands=[['regenerate_f4.py','pivots'],['regenerate_f4.py','tests']]
  for script_args in commands:
   with (out/('regeneration-'+script_args[1]+'.out')).open('w') as log:
    subprocess.run([sys.executable,'-u']+script_args,cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 else:
  if not args.archive:raise ValueError('Pass --archive, or use --regenerate')
  archived=json.loads((HERE/'ARCHIVES-F4.json').read_text())
  expected={r['file']:r for r in archived['archives']}
  assert len(args.archive)==len(expected) and {p.name for p in args.archive}==set(expected), 'Supply each archive volume exactly once'
  extracted=set()
  for archive in args.archive:
   row=expected[archive.name]
   assert archive.stat().st_size==row['bytes'],('Archive size differs',archive)
   with archive.open('rb') as stream:
    assert hashlib.file_digest(stream,'sha256').hexdigest()==row['sha256'],('Archive hash differs',archive)
   with zipfile.ZipFile(archive) as z:
    for info in z.infolist():
     p=Path(info.filename)
     if p.is_absolute() or '..' in p.parts:raise ValueError('Unsafe archive member')
     if p.parts and p.parts[0]=='certificats-poids-quatre-v03':p=Path(*p.parts[1:])
     if p.parts and (p.parts[0] in ('f4_pivots','f4_families') or p.parts[0].startswith('f4_tests_')):
      target=out/p
      if info.is_dir():target.mkdir(parents=True,exist_ok=True)
      else:
       assert str(p) not in extracted,('Duplicate archive member',str(p))
       extracted.add(str(p))
       target.parent.mkdir(parents=True,exist_ok=True)
       with z.open(info) as src,target.open('wb') as dst:shutil.copyfileobj(src,dst)
  assert len(extracted)==archived['total_files'],('Incomplete archive extraction',len(extracted))
 with (out/'verification.out').open('w') as log:
  subprocess.run([sys.executable,'-u','verify_f4_strict.py'],cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 with (out/'zero-check.out').open('w') as log:
  subprocess.run([sys.executable,'-u','check_zero_products.py','--work',str(out),'--output',str(out/'zero-products-independent.json'),'--seconds','3600'],cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 with (out/'aggregation.out').open('w') as log:
  subprocess.run([sys.executable,'-u','aggregate_f4.py','--zeros',str(out/'zero-products-independent.json')],cwd=out,env=env,stdout=log,stderr=subprocess.STDOUT,check=True)
 print('Exact F4 verification completed; see',out/'COMPLETE-F4.json')

if __name__=='__main__':main()
