"""Verify the complete repository-side Article III delivery with the standard library."""
from pathlib import Path
import hashlib,json
ROOT=Path(__file__).resolve().parent
def main():
 data=json.loads((ROOT/'MANIFEST-III.json').read_text(encoding='utf-8'));rows=data['repository_files'];errors=[]
 for n,row in rows.items():
  rel=Path(n)
  if rel.is_absolute() or '..' in rel.parts:raise ValueError('Unsafe path')
  p=ROOT/rel
  if not p.is_file() or p.stat().st_size!=row['bytes'] or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:errors.append(n)
 actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and not {'.git','__pycache__','.venv'}.intersection(p.relative_to(ROOT).parts) and p.suffix!='.pyc' and p.name!='MANIFEST-III.json'}
 errors+=sorted(actual-set(rows))
 if errors:raise SystemExit('Manifest mismatch: '+', '.join(errors))
 print('PASS:',len(rows),'repository-side Article III files. Complete payloads are the fixed release archives.')
if __name__=='__main__':main()
