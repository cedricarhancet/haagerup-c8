"""Verify this standalone Article II delivery with the standard library."""
from pathlib import Path
import hashlib,json
def main():
 root=Path(__file__).resolve().parent
 data=json.loads((root/'MANIFEST-II.json').read_text(encoding='utf-8'));errors=[]
 for name,row in data.items():
  rel=Path(name)
  if rel.is_absolute() or '..' in rel.parts:errors.append('Unsafe: '+name);continue
  p=root/rel
  if not p.is_file():errors.append('Missing: '+name);continue
  if p.stat().st_size!=row['bytes'] or hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:errors.append('Changed: '+name)
 actual={p.relative_to(root).as_posix() for p in root.rglob('*') if p.is_file() and not {'.git','.venv','__pycache__'}.intersection(p.relative_to(root).parts) and p.suffix!='.pyc' and p.name!='MANIFEST-II.json'}
 errors+=['Unexpected: '+n for n in sorted(actual-set(data))]
 if errors:raise SystemExit('\n'.join(errors))
 print('PASS:',len(data),'files match Article II delivery')
if __name__=='__main__':main()
