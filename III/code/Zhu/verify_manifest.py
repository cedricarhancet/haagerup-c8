"""Check the complete extracted supplementary archive using the standard library."""
from pathlib import Path
import hashlib,sys
ROOT=Path(__file__).resolve().parent
def main():
 records={};errors=[]
 for line in (ROOT/'SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
  h,n=line.split(maxsplit=1);n=n.lstrip('*');rel=Path(n)
  if rel.is_absolute() or '..' in rel.parts:raise ValueError('Unsafe path')
  if n in records:raise ValueError('Duplicate path')
  records[n]=h;p=ROOT/rel
  if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=h:errors.append(n)
 actual={p.relative_to(ROOT).as_posix() for p in ROOT.rglob('*') if p.is_file() and '__pycache__' not in p.parts and p.suffix!='.pyc' and p.name!='SHA256SUMS.txt'}
 errors+=sorted(actual-set(records))
 if errors:raise SystemExit('Archive manifest mismatch: '+', '.join(errors))
 print('PASS:',len(records),'files match the complete supplementary archive')
if __name__=='__main__':main()
