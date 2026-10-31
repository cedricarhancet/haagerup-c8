"""Package immutable completed F4 coordinate streams, excluding working data."""
from pathlib import Path
import argparse,hashlib,json,zipfile
ap=argparse.ArgumentParser();ap.add_argument('--source',type=Path,required=True);ap.add_argument('--output',type=Path,required=True);a=ap.parse_args();source=a.source.resolve();out=a.output.resolve();out.mkdir(exist_ok=True)
dirs=['f4_pivots','f4_families']+sorted(p.name for p in source.glob('f4_tests_*') if p.is_dir())
assert len(dirs)==12
files=[]
for name in dirs:
 d=source/name;m=json.loads((d/'COMPLETE.json').read_text());assert m['post_write_full_reread']
 actual={p.name for p in d.iterdir() if p.is_file()};expected={'COMPLETE.json'}|{f"blk_{b['block']:03}.pkl" for b in m['blocks']};assert actual==expected
 files+=[d/'COMPLETE.json']+[d/f"blk_{b['block']:03}.pkl" for b in sorted(m['blocks'],key=lambda b:b['block'])]
archives=[];z=None;part=0;current_files=0;raw=0
for p in files:
 if z is None:
  part+=1;path=out/f'Certificats-poids-quatre-v03-partie-{part}.zip';assert not path.exists();z=zipfile.ZipFile(path,'w',zipfile.ZIP_DEFLATED,compresslevel=6);current_files=0;raw=0
 z.write(p,'certificats-poids-quatre-v03/'+p.relative_to(source).as_posix());current_files+=1;raw+=p.stat().st_size
 if z.fp.tell()>400*1024*1024:
  z.close();assert path.stat().st_size<500*1024*1024;archives.append({'file':path.name,'files':current_files,'raw_bytes':raw});print('part closed',path.name,path.stat().st_size,flush=True);z=None
if z is not None:z.close();archives.append({'file':path.name,'files':current_files,'raw_bytes':raw})
if len(archives)==1:
 old=out/archives[0]['file'];new=out/'Certificats-poids-quatre-v03.zip';assert not new.exists();old.rename(new);archives[0]['file']=new.name
for r in archives:
 p=out/r['file'];assert p.stat().st_size<500*1024*1024
 with zipfile.ZipFile(p) as z:assert z.testzip() is None
 with p.open('rb') as f:r['sha256']=hashlib.file_digest(f,'sha256').hexdigest()
 r['bytes']=p.stat().st_size;print('archive checked',r,flush=True)
result={'status':'PASS','archive_root':'certificats-poids-quatre-v03/','coordinate_domains':dirs,'total_files':len(files),'archives':archives,'scope':'Complete coordinate streams; mathematical verification receipt and source code are supplied in the main dossier.'}
(out/'ARCHIVES-F4.json').write_text(json.dumps(result,indent=2))
