"""Fetch fixed public supplementary archives and verify them before use."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,stat,subprocess,urllib.request,zipfile
ROOT=Path(__file__).resolve().parent
def digest(p):
 with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def load(p):return json.loads(p.read_text(encoding='utf-8'))
def verify_payload(root):
 manifest=root/'SHA256SUMS.txt'
 if manifest.is_file():
  rows={}
  for line in manifest.read_text(encoding='utf-8').splitlines():
   h,n=line.split(maxsplit=1);rows[n.lstrip('*')]={'sha256':h}
 else:rows=load(root/'MANIFEST-II.json')
 for n,row in rows.items():
  rel=Path(n)
  if rel.is_absolute() or '..' in rel.parts:raise ValueError('Unsafe manifest path')
  p=root/rel
  if not p.is_file() or digest(p)!=row['sha256'] or ('bytes' in row and p.stat().st_size!=row['bytes']):raise RuntimeError('Manifest mismatch: '+n)
 return len(rows)
def fetch(keys,destination,local_assets=None):
 destination=Path(destination).resolve();repo=ROOT.parent if (ROOT.parent/'MANIFEST-PUBLIC.json').is_file() else ROOT
 if destination==repo or destination.is_relative_to(repo):raise ValueError('Download outside the repository')
 assets=load(ROOT/'ASSETS.json');destination.mkdir(parents=True,exist_ok=True);result={}
 for key in keys:
  spec=assets[key];p=destination/spec['filename'];part=p.with_name(p.name+'.part')
  if not p.exists():
   if part.exists():raise FileExistsError('Use a fresh download cache')
   curl=shutil.which('curl.exe' if os.name=='nt' else 'curl')
   local=Path(local_assets)/spec['filename'] if local_assets is not None else None
   if local is not None:
    if not local.is_file():raise FileNotFoundError('Missing offline asset: '+str(local))
    shutil.copyfile(local,part)
   elif curl:
    r=subprocess.run([curl,'-L','--fail','--silent','--show-error','--max-time','600',spec['url'],'-o',str(part)],capture_output=True,text=True)
    if r.returncode:raise RuntimeError(r.stderr)
   else:
    request=urllib.request.Request(spec['url'],headers={'User-Agent':'Haagerup-c8-supplements'})
    with urllib.request.urlopen(request,timeout=600) as stream,part.open('xb') as out:shutil.copyfileobj(stream,out)
   if part.stat().st_size!=spec['bytes'] or digest(part)!=spec['sha256']:raise RuntimeError('Downloaded archive hash mismatch: '+key)
   part.rename(p)
  if p.stat().st_size!=spec['bytes'] or digest(p)!=spec['sha256']:raise RuntimeError('Cached archive hash mismatch: '+key)
  if spec.get('kind')=='file':result[key]=str(p);print('PASS',key,spec['sha256'],flush=True);continue
  extracted=destination/('extracted-'+key);root=extracted/spec['root']
  if not extracted.exists():
   with zipfile.ZipFile(p) as z:
    if z.testzip() is not None:raise RuntimeError('ZIP CRC error')
    for info in z.infolist():
     if not (extracted/info.filename).resolve().is_relative_to(extracted.resolve()) or stat.S_ISLNK(info.external_attr>>16):raise ValueError('Unsafe ZIP entry')
    z.extractall(extracted)
  count=verify_payload(root);result[key]=str(root);print('PASS',key,spec['sha256'],count,'manifest entries',flush=True)
 (destination/'EXTRACTED.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8');return result
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--local-assets',type=Path);p.add_argument('--packages',nargs='+',choices=['calculations','zhu','II','I_source','III_source'],default=['calculations','zhu','II','I_source','III_source']);a=p.parse_args();fetch(a.packages,a.output,a.local_assets)
if __name__=='__main__':main()
