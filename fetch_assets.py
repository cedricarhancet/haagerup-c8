#!/usr/bin/env python3
"""Telecharger les deux volumes publics F4 et verifier leurs SHA-256."""
import argparse, hashlib, json, urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parent
def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args();out=a.output.resolve();out.mkdir(parents=True,exist_ok=False)
    records=[]
    for r in json.loads((ROOT/'EXTERNAL-ASSETS.json').read_text(encoding='utf-8'))['assets']:
        if Path(r['file']).name!=r['file']:raise ValueError('Nom de fichier invalide')
        path=out/r['file'];h=hashlib.sha256();n=0
        print('DOWNLOAD',r['file'],flush=True)
        with urllib.request.urlopen(r['url'],timeout=180) as src,path.open('xb') as dst:
            while True:
                b=src.read(4*1024*1024)
                if not b:break
                dst.write(b);h.update(b);n+=len(b)
        if n!=r['bytes'] or h.hexdigest()!=r['sha256']:
            raise RuntimeError('Empreinte ou taille incorrecte : '+r['file'])
        records.append(dict(file=r['file'],bytes=n,sha256=h.hexdigest(),status='PASS'))
        print('PASS',r['file'],flush=True)
    (out/'DOWNLOAD-RECEIPT.json').write_text(json.dumps(records,indent=2)+'\n',encoding='utf-8')
if __name__=='__main__':main()
