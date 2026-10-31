from pathlib import Path
import hashlib,json,pickle,runpy
old=pickle.loads(Path('exact_S3_basis.pkl').read_bytes())
runpy.run_path('certify_w3.py',run_name='__main__')
r=json.loads(Path('certify_w3.json').read_text())
assert r['U2_rank']==120 and r['S3_rank_lower_bound']==956
for k,v in r.items():
 if isinstance(v,list):assert len(v)==2 and v[1]==0,(k,v)
new=pickle.loads(Path('exact_S3_basis.pkl').read_bytes())
assert new==old,'Regenerated S3 basis differs from archived input'
r['exact_basis_equal_to_v02']=True
r['sha256']=hashlib.sha256(Path('exact_S3_basis.pkl').read_bytes()).hexdigest()
Path('certify_w3_validated.json').write_text(json.dumps(r,indent=2))
print('ALL CHECKS PASSED; regenerated basis equals archived input',flush=True)
