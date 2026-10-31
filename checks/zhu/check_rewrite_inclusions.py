import argparse,pickle,json
from collections import Counter
from functools import lru_cache
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);args=ap.parse_args();d=args.directory
D=pickle.load(open(d/'z4_upper_algebra_modp.pkl','rb'));rules=D['rules'];P=8388571
lengths=sorted({len(q) for q in rules},reverse=True)
def add(out,row,c=1):
 for q,a in row.items():
  t=(out.get(q,0)+c*a)%P
  if t:out[q]=t
  else:out.pop(q,None)
@lru_cache(None)
def nf(q):
 for L in lengths:
  for i in range(len(q)-L+1):
   p=q[i:i+L]
   if p in rules:
    z={}
    for w,c in rules[p].items():add(z,nf(q[:i]+w+q[i+L:]),c)
    return z
 return {q:1}
count=0;bad=[]
for q,r in rules.items():
 for L in lengths:
  if L>=len(q):continue
  for i in range(len(q)-L+1):
   p=q[i:i+L]
   if p not in rules:continue
   count+=1;diff={}
   for w,c in r.items():add(diff,nf(w),c)
   for w,c in rules[p].items():add(diff,nf(q[:i]+w+q[i+L:]),-c)
   if diff:bad.append({'word':q,'subword':p,'position':i,'residual':[[list(w),c] for w,c in diff.items()]})
result={'rules_by_length':dict(Counter(map(len,rules))),'inclusions':count,'nonresolving':len(bad),'examples':bad[:5]}
(d/'inclusion-check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
