"""Check termination, inclusion and overlap ambiguities of the stored Zhu rules.
Uses a leftmost/shortest reduction order, different from the discovery script.
No VOA engine or numerical linear algebra is needed for this verification.
"""
import argparse,pickle,json,time
from collections import Counter,defaultdict
from functools import lru_cache
from pathlib import Path
ap=argparse.ArgumentParser();ap.add_argument('directory',type=Path);args=ap.parse_args();d=args.directory
D=pickle.load(open(d/'z4_upper_algebra_modp.pkl','rb'));rules=D['rules']
source=pickle.load(open(d/'corner_rows_modp.pkl','rb'));P=source['P'];AR=source['ARROWS']
lengths=sorted({len(q) for q in rules});start=time.monotonic()
def source_of(q):return AR[q[-1]]['source']
def target_of(q):return AR[q[0]]['target']
def add(out,row,c=1):
 for q,a in row.items():
  t=(out.get(q,0)+c*a)%P
  if t:out[q]=t
  else:out.pop(q,None)
for q,row in rules.items():
 assert len(q)>=2
 assert all(AR[q[i]]['source']==AR[q[i+1]]['target'] for i in range(len(q)-1))
 for w,c in row.items():
  assert 0<c<P
  assert len(w)<len(q) or (len(w)==len(q) and w>q), ('nondecreasing',q,w)
  assert all(AR[w[i]]['source']==AR[w[i+1]]['target'] for i in range(len(w)-1))
  assert (source_of(w),target_of(w))==(source_of(q),target_of(q)) if w else source_of(q)==target_of(q)
@lru_cache(None)
def nf(q):
 for i in range(len(q)):
  for L in lengths:
   p=q[i:i+L]
   if len(p)==L and p in rules:
    z={}
    for w,c in rules[p].items():add(z,nf(q[:i]+w+q[i+L:]),c)
    return z
 return {q:1}
def normal(poly):
 z={}
 for q,c in poly.items():add(z,nf(q),c)
 return z
inclusions=0
for q,row in rules.items():
 for L in lengths:
  if L>=len(q):continue
  for i in range(len(q)-L+1):
   sub=q[i:i+L]
   if sub not in rules:continue
   inclusions+=1;other={}
   for w,c in rules[sub].items():add(other,nf(q[:i]+w+q[i+L:]),c)
   assert normal(row)==other, ('inclusion',q,sub,i)
print('inclusions checked:',inclusions,flush=True)
starts=defaultdict(list)
for q in rules:starts[q[0]].append(q)
overlaps=0
for q1,row1 in rules.items():
 for k in range(1,len(q1)):
  suffix=q1[k:]
  for q2 in starts[suffix[0]]:
   if len(q2)<=len(suffix) or q2[:len(suffix)]!=suffix:continue
   overlaps+=1;left={};right={}
   for w,c in row1.items():add(left,nf(w+q2[len(suffix):]),c)
   for w,c in rules[q2].items():add(right,nf(q1[:k]+w),c)
   assert left==right, ('overlap',q1,q2,k)
active=D['active'];words=[[(a,) for a in active]]
while words[-1] and len(words)<8:
 candidates=[q+(a,) for q in words[-1] for a in active if AR[q[-1]]['source']==AR[a]['target']]
 words.append([q for q in candidates if all(q[i:i+L] not in rules for L in lengths for i in range(len(q)-L+1))])
assert len(words)==4 and not words[-1]
assert [len(q) for q in words]==[239,506,5,0]
assert len(D['verts'])+sum(map(len,words))==762
result={'status':'PASS','prime':P,'rules_by_length':dict(Counter(map(len,rules))),'inclusion_ambiguities':inclusions,'overlap_ambiguities':overlaps,'nonresolving':0,'normal_paths_by_length':[len(q) for q in words],'vertices':len(D['verts']),'dimension':762,'elapsed_seconds':round(time.monotonic()-start,3)}
(d/'rewriting-system-check.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
