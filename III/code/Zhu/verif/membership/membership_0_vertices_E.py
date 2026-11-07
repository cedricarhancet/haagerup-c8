import pickle, sys, time
from collections import defaultdict
sys.path.insert(0,'/home/claude/hj/work/Calculs/zhu-exact')
import kexact as KX, Kfield as K
from flint import nmod_mat
from gmpy2 import mpq
t0=time.time()
ZE='/home/claude/hj/work/Calculs/zhu-exact'
U=pickle.load(open(ZE+'/upper_exact.pkl','rb')); D=pickle.load(open(ZE+'/rows_exact.pkl','rb'))
AR=D['ARROWS']; P=KX.P; killed=set(U['killed'])
byep=defaultdict(list)
for ep,o,r in D['ROWS']: byep[ep].append(r)
del D
def piv(A, nr, nc):
    R, r = A.rref(); J=[]; j=0
    for i in range(r):
        while int(R[i,j])==0: j+=1
        J.append(j); j+=1
    return r, J
def inspan(rows, target):
    keys=sorted({k for v in rows+[target] for k in v}); idx={k:i for i,k in enumerate(keys)}; N=len(keys)
    red=[KX.kv_red(v,idx,N) for v in rows]
    A=nmod_mat(len(rows),N,[x for v in red for x in v],P)
    r,J=piv(A,len(rows),N)
    B=nmod_mat(len(rows)+1,N,[x for v in red+[KX.kv_red(target,idx,N)] for x in v],P)
    if B.rank()>r: return 'not in span mod p', r
    AJt=nmod_mat(r,len(rows),[int(A[i,c]) for c in J for i in range(len(rows))],P)
    r2,I=piv(AJt,r,len(rows)); assert r2==r
    At=[[rows[I[a]].get(keys[J[b]],K.Z) for a in range(r)] for b in range(r)]
    Bt=[[target.get(keys[J[b]],K.Z)] for b in range(r)]
    X=KX.solve_K(At,Bt)
    acc={}
    for a in range(r):
        c=X[a][0]
        if K.iszero(c): continue
        for k,v in rows[I[a]].items(): acc[k]=K.add(acc.get(k,K.Z),K.mul(c,v))
    for k,v in target.items(): acc[k]=K.add(acc.get(k,K.Z),K.scal(mpq(-1),v))
    ok=all(K.iszero(v) for v in acc.values())
    return ('IN SPAN, exact identity verified' if ok else 'FAILED exact verification'), r
for v in (12,13):
    s,r=inspan(byep[(v,v)],{():K.ONE})
    print('vertex',v,'rows',len(byep[(v,v)]),'rank',r,':',s,'%.0fs'%(time.time()-t0),flush=True)
cnt=defaultdict(int)
for a,l in sorted(U['alllin'].items()):
    ep=(AR[a]['target'],AR[a]['source'])
    if ep[0] in killed or ep[1] in killed: continue
    tgt={(a,):K.ONE}
    for w,c in l.items(): tgt[w]=K.add(tgt.get(w,K.Z),K.scal(mpq(-1),c))
    s,r=inspan(byep[ep],tgt); cnt[s]+=1
    print('E arrow',a,'ep',ep,'rows',len(byep[ep]),'rank',r,':',s,'%.0fs'%(time.time()-t0),flush=True)
print(dict(cnt),flush=True)
