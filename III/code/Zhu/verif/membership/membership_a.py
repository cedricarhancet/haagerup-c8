# certify in I_K, independently of z4x's bookkeeping: (a) the 21 stage-1 eliminations, (b) the final rules lying in
# span{ nf_full(generators at ep) }, where nf_full kills paths through 12,13 (certified: e_12, e_13 in I_K) and
# substitutes eliminations already certified in I_K.  Everything exact over K; pivots mod p; exact verification.
import pickle, sys, time
from collections import defaultdict
from functools import lru_cache
sys.path.insert(0,'/home/claude/hj/work/Calculs/zhu-exact')
import kexact as KX, Kfield as K
from flint import nmod_mat
from gmpy2 import mpq
t0=time.time()
def log(*a): print(*a,'%.0fs'%(time.time()-t0),flush=True)
ZE='/home/claude/hj/work/Calculs/zhu-exact'
U=pickle.load(open(ZE+'/upper_exact.pkl','rb')); D=pickle.load(open(ZE+'/rows_exact.pkl','rb'))
AR=D['ARROWS']; P=KX.P; killed=set(U['killed']); ALL=U['alllin']; RULES=U['rules']
byep=defaultdict(list)
for ep,o,r in D['ROWS']:
    if not (set(ep)&killed): byep[ep].append(r)
del D
s_=lambda a:AR[a]['source']; t_=lambda a:AR[a]['target']
def piv(A):
    R, r = A.rref(); J=[]; j=0
    for i in range(r):
        while int(R[i,j])==0: j+=1
        J.append(j); j+=1
    return r, J
def certify(rows, targets):
    """returns list of booleans: target in K-span(rows), with exact verification"""
    keys=sorted({k for v in rows+targets for k in v}); idx={k:i for i,k in enumerate(keys)}; N=len(keys)
    red=[KX.kv_red(v,idx,N) for v in rows]
    A=nmod_mat(len(rows),N,[x for v in red for x in v],P); r,J=piv(A)
    out=[]; cand=[]
    for t in targets:
        B=nmod_mat(len(rows)+1,N,[x for v in red+[KX.kv_red(t,idx,N)] for x in v],P)
        cand.append(B.rank()==r)
    if not any(cand): return [False]*len(targets)
    AJt=nmod_mat(r,len(rows),[int(A[i,c]) for c in J for i in range(len(rows))],P)
    r2,I=piv(AJt)
    At=[[rows[I[a]].get(keys[J[b]],K.Z) for a in range(r)] for b in range(r)]
    TT=[t for t,c in zip(targets,cand) if c]
    Bt=[[t.get(keys[J[b]],K.Z) for t in TT] for b in range(r)]
    X=KX.solve_K(At,Bt)
    res=[]; m=0
    for t,c in zip(targets,cand):
        if not c: res.append(False); continue
        acc={}
        for a in range(r):
            x=X[a][m]
            if K.iszero(x): continue
            for k,v in rows[I[a]].items(): acc[k]=K.add(acc.get(k,K.Z),K.mul(x,v))
        for k,v in t.items(): acc[k]=K.add(acc.get(k,K.Z),K.scal(mpq(-1),v))
        res.append(all(K.iszero(v) for v in acc.values())); m+=1
    return res
def elimrow(a):
    t={(a,):K.ONE}
    for w,c in ALL[a].items(): t[w]=K.add(t.get(w,K.Z),K.scal(mpq(-1),c))
    return t
E=[a for a in ALL if not ({s_(a),t_(a)}&killed)]
# certified set: arrows touching killed vertices are in the ideal of e_12,e_13 (certified separately), so their
# eliminations are harmless; stage-0 E eliminations are certified directly below
cert=set()
def make_nf(S):
    @lru_cache(None)
    def nf(p):
        if any(s_(a) in killed or t_(a) in killed for a in p): return {}
        for i,a in enumerate(p):
            if a in S:
                o={}
                for w,c in ALL[a].items():
                    for v,d in nf(p[:i]+w+p[i+1:]).items(): o[v]=K.add(o.get(v,K.Z),K.mul(c,d))
                return {k:v for k,v in o.items() if not K.iszero(v)}
        return {p:K.ONE}
    return nf
def apply(nf,row):
    o={}
    for p,c in row.items():
        for v,d in nf(p).items(): o[v]=K.add(o.get(v,K.Z),K.mul(c,d))
    return {k:v for k,v in o.items() if not K.iszero(v)}
for rnd in range(4):
    nf=make_nf(frozenset(cert))
    eps=defaultdict(list)
    for a in E:
        if a not in cert: eps[(t_(a),s_(a))].append(a)
    new=set()
    for ep,arrs in eps.items():
        rows=[x for x in (apply(nf,r) for r in byep[ep]) if x]
        ok=certify(rows,[elimrow(a) for a in arrs])
        new|={a for a,o in zip(arrs,ok) if o}
    cert|=new
    log('round',rnd,'certified E eliminations',len(cert),'of',len(E))
    if not new: break
# (b) rules in span of nf_full(generators)
nf=make_nf(frozenset(cert))
rep=defaultdict(list)
for q in RULES: rep[(t_(q[0]),s_(q[-1]))].append(q)
G0=[]
for ep,qs in sorted(rep.items()):
    rows=[x for x in (apply(nf,r) for r in byep[ep]) if x]
    tg=[]
    for q in qs:
        t={q:K.ONE}
        for w,c in RULES[q].items(): t[w]=K.add(t.get(w,K.Z),K.scal(mpq(-1),c))
        tg.append(t)
    res=certify(rows,tg) if rows else [False]*len(qs)
    G0+=[q for q,o in zip(qs,res) if o]
log('certified rules',len(G0))
pickle.dump({'G0':G0,'cert_elim':sorted(cert)},open('/home/claude/hj/work/verif/membership/G0.pkl','wb'))
