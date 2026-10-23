import sys, itertools, random
from fractions import Fraction as Fr
from orient import setup, csort_w0, inversion_order
import networkx as nx
def count_antichains(N,covers):
    G=nx.DiGraph(); G.add_nodes_from(range(N)); G.add_edges_from(covers)
    TC=nx.transitive_closure_dag(G); comp=[[TC.has_edge(i,j) or TC.has_edge(j,i) for j in range(N)] for i in range(N)]
    cnt=0
    def rec(i,chosen):
        nonlocal cnt
        if i==N: cnt+=1; return
        rec(i+1,chosen)
        if all(not comp[i][j] for j in chosen): rec(i+1,chosen+[i])
    rec(0,[]); return cnt
def make_toggles(N,covers):
    down={i:[] for i in range(N)}; up={i:[] for i in range(N)}
    for a,b in covers: up[a].append(b); down[b].append(a)
    def T(f,v):
        s=sum((f[u] for u in down[v]),Fr(0)) if down[v] else Fr(1)   # 0-hat label 1
        t=sum((1/f[u] for u in up[v]),Fr(0)) if up[v] else Fr(1)     # 1-hat label 1
        g=dict(f); g[v]=s/(f[v]*t); return g
    return T
def order_of(seq,T,N,maxit=200,trials=2,seed=0):
    rnd=random.Random(seed); ords=set()
    for _ in range(trials):
        f0={i:Fr(rnd.randint(1,97),rnd.randint(1,89)) for i in range(N)}
        f=f0
        for k in range(1,maxit+1):
            for v in seq: f=T(f,v)
            if f==f0: ords.add(k); break
        else: ords.add(None)
    return ords
if __name__=='__main__':
    typ=sys.argv[1]; n=int(sys.argv[2])
    A,R,idx,covers=setup(typ,n); N=len(R)
    print(typ,n,'N',N,'antichains',count_antichains(N,covers))
    T=make_toggles(N,covers)
    lin=sorted(range(N),key=lambda i:(float(sum(R[i])),i))
    print(' rowmotion order:',order_of(list(reversed(lin)),T,N))
    for perm in itertools.permutations(range(len(A))):
        word=csort_w0(A,list(perm),N); seq=inversion_order(A,word,idx)
        print(' c=',''.join(str(i+1) for i in perm),'K_c order:',order_of(seq,T,N),flush=True)
