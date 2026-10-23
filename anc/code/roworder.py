import sys, random, networkx as nx, sympy as sp
from roots import cartan, positive_roots
def poset(typ,n):
    A=cartan(typ,n); R=positive_roots(A); N=len(R)
    Rf=[[float(x) for x in r] for r in R]
    G=nx.DiGraph(); G.add_nodes_from(range(N))
    G.add_edges_from([(i,j) for i in range(N) for j in range(N) if i!=j and all(Rf[j][k]-Rf[i][k]>=-1e-12 for k in range(n))])
    return N,list(nx.transitive_reduction(G).edges())
def order(N,cov,maxit=80):
    down={i:[] for i in range(N)}; up={i:[] for i in range(N)}
    for a,b in cov: up[a].append(b); down[b].append(a)
    G=nx.DiGraph(cov); G.add_nodes_from(range(N)); top=list(nx.topological_sort(G))[::-1]
    rnd=random.Random(1); f0=[rnd.uniform(0.5,2) for _ in range(N)]; A,B=1.7,0.6; f=list(f0)
    for k in range(1,maxit+1):
        for v in top:
            s=sum(f[u] for u in down[v]) if down[v] else A
            t=sum(1/f[u] for u in up[v]) if up[v] else 1/B
            f[v]=s/(f[v]*t)
        if max(abs(x-y)/abs(y) for x,y in zip(f,f0))<1e-9: return k
    return None
for t,n in [('A',2),('A',3),('A',4),('A',5),('B',2),('B',3),('B',4),('B',5),('B',6),('D',4)]:
    N,cov=poset(t,n); print(t,n,'birational rowmotion order (float detection):',order(N,cov))
