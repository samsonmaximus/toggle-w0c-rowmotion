# Randomized identity test of birational rowmotion periodicity over F_p, with a no-cancellation degree bound.
import json, random, sys, networkx as nx
def setup(cov,N):
    down={i:[] for i in range(N)}; up={i:[] for i in range(N)}
    for a,b in cov: up[a].append(b); down[b].append(a)
    G=nx.DiGraph(); G.add_nodes_from(range(N)); G.add_edges_from(cov)
    return down,up,list(nx.topological_sort(G))[::-1]
def degree_bound(cov,N,period):
    down,up,order=setup(cov,N)
    # each label is a fraction num/den with (dn, dd) total degree bounds (in all 17 variables)
    f=[(1,0)]*N; A=(1,0); B=(1,0)
    def add(a,b): return (max(a[0]+b[1],b[0]+a[1]), a[1]+b[1])
    def inv(a): return (a[1],a[0])
    def mul(a,b): return (a[0]+b[0],a[1]+b[1])
    for it in range(period):
        for v in order:
            s=A if not down[v] else None
            for u in down[v]: s=f[u] if s is None else add(s,f[u])
            t=inv(B) if not up[v] else None
            for u in up[v]: t=inv(f[u]) if t is None else add(t,inv(f[u]))
            f[v]=mul(s,inv(mul(f[v],t)))
    return max(max(x) for x in f)
def test(cov,N,period,p,trials,seed):
    down,up,order=setup(cov,N); rnd=random.Random(seed); ok=0
    for _ in range(trials):
        while True:
            x=[rnd.randrange(1,p) for _ in range(N)]; A=rnd.randrange(1,p); B=rnd.randrange(1,p)
            f=list(x); bad=False
            try:
                for it in range(period):
                    for v in order:
                        s=A if not down[v] else sum(f[u] for u in down[v])%p
                        t=pow(B,p-2,p) if not up[v] else sum(pow(f[u],p-2,p) for u in up[v])%p
                        den=f[v]*t%p
                        if den==0 or s==0: raise ZeroDivisionError
                        f[v]=s*pow(den,p-2,p)%p
            except ZeroDivisionError: continue
            break
        if f==x: ok+=1
        else: return False
    return ok
if __name__=='__main__':
    d=json.load(open('h3_map.json')); cov=[tuple(e) for e in d['covers']]
    per=int(sys.argv[1])
    D=degree_bound(cov,15,per)
    print('no-cancellation total degree bound for numerators/denominators after',per,'iterations:',D)
    primes=[(1<<61)-1, (1<<89)-1, (1<<107)-1, (1<<127)-1]
    for p in primes:
        r=test(cov,15,per,p,25,p%1000)
        print('p=2^k-1 bits',p.bit_length(),'all 25 random points fixed:',r)
