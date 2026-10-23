# Deterministic check Row^5 == Row^{-5} on H3 Cuntz-Stump poset, labels 1 at 0hat,1hat (or symbolic with ab=1 flag)
import json, time, flint, sys
symb = len(sys.argv)>1 and sys.argv[1]=='ab'
iso=json.load(open('h3_iso.json')); cov=[tuple(c) for c in iso["('lexgreater', 'N^6 (D6 order)')"][1]]
n=15
names=tuple([f"x{k}" for k in range(n)]+(['A','B'] if symb else []))
ctx=flint.fmpz_mpoly_ctx.get(names,'degrevlex'); gens=ctx.gens(); xs=gens[:n]
one=ctx.from_dict({})+1
down={k:[] for k in range(n)}; up={k:[] for k in range(n)}
for i,j in cov: down[j-1].append(i-1); up[i-1].append(j-1)
mins=[k for k in range(n) if not down[k]]; rank={k:0 for k in mins}; ch=True
while ch:
    ch=False
    for k in range(n):
        if k not in rank and down[k] and all(x in rank for x in down[k]): rank[k]=rank[down[k][0]]+1; ch=True
row=sorted(range(n),key=lambda k:-rank[k])
def red(p,q):
    g=p.gcd(q)
    if g!=one: p=p//g; q=q//g
    if q.leading_coefficient()<0: p=-p; q=-q
    return (p,q)
def add(F,G): return red(F[0]*G[1]+G[0]*F[1], F[1]*G[1])
def mul(F,G): return red(F[0]*G[0], F[1]*G[1])
def invf(F): return red(F[1],F[0])
A=(gens[n],one) if symb else (one,one); B=(gens[n+1],one) if symb else (one,one)
def tog(f,v):
    s=None
    for u in down[v]: s=f[u] if s is None else add(s,f[u])
    if s is None: s=A
    tt=None
    for u in up[v]:
        iu=invf(f[u]); tt=iu if tt is None else add(tt,iu)
    if tt is None: tt=invf(B)
    f[v]=mul(s, invf(mul(f[v],tt)))
t0=time.time()
F={k:(xs[k],one) for k in range(n)}; G=dict(F)
for it in range(1,6):
    for v in row: tog(F,v)
    for v in reversed(row): tog(G,v)
    print(f"step {it}: terms F {max(len(p) for P in F.values() for p in P)} G {max(len(p) for P in G.values() for p in P)} t={time.time()-t0:.0f}s",flush=True)
eq=all(F[k][0]==G[k][0] and F[k][1]==G[k][1] for k in range(n))
print("Row^5 == Row^-5 (hence Row^10 = id):",eq, flush=True)
# also check Row^5 != id and Row^2 != id
print("Row^5 == id:", all(F[k][0]==xs[k] and F[k][1]==one for k in range(n)))
