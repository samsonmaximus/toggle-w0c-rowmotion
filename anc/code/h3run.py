import json, itertools
from fractions import Fraction as Fr
from roots import cartan, positive_roots
from orient import csort_w0, inversion_order, circulations
from birat import make_toggles, order_of
d=json.load(open('h3_map.json')); p2r={int(k):v for k,v in d['poset_to_root'].items()}; r2p={v:k for k,v in p2r.items()}
cov=[tuple(e) for e in d['covers']]; N=15
A=cartan('H3',3); R=positive_roots(A); idx={r:i for i,r in enumerate(R)}
T=make_toggles(N,cov)
import networkx as nx
G=nx.DiGraph(cov); rk={v:0 for v in range(N)}
for v in nx.topological_sort(G):
    for u in G.predecessors(v): rk[v]=max(rk[v],rk[u]+1)
lin=sorted(range(N),key=lambda i:(rk[i],i))
pos_row={x:N-1-k for k,x in enumerate(lin)}; pos_inv={x:k for k,x in enumerate(lin)}
cr=circulations(N,cov,pos_row); ci=circulations(N,cov,pos_inv)
print('H3 rowmotion order:',order_of(list(reversed(lin)),T,N,maxit=60,trials=3))
for perm in itertools.permutations(range(3)):
    word=csort_w0(A,list(perm),N); seq=[r2p[i] for i in inversion_order(A,word,idx)]
    pos={x:k for k,x in enumerate(seq)}; cc=circulations(N,cov,pos)
    tag='ROW' if cc==cr else ('ROW^-1' if cc==ci else 'OTHER')
    print('c=',''.join(str(i+1) for i in perm),'word',''.join(str(i+1) for i in word),'orientation',tag,'K_c order:',order_of(seq,T,N,maxit=60,trials=3),flush=True)
