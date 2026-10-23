# H3 root poset (Cuntz–Stump): realised inside the D6 root poset. From each of the 15 pairs of Cuntz–Stump Fig. 2 take the
# lexicographically greater vector (v-coordinates); this is the poset used in the paper (h3_poset.json).
import itertools, networkx as nx, sympy as sp
from birat import count_antichains
pairs_flat=[(0,0,0,0,1,1),(1,-1,0,0,0,0),(0,0,0,1,-1,0),(0,1,-1,0,0,0),(0,0,0,1,0,1),(1,0,-1,0,0,0),(0,0,1,-1,0,0),(0,0,0,0,1,-1),
 (0,0,1,0,0,-1),(0,1,0,-1,0,0),(0,0,1,0,1,0),(1,0,0,-1,0,0),(0,1,0,0,-1,0),(0,0,0,1,0,-1),(0,1,0,0,0,-1),(0,0,1,0,-1,0),
 (0,1,0,1,0,0),(1,0,0,0,-1,0),(0,1,1,0,0,0),(1,0,0,0,0,-1),(1,0,0,0,0,1),(0,0,0,1,1,0),(1,0,0,0,1,0),(0,0,1,0,0,1),
 (1,0,0,1,0,0),(0,1,0,0,0,1),(1,0,1,0,0,0),(0,1,0,0,1,0),(1,1,0,0,0,0),(0,0,1,1,0,0)]
# simple roots of D6: a_i = v_i - v_{i+1} (i=1..5), a_6 = v5+v6. Convert v-coords -> simple-root coords.
M=sp.Matrix([[1,-1,0,0,0,0],[0,1,-1,0,0,0],[0,0,1,-1,0,0],[0,0,0,1,-1,0],[0,0,0,0,1,-1],[0,0,0,0,1,1]]).T  # columns = simple roots in v-coords
Minv=M.inv()
def eps(v): return tuple(int(x) for x in Minv*sp.Matrix(v))
tau=(1+sp.sqrt(5))/2
def gam(a): return tuple(sp.simplify(a[2*i]+a[2*i+1]*tau) for i in range(3))
left=[pairs_flat[2*i] for i in range(15)]; right=[pairs_flat[2*i+1] for i in range(15)]
L=[eps(v) for v in left]; Rr=[eps(v) for v in right]
for a,b in zip(L,Rr):
    assert all(x>=0 for x in a) and all(x>=0 for x in b)
    # (The tau-relation between the two members of a pair holds in Cuntz-Stump's gamma-coordinates, which use a
    #  different ordering of the simple roots; it is not needed here. The H3 labelling used is the table in h3_map.json.)
def le(a,b): return all(y-x>=0 for x,y in zip(a,b))
def poset_from(V):
    N=len(V); G=nx.DiGraph(); G.add_nodes_from(range(N))
    G.add_edges_from([(i,j) for i in range(N) for j in range(N) if i!=j and le(V[i],V[j])])
    return list(nx.transitive_reduction(G).edges())
for name,V in [('left',L),('right',Rr)]:
    cov=poset_from(V); N=len(V)
    rk={}
    G=nx.DiGraph(cov); 
    print(name,'covers',len(cov),'antichains',count_antichains(N,cov), 'D6-heights',sorted(sum(v) for v in V))
if __name__=='__main__':
    import json
    LEX=[eps(max(a,b)) for a,b in zip(left,right)]          # lexicographically greater representative of each pair
    print('lex-greater: covers',len(poset_from(LEX)),'antichains',count_antichains(15,poset_from(LEX)))
    d=json.load(open('h3_poset.json')); D=[tuple(v) for v in d['D6coords']]
    assert set(LEX)==set(D), 'lex-greater set differs from h3_poset.json'
    assert {tuple(e) for e in poset_from(D)}=={tuple(e) for e in d['covers']}, 'covers differ'
    print('OK: lex-greater representatives reproduce h3_poset.json (same 15 elements, same covers); file not modified')
