import sys, itertools, networkx as nx, sympy as sp
from roots import cartan, positive_roots, reflect
def setup(typ,n):
    A=cartan(typ,n); R=positive_roots(A); idx={r:i for i,r in enumerate(R)}
    N=len(R); nn=len(A)
    def le(a,b): return all(sp.simplify(b[k]-a[k])>=0 for k in range(nn))
    rel=[[le(R[i],R[j]) and i!=j for j in range(N)] for i in range(N)]
    covers=[(i,j) for i in range(N) for j in range(N) if rel[i][j] and not any(rel[i][k] and rel[k][j] for k in range(N))]
    return A,R,idx,covers
def apply_word(A,word,v):
    for i in reversed(word): v=reflect(A,i,v)
    return v
def csort_w0(A,c,N):
    word=[]; k=0
    while len(word)<N:
        s=c[k%len(c)]; k+=1
        # length increases iff u(alpha_s) > 0 where u = product of word
        alpha=tuple(1 if j==s else 0 for j in range(len(A)))
        v=apply_word(A,word,alpha)
        if all(sp.simplify(x)>=0 for x in v): word.append(s)
    return word
def inversion_order(A,word,idx):
    seq=[]
    for k,s in enumerate(word):
        alpha=tuple(1 if j==s else 0 for j in range(len(A)))
        seq.append(idx[apply_word(A,word[:k],alpha)])
    return seq
def circulations(N,covers,pos):
    # pos[i] = time element i is toggled. Orientation: edge from earlier to later.
    G=nx.Graph(); G.add_nodes_from(range(N)); G.add_edges_from(covers)
    basis=nx.cycle_basis(G); out=[]
    for cyc in basis:
        tot=0
        for a,b in zip(cyc,cyc[1:]+cyc[:1]):
            tot+= 1 if pos[a]<pos[b] else -1
        out.append(tot)
    return out
if __name__=='__main__':
    typ=sys.argv[1]; n=int(sys.argv[2])
    A,R,idx,covers=setup(typ,n); N=len(R)
    rank=lambda r: float(sum(r))
    # rowmotion: toggle top first -> pos = reverse of a linear extension
    lin=sorted(range(N),key=lambda i:(rank(R[i]),i))
    pos_row={x:N-1-t for t,x in enumerate(lin)}; pos_inv={x:t for t,x in enumerate(lin)}
    cr=circulations(N,covers,pos_row); ci=circulations(N,covers,pos_inv)
    nn=len(A); res=set()
    for perm in itertools.permutations(range(nn)):
        c=list(perm); word=csort_w0(A,c,N); seq=inversion_order(A,word,idx)
        assert sorted(seq)==list(range(N))
        pos={x:t for t,x in enumerate(seq)}
        cc=circulations(N,covers,pos)
        tag='ROW' if cc==cr else ('ROW^-1' if cc==ci else 'OTHER')
        res.add(tag); print(typ,n,'c=',''.join(str(i+1) for i in c),'->',tag, '' if tag!='OTHER' else (cc,cr))
    print(typ,n,'N',N,'covers',len(covers),'summary',res)
