# Root systems (simple roots in a basis where reflections are exact), positive roots as integer/golden coords in simple-root basis.
import itertools, sympy as sp
phi=(1+sp.sqrt(5))/2
def cartan(typ,n):
    # returns symmetrizable Cartan matrix A with s_i(alpha_j) = alpha_j - A[i][j] alpha_i  (A[i][j] = <alpha_i^vee, alpha_j>)
    A=[[0]*n for _ in range(n)]
    for i in range(n): A[i][i]=2
    if typ=='A':
        for i in range(n-1): A[i][i+1]=A[i+1][i]=-1
    elif typ=='B':  # node n-1 short root? B_n: alpha_{n-1} short
        for i in range(n-1): A[i][i+1]=A[i+1][i]=-1
        A[n-2][n-1]=-1; A[n-1][n-2]=-2
    elif typ=='D':
        for i in range(n-2): A[i][i+1]=A[i+1][i]=-1
        A[n-3][n-1]=A[n-1][n-3]=-1
    elif typ=='H3':
        A=[[2,-phi,0],[-phi,2,-1],[0,-1,2]]
    elif typ=='I':  # dihedral I2(m), n=m
        m=n; c=2*sp.cos(sp.pi/m); A=[[2,-c],[-c,2]]
    return [[sp.nsimplify(x) for x in r] for r in A]
def reflect(A,i,v):
    n=len(A); v=list(v); coef=sum(A[i][j]*v[j] for j in range(n))  # <alpha_i^vee, v>
    v[i]=sp.simplify(v[i]-coef); return tuple(v)
def positive_roots(A):
    n=len(A); simple=[tuple(1 if j==i else 0 for j in range(n)) for i in range(n)]
    roots=set(simple); frontier=list(simple)
    while frontier:
        new=[]
        for r in frontier:
            for i in range(n):
                s=reflect(A,i,r)
                if all(sp.simplify(x)>=0 for x in s) and s not in roots and any(sp.simplify(x)!=0 for x in s):
                    roots.add(s); new.append(s)
        frontier=new
    return sorted(roots,key=lambda r:(float(sum(r)),[float(x) for x in r]))
