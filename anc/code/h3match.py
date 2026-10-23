import json, itertools, sympy as sp
from roots import cartan, positive_roots
d=json.load(open('h3_poset.json')); V=[tuple(v) for v in d['D6coords']]; cov=[tuple(e) for e in d['covers']]
t=(1+5**0.5)/2
gam=lambda a:(a[0]+a[5]*t, a[1]+a[3]*t, a[4]+a[2]*t)
G=[gam(a) for a in V]
A=cartan('H3',3); R=positive_roots(A); Rf=[tuple(float(x) for x in r) for r in R]
def close(u,v): return all(abs(x-y)<1e-9 for x,y in zip(u,v))
for perm in itertools.permutations(range(3)):
    m={}
    for i,g in enumerate(G):
        gp=tuple(g[perm[k]] for k in range(3))
        for j,r in enumerate(Rf):
            for s in (1,t,1/t):
                if close(gp,tuple(s*x for x in r)): m[i]=j
    if len(m)==15 and len(set(m.values()))==15:
        print('perm',perm,'bijection found'); json.dump({'poset_to_root':m,'covers':cov},open('h3_map.json','w')); break
else: print('no bijection')
