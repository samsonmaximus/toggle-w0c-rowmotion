# I2(m) root poset (Armstrong/Cuntz-Stump): minima a,b; chain c1<...<c_{m-2}; a,b<c1. Birational rowmotion order.
from birat import make_toggles, order_of
for m in range(3,15):
    N=m; a,b=0,1; chain=list(range(2,m))
    cov=[(a,2),(b,2)]+[(chain[k],chain[k+1]) for k in range(len(chain)-1)]
    T=make_toggles(N,cov); order=list(reversed([a,b]+chain))
    o=order_of(order,T,N,maxit=100,trials=3)
    import math
    print('m',m,'rowmotion order',o,'lcm(2,m)=',m*2//math.gcd(2,m))
