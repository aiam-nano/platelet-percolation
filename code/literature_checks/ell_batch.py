import numpy as np, math, json, sys
from ell_test import run
G={2:0.2629,3:0.2289,4:0.2003,5:0.1757,8:0.1262,10:0.1058,100:0.01248,1000:0.001275,2000:0.000637}
out=[]
for L,nr in ((24,16),(36,10)):
    for u in (2,3,5,10,20,50,100,300,1000,2000):
        nc=run(3,u,L,nr,5000+u+L,0.5 if u>=5 else 0.4)
        n=nc.mean(); e=nc.std()/math.sqrt(len(nc)); vol=4/3*math.pi/u
        phi=n*vol; row=dict(L=L,u=u,nreal=len(nc),n_c=n,err=e,phi_c=phi,phi_u=phi*u,garboczi=G.get(u))
        if G.get(u): row['ratio_sim_over_garboczi']=phi/G[u]
        out.append(row); print(row,flush=True)
json.dump(out,open('ell_batch.json','w'),indent=1)
