import numpy as np, math, time, sys
from sim2e import *
FR=0.6
def run(eps,L,nreal,seed):
    rng=np.random.default_rng(seed); out=[]
    for _ in range(nreal):
        vol=4/3*math.pi*eps; n=int(FR/vol*L**3)   # enough particles to pass threshold (nominal phi ~0.12)
        pos=rng.random((n,3))*L; rot=random_rotations(n,rng)
        rad=np.ones(n); hth=np.full(n,float(eps)); brad=np.full(n,float(eps))
        Nc=run_once(3,pos,rot,rad,hth,L,brad)
        if Nc>0: out.append(Nc/L**3*vol)
    return np.array(out)
for eps,L,nr in ((2,40,12),(5,100,12),(10,120,10),(10,240,8),(20,200,8),(20,400,6)):
    t0=time.time(); ph=run(eps,L,nr,eps*7+L)
    occ=1-np.exp(-ph)
    print(f"prolate eps={eps} L={L} n={len(ph)}: nominal phi_c={ph.mean():.4f}+-{ph.std()/math.sqrt(len(ph)):.4f}  occupied={occ.mean():.4f}  [{time.time()-t0:.0f}s]",flush=True)
