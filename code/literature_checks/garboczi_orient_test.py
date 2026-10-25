"""Hypothesis (Yi and Sastry 2004): Garboczi et al. 1995 sampled particle orientations with a density that was not uniform on the sphere.
Test: thin oblate spheroids whose symmetry-axis polar angle is drawn uniformly in theta (not in cos theta)."""
import numpy as np, math, sys
from sim2e import *
def run(u,L,nreal,seed,mode):
    rng=np.random.default_rng(seed); out=[]
    for _ in range(nreal):
        n=int(0.5*L**3)
        pos=rng.random((n,3))*L
        if mode=='uniform': rot=random_rotations(n,rng)
        else:
            th=rng.random(n)*math.pi if mode=='theta' else None
            ph=rng.random(n)*2*math.pi
            nrm=np.stack([np.sin(th)*np.cos(ph),np.sin(th)*np.sin(ph),np.cos(th)],1)
            rot=rotation_with_normal(nrm,rng)
        rad=np.ones(n); hth=np.full(n,1.0/u); brad=np.sqrt(rad**2+hth**2)
        Nc=run_once(3,pos,rot,rad,hth,L,brad)
        if Nc>0: out.append(Nc/L**3)
    return np.array(out)
for u in (10,100,1000):
    for mode in ('uniform','theta'):
        nc=run(u,30,12,321+u,mode); vol=4/3*math.pi/u
        phi=nc*vol; occ=1-np.exp(-phi)
        print(f"u={u} orientation={mode}: occupied phi_c={occ.mean():.5g}+-{occ.std()/math.sqrt(len(occ)):.2g}  phi*u={ (phi*u).mean():.4f}",flush=True)
