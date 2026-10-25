import numpy as np, math, time, sys, pickle
from numba import njit
from sim2 import *

@njit(cache=True)
def shell_P2(pos, nrm, L, rmax):
    N=pos.shape[0]; s=0.0; c=0
    for i in range(N):
        for j in range(i+1,N):
            dx=pos[j,0]-pos[i,0]; dy=pos[j,1]-pos[i,1]; dz=pos[j,2]-pos[i,2]
            dx-=L*round(dx/L); dy-=L*round(dy/L); dz-=L*round(dz/L)
            if dx*dx+dy*dy+dz*dz<rmax*rmax:
                cg=nrm[i,0]*nrm[j,0]+nrm[i,1]*nrm[j,1]+nrm[i,2]*nrm[j,2]
                s+=1.5*cg*cg-0.5; c+=1
    return s/max(c,1), c

def fisher(kappa, size, rng):
    u=rng.random(size)
    if kappa<1e-9: return 1-2*u
    return 1+np.log(u+(1-u)*np.exp(-2*kappa))/kappa

def aligned_normals(pos, L, ell, kappa, rng):
    n=len(pos)
    if kappa<1e-9:
        v=rng.normal(size=(n,3)); return v/np.linalg.norm(v,axis=1)[:,None]
    nseed=max(1,int(round(L**3/ell**3)))
    seeds=rng.random((nseed,3))*L
    dirs=rng.normal(size=(nseed,3)); dirs/=np.linalg.norm(dirs,axis=1)[:,None]
    idx=np.empty(n,int)
    for a in range(0,n,4000):
        d=pos[a:a+4000,None,:]-seeds[None,:,:]; d-=L*np.round(d/L)
        idx[a:a+4000]=np.argmin((d**2).sum(2),axis=1)
    c0=dirs[idx]
    ct=fisher(kappa,n,rng); st=np.sqrt(np.maximum(0,1-ct**2)); ph=rng.random(n)*2*np.pi
    a=np.where(np.abs(c0[:,[0]])<0.9,np.array([[1.0,0,0]]),np.array([[0,1.0,0]]))
    e1=np.cross(c0,a); e1/=np.linalg.norm(e1,axis=1)[:,None]; e2=np.cross(c0,e1)
    return ct[:,None]*c0+st[:,None]*(np.cos(ph)[:,None]*e1+np.sin(ph)[:,None]*e2)

def one(L, ell, kappa, rng):
    n=int(2.5*L**3)
    pos=rng.random((n,3))*L
    nrm=aligned_normals(pos,L,ell,kappa,rng)
    rot=rotation_with_normal(nrm,rng)
    rad=np.ones(n)
    Nc=run_once(0,pos,rot,rad,np.zeros(n),L,rad)
    if Nc<0: return None
    S,cnt=shell_P2(pos[:Nc],nrm[:Nc],L,2.0)
    Q=1.5*(nrm[:Nc,:,None]*nrm[:Nc,None,:]).mean(0)-0.5*np.eye(3)
    Sm=np.linalg.eigvalsh(Q).max()
    return Nc/L**3, S, Sm

if __name__=="__main__":
    L=24.0; out={}
    t0=time.time()
    for ell in [3.0,6.0]:
        for kappa in [0.0,1.0,3.0,10.0,30.0]:
            rng=np.random.default_rng(int(ell*100+kappa*7))
            res=[one(L,ell,kappa,rng) for _ in range(20 if kappa>0 else 20)]
            res=[r for r in res if r]
            a=np.array(res)
            out[(ell,kappa)]=a
            print(f"ell={ell} kappa={kappa}: n_c R^3={a[:,0].mean():.4f}  S_loc={a[:,1].mean():.3f}  S_macro={a[:,2].mean():.3f}  [{time.time()-t0:.0f}s]",flush=True)
            pickle.dump(out,open("align.pkl","wb"))
