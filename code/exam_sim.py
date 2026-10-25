# Sealed ground truth for the baseline-arms forecasting exam (cases F1-F12). Run before any chat output exists.
import numpy as np, math, time, sys, pickle, json
from sim2 import *
from align import aligned_normals
import gap_rip
KAP={0.5:4.733,0.6:6.312}
def lognormal_R(n,s,rng,rcap=6.0):
    R=rng.lognormal(-0.5*s*s,s,size=n); return R[R<=rcap]
def build(case,n,rng):
    """returns R, t (full thickness), kappa, ell"""
    kap=0.0; ell=6.0
    if case=="F1": R=np.ones(n); t=np.full(n,2/50)
    elif case=="F2": R=np.ones(n); t=np.full(n,2/8)
    elif case=="F3": R=np.ones(n); t=np.full(n,4*2/60)
    elif case=="F4": R=lognormal_R(n,.4,rng); t=np.full(len(R),2/80)
    elif case=="F5":
        R=lognormal_R(n,.5,rng); th=R**0.5; t=(2/40)*th/np.mean(th)
    elif case=="F6":
        R=lognormal_R(n,.5,rng); z=np.log(R); beta=0.8*0.7/0.5
        lt=beta*z+math.sqrt(1-0.64)*0.7*rng.normal(size=len(R)); th=np.exp(lt); t=(2/60)*th/np.mean(th)
    elif case=="F7":
        R=np.where(rng.random(n)<0.7,0.5,1.5); t=np.full(n,0.02)
    elif case=="F8": R=np.ones(n); t=np.full(n,2/40); kap=KAP[0.6]
    elif case=="F9": R=np.ones(n); t=np.full(n,3*2/60); kap=KAP[0.6]
    elif case=="F12":
        R=lognormal_R(n,.3,rng); t=np.full(len(R),3*2/60); kap=KAP[0.5]
    return R,t,kap,ell
def one(case,L,nfac,rng):
    n=int(nfac*L**3); R,t,kap,ell=build(case,n,rng); n=len(R)
    pos=rng.random((n,3))*L
    if kap>0:
        nrm=aligned_normals(pos,L,ell,kap,rng); rot=rotation_with_normal(nrm,rng)
    else: rot=random_rotations(n,rng)
    hth=t/2; brad=np.sqrt(R**2+hth**2)
    Nc=run_once(1,pos,rot,R,hth,L,brad)
    if Nc<0: return None
    return Nc/L**3*np.mean(math.pi*R[:Nc]**2*t[:Nc]), Nc/L**3
if __name__=="__main__":
    tag=sys.argv[1]; cases=sys.argv[2:]; res={}; t0=time.time()
    for c in cases:
        rng=np.random.default_rng(31000+sum(map(ord,c)))
        if c in ("F10","F11"):
            gap_rip.RCAP=6.0
            sig=0.0 if c=="F10" else 0.4
            o=gap_rip.paired(32.0,sig,1.2,1.0,48,31000+sum(map(ord,c)))
            res[c]=o; r=o[:,1].mean()/o[:,0].mean()
            print(f"{c} ratio rippled/flat n_c = {r:.4f} (n={len(o)}) [{time.time()-t0:.0f}s]",flush=True)
        else:
            L=32.0; nfac=2.5; out=[]
            for _ in range(24):
                r=one(c,L,nfac,rng)
                if r: out.append(r)
            a=np.array(out); res[c]=a; x=a[:,0].mean()/1.015
            print(f"{c} n={len(a)} nV={a[:,0].mean():.5f} phi={1-math.exp(-x):.5f} [{time.time()-t0:.0f}s]",flush=True)
        pickle.dump(res,open(f"exam/truth_{tag}.pkl","wb"))
