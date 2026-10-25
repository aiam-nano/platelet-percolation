import numpy as np, math, time, sys, pickle
from sim2 import *
from gap_lib import *
def one(L,cid,P,rng):
    kind,kw,A,beta,nrm,mb=P
    n=int(2.5*L**3); R=sample(kind,n,rng,**kw)*nrm; n=len(R)
    t=(2.0/A)*(R**beta/mb if beta>0 else np.ones(n)); hth=t/2
    pos=rng.random((n,3))*L; rot=random_rotations(n,rng)
    brad=np.sqrt(R**2+hth**2)
    Nc=run_once(1,pos,rot,R,hth,L,brad)
    if Nc<0: return None
    return Nc/L**3*math.pi*np.mean(R[:Nc]**2*t[:Nc]), Nc/L**3, np.mean(R[:Nc])
if __name__=="__main__":
    tag=sys.argv[1]; L=float(sys.argv[2]); nreal=int(sys.argv[3]); seed0=int(sys.argv[4]); ids=sys.argv[5:]; res={}; t0=time.time()
    for cid in ids:
        P=cyl_params(cid); rng=np.random.default_rng(seed0+sum(map(ord,cid))); out=[]
        for _ in range(nreal):
            r=one(L,cid,P,rng)
            if r: out.append(r)
            if len(out)==1 and _==0: print(f"  {cid} first realisation done [{time.time()-t0:.0f}s]",flush=True)
        a=np.array(out); res[cid]=a
        print(f"{cid} L={L} n={len(a)}: phi_c={100*a[:,0].mean():.4f}% +- {100*a[:,0].std()/math.sqrt(len(a)):.4f} [{time.time()-t0:.0f}s]",flush=True)
        pickle.dump(res,open(f"gapcyl_{tag}.pkl","wb"))
