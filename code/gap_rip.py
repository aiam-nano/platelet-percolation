import numpy as np, math, time, sys, pickle
from sim4 import *
RCAP=3.5
CASES={"P1":(0.3,2.0,1.0),"P2":(0.5,2.0,1.0),"P3":(0.4,1.5,0.7),"P4":(0.5,3.0,1.5)}
def paired(L,sig,a,lam,nreal,seed,nfac=1.5,validate=False):
    shF=make_sheet('flat',1.0,h=1/8); shR=make_sheet('ripple',1.0,h=1/8,amp=a*lam/(2*math.pi),lam=lam)
    rng=np.random.default_rng(seed); out=[]
    for _ in range(nreal):
        n=int(nfac*L**3); R=rng.lognormal(-0.5*sig**2,sig,size=n); R=R[R<=RCAP]; n=len(R); R=R/np.mean(R)*1.0 if False else R
        pos=rng.random((n,3))*L; rot=random_rotations(n,rng); z=np.zeros(n)
        bF=shF['brad']*R; bR=shR['brad']*R
        nF=run_once(2,pos,rot,R,z,L,bF,sheet=shF); nR=run_once(2,pos,rot,R,z,L,bR,sheet=shR)
        if validate:
            n0=run_once(0,pos,rot,R,z,L,R*1.0001)
            out.append((n0/L**3 if n0>0 else np.nan, nF/L**3 if nF>0 else np.nan))
        elif nF>0 and nR>0: out.append((nF/L**3,nR/L**3))
    return np.array(out)
if __name__=="__main__":
    tag=sys.argv[1]; L=float(sys.argv[2]); nreal=int(sys.argv[3]); ids=sys.argv[4:]; res={}; t0=time.time()
    for i in ids:
        sig,a,lam=CASES[i]; o=paired(L,sig,a,lam,nreal,9000+int(i[1:]))
        r=o[:,1]/o[:,0]; res[i]=o
        print(f"{i} sig={sig} a={a} lam={lam} L={L} n={len(o)}: flat n_c={o[:,0].mean():.4f} ripple n_c={o[:,1].mean():.4f} ratio(of means)={o[:,1].mean()/o[:,0].mean():.4f} paired-mean {r.mean():.4f}+-{r.std()/math.sqrt(len(r)):.4f} [{time.time()-t0:.0f}s]",flush=True)
        pickle.dump(res,open(f"gaprip_{tag}.pkl","wb"))
