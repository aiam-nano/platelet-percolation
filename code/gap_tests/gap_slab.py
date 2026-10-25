import numpy as np, math, time, sys, pickle
from sim4 import *
CASES={"TB1":(1.5,1.0,40,5),"TB2":(2.0,0.8,25,6),"TB3":(1.0,0.7,60,4),"TB4":(2.5,1.2,30,5),"TB5":(3.0,0.7,50,4)}
def paired(L,a,lam,A,M,nreal,seed,nfac=1.5):
    t=2.0/A
    shF=make_sheet('flat',1.0,h=1/8); shE=make_sheet('egg',1.0,h=1/8,amp=a*lam/(2*math.pi),lam=lam); rho=shE['rho']
    tvF=t; tvE=t/rho**2
    oF=np.linspace(-tvF/2,tvF/2,M); oE=np.linspace(-tvE/2,tvE/2,M)
    bF=np.full(1,0.0); rng=np.random.default_rng(seed); out=[]
    for _ in range(nreal):
        n=int(nfac*L**3); pos=rng.random((n,3))*L; rot=random_rotations(n,rng); rad=np.ones(n); z=np.zeros(n)
        nF=run_once(2,pos,rot,rad,z,L,np.full(n,shF['brad']+tvF/2+1e-6),sheet=shF,offs=oF)
        nE=run_once(2,pos,rot,rad,z,L,np.full(n,shE['brad']+tvE/2+1e-6),sheet=shE,offs=oE)
        if nF>0 and nE>0: out.append((nF/L**3,nE/L**3))
    return np.array(out),rho
if __name__=="__main__":
    tag=sys.argv[1]; L=float(sys.argv[2]); nreal=int(sys.argv[3]); ids=sys.argv[4:]; res={}; t0=time.time()
    for i in ids:
        a,lam,A,M=CASES[i]; o,rho=paired(L,a,lam,A,M,nreal,8000+int(i[2:]))
        r=o[:,1]/o[:,0]; res[i]=o
        print(f"{i} a={a} lam={lam} A={A} M={M} rho={rho:.4f} L={L} n={len(o)}: flat n_c={o[:,0].mean():.4f} egg n_c={o[:,1].mean():.4f} ratio(of means)={o[:,1].mean()/o[:,0].mean():.4f} paired-mean {r.mean():.4f}+-{r.std()/math.sqrt(len(r)):.4f} [{time.time()-t0:.0f}s]",flush=True)
        pickle.dump(res,open(f"gapslab_{tag}.pkl","wb"))
