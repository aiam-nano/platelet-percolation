import numpy as np, math, time, sys, pickle
from sim4 import *
from fieldgen import make_field, field_stats, k_rf, HG, GH
CASES={"R1":(4,0.4),"R2":(4,0.8),"R3":(4,1.2),"R4":(2,0.8)}
def reals(sh,L,nreal,rng,nfac=1.5):
    out=[]
    for _ in range(nreal):
        n=int(nfac*L**3); pos=rng.random((n,3))*L; rot=random_rotations(n,rng)
        Nc=run_once(2,pos,rot,np.ones(n),np.zeros(n),L,np.full(n,sh['brad']),sheet=sh)
        if Nc>0: out.append(Nc/L**3)
    return out
if __name__=="__main__":
    tag=sys.argv[1]; L=float(sys.argv[2]); nf=int(sys.argv[3]); nr=int(sys.argv[4]); ids=sys.argv[5:]; res={}; t0=time.time()
    if 'flat' in ids:
        fl=np.array(reals(make_sheet('flat',1.0,h=1/8),L,nf*nr*2,np.random.default_rng(4242))); res['flat']=fl
        print(f"flat n_c={fl.mean():.4f}+-{fl.std()/math.sqrt(len(fl)):.4f} n={len(fl)} [{time.time()-t0:.0f}s]",flush=True)
    for i in ids:
        if i=='flat': continue
        p,s=CASES[i]; per=[]; ks=[]
        for sd in range(1,nf+1):
            g=make_field(1000*int(i[1:])+sd,s,p); rho,h=field_stats(g); ks.append(k_rf(rho,h))
            sh=make_sheet('field',1.0,h=1/8,fgrid=g,hg=HG,Gh=GH)
            per.append(np.mean(reals(sh,L,nr,np.random.default_rng(7000+1000*int(i[1:])+sd))))
        res[i]=dict(per=np.array(per),kpred=np.array(ks))
        print(f"{i} p={p} s={s}: mean n_c {np.mean(per):.4f} fields={len(per)} [{time.time()-t0:.0f}s]",flush=True)
        pickle.dump(res,open(f"gapfield_{tag}.pkl","wb"))
