import numpy as np, math, time, sys
from sim2e import *
def run(mode,u,L,nreal,seed,Nfac=0.6):
    rng=np.random.default_rng(seed); out=[]
    for _ in range(nreal):
        # polar semi-axis / half thickness = 1/u  (u = R/b for ellipsoid, = 2R/t for cylinder)
        n=int(Nfac*L**3)
        pos=rng.random((n,3))*L; rot=random_rotations(n,rng)
        rad=np.ones(n); hth=np.full(n,1.0/u); brad=np.sqrt(rad**2+hth**2)
        Nc=run_once(mode,pos,rot,rad,hth,L,brad)
        if Nc>0: out.append(Nc/L**3)
    return np.array(out)
if __name__=="__main__":
    mode=int(sys.argv[1]); u=float(sys.argv[2]); L=float(sys.argv[3]); nr=int(sys.argv[4])
    t0=time.time()
    nc=run(mode,u,L,nr,1000+int(u)+int(L),Nfac=float(sys.argv[5]) if len(sys.argv)>5 else 0.6)
    n=nc.mean(); e=nc.std()/math.sqrt(len(nc))
    if mode==3: vol=4/3*math.pi/u      # spheroid volume with a=1, b=1/u
    else: vol=math.pi*2/u              # cylinder radius 1, thickness 2/u
    phi=n*vol; eta=n*4/3*math.pi       # eta in Yi-Tawerghi convention (4/3 pi n r^3)
    print(f"mode={mode} u={u} L={L} n={len(nc)}  n_c={n:.4f}+-{e:.4f}  phi_c={phi:.5f}  phi_c*u={phi*u:.4f}  eta_c={eta:.4f}  [{time.time()-t0:.0f}s]",flush=True)
