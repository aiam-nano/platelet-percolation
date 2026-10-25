import numpy as np, math
from sim4 import field_rho
HG=1/16; GH=21/16
def make_field(seed, s_rms, p=4, Lbox=4.0, lam_min=0.6, lam_max=2.0):
    N=int(round(Lbox/HG)); rng=np.random.default_rng(seed)
    k=2*math.pi*np.fft.fftfreq(N,d=HG); KX,KY=np.meshgrid(k,k); q=np.hypot(KX,KY)
    mask=(q>=2*math.pi/lam_max)&(q<=2*math.pi/lam_min)
    amp=np.where(mask,q**(-p/2.0),0.0)
    coef=amp*(rng.normal(size=(N,N))+1j*rng.normal(size=(N,N)))
    f=np.real(np.fft.ifft2(coef))
    f=np.fft.fftshift(f)            # index N/2 <-> x=0
    c=N//2; m=21
    g=f[c-m:c+m+1,c-m:c+m+1].copy()   # rows = y, cols = x ; (2m+1)=43 points spanning [-GH,GH]
    n=g.shape[0]; ax=-GH+HG*np.arange(n); X,Y=np.meshgrid(ax,ax); disc=(X*X+Y*Y<=1.0)
    gy,gx=np.gradient(g,HG); s=math.sqrt(np.mean((gx*gx+gy*gy)[disc]))
    g*=s_rms/s; g-=g[disc].mean()
    return g
def field_stats(g):
    n=g.shape[0]; ax=-GH+HG*np.arange(n); X,Y=np.meshgrid(ax,ax)
    rho=field_rho(g,HG,GH,1.0); foot=(X*X+Y*Y<=rho*rho)
    h=g[foot]-g[foot].mean(); hrms=float(np.sqrt(np.mean(h*h)))
    return rho,hrms
def k_rf(rho,hrms): return 1+(rho**-3-1)*math.exp(-5.75*hrms)
if __name__=="__main__":
    for nm,p,s in [("R1",4,0.4),("R2",4,0.8),("R3",4,1.2),("R4",2,0.8)]:
        ks=[];rs=[];hs=[]
        for sd in range(1,11):
            g=make_field(1000*int(nm[1:])+sd,s,p); rho,h=field_stats(g); ks.append(k_rf(rho,h)); rs.append(rho); hs.append(h)
        print(nm,f"p={p} s_rms={s}: rho mean {np.mean(rs):.3f} h_rms mean {np.mean(hs):.3f} | K_rf mean {np.mean(ks):.3f} (range {min(ks):.3f}-{max(ks):.3f})")
