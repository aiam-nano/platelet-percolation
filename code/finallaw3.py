# finallaw3: revised alignment factor (prefactor 0.588, gamma 0.576, exponent 2) refitted to the 24 alignment points; all other constants as finallaw2. finallaw2.py stays as the record of the law scored on 24 Oct.
import math, numpy as np
from scipy.special import ellipe
from mylaw import c0,c1,c2,P
B_,KAPPA,GAMMA,LAM0=0.663,0.625,0.576,0.660   # GAMMA refitted with the alignment prefactor free (see calibrate_all.py)
def kx(a):
    th=np.linspace(0,np.pi/2,20001); return 1/np.mean(np.sqrt(1+a*a*np.cos(th)**2))
def Kell(a):
    k=kx(a); return math.pi/(2*k*ellipe(1-k*k))
def Kal(S,d): return 1/(1-0.588*S*(1+GAMMA/(d/2)**2)) if S>0 else 1.0
CRIP=0.647
def Krip(a,lam): return 1+(Kell(a)-1)*math.exp(-CRIP*a*lam) if a>0 else 1.0
def phi(A,N=1,sig=0.0,S=0.0,d=4.0,a=0.0,lam=1.0):
    Ae=A/N; K=Kal(S,d)*Krip(a,lam)
    base=c0*(1+B_/(Ae+1.5))/(Ae+c1+c2/Ae)
    dil=(Ae+c1+c2/Ae)/(Ae/K+c1+c2/Ae)
    return base*dil*P(sig,Ae,KAPPA)      # volume fraction
if __name__=="__main__":
    print("Earlier sets (law vs sim/1.015):")
    for nm,args,sim in [("Set1",dict(A=100,N=3,sig=.3),3.391),("Set2",dict(A=100,N=3,S=.3,d=4),4.863),("Set3",dict(A=30,sig=.4,S=.25,d=4),3.823),("Set4",dict(A=20,N=2,sig=.3,S=.2,d=4),10.796),("Set5",dict(A=20,a=1,lam=1),6.98)]:
        p=100*phi(**args); s=sim/1.015; print(f" {nm}: law {p:.3f} sim/FS {s:.3f} law/sim-1 {100*(p/s-1):+.1f}%")
    F=[("F1",dict(A=200,N=5,sig=.5)),("F2",dict(A=15,sig=.2,S=.30,d=3)),("F3",dict(A=60,N=2,S=.35,d=5)),
       ("F4",dict(A=40,a=1.2,lam=1.2,S=.20,d=4)),("F5",dict(A=30,a=3.0,lam=1.0)),("F6",dict(A=100,sig=.45,S=.25,d=6)),
       ("F7",dict(A=24,N=3,sig=.3)),("F8",dict(A=20,a=1.5,lam=0.8))]
    print("Fresh held-out predictions (vol%):")
    for nm,args in F: print(f" {nm} {args}: {100*phi(**args):.3f}")
