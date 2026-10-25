import math, numpy as np
from scipy.optimize import brentq
import finallaw2 as L2
pi=math.pi
# ---------- tactoid-number polydispersity ----------
def dist(kind):
    if kind=='geo3':  p=1/3; ks=np.arange(1,41); w=p*(1-p)**(ks-1)
    elif kind=='bi19': ks=np.array([1,9]); w=np.array([.5,.5])
    elif kind=='uni7': ks=np.arange(1,8); w=np.ones(7)
    elif kind=='geo5': p=0.2; ks=np.arange(1,61); w=p*(1-p)**(ks-1)
    elif kind=='bi14': ks=np.array([1,4]); w=np.array([.7,.3])
    return ks, w/w.sum()
TACT=[("T1",100,'geo3'),("T2",100,'bi19'),("T3",60,'uni7'),("T4",200,'geo5'),("T5",30,'bi14')]
def eig_factor(A,ks,w):
    D=2.0; t0=2.0/A; t=ks*t0
    a0=pi**2*D**3/8; a1=pi*(pi+3)*D**2/8; a2=pi*D/2
    C=np.array([[a0,a1],[a1,a2]]); m1=(w*t).sum(); m2=(w*t*t).sum(); G=np.array([[1,m1],[m1,m2]])
    lam=max(np.linalg.eigvals(C@G).real); tb=m1
    return (a0+2*a1*tb+a2*tb*tb)/lam
if __name__=="__main__":
    print("tactoid cases: id A dist Nbar Nw | eig factor | pred (number-average rule) | rival (weight-average rule)")
    for i,A,k in TACT:
        ks,w=dist(k); Nb=(w*ks).sum(); Nw=(w*ks**2).sum()/Nb; f=eig_factor(A,ks,w)
        p=100*L2.phi(A=A,N=Nb)*f; r=100*L2.phi(A=A,N=Nw)*f
        print(i,A,k,f"Nbar {Nb:.3f} Nw {Nw:.3f} | f_eig {f:.4f} | pred {p:.3f}% | rival {r:.3f}%")
    # ---------- biaxial egg-box ----------
    def rho_of(a,lam):
        A=a*lam/(2*pi); k=2*pi/lam; n=1201; g=np.linspace(-1.2,1.2,n); X,Y=np.meshgrid(g,g); dA=(g[1]-g[0])**2
        s2=(A/math.sqrt(2)*k)**2*(np.cos(k*X)**2+np.cos(k*Y)**2); dens=np.sqrt(1+s2); r2=X*X+Y*Y
        f=lambda rho: (dens*(r2<=rho*rho)).sum()*dA-pi
        return brentq(f,0.5,1.0)
    print("\nbiaxial cases: id a lam | rho | H1 uniaxial K_rip | Pb = 1+(rho^-3-1)exp(-0.647 a lam)")
    for i,a,lam in [("W1",1.0,1.0),("W2",2.0,1.0),("W3",1.5,0.7),("W4",3.0,1.5),("W5",0.7,2.0),("W6",2.5,0.6)]:
        rho=rho_of(a,lam); h1=L2.Krip(a,lam); pb=1+(rho**-3-1)*math.exp(-0.647*a*lam)
        print(i,a,lam,f"| rho {rho:.4f} | H1 {h1:.3f} | Pb {pb:.3f}")
