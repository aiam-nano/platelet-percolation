import math, json, numpy as np
import finallaw2 as L2
from finallaw2 import c0,c1,c2,B_,KAPPA
from mylaw import P as Pmy
pi=math.pi
# ---------- samplers (deterministic big samples for moments) ----------
def sample(kind, n, rng, **kw):
    if kind=='lognormal':
        s=kw['sig']; R=rng.lognormal(-0.5*s*s,s,size=n); R=R[R<=kw.get('rcap',6.0)]
    elif kind=='bimodal':
        lo,hi,p=kw['lo'],kw['hi'],kw['p']; R=np.where(rng.random(n)<p,lo,hi)
    elif kind=='pareto':
        a=kw['alpha']; R=(1-rng.random(n))**(-1/a); R=R[R<=kw['rcap']]
    return R
# normalisation constants so that the number mean is 1 for continuous cases
def normaliser(kind,**kw):
    r=np.random.default_rng(12345); R=sample(kind,2_000_000,r,**kw); return 1.0/R.mean()
def Vex(Ri,ti,Rj,tj):
    Vi=pi*Ri**2*ti; Vj=pi*Rj**2*tj
    Si=2*pi*Ri**2+2*pi*Ri*ti; Sj=2*pi*Rj**2+2*pi*Rj*tj
    Mi=pi*ti+pi**2*Ri; Mj=pi*tj+pi**2*Rj
    return Vi+Vj+(Si*Mj+Sj*Mi)/(4*pi)
def eig_ratio(R,t):
    """(<Vp>/lam_max)/(Vp0/lam0) relative to a monodisperse cylinder with R=1 and t=<t>."""
    n=len(R); m=min(n,700); idx=np.minimum(((np.arange(m)+0.5)/m*n).astype(int),n-1)
    order=np.argsort(R); Rs=R[order][idx]; ts=t[order][idx]    # stratified representatives (sorted by R)
    K=Vex(Rs[:,None],ts[:,None],Rs[None,:],ts[None,:])/len(Rs)
    lam=max(np.linalg.eigvals(K).real); Vp=np.mean(pi*Rs**2*ts)
    tb=ts.mean(); Rb=R.mean(); lam0=Vex(Rb,tb,Rb,tb); Vp0=pi*Rb**2*tb
    return (Vp/lam)/(Vp0/lam0)
def pred(R,t,drift=True):
    sl=np.std(np.log(R)); tb=t.mean(); A_ref=2.0/tb
    r=eig_ratio(R,t)*((1+KAPPA*sl**2) if drift else 1.0)
    return 100*L2.phi(A=A_ref)*r, r, sl
out={}
# sanity: beta=0 lognormal against the existing law
rng=np.random.default_rng(5)
for sig,A in [(0.3,100),(0.5,60),(0.45,100)]:
    R=sample('lognormal',400000,rng,sig=sig); R=R/R.mean(); t=np.full_like(R,2.0/A)
    p,r,sl=pred(R,t); print(f"sanity lognormal sig {sig} A {A}: general {p:.3f}%  finallaw2 {100*L2.phi(A=A,sig=sig):.3f}%")
# ---------- G1 averaging conversion ----------
print("\nG1 apparent excess  phi_c*A_x/1.442 for number / area-weighted D[3,2] / volume-weighted D[4,3]")
for A,s in [(100,0.3),(100,0.5),(30,0.5)]:
    ph=L2.phi(A=A,sig=s)
    row=[100*ph]+[ph*A*math.exp(k*s*s)/c0 for k in (0,2,3)]
    print(f"G1 A_n={A} sig={s}: law phi {row[0]:.3f}% | excess vs 1.442/A_x: number {row[1]:.3f}, area-wt {row[2]:.3f}, vol-wt {row[3]:.3f}"); out[f"G1_{A}_{s}"]=row
# ---------- G2 non-lognormal ----------
print("\nG2 non-lognormal distributions, A=60 common thickness")
G2={"B1":("bimodal",dict(lo=2/3,hi=2.0,p=0.75)),"B2":("bimodal",dict(lo=0.588235,hi=3.333333,p=0.85)),
    "H1":("lognormal",dict(sig=0.8,rcap=5.0)),"Q1":("pareto",dict(alpha=3.5,rcap=5.0))}
for nm,(kind,kw) in G2.items():
    rng=np.random.default_rng(99); R=sample(kind,600000,rng,**kw)
    if kind!='bimodal': R=R*normaliser(kind,**kw)
    R=R/R.mean(); t=np.full_like(R,2.0/60)
    p,r,sl=pred(R,t); pn,_,_=pred(R,t,drift=False)
    print(f"{nm} {kind} {kw}: mean R {R.mean():.3f} sigma_ln {sl:.3f} m2 {np.mean(R**2):.3f} | pred {p:.3f}% (no-drift {pn:.3f}%) eigen ratio {r:.4f}"); out[nm]=dict(pred=p,nodrift=pn,sl=sl)
# ---------- G3 thickness-size covariance ----------
print("\nG3 allometric covariance t ~ R^beta (number-mean thickness 2/A)")
for nm,sig,beta,A in [("C1",0.4,1.0,60),("C2",0.5,0.5,40),("C3",0.3,1.0,20),("C4",0.3,2.0,100)]:
    rng=np.random.default_rng(7); R=sample('lognormal',600000,rng,sig=sig); R=R/R.mean()
    t=R**beta; t=t*(2.0/A)/t.mean()
    p,r,sl=pred(R,t); pn,_,_=pred(R,t,drift=False)
    # comparison: same thickness mean but no covariance
    p0,_,_=pred(R,np.full_like(R,2.0/A))
    print(f"{nm} sig {sig} beta {beta} A {A}: pred {p:.3f}% (no-drift {pn:.3f}%) | uncorrelated-thickness law {p0:.3f}% | law finallaw2 sig {100*L2.phi(A=A,sig=sig):.3f}%"); out[nm]=dict(pred=p,nodrift=pn,p0=p0)
# ---------- G5 thick biaxial ----------
from new_preds import L2 as _L2
from scipy.optimize import brentq
def rho_of(a,lam):
    A_=a*lam/(2*pi); k=2*pi/lam; n=1201; g=np.linspace(-1.3,1.3,n); X,Y=np.meshgrid(g,g); dA=(g[1]-g[0])**2
    s2=(A_/math.sqrt(2)*k)**2*(np.cos(k*X)**2+np.cos(k*Y)**2); dens=np.sqrt(1+s2); r2=X*X+Y*Y
    return brentq(lambda rr:(dens*(r2<=rr*rr)).sum()*dA-pi,0.5,1.08)
def phiK(A,K):
    base=c0*(1+B_/(A+1.5))/(A+c1+c2/A); dil=(A+c1+c2/A)/(A/K+c1+c2/A); return base*dil
print("\nG5 thick biaxial: ratio to flat slab")
for nm,a,lam,A,M in [("TB1",1.5,1.0,40,5),("TB2",2.0,0.8,25,6),("TB3",1.0,0.7,60,4),("TB4",2.5,1.2,30,5),("TB5",3.0,0.7,50,4)]:
    rho=rho_of(a,lam); K=1+(rho**-3-1)*math.exp(-0.647*a*lam); r=phiK(A,K)/phiK(A,1.0)
    print(f"{nm} a {a} lam {lam} A {A} M {M}: rho {rho:.4f} K_thin {K:.3f} -> thick ratio pred {r:.3f}"); out[nm]=dict(rho=rho,K=K,ratio=r)
# ---------- G4 polydispersity x ripple (multiplicative hypothesis) ----------
print("\nG4 ripple x polydispersity (self-similar sheets): predicted ratio = K_rip(a,lam)")
for nm,s,a,lam in [("P1",0.3,2.0,1.0),("P2",0.5,2.0,1.0),("P3",0.4,1.5,0.7),("P4",0.5,3.0,1.5)]:
    print(f"{nm} sig {s} a {a} lam {lam}: K_rip {L2.Krip(a,lam):.3f}"); out[nm]=L2.Krip(a,lam)
json.dump(out,open('gap_preds.json','w'),indent=1,default=float)
