import numpy as np, math
from scipy.optimize import least_squares
from scipy.special import ellipe
from mylaw import *
FS=1.015
# ---- baseline
base=[(5,.88),(10,1.11),(20,1.27),(50,1.38),(100,1.43),(200,1.45),(500,1.47),(33.333,1.340),(30,1.310),(10,1.093),(50,1.3645),(100,1.400)]
def phiA(Aeff,b): return c0*(1+b/(Aeff+1.5))*Aeff/(Aeff+c1+c2/Aeff)
def res_base(p):
    return [ (phiA(A,p[0])/(y/FS)-1)/0.015 for A,y in base]
rb=least_squares(res_base,[0.5]); b=rb.x[0]
print("baseline: b = %.3f"%b)
for A,y in base: print(f"  A={A}: sim/FS {y/FS:.3f}  law {phiA(A,b):.3f}  ({100*(phiA(A,b)/(y/FS)-1):+.1f}%)")
def base_phi(Aeff): return c0*(1+b/(Aeff+1.5))/(Aeff+c1+c2/Aeff)
def dil(Aeff,K): return (Aeff+c1+c2/Aeff)/(Aeff/K+c1+c2/Aeff)   # multiplicative effect of anisotropy K on threshold
# ---- polydispersity
pol=[(1e6,.1,1.007,.006),(1e6,.2,1.021,.006),(1e6,.3,1.053,.008),(1e6,.4,1.099,.008),(1e6,.45,1.136,.012)]  # B_c ratios (thin, eigen)
fin=[(33.333,.3,0.843,.02),(30,.4,0.797,.02),(10,.3,0.925,.02)]
def thinP0(s): return P(s,1e6,0.0)
def res_pol(p):
    r=[]
    for A,s,Br,e in pol: r.append((1+p[0]*s*s)/Br-1)
    for A,s,y,e in fin: r.append((P(s,A,p[0])/y-1)/e*0.01*0+ (P(s,A,p[0])/y-1)/0.02)
    return [x/ (0.01 if i<len(pol) else 1) if False else x for i,x in enumerate(r)]
def res_pol2(p):
    r=[((1+p[0]*s*s)/Br-1)/(e/Br) for A,s,Br,e in pol]
    r+=[(P(s,A,p[0])/y-1)/0.02 for A,s,y,e in fin]
    return r
rp=least_squares(res_pol2,[0.5]); kappa=rp.x[0]; print("kappa = %.3f"%kappa)
for A,s,y,e in fin: print(f"  poly A={A} s={s}: sim {y:.3f} law {P(s,A,kappa):.3f} ({100*(P(s,A,kappa)/y-1):+.1f}%)")
for A,s,Br,e in pol: print(f"  thin s={s}: B ratio {Br:.3f} law {1+kappa*s*s:.3f}")
# ---- alignment data: (Aeff, d, S, ratio, err)
al=[(33.333,2,.119,1.0995,.012),(33.333,2,.146,1.1358,.012),(33.333,2,.166,1.141,.012),(33.333,2,.193,1.177,.012),
    (33.333,4,.282,1.206,.012),(33.333,4,.314,1.2124,.012),
    (33.333,6,.180,1.0858,.012),(33.333,6,.268,1.1761,.012),(33.333,6,.355,1.236,.012),
    (33.333,8,.252,1.1431,.012),(33.333,8,.398,1.268,.012),(33.333,8,.478,1.3736,.012),
    (30,4,.239,1.1672,.012),(30,4,.283,1.1932,.012),(30,8,.339,1.2216,.012),(30,8,.397,1.263,.012),(30,8,.476,1.3587,.012),
    (10,4,.181,1.083,.012),(10,4,.236,1.1076,.012),
    (1e6,3,.044,1.0425,.01),(1e6,3,.217,1.188,.01),(1e6,3,.330,1.361,.01),(1e6,6,.073,1.0588,.01),(1e6,6,.356,1.276,.01)]
def Kal(S,d,g,pw): return 1/(1-0.625*S*(1+g/(d/2)**pw))
def res_al(p): return [(dil(A,Kal(S,d,p[0],p[1]))/r-1)/(e/r) for A,d,S,r,e in al]
ra=least_squares(res_al,[0.2,1.0],bounds=([0,0.2],[3,4])); gam,pw=ra.x; print("alignment: gamma=%.3f p=%.3f"%(gam,pw))
for A,d,S,r,e in al: print(f"  align A={A:g} d={d} S={S}: sim {r:.3f} law {dil(A,Kal(S,d,gam,pw)):.3f} ({100*(dil(A,Kal(S,d,gam,pw))/r-1):+.1f}%)")
# ---- ripple
def kx(a):
    th=np.linspace(0,np.pi/2,20001); return 1/np.mean(np.sqrt(1+a*a*np.cos(th)**2))
def Kell(a):
    k=kx(a); return math.pi/(2*k*ellipe(1-k*k))
rip=[(1e6,.5,1,1.054,.012),(1e6,1,1,1.132,.012),(1e6,1.5,1,1.245,.02),(1e6,2,1,1.314,.012),(1e6,1.5,.5,1.405,.04),(1e6,1.5,2,1.073,.012),
     (20,1,1,1.099,.015),(50,1.5,1,1.231,.012),(20,1.5,2,1.010,.012),(100,1,1.5,1.100,.012),(20,2,1,1.205,.012)]
def Krip(a,lam,l0,pw2,amp): return 1+amp*(Kell(a)-1)/(1+(lam/l0)**pw2)
def res_rip(p): return [(dil(A,Krip(a,l,p[0],p[1],p[2]))/r-1)/(e/r) for A,a,l,r,e in rip]
rr=least_squares(res_rip,[0.7,2.0,1.0],bounds=([0.1,0.5,0.3],[3,6,2]))
l0,p2,amp=rr.x; print("ripple: lam0=%.3f p=%.3f amp=%.3f"%(l0,p2,amp))
for A,a,l,r,e in rip: print(f"  ripple A={A:g} a={a} lam={l}: sim {r:.3f} law {dil(A,Krip(a,l,l0,p2,amp)):.3f} ({100*(dil(A,Krip(a,l,l0,p2,amp))/r-1):+.1f}%)")
np.save('fitpars.npy',np.array([b,kappa,gam,pw,l0,p2,amp]))
