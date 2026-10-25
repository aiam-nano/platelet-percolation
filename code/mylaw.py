import numpy as np, math
pi=math.pi
c0=1.4421; c1=2*(pi+3)/pi; c2=4/pi
def moments(sig,kmax=6):
    return [math.exp(k*(k-1)*sig**2/2) for k in range(kmax+1)]
def lam_max(sig,A):
    t=1.0/A; m=moments(sig)
    C={(2,1):pi**2/16,(1,2):pi**2/16,(2,0):3*pi*t/8,(0,2):3*pi*t/8,(1,1):pi**2*t/4,(1,0):pi*t*t/4,(0,1):pi*t*t/4}
    M=np.zeros((3,3))
    for (k,l),c in C.items():
        for a in range(3): M[k,a]+=c*m[l+a]
    return max(np.linalg.eigvals(M).real)
def P(sig,A,kappa):
    """polydispersity factor: phi_c(sigma)/phi_c(0) at fixed number-mean diameter/thickness"""
    if sig==0: return 1.0
    t=1.0/A; lm=lam_max(0.0,A); ls=lam_max(sig,A)
    return (lm*math.exp(sig**2)/ls)*(1+kappa*sig**2)
if __name__=="__main__":
    data=[("thin sigma=.2",1e6,.2,None),("thin .3",1e6,.3,None),("thin .4",1e6,.4,None),
          ("s1 A33 s.3",33.3333,.3,0.843),("s3p A30 s.4",30,.4,0.797),("s4p A10 s.3",10,.3,0.925)]
    for nm,A,s,sim in data:
        print(nm, "kappa=0: %.3f  kappa=.57: %.3f"%(P(s,A,0),P(s,A,.57)), "sim",sim)
