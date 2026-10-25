import numpy as np, math
def sample(kind, n, rng, **kw):
    if kind=='lognormal':
        s=kw['sig']; R=rng.lognormal(-0.5*s*s,s,size=n); R=R[R<=kw.get('rcap',6.0)]
    elif kind=='bimodal':
        lo,hi,p=kw['lo'],kw['hi'],kw['p']; R=np.where(rng.random(n)<p,lo,hi)
    elif kind=='pareto':
        a=kw['alpha']; R=(1-rng.random(n))**(-1/a); R=R[R<=kw['rcap']]
    return R
def normaliser(kind,**kw):
    r=np.random.default_rng(12345); R=sample(kind,2_000_000,r,**kw); return 1.0/R.mean()
# cylinder cases: id -> (kind, kw, A_number_mean, beta)   thickness t_i=(2/A) R_i^beta/<R^beta>
CYL={
 "G1a":("lognormal",dict(sig=0.3,rcap=6.0),100,0.0),
 "G1b":("lognormal",dict(sig=0.5,rcap=6.0),100,0.0),
 "G1c":("lognormal",dict(sig=0.5,rcap=6.0),30,0.0),
 "B1":("bimodal",dict(lo=2/3,hi=2.0,p=0.75),60,0.0),
 "B2":("bimodal",dict(lo=0.588235,hi=3.333333,p=0.85),60,0.0),
 "H1":("lognormal",dict(sig=0.8,rcap=5.0),60,0.0),
 "Q1":("pareto",dict(alpha=3.5,rcap=5.0),60,0.0),
 "C1":("lognormal",dict(sig=0.4,rcap=6.0),60,1.0),
 "C2":("lognormal",dict(sig=0.5,rcap=6.0),40,0.5),
 "C3":("lognormal",dict(sig=0.3,rcap=6.0),20,1.0),
 "C4":("lognormal",dict(sig=0.3,rcap=6.0),100,2.0),
}
def cyl_params(cid):
    kind,kw,A,beta=CYL[cid]
    nrm=1.0 if kind=='bimodal' else normaliser(kind,**kw)
    if kind=='lognormal' and cid.startswith('G1'): nrm=1.0   # as in comb.py
    r=np.random.default_rng(54321); Rr=sample(kind,2_000_000,r,**kw)*nrm
    mb=np.mean(Rr**beta) if beta>0 else 1.0
    return kind,kw,A,beta,nrm,mb
