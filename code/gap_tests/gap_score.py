import pickle, math, numpy as np, glob
FS=1.015
def load(pat):
    d={}
    for f in sorted(glob.glob(pat)): d.update(pickle.load(open(f,'rb')))
    return d
cyl=load('gapcyl_*.pkl')
PRED={"G1a":(1.212,None,6),"G1b":(0.932,None,6),"G1c":(2.985,None,6),
 "B1":(1.554,1.361,8),"B2":(0.979,0.790,8),"H1":(1.276,0.926,8),"Q1":(1.742,1.667,8),
 "C1":(2.396,2.178,6,1.789),"C2":(2.895,2.505,6,2.315),"C3":(6.340,6.003,6,5.506),"C4":(1.699,1.609,6,1.214)}
print("Cylinder tests (sim/1.015 vs pre-registered prediction)")
for k,v in PRED.items():
    if k not in cyl: continue
    a=cyl[k][:,0]*100; m=a.mean()/FS; se=a.std(ddof=1)/math.sqrt(len(a))/FS
    err=100*(v[0]/m-1)  # (pred/sim-1)
    s=f"{k}: sim {m:.3f}+-{se:.3f}% (raw {a.mean():.3f}) pred {v[0]:.3f} err(pred/sim-1) {err:+.1f}% window ±{v[2]} {'PASS' if abs(err)<=v[2] else 'FAIL'}"
    if v[1]: s+=f" | no-drift {v[1]:.3f} ({100*(v[1]/m-1):+.1f}%)"
    if len(v)>3: s+=f" | rival {v[3]:.3f} ({100*(v[3]/m-1):+.1f}%)"
    print(s)
if all(k in cyl for k in ("G1a","G1b","G1c")):
    print("G1 apparent excess (sim/FS): number / area-wt / vol-wt")
    for k,A,s in [("G1a",100,.3),("G1b",100,.5),("G1c",30,.5)]:
        ph=cyl[k][:,0].mean()/FS
        print(k,[round(ph*A*math.exp(j*s*s)/1.4421,3) for j in (0,2,3)])
rip=load('gaprip_*.pkl'); P={"P1":1.293,"P2":1.293,"P3":1.340,"P4":1.107}
print("\nG4 ripple x polydispersity")
for k,v in P.items():
    if k in rip:
        o=rip[k]; r=o[:,1]/o[:,0]; m=r.mean(); se=r.std(ddof=1)/math.sqrt(len(r)); rm=o[:,1].mean()/o[:,0].mean()
        print(f"{k}: paired ratio {m:.3f}+-{se:.3f} (ratio of means {rm:.3f}) pred {v} err(pred/sim-1) {100*(v/rm-1):+.1f}% {'PASS' if abs(v/rm-1)<=0.08 else 'FAIL'}")
sl=load('gapslab_*.pkl'); S={"TB1":(1.244,8),"TB2":(1.344,8),"TB3":(1.206,8),"TB4":(1.228,8),"TB5":(1.543,15)}
print("\nG5 thick biaxial slabs")
for k,(v,w) in S.items():
    if k in sl:
        o=sl[k]; r=o[:,1]/o[:,0]; rm=o[:,1].mean()/o[:,0].mean(); se=r.std(ddof=1)/math.sqrt(len(r))
        print(f"{k}: ratio of means {rm:.3f} (paired {r.mean():.3f}+-{se:.3f}) pred {v} err(pred/sim-1) {100*(v/rm-1):+.1f}% window ±{w} {'PASS' if abs(v/rm-1)<=w/100 else 'FAIL'}")
fd=load('gapfield_*.pkl'); K={"R1":1.080,"R2":1.146,"R3":1.199,"R4":1.200}
print("\nG6 random wrinkle fields")
if 'flat' in fd:
    fl=fd['flat'].mean()
    for k,v in K.items():
        if k in fd:
            per=fd[k]['per']/fl; kp=fd[k]['kpred']; m=per.mean()
            c=np.corrcoef(per,kp)[0,1]
            print(f"{k}: sim ratio {m:.3f}+-{per.std(ddof=1)/math.sqrt(len(per)):.3f} pred (registered) {v}, mean per-field K {kp.mean():.3f} err {100*(v/m-1):+.1f}% {'PASS' if abs(v/m-1)<=0.08 else 'FAIL'} | per-field corr {c:.2f}")
