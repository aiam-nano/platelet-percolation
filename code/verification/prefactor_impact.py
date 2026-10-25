"""Impact of replacing the alignment factor (prefactor 0.625, gamma 0.392) by the freely refitted (0.588, 0.576), exponent 2.
Run from the workdir: python3 verification/prefactor_impact.py"""
import sys, math, json, numpy as np
sys.path.insert(0, ".")
import finallaw2 as F
from mylaw import c0, c1, c2
def kal_old(S, d): return 1/(1-0.625*S*(1+0.392/(d/2)**2))
def kal_new(S, d): return 1/(1-0.588*S*(1+0.576/(d/2)**2))
def phi_with(kal, A, N=1, sig=0.0, S=0.0, d=4.0):
    Ae = A/N; K = kal(S, d) if S > 0 else 1.0
    base = F.c0*(1+F.B_/(Ae+1.5))/(Ae+c1+c2/Ae); dil = (Ae+c1+c2/Ae)/(Ae/K+c1+c2/Ae)
    return base*dil*F.P(sig, Ae, F.KAPPA)
print("Ratio phi_new/phi_old - 1 (per cent) over a grid; rows S_local, columns d (radii); aspect ratio A=40 and A=1e6")
for A in (40, 1e6):
    print(f"A={A:g}")
    print("   S\\d " + " ".join(f"{d:6g}" for d in (2, 3, 4, 6, 8)))
    for S in (0.05, 0.1, 0.2, 0.3, 0.4, 0.5):
        print(f"  {S:4.2f} " + " ".join(f"{100*(phi_with(kal_new,A,S=S,d=d)/phi_with(kal_old,A,S=S,d=d)-1):+6.2f}" for d in (2, 3, 4, 6, 8)))
# exam
t = json.load(open("exam/truth_sealed.json"))
Sloc = {0.5: 0.1662, 0.6: 0.2394}
print("\nExam aligned cases, law error vs truth with old and new alignment factor (S_local from generator, d=6):")
for k, (A, N, sig, Sd) in {"F8": (40, 1, 0, 0.6), "F9": (60, 3, 0, 0.6), "F12": (60, 3, 0.3, 0.5)}.items():
    r = []
    for kal in (kal_old, kal_new):
        occ = 1-math.exp(-phi_with(kal, A, N, sig, Sloc[Sd], 6)); r.append(100*(occ/t[k]-1))
    print(f"  {k}: old {r[0]:+.1f}%  new {r[1]:+.1f}%")
# calibration fit residuals for alignment points (24) old vs new
exec(open("verification/refit_alignment_prefactor.py").read().split("def dil")[0].split("from mylaw")[0].replace("sys.path.insert(0, \"..\")",""))
