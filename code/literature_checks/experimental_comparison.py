"""Compare the frozen-candidate law (finallaw3) with experimental thresholds in the papers folder.
Law = geometric, contact-defined, random orientation, flat discs, nominal volume fraction nV (occupied = 1-exp(-nV)).
All thresholds converted to volume per cent. Ratio = measured / law.  >1 means measured above geometric prediction."""
import sys, csv, math
sys.path.insert(0,'..'); sys.path.insert(0,'.')
from finallaw3 import phi
from scipy.optimize import brentq

def occ(x): return 1-math.exp(-x)
def wt2vol(w, rho_f=2.2, rho_m=1.2):
    return 100*(w/rho_f)/((w/rho_f)+(100-w)/rho_m)

# name, paper, A_lo, A_hi, sigma, thr_lo(vol%), thr_hi(vol%), kind, note
cases=[
("Few-layer graphene in NMP","Yielding and flow (Soft Matter), D=1+-0.5 um, h~1 nm",1000,1000,0.47,0.195,0.205,"rheological (rigidity)","no polymer, no tunnelling; cleanest contact-defined case"),
("xGn-25 in PLA","Sabzi 2013 (manufacturer D=25 um, t=7 nm)",3571,3571,0.0,0.037,0.045,"rheological and electrical","nominal sizes, no distribution"),
("N02 in PLA","Sabzi 2013 (D<=10 um, t<1 nm)",10000,10000,0.0,0.0031,0.0040,"rheological and electrical","upper-bound A; thickness is 'less than 1 nm'"),
("GNP in UPR (12.1 um / 51.5 nm)","Lu 2006",235,235,0.0,0.64,0.64,"electrical","authors' excluded-volume estimate 0.67 vol%"),
("Foliated graphite in nylon 6","Weng 2005 (same filler as Lu)",235,235,0.0,0.74,0.74,"electrical","tunnelling invoked by authors"),
("Graphite sheets in PMMA (UG)","Zheng 2002 (0.4 mm?/5 um? form factor ~80)",80,100,0.0,2.0,2.0,"electrical","units in the extraction look garbled; form factor ~80-100 as stated"),
("Graphite nanosheets in PMMA","Chen 2003 (A 100-500)",100,500,0.0,0.31,0.31,"electrical","size range only"),
("GNS in TPU","Jan 2016 (A=1272, L=1.4 um)",1272,1272,0.0,0.55,0.55,"electrical","3-4 layers; aggregation noted"),
("GNP-L in PLA","Gao 2017 (8 um/42 nm in composite, 14 um/40 nm pristine)",190,350,0.0,2.0,4.1,"electrical","7 wt% already conductive; threshold between 5 and 7 wt%"),
("GNP-S in PLA","Gao 2017 (0.7 um/13 nm in composite, 1.2 um/10 nm pristine)",54,120,0.0,5.0,7.5,"electrical","no conduction up to 13 wt%"),
("xGnP-1 in PP (powder coated)","Kalaitzidou 2010 (1 um, ~5 nm)",200,200,0.0,0.1,0.1,"electrical","segregated network"),
("xGnP-1 in PP (melt mixed)","Kalaitzidou 2010",200,200,0.0,8.0,9.0,"electrical","random dispersion"),
("xGnP-15 in PP (powder coated)","Kalaitzidou 2010 (15 um, ~5 nm)",3000,3000,0.0,0.3,0.5,"electrical","segregated network"),
("xGnP-15 in PP (melt mixed)","Kalaitzidou 2010",3000,3000,0.0,8.0,10.0,"electrical","random dispersion"),
("GNP in PLA (AFM)","Lashgari 2025 (150 nm / 4.5 nm)",33,33,0.0,2.5,3.4,"rheological and electrical","authors' excluded-volume estimate 1.7 vol%"),
("xGnP-M-15 in epoxy","King 2012 (15 um / 7 nm)",2143,2143,0.0,0.99,0.99,"electrical","10 nm tunnelling distance assumed by authors"),
]
rows=[]
for n,p,Alo,Ahi,sig,tlo,thi,kind,note in cases:
    Am=math.sqrt(Alo*Ahi)
    pm=100*phi(Am,sig=sig)
    plo=100*phi(Ahi,sig=sig); phi_=100*phi(Alo,sig=sig)
    tm=math.sqrt(tlo*thi)
    r=tm/pm
    rlo=tlo/phi_; rhi=thi/plo
    rows.append((n,p,Alo,Ahi,sig,round(plo,4),round(phi_,4),tlo,thi,round(r,2),round(rlo,2),round(rhi,2),kind,note))
with open('experimental_comparison.csv','w',newline='') as f:
    w=csv.writer(f); w.writerow(["system","source","A_low","A_high","sigma_lnR","law_vol%_at_A_high","law_vol%_at_A_low","measured_low_vol%","measured_high_vol%","ratio_mid","ratio_min","ratio_max","kind","note"]); w.writerows(rows)
for r in rows: print(f"{r[0][:34]:34s} A={r[2]}-{r[3]} law={r[5]}-{r[6]} meas={r[7]}-{r[8]} ratio={r[9]} [{r[10]},{r[11]}] {r[12]}")
# implied A (isotropic, flat, monodisperse) for the cases where thickness/lateral size is not independent
for name,t in [("Majidian 2021 rheological 0.64 vol%",0.64),("Nassira 2016 electrical 0.033 vol%",0.033),("Yuan 2015 predicted 0.1 vol%",0.1),("Yuan 2015 measured 4.7 wt% in PDMS",wt2vol(4.7,2.0,0.97)),("Gao 2014 ABS 0.13 vol% (authors infer A=272)",0.13),("Kim 2010 TPU <0.5 vol% (authors infer 220-800)",0.5)]:
    try:
        A=brentq(lambda a:100*phi(a)-t,5,1e6)
    except Exception as e: A=float('nan')
    print(f"implied A for {name}: {A:.0f}")
print("Yuan2015 4.7 wt% =",round(wt2vol(4.7,2.0,0.97),2),"vol%")
