"""Compare the model (finallaw3) with experimental thresholds from the literature (16 cases from 11 papers).
Model = geometric, contact-defined, random orientation, flat discs, nominal volume fraction nV. Thresholds are in
volume per cent; ratio = measured / model, so a ratio above 1 means the measured threshold lies above the geometric prediction.
The values were transcribed from the source papers with AI assistance (Supplementary Table 4 of the paper). Run from this folder (needs numpy and scipy;
finallaw3.py and mylaw.py are in ../)."""
import sys, csv, math, os
sys.path.insert(0,'..'); sys.path.insert(0,'.')
from finallaw3 import phi
# name, paper, A_lo, A_hi, sigma, thr_lo, thr_hi (vol%), kind, note, (spare)
cases = [
("Few-layer graphene in NMP","Barwich 2015, D=1+-0.5 um, h~1 nm",1000,1000,0.47,0.195,0.205,"rheological (rigidity)","no polymer, no tunnelling",""),
("xGn-25 in PLA","Sabzi 2013 (manufacturer D=25 um, t=7 nm)",3571,3571,0.0,3.7,4.5,"rheological and electrical","Table 3 values are volume fractions printed as vol%",""),
("N02 in PLA","Sabzi 2013 (D<=10 um, t<1 nm)",10000,10000,0.0,0.31,0.40,"rheological and electrical","as above; A is for t = 1 nm",""),
("GNP in UPR (12.1 um / 51.5 nm)","Lu 2006 (paper also prints D = 12.5 um)",235,235,0.0,0.64,0.64,"electrical","fitted 0.64; 0.67 is an empirical estimate",""),
("Foliated graphite in nylon 6","Weng 2005 (same filler as Lu)",235,235,0.0,0.74,0.74,"electrical","",""),
("Graphite sheets in PMMA (UG)","Zheng 2002 (0.4 mm / 5 um, form factor ~80)",80,100,0.0,2.0,2.0,"electrical","A = 80 measured, 100 from the authors' Eq. 6",""),
("Graphite nanosheets in PMMA","Chen 2003 (A 100-500)",100,500,0.0,0.31,0.31,"electrical","size range only",""),
("GNS in TPU","Jan 2016 (L=1.4 um, 3-4 layers)",1272,1272,0.0,0.55,0.55,"electrical","t = 3.3 x 0.335 nm assumed; paper prints 0.035 nm per layer",""),
("GNP-L in PLA","Gao 2017",190,350,0.0,2.9,4.1,"electrical","5-7 wt%, rho_PLA 1.24, rho_GNP 2.2",""),
("GNP-S in PLA","Gao 2017",54,120,0.0,5.9,7.8,"electrical","10-13 wt% (unannealed)",""),
("xGnP-1 in PP (powder coated)","Kalaitzidou 2010 (<1 um, ~5 nm)",200,200,0.0,0.1,0.1,"electrical","segregated network",""),
("xGnP-1 in PP (melt mixed)","Kalaitzidou 2010",200,200,0.0,8.0,9.0,"electrical","",""),
("xGnP-15 in PP (powder coated)","Kalaitzidou 2010 (15 um, ~5 nm)",3000,3000,0.0,0.3,0.3,"electrical","segregated network; 0.3 to 0.5 applies to slow cooling",""),
("xGnP-15 in PP (melt mixed)","Kalaitzidou 2010",3000,3000,0.0,8.0,9.0,"electrical","8-10 vol% in the paper is PAN carbon fibre",""),
("GNP in PLA (AFM)","Lashgari 2025 (150 nm / 4.5 nm)",33,33,0.0,2.5,3.4,"rheological and electrical","",""),
("xGnP-M-15 in epoxy","King 2012 (15 um / 7 nm)",2143,2143,0.0,0.99,0.99,"electrical","tunnelling assumed by authors",""),
]
rows=[]
for n,p,Alo,Ahi,sig,tlo,thi,kind,note,_ in cases:
    Am=math.sqrt(Alo*Ahi); pm=100*phi(Am,sig=sig); plo=100*phi(Ahi,sig=sig); phi_=100*phi(Alo,sig=sig)
    tm=math.sqrt(tlo*thi)
    rows.append((n,p,Alo,Ahi,sig,round(plo,4),round(phi_,4),tlo,thi,round(tm/pm,2),round(tlo/phi_,2),round(thi/plo,2),kind,note))
out='experimental_comparison.csv'
with open(out,'w',newline='') as f:
    w=csv.writer(f); w.writerow(["system","source","A_low","A_high","sigma_lnR","law_vol%_at_A_high","law_vol%_at_A_low","measured_low_vol%","measured_high_vol%","ratio_mid","ratio_min","ratio_max","kind","note"]); w.writerows(rows)
for r in rows: print(f"{r[0][:34]:34s} ratio={r[9]:7.2f} [{r[10]},{r[11]}]")
