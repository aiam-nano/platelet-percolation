"""T9: the mesh-versus-implicit-surface sheet contact (used for the ripple cases F10/F11) must reproduce the exact disc contact when the sheet is flat.
Run from the workdir: python3 verification/sheet_flat_check.py"""
import sys, numpy as np
sys.path.insert(0, ".")
import gap_rip
gap_rip.RCAP = 6.0
for sig in (0.0, 0.4):
    o = gap_rip.paired(24.0, sig, 1.2, 1.0, 30, 5151, validate=True); d = o[:, 1]-o[:, 0]
    print(f"sigma={sig}: n={len(o)} disc {np.nanmean(o[:,0]):.4f} flat sheet {np.nanmean(o[:,1]):.4f} ratio {np.nanmean(o[:,1])/np.nanmean(o[:,0]):.4f}; identical N_c in {(np.abs(d)<1e-12).sum()} of {len(o)}")
