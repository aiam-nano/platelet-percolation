from fig_common import *
D = [('G1a','Averaging conversions',-0.7,6,'pass'),('G1b','Averaging conversions',0.6,6,'pass'),('G1c','Averaging conversions',0.3,6,'pass'),
('B1 bimodal','Non-log-normal distributions',7.6,8,'pass'),('B2 bimodal','Non-log-normal distributions',18.7,8,'fail'),('H1 log-normal 0.8','Non-log-normal distributions',23.6,8,'fail'),('Q1 Pareto','Non-log-normal distributions',-6.7,8,'pass'),
('C1','Thickness tied to radius',-1.5,6,'pass'),('C2','Thickness tied to radius',0.5,6,'pass'),('C3','Thickness tied to radius',-0.2,6,'pass'),('C4','Thickness tied to radius',0.2,6,'pass'),
('P1','Ripple with dispersity',-3.5,8,'pass'),('P2','Ripple with dispersity',-4.4,8,'pass'),('P3','Ripple with dispersity',-0.7,8,'pass'),('P4','Ripple with dispersity',4.8,8,'pass'),
('TB1','Thick biaxial slabs',0.2,8,'pass'),('TB2','Thick biaxial slabs',-1.2,8,'pass'),('TB3','Thick biaxial slabs',1.9,8,'pass'),('TB4','Thick biaxial slabs',-1.5,8,'pass'),('TB5','Thick biaxial slabs',-3.2,15,'pass'),
('R1','Random wrinkle fields',1.5,8,'pass'),('R2','Random wrinkle fields',0.2,8,'pass'),('R3','Random wrinkle fields',-4.5,8,'pass'),('R4','Random wrinkle fields',4.4,8,'pass'),
('F1','Sealed exam (model with refitted prefactor)',0.3,0,'exam'),('F2','Sealed exam (model with refitted prefactor)',0.4,0,'exam'),('F3','Sealed exam (model with refitted prefactor)',-0.4,0,'exam'),('F4','Sealed exam (model with refitted prefactor)',0.0,0,'exam'),
('F10 ripple','Sealed exam (model with refitted prefactor)',1.9,0,'exam'),('F11 ripple','Sealed exam (model with refitted prefactor)',0.8,0,'exam'),
('F8 aligned','Sealed exam (model with refitted prefactor)',1.1,0,'post'),('F9 aligned','Sealed exam (model with refitted prefactor)',1.0,0,'post'),('F12 aligned','Sealed exam (model with refitted prefactor)',3.8,0,'post')]
STEP = 0.15; GAP = 0.12
y = 0.0; ys = []; heads = []; prev = None
for d in D:
    if d[1] != prev:
        if prev is not None: y += GAP
        heads.append((y, d[1])); y += STEP; prev = d[1]
    ys.append(y); y += STEP
Ht = y + 0.62
fig = plt.figure(figsize=(7.2, Ht)); ax = fig.add_axes([0.30, 0.55 / Ht, 0.66, (y - 0.0) / Ht]); ax.set_ylim(y - 0.0 + 0.0, -0.1)
ax.set_xlim(-26, 26)
for t in (-20, -10, 10, 20): ax.axvline(t, color='#d9d9d9', lw=1.0, zorder=0)
ax.axvline(0, color=COL['base'], lw=1.0, zorder=1)
for d, yy in zip(D, ys):
    if d[3] > 0: ax.add_patch(plt.Rectangle((-d[3], yy - STEP * 0.38), 2 * d[3], STEP * 0.76, fc='#e6e6e6', ec='none', zorder=0.5))
    kind = d[4]
    if kind == 'pass': ax.plot(d[2], yy, 'o', ms=4.5, color='#6e6e6e', mec='#6e6e6e', zorder=3)
    elif kind == 'fail': ax.plot(d[2], yy, 'D', ms=4.5, color=COL['bad'], mec=COL['bad'], zorder=3)
    elif kind == 'exam': ax.plot(d[2], yy, 'o', ms=4.5, color=COL['stack'], mec=COL['stack'], zorder=3)
    else: ax.plot(d[2], yy, 'o', ms=4.5, mfc='white', mec=COL['stack'], mew=1.2, zorder=3)
    ax.text(-27.2, yy, d[0], ha='right', va='center', transform=ax.transData, clip_on=False)
for hy, h in heads: ax.text(-33.0, hy, h, ha='left', va='center', fontweight='bold', clip_on=False)
ax.spines['left'].set_visible(False); ax.set_yticks([]); ax.set_xticks([-20, -10, 0, 10, 20]); ax.set_xticklabels(['−20', '−10', '0', '+10', '+20'])
ax.set_xlabel('Fractional difference of model from simulation (%)')
write_csv('fig5_prediction_errors', ['test', 'family', 'error_percent', 'window_half_width_percent', 'status'], [(d[0], d[1], d[2], d[3] if d[3] else '', d[4]) for d in D])
savefig(fig, 'Fig5_errors'); print(Ht)
