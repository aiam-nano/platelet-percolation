from fig_common import *
fig = plt.figure(figsize=(7.2, 6.2))
gs = fig.add_gridspec(2, 6, hspace=0.5, wspace=2.6)
axs = [fig.add_subplot(gs[0, 0:2]), fig.add_subplot(gs[0, 2:4]), fig.add_subplot(gs[0, 4:6]), fig.add_subplot(gs[1, 0:3]), fig.add_subplot(gs[1, 3:6])]
# (a) aspect ratio -------------------------------------------------------------------------
ax = axs[0]; A = np.logspace(0.5, 3.1, 80)
ax.plot(A, 100 * c0 / A, color=COL['base'], ls=(0, (4, 3)), lw=1.0, label='thin limit $1.44/A$')
ax.plot(A, [law(a) for a in A], color=COL['base'], label='model')
pts = [(a, 100 * y / FS / a) for a, y in FL.base]   # fitlaw stores phi_c*A (fraction x A)
ax.errorbar(*zip(*pts), yerr=[0.015 * y for _, y in pts], fmt='o', ms=4.5, color=COL['sim'], mfc='white', mew=1.0, lw=1.0, label='simulation')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('aspect ratio $A=2R/t$'); ax.set_ylabel('threshold $\\phi_c$ (vol. %)')
ax.legend(frameon=False, loc='upper right', handlelength=1.4); ax.set_ylim(0.05, 2e3); ax.set_yticks([0.1, 1, 10]); ax.minorticks_off()
write_csv('fig3a_aspect', ['A', 'sim_vol_pct_div_1.015', 'model_vol_pct'], [(a, round(y, 4), round(law(a), 4)) for a, y in pts])
# (b) restacking ----------------------------------------------------------------------------
ax = axs[1]; Nn = np.logspace(0, math.log10(25), 60)
ax.plot(Nn, [law(100, N=n) for n in Nn], color=COL['stack'], label='model, $A=100$')
bs = {}
for a, y in FL.base:
    for n in (1, 2, 3, 5, 10, 20):
        if abs(a - 100 / n) < 0.02: bs.setdefault(n, []).append(100 * y / FS / a)
bx = [n for n in sorted(bs) for _ in bs[n]]; by = [v for n in sorted(bs) for v in bs[n]]
ax.plot(bx, by, 'o', ms=4.5, color=COL['sim'], mfc='white', mew=0.9, label='simulation')
ax.set_xscale('log'); ax.set_yscale('log'); ax.set_xlabel('layers per tactoid $N$'); ax.set_ylabel('threshold $\\phi_c$ (vol. %)')
ax.set_xticks([1, 2, 5, 10, 20]); ax.set_xticklabels(['1', '2', '5', '10', '20']); ax.set_yticks([1, 2, 5, 10, 20]); ax.set_yticklabels(['1', '2', '5', '10', '20']); ax.minorticks_off(); ax.legend(frameon=False, loc='upper left', handlelength=1.6)
write_csv('fig3b_layers', ['N', 'sim_vol_pct_div_1.015', 'model_vol_pct_A100'], [(n, round(v, 4), round(law(100, N=n), 4)) for n, v in zip(bx, by)])
# (c) dispersity ----------------------------------------------------------------------------
ax = axs[2]; sg = np.linspace(0, 0.95, 60)
for A_, ls in ((10, ':'), (30, '-'), (100, '--')):
    ax.plot(sg, [law(A_, sig=s) / law(A_) for s in sg], color=COL['disp'], ls=ls, label=f'model, $A={A_}$')
ax.axvspan(0.55, 1.0, color='#e6e6e6', lw=0, zorder=0)
ax.text(0.77, 0.97, 'not claimed', rotation=90, va='top', ha='center', color=COL['base'])
# simulated points: finite-A fits and pre-registered G1, H1 (number-mean radius fixed), relative to the model's monodisperse value at the same A
rows = []
for A_, s, y, e in FL.fin: rows.append((A_, s, y, 0.02, 'calibration'))
gc = pickle.load(open(os.path.join(SP, 'gapcyl_a.pkl'), 'rb')); gb = pickle.load(open(os.path.join(SP, 'gapcyl_b.pkl'), 'rb'))
for k, A_, s in (('G1a', 100, .3), ('G1b', 100, .5), ('G1c', 30, .5)):
    a = gc[k][:, 0] * 100 / FS; rows.append((A_, s, a.mean() / law(A_), a.std(ddof=1) / math.sqrt(len(a)) / law(A_), 'pre-registered'))
a = gb['H1'][:, 0] * 100 / FS; rows.append((60, 0.8, a.mean() / law(60), a.std(ddof=1) / math.sqrt(len(a)) / law(60), 'pre-registered, failed'))
mk = {'calibration': ('o', COL['sim'], 'white'), 'pre-registered': ('s', COL['sim'], COL['sim']), 'pre-registered, failed': ('X', COL['sim'], COL['sim'])}
for lab in mk:
    r = [x for x in rows if x[4] == lab]
    ax.errorbar([x[1] for x in r], [x[2] for x in r], yerr=[x[3] for x in r], fmt=mk[lab][0], ms=4.5, color=mk[lab][1], mfc=mk[lab][2], mew=1.0, lw=1.0, label={'calibration': 'calibration', 'pre-registered': 'pre-registered', 'pre-registered, failed': 'pre-registered, failed'}[lab])
ax.set_xlabel('$\\sigma$, standard deviation of ln $R$'); ax.set_ylabel('Ratio of $\\phi_c$ to that of\nthe monodisperse system'); ax.set_ylim(0.3, 1.75); ax.set_yticks([0.4, 0.6, 0.8, 1.0]); ax.set_xlim(0, 0.95)
ax.legend(frameon=False, loc='upper left', handlelength=1.6, bbox_to_anchor=(0.0, 1.02), labelspacing=0.25)
write_csv('fig3c_dispersity', ['A', 'sigma', 'sim_over_law_monodisperse', 'se', 'set'], [(a, s, round(y, 4), round(e, 4), k) for a, s, y, e, k in rows])
# (d) alignment -----------------------------------------------------------------------------
ax = axs[3]; Sv = np.linspace(0, 0.52, 40); cols = {2: '#E7BBD3', 4: '#CC79A7', 6: '#9C4A78', 8: '#5E2146'}
for d in (2, 4, 6, 8):
    ax.plot(Sv, [LAW.phi(33.333, S=s, d=d) / LAW.phi(33.333) for s in Sv], color=cols[d], label=f'model, $d={d}$')
al = [x for x in FL.al if x[0] < 100]
for d in (2, 4, 6, 8):
    r = [x for x in al if x[1] == d]
    if r: ax.errorbar([x[2] for x in r], [x[3] for x in r], yerr=[x[4] for x in r], fmt='o', ms=4.5, color=cols[d], mfc='white', mew=1.0, lw=1.0)
ax.set_xlabel('local pair order $S$'); ax.set_ylabel('Ratio of $\\phi_c$ to that of\nthe isotropic system'); ax.set_xlim(0, 0.52)
ax.legend(frameon=False, loc='upper left', ncol=2, handlelength=1.6)
write_csv('fig3d_alignment', ['A_eff', 'domain_diameter_d', 'S', 'sim_ratio', 'se'], [(a, d, S, r, e) for a, d, S, r, e in al])
# (e) ripple (thin limit) --------------------------------------------------------------------
ax = axs[4]; av = np.linspace(0, 6.2, 80); lc = {0.5: '#7FD1B9', 1.0: '#009E73', 2.0: '#00553E'}
r = pickle.load(open(os.path.join(SP, 'ripthin2_G.pkl'), 'rb')); fl = r['flat']; m0 = fl.mean(); e0 = fl.std() / math.sqrt(len(fl)); D = []
for k, v in r.items():
    if k == 'flat': continue
    a_, lam = k; R = v.mean() / m0; e = R * math.sqrt((v.std() / math.sqrt(len(v)) / v.mean()) ** 2 + (e0 / m0) ** 2); D.append((a_, lam, R, e))
for lam in (0.5, 1.0, 2.0):
    ax.plot(av, [LAW.Krip(a_, lam) for a_ in av], color=lc[lam], label=f'model, $\\lambda={lam:g}$')
    q = [x for x in D if x[1] == lam]
    ax.errorbar([x[0] for x in q], [x[2] for x in q], yerr=[x[3] for x in q], fmt='o', ms=4.5, color=lc[lam], mfc='white', mew=1.0, lw=1.0)
ax.axvspan(3.75, 6.3, color='#e6e6e6', lw=0, zorder=0); ax.text(5.0, 1.47, 'not claimed\n($a\\geq4$)', ha='center', va='center', color=COL['base'])
ax.set_xlim(0, 6.2); ax.set_xlabel('peak slope $a$'); ax.set_ylabel('Ratio of $\\phi_c$ to that of\na flat sheet'); ax.legend(frameon=False, loc='upper left', handlelength=1.6)

write_csv('fig3e_ripple', ['peak_slope', 'wavelength', 'sim_ratio', 'se'], [(a, l, round(x, 4), round(e, 4)) for a, l, x, e in D])
for lab, ax in zip('abcde', axs): panel_label(ax, lab)
savefig(fig, 'Fig4_dose_response')
print('ok')
