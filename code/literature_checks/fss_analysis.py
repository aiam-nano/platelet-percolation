"""Summarise shape_fss.jsonl (periodic, first wrapping cycle), span_open.jsonl (open faces, first spanning cluster, Yi-Tawerghi style)
and ell_fss.jsonl (spheroids).  Writes fss_analysis.txt."""
import json, math, sys
import numpy as np
NU = 0.8765
YT = {'disc': 0.9614, 'square': 0.8647, 'triangle': 0.7295}
def load(f):
    try: return [json.loads(l) for l in open(f)]
    except FileNotFoundError: return []
def wfit(x, y, e, deg=1):
    A = np.vander(x, deg + 1, increasing=True); W = np.diag(1 / np.array(e) ** 2)
    cov = np.linalg.inv(A.T @ W @ A); b = cov @ A.T @ W @ y
    chi2 = float(((A @ b - y) / np.array(e)) @ ((A @ b - y) / np.array(e)))
    return b, np.sqrt(np.diag(cov)), chi2
out = []
P = lambda *a: (print(*a), out.append(' '.join(str(x) for x in a)))
rows = load('shape_fss.jsonl')
P("== periodic box, first wrapping cycle (any direction), zero-thickness plates ==")
res = {}
for sh in ('disc', 'square', 'triangle'):
    r = sorted([x for x in rows if x['shape'] == sh], key=lambda x: x['L'])
    if not r: continue
    P(sh, ' '.join(f"L{x['L']}:{x['eta']:.4f}+-{x['err']:.4f}(n={x['n']})" for x in r))
    big = [x for x in r if x['L'] >= 32]
    w = np.array([1 / x['err'] ** 2 for x in big]); m = float((w * np.array([x['eta'] for x in big])).sum() / w.sum()); se = float(1 / math.sqrt(w.sum()))
    L = np.array([x['L'] for x in r], float); y = np.array([x['eta'] for x in r]); e = np.array([x['err'] for x in r])
    b, s, chi2 = wfit(L ** (-1 / NU), y, e)
    P(f"  plateau (weighted mean, L>=32): {m:.4f} +- {se:.4f}; fit eta_inf + a L^(-1/nu) over all L: {b[0]:.4f} +- {s[0]:.4f}, a = {b[1]:.3f}+-{s[1]:.3f}, chi2 = {chi2:.1f} for {len(r) - 2} dof")
    P(f"  Yi and Tawerghi {YT[sh]:.4f}: plateau/YT = {m / YT[sh]:.4f} +- {se / YT[sh]:.4f}")
    res[sh] = (m, se)
if 'disc' in res:
    for sh in ('square', 'triangle'):
        if sh in res:
            ratio = res[sh][0] / res['disc'][0]; er = ratio * math.hypot(res[sh][1] / res[sh][0], res['disc'][1] / res['disc'][0])
            P(f"  ratio {sh}/disc = {ratio:.4f} +- {er:.4f}   Yi and Tawerghi: {YT[sh] / YT['disc']:.4f}")
rows = load('span_open.jsonl')
P("\n== open faces, first spanning cluster (Yi and Tawerghi style), single-axis mean and any-axis ==")
for sh in ('disc', 'square', 'triangle'):
    r = sorted([x for x in rows if x['shape'] == sh], key=lambda x: x['L'])
    if not r: continue
    P(sh, ' '.join(f"L{x['L']}:{x['eta_axis_mean']:.4f}+-{x['eta_axis_err']:.4f}" for x in r))
    if len(r) < 3: continue
    L = np.array([x['L'] for x in r], float); y = np.array([x['eta_axis_mean'] for x in r]); e = np.array([x['eta_axis_err'] for x in r])
    b1, s1, c1 = wfit(1 / L, y, e); b2, s2, c2 = wfit(L ** (-1 / NU), y, e)
    sel = L >= 24
    P(f"  linear in r = 1/L (their extrapolation): eta_c = {b1[0]:.4f} +- {s1[0]:.4f} (chi2 {c1:.1f}/{len(r) - 2});  L^(-1/nu) fit: {b2[0]:.4f} +- {s2[0]:.4f} (chi2 {c2:.1f}/{len(r) - 2});  Yi and Tawerghi {YT[sh]:.4f}")
    P(f"  any-axis first-span mean: " + ' '.join(f"L{x['L']}:{x['eta_any']:.4f}" for x in r))
rows = load('ell_fss.jsonl')
P("\n== oblate spheroids, periodic box, first wrapping cycle ==")
for u in sorted(set(x['u'] for x in rows)):
    r = sorted([x for x in rows if x['u'] == u], key=lambda x: x['L'])
    P(f"u = {u:g}", ' '.join(f"L{x['L']}:{x['eta']:.4f}+-{x['err']:.4f}(n={x['n']})" for x in r))
    big = [x for x in r if x['L'] >= 32]
    if big:
        w = np.array([1 / x['err'] ** 2 for x in big]); m = float((w * np.array([x['eta'] for x in big])).sum() / w.sum()); se = float(1 / math.sqrt(w.sum()))
        L = np.array([x['L'] for x in r], float); y = np.array([x['eta'] for x in r]); e = np.array([x['err'] for x in r])
        b, s, chi2 = wfit(L ** (-1 / NU), y, e)
        P(f"  plateau (L>=32) {m:.4f} +- {se:.4f};  fit eta_inf = {b[0]:.4f} +- {s[0]:.4f} (chi2 {chi2:.1f}/{len(r) - 2});  vs Yi-Tawerghi 0.9614: {m / 0.9614:.4f};  vs Garboczi 1.27: {m / 1.27:.4f}")
open('fss_analysis.txt', 'w').write('\n'.join(out) + '\n')
