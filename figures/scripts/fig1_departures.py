from scenes import *
from matplotlib.patches import Ellipse, Arc, Circle, Polygon
BL = COL['stack']; TINT = '#d3e6f3'; GR = '#8c8c8c'; DK = '#262626'
LWS = 1.25
def sk_a(ax):
    cx, cy, w, h, t = 4.4, 4.3, 6.4, 2.2, 0.9
    ax.add_patch(Arc((cx, cy - t), w, h, theta1=180, theta2=360, ec=GR, lw=LWS))
    ax.plot([cx - w / 2] * 2, [cy, cy - t], color=GR, lw=LWS); ax.plot([cx + w / 2] * 2, [cy, cy - t], color=GR, lw=LWS)
    ax.add_patch(Ellipse((cx, cy), w, h, fc=TINT, ec=BL, lw=1.6))
    ax.plot([cx, cx + w / 2], [cy, cy], color=GR, lw=1.0, ls=(0, (3, 2))); ax.text(cx + w / 4, cy + 0.25, '$R$', ha='center', va='bottom')
    xb = cx + w / 2 + 0.9; ax.plot([xb, xb], [cy - t, cy], color=DK, lw=1.0)
    for y in (cy - t, cy): ax.plot([xb - 0.3, xb + 0.3], [y, y], color=DK, lw=1.0)
    ax.text(xb + 0.55, cy - t / 2, '$t$', va='center')
def sk_b(ax):
    cx, cy, w, h, t, N = 4.0, 5.2, 5.6, 1.9, 0.95, 3
    for k in range(N, 0, -1):
        ax.add_patch(Arc((cx, cy - k * t), w, h, theta1=180, theta2=360, ec=GR, lw=LWS))
    for s in (-1, 1): ax.plot([cx + s * w / 2] * 2, [cy, cy - N * t], color=GR, lw=LWS)
    ax.add_patch(Ellipse((cx, cy), w, h, fc=TINT, ec=BL, lw=1.6))
    xb = cx + w / 2 + 0.6; yt, yb = cy - h / 2 * 0 , cy - N * t
    ax.plot([xb, xb], [yb, yt], color=DK, lw=1.0)
    for y in (yb, yt): ax.plot([xb, xb - 0.3], [y, y], color=DK, lw=1.0)
    ax.text(xb + 0.4, (yb + yt) / 2, '$N$ layers', va='center')
def sk_c(ax):
    E = [(3.2, 5.4, 4.6, 2.1, 0), (9.2, 6.5, 3.6, 1.6, 0), (9.6, 3.9, 3.4, 1.5, 0), (6.0, 2.9, 3.8, 1.7, 1), (2.0, 2.2, 1.6, 0.9, 0)]
    for cx, cy, w, h, hl in E:
        ax.add_patch(Ellipse((cx, cy), w, h, fc=TINT if hl else 'white', ec=BL if hl else GR, lw=1.6 if hl else LWS))
def sk_d(ax):
    cx, cy, r = 6.0, 4.4, 2.9
    ax.add_patch(Circle((cx, cy), r, fc=TINT, ec=BL, lw=1.6, ls=(0, (3, 2))))
    def pin(x, y, ang, L=1.15, c=DK):
        dx, dy = L * math.cos(ang), L * math.sin(ang); ax.plot([x, x + dx], [y, y + dy], color=c, lw=1.25, solid_capstyle='round'); ax.plot([x + dx], [y + dy], 'o', ms=3.2, color=c, mew=0)
    rg = np.random.default_rng(3)
    for (x, y) in [(4.9, 3.5), (6.0, 3.1), (7.2, 3.6), (4.6, 4.7), (6.7, 4.5), (5.6, 5.4), (7.0, 5.6), (5.3, 4.0)]:
        pin(x, y, math.radians(90 + rg.normal(0, 7)))
    for (x, y, a) in [(1.6, 5.8, 15), (1.2, 2.6, 160), (10.6, 5.8, 35), (9.6, 1.9, 60)]:
        pin(x, y, math.radians(a), c=GR)
    yb = 0.85; ax.plot([cx - r, cx + r], [yb, yb], color=DK, lw=1.0)
    for s in (-1, 1): ax.plot([cx + s * r] * 2, [yb - .25, yb + .25], color=DK, lw=1.0)
def sk_e(ax):
    x = np.linspace(1, 11, 300); y = 4.2 + 1.3 * np.sin(2 * np.pi * (x - 1) / 5)
    ax.plot(x, y, color=BL, lw=1.6)
    k = 1.3 * 2 * np.pi / 5; xs = np.array([5.0, 7.2]); ax.plot(xs, 4.2 + k * (xs - 6), color=GR, lw=1.0, ls=(0, (3, 2)))
    ax.text(7.4, 6.2, 'slope $a$', va='bottom')
    yb = 1.5; ax.plot([3.75, 8.75], [yb, yb], color=DK, lw=1.0)
    for xx in (3.75, 8.75): ax.plot([xx, xx], [yb - .25, yb + .25], color=DK, lw=1.0)
    ax.text(6.25, yb + 0.35, '$\\lambda$', ha='center', va='bottom')
SK = [sk_a, sk_b, sk_c, sk_d, sk_e]
TT = ['Finite thickness', 'Restacked tactoid', 'Log-normal radii', 'Aligned domains', 'Sinusoidal ripple']
CAP = ['aspect ratio\n$A=2R/t$', 'effective aspect ratio\n$A_e=A/N$', 'log-normal $\\sigma$,\nstandard deviation\nof ln $R$', 'domain diameter $d$,\nlocal order $S$', 'peak slope $a$,\nwavelength $\\lambda$']
keys = ['cyl10', 'stack3', 'ln5', 'al_b', 'rip3']
PAR = ['$A=10$', '$A=30$, $N=3$', '$\\sigma=0.5$', '$d=4$, $S\\approx0.4$', '$a=3$, $\\lambda=1$']
W, H = 7.2, 4.2
fig = plt.figure(figsize=(W, H)); cw = W / 5
def inch(x, y): return (x / W, 1 - y / H)
info = []
for j in range(5):
    x0 = j * cw
    ax = fig.add_axes([x0 / W + 0.004, 1 - 1.3 / H, cw / W - 0.008, 1.0 / H]); ax.set_xlim(0, 12); ax.set_ylim(0, 8); ax.set_aspect('equal'); ax.axis('off'); SK[j](ax)
    fig.text(*inch(x0 + 0.04, 0.05), 'abcde'[j], fontweight='bold', va='top', ha='left')
    fig.text(*inch(x0 + 0.2, 0.05), TT[j], va='top', ha='left')
    fig.text(*inch(x0 + cw / 2, 1.34), CAP[j], va='top', ha='center', linespacing=1.25)
    fig.text(*inch(x0 + cw / 2, 1.84), PAR[j], va='top', ha='center')
    ax = fig.add_axes([x0 / W + 0.002, 0.0, cw / W - 0.004, 1 - 2.07 / H]); inf = draw_scene(ax, keys[j], title=False)
    info.append((keys[j], inf['count'], inf['span'], inf['S']))
write_csv('fig1_configurations', ['scene', 'flakes_at_first_wrap', 'flakes_in_spanning_cluster', 'S_local_pair'], info)
savefig(fig, 'Fig1_departures'); print(info)
