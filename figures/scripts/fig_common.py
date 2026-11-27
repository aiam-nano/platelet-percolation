import sys, os, math, io, contextlib, csv, pickle, json
HERE = os.path.dirname(os.path.abspath(__file__)); FIG = os.path.abspath(os.path.join(HERE, '..'))
CODE = os.path.abspath(os.path.join(HERE, '..', '..', 'code')); LIT = os.path.join(CODE, 'literature_checks')
# DATA: the result files (pickles) behind the figures; they are in the Zenodo data record. Set PERC_DATA to their folder.
SP = DATA = os.environ.get('PERC_DATA', os.path.join(FIG, 'data'))
for _p in (DATA, CODE, os.path.join(CODE, 'gap_tests'), LIT): sys.path.insert(0, _p)
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
with contextlib.redirect_stdout(io.StringIO()):
    import fitlaw as FL
    import finallaw3 as LAW
from mylaw import c0, c1, c2, P as Pdisp
FS = 1.015
# colour language (Okabe-Ito, colour-blind safe): one colour per departure, used in every figure
COL = dict(base='#4d4d4d', stack='#0072B2', disp='#E69F00', align='#CC79A7', ripple='#009E73', sim='#000000', grey='#8c8c8c', bad='#D55E00', sky='#56B4E9')
from matplotlib import font_manager as _fm
FONT = 'Arial' if any('Arial' == f.name for f in _fm.fontManager.ttflist) else 'Liberation Sans'
plt.rcParams.update({'font.family': FONT, 'font.size': 8, 'axes.titlesize': 8, 'axes.labelsize': 8, 'xtick.labelsize': 8, 'ytick.labelsize': 8,
    'legend.fontsize': 8, 'axes.spines.top': False, 'axes.spines.right': False, 'axes.linewidth': 1.0, 'xtick.major.width': 1.0, 'ytick.major.width': 1.0,
    'xtick.major.size': 3.5, 'ytick.major.size': 3.5, 'lines.linewidth': 1.25, 'lines.markeredgewidth': 1.0, 'patch.linewidth': 1.0, 'figure.dpi': 150, 'savefig.dpi': 600,
    'axes.titleweight': 'normal', 'axes.titlelocation': 'left', 'figure.facecolor': 'white', 'axes.facecolor': 'white', 'savefig.facecolor': 'white',
    'mathtext.fontset': 'custom', 'mathtext.rm': FONT, 'mathtext.it': FONT + ':italic', 'mathtext.bf': FONT + ':bold', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'svg.fonttype': 'none'})
def panel_label(ax, s, dx=-0.16, dy=1.04):
    ax.text(dx, dy, s, transform=ax.transAxes, fontsize=8, fontweight='bold', va='bottom', ha='left')
def law(A, **kw): return 100 * LAW.phi(A, **kw)   # volume per cent
def savefig(fig, name, pdf=True):
    from PIL import Image
    p = os.path.join(FIG, name + '.png')
    fig.savefig(p, bbox_inches='tight', pad_inches=0.04)
    Image.open(p).convert('RGB').save(p, dpi=(600, 600))
    if pdf: fig.savefig(os.path.join(FIG, name + '.pdf'), bbox_inches='tight', pad_inches=0.04)
    plt.close(fig)
SRC = os.path.join(FIG, 'source_data'); os.makedirs(SRC, exist_ok=True)
def write_csv(name, header, rows):
    with open(os.path.join(SRC, name + '.csv'), 'w', newline='') as f:
        w = csv.writer(f); w.writerow(header); w.writerows(rows)
