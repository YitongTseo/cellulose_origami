import os, csv
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

HERE = os.path.dirname(os.path.abspath(__file__))
D = os.path.join(HERE, '..', 'raw')
OUT = HERE
FIG = os.path.join(HERE, '..', 'figures')
rcParams.update({'font.family': 'Arial', 'font.size': 9.75, 'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold', 'mathtext.default': 'regular', 'axes.linewidth': 0.6, 'pdf.fonttype': 42})

# Sample assignment per Yitong, 2026-10-08 (same names / colours as the XPS figure). Figure order: bottom -> top.
FIGS = [('WT_untreated', 'WT untreated', '#F5B041'),
        ('WT_NaOH_treated', 'WT NaOH treated', '#dd9789'),
        ('BslA_NaOH_treated', 'BslA NaOH treated', '#5DADE2')]   # yellow / salmon / blue, as in Figures 3 and 6
ORDER = FIGS

def load(n):
    a = np.loadtxt(os.path.join(D, n + '.csv'), delimiter=','); o = np.argsort(a[:, 0]); return a[o, 0], a[o, 1]
def peak(x, y, lo, hi):
    m = (x >= lo) & (x <= hi); i = np.argmax(y[m]); return y[m][i], x[m][i]
def height(x, y, lo, hi, b1, b2):
    def anchor(w): m = (x >= w - 10) & (x <= w + 10); return x[m][np.argmin(y[m])], y[m].min()
    xa, ya = anchor(b1); xb, yb = anchor(b2); h, pos = peak(x, y, lo, hi)
    return h - (ya + (yb - ya) * (pos - xa) / (xb - xa)), pos
def val(x, y, w): return y[np.argmin(abs(x - w))]

spec = {n: load(n) for n, _, _ in ORDER}

# ---------------- figure: (A) full range, (B) fingerprint region, (C) amide region overlaid
BANDS = [(3340, 'O–H'), (2900, 'C–H'), (1740, 'C=O (ester)'), (1650, 'amide I'), (1540, 'amide II'),
         (1430, 'CH$_2$'), (1160, 'C–O–C'), (1042, '1055 / 1030\nC–O'), (898, 'β-glycosidic')]
fig = plt.figure(figsize=(7.2, 6.6))
gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.02], width_ratios=[1.75, 1], hspace=0.62, wspace=0.22,
                      left=0.065, right=0.985, top=0.9, bottom=0.07)
a1 = fig.add_subplot(gs[0, :]); a2 = fig.add_subplot(gs[1, 0]); a3 = fig.add_subplot(gs[1, 1])
OFF = 0.85
norm = {}
for i, (n, name, col) in enumerate(FIGS):
    x, y = spec[n]; yn = y / peak(x, y, 1000, 1070)[0]; norm[n] = (x, yn); k = x >= 420
    a1.plot(x[k], yn[k] + OFF * i, color=col, lw=0.9)
    a1.text(2500, OFF * i + yn[(x > 2250) & (x < 2750)].max() + 0.07, name, color=col, fontsize=10.4, va='bottom', ha='center')
    m = (x >= 800) & (x <= 1800); a2.plot(x[m], yn[m] + OFF * i, color=col, lw=0.9)
top = OFF * (len(FIGS) - 1) + 1.2
for ax, lo, hi in ((a1, 4000, 400), (a2, 1800, 800)):
    ax.set_xlim(lo, hi); ax.set_ylim(-0.05, top); ax.set_yticks([])
    ax.set_xlabel('Wavenumber (cm$^{-1}$)')
    for w, lab in BANDS:
        if not (hi <= w <= lo): continue
        if ax is a1 and w < 1800: continue          # in the full-range panel, label only the high-wavenumber bands
        for ww in ((1055, 1030) if w == 1042 else (w,)): ax.axvline(ww, color='#b5b5b5', lw=0.5, ls=(0, (3, 3)), zorder=0)
        ax.text(w, top + 0.04, (lab if w == 1042 else f'{w}\n{lab}'), rotation=90, fontsize=7.15, color='#444', ha='center', va='bottom')
a1.axvspan(1800, 800, color='#000', alpha=0.04, lw=0)
a1.text(1790, top - 0.06, 'region shown in (b)', fontsize=7.15, color='#666', ha='left', va='top')
for ax in (a1, a2): ax.set_ylabel('Absorbance (norm., offset)')

# C: amide I / II, same normalisation, no offset, linear baseline through 1800 and 1480 cm-1 removed.
# Each amide I maximum is marked and labelled with its position (the BslA band sits ~17 cm-1 lower).
pkI = {}
for n, name, col in FIGS:
    x, yn = norm[n]; m = (x >= 1480) & (x <= 1800)
    y0, y1 = val(x, yn, 1480), val(x, yn, 1800)
    yb = yn[m] - (y0 + (y1 - y0) * (x[m] - 1480) / 320)
    a3.plot(x[m], yb, color=col, lw=1.1, label=name)
    k = (x[m] >= 1625) & (x[m] <= 1670); i = np.argmax(yb[k]); px, py = x[m][k][i], yb[k][i]; pkI[n] = px
    a3.plot([px, px], [0, py], color=col, lw=0.7, ls=(0, (2, 2)))
    a3.plot(px, py, 'o', color=col, ms=3.2, mec='white', mew=0.5, zorder=5)
    a3.text(px, py + 0.022, f'{px:.0f}', color=col, fontsize=8.45, fontweight='bold', ha='center', va='bottom',
            bbox=dict(fc='white', ec='none', pad=0.5, alpha=0.85), zorder=6)
xw = np.mean([pkI[n] for n, _, _ in FIGS[:2]]); xb = pkI[FIGS[2][0]]; ya = 0.79
a3.annotate('', xy=(xb, ya), xytext=(xw, ya), arrowprops=dict(arrowstyle='->', color='#222', lw=0.8, shrinkA=0, shrinkB=0))
a3.text(1490, ya, f'amide I maximum\n{xw:.0f} → {xb:.0f} cm$^{{-1}}$', fontsize=7.8, color='#222', ha='right', va='center', linespacing=1.25)
for w, lab in ((1650, 'amide I'), (1540, 'amide II')):
    a3.text(w - (10 if w == 1650 else 0), 1.02, lab, transform=a3.get_xaxis_transform(), fontsize=8.45, color='#444', ha='center', va='bottom')
a3.axhline(0, color='#999', lw=0.4)
a3.set_xlim(1800, 1480); a3.set_ylim(-0.03, 1.15)
a3.set_xlabel('Wavenumber (cm$^{-1}$)'); a3.set_ylabel('Absorbance / A$_{1030}$')
a3.tick_params(labelsize=8.45)
a3.legend(loc='upper left', frameon=False, fontsize=7.8, handlelength=1.4, labelspacing=0.3)
for ax, t in ((a1, '(a)'), (a2, '(b)'), (a3, '(c)')):
    ax.text(-0.16 if ax is a3 else -0.05, 1.03, t, transform=ax.transAxes, fontsize=13, va='bottom', ha='left')
fig.savefig(os.path.join(FIG, 'FTIR_figure.png'), dpi=300); fig.savefig(os.path.join(FIG, 'FTIR_figure.pdf'))

# ---------------- indices table
rows = []
for n, name, _ in ORDER:
    x, y = spec[n]; yn = y / peak(x, y, 1000, 1070)[0]
    h1372, _ = height(x, yn, 1360, 1385, 1345, 1400); h2900, _ = height(x, yn, 2880, 2920, 2800, 3000)
    h1430, p1430 = height(x, yn, 1415, 1440, 1400, 1450); h898, p898 = height(x, yn, 888, 905, 870, 925)
    hAmI, pAmI = height(x, yn, 1625, 1670, 1590, 1720); hAmII, _ = height(x, yn, 1520, 1560, 1490, 1590)
    _, pOH = peak(x, yn, 3200, 3500)
    smeared = n == 'BslA_NaOH_treated'
    rows.append(dict(sample=name, file=n, OH_pos=pOH, CH2_pos=p1430,
                     TCI=h1372 / h2900, LOI=(h1430 / h898 if not smeared and h898 > 0.01 else np.nan),
                     HBI=val(x, yn, 3400) / val(x, yn, 1320), amideI_1030=hAmI, amideII_1030=hAmII, amideI_pos=pkI[n]))
with open(os.path.join(OUT, 'FTIR_table.csv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader()
    for r in rows: w.writerow({k: (f'{v:.2f}' if isinstance(v, float) else v) for k, v in r.items()})
for r in rows: print({k: (round(v, 2) if isinstance(v, float) else v) for k, v in r.items()})
print('saved')
