"""Publication figure: (A) surveys, (B) hi-res C 1s / O 1s / N 1s fits. Plus a CSV summary table."""
import os, io, contextlib, csv
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import rcParams

with contextlib.redirect_stdout(io.StringIO()):
    import fit_hires as F          # runs the fits on import
HERE = os.path.dirname(os.path.abspath(__file__)); OUT = HERE; FIG = os.path.join(HERE, '..', 'figures')

rcParams.update({'font.family': 'Arial', 'font.size': 9.75, 'mathtext.fontset': 'custom', 'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic', 'mathtext.bf': 'Arial:bold', 'mathtext.default': 'regular', 'axes.linewidth': 0.6,
                 'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'pdf.fonttype': 42})

# (fit key, survey file, hi-res file, display name). Sample assignment per Yitong, 2026-10-08.
SAMPLES = [('WT untreated', 'raw/WT_untreated_survey.csv', 'raw/WT_untreated_hires.csv', 'WT untreated'),
           ('WT NaOH treated', 'raw/WT_NaOH_treated_survey.spe', 'raw/WT_NaOH_treated_hires.spe', 'WT NaOH treated'),
           ('BslA NaOH treated', 'raw/BslA_NaOH_treated_survey.csv', 'raw/BslA_NaOH_treated_hires.csv', 'BslA NaOH treated')]
EXTRA = []
COL = {'WT untreated': '#F5B041', 'WT NaOH treated': '#dd9789', 'BslA NaOH treated': '#5DADE2'}   # yellow / salmon / blue, as in Figures 3 and 6
CFILL = ['#e8504f', '#4fd34f', '#4f4fe8']    # component colours, fixed order (KherveFitting-style red / green / blue)
FITC = '#3b3bd6'; RESC = '#7fdc7f'

# facility (PHI MultiPak) atomic concentrations from the 260904 session PPTX: C, N, O
PHI = {'WT untreated': (72.67, 3.80, 23.53), 'BslA NaOH treated': (58.55, 7.17, 34.27)}
RSF = {'C1s': 4.760, 'N1s': 7.918, 'O1s': 12.211}      # MultiPak corrected RSFs, same PPTX

def dwell(rel):
    """Seconds per point for each region, from the .spe header (SpectralRegDef); .spe/csv values are counts/s."""
    spe = rel if rel.endswith('.spe') else rel.replace('_hires.csv', '_hires.spe')
    head = open(os.path.join(F.XPS, spe), 'rb').read().split(b'EOFH')[0].decode('latin-1')
    return {l.split()[3]: float(l.split()[11]) for l in head.splitlines() if l.startswith('SpectralRegDef:')}

def quant(rel):
    """C, N, O at.% from Shirley areas of the hi-res regions / MultiPak corrected RSFs. Used where no
    facility table exists (261006); reproduces the 260904 facility values to within 1 at.%."""
    b = F.load(rel); pk = {k: v[np.argmax(v[:, 1]), 0] for k, v in b.items()}; A = {}
    for k, (lo, hi) in {'C1s': (pk['C1s'] - 5, pk['C1s'] + 7.5), 'O1s': (pk['O1s'] - 5, pk['O1s'] + 5),
                        'N1s': (pk['C1s'] + 110.2, pk['C1s'] + 120.2)}.items():
        x, y = F.crop(b[k], lo, hi)
        A[k] = (y - F.shirley(x, y, n_end=3)).sum() * abs(np.median(np.diff(x))) / RSF[k]
    t = sum(A.values())
    return tuple(100 * A[k] / t for k in ('C1s', 'N1s', 'O1s'))

AT, AT_SRC = {}, {}
for key, sur, hires, name in SAMPLES + EXTRA:
    AT[key], AT_SRC[key] = (PHI[key], 'MultiPak') if key in PHI else (quant(hires), 'this script')

# ---------------- surveys
def load_survey(rel):
    b = F.load(rel)['Su1s']
    x, y = b[:, 0], b[:, 1]
    # drop isolated doubled channels (261006 Control strip survey: every other point x2 over 1089-1065 eV)
    nb = np.r_[y[1], (y[:-2] + y[2:]) / 2, y[-2]]
    keep = ~((y > 1.5 * nb) & (x > 20))
    if (~keep).any(): print(f'{rel}: dropped {(~keep).sum()} doubled survey points, {x[~keep].min():.1f}-{x[~keep].max():.1f} eV')
    return x[keep], y[keep]

n = len(SAMPLES); DY = 0.7
fig = plt.figure(figsize=(7.2, 9.0))
gs = fig.add_gridspec(5, 3, height_ratios=[1.1, 0.03, 1, 1, 1], hspace=0.58, wspace=0.2,   # row 1 is a spacer
                      left=0.08, right=0.985, top=0.935, bottom=0.045)
axA = fig.add_subplot(gs[0, :])
ytop = DY * (n - 1) + 1.15
for i, (key, sur, hires, name) in enumerate(SAMPLES):
    x, y = load_survey(sur)
    x = x + F.results[key]['C1s']['shift']           # same charge correction as the hi-res fits
    m = (x > 278) & (x < 292); hC = y[m].max() - np.median(y[(x > 292) & (x < 296)])
    axA.plot(x, np.minimum(y / hC, ytop - DY * i) + DY * i, color=COL[key], lw=0.7)
    axA.text(262, DY * i + 0.13, name, color=COL[key], fontsize=10.4, ha='left', va='bottom')
for el, be, dy in [('O KLL', 978, 0), ('O 1s', 532, 0), ('N 1s', 400, 0), ('C 1s', 285, 0),
                   ('Si 2s', 153, 0), ('P 2p', 133, 0.12), ('Si 2p', 102, 0), ('Na 1s', 1072, 0)]:
    axA.axvline(be, color='#bbb', lw=0.4, ls=(0, (2, 3)), zorder=0)
    axA.text(be, ytop + 0.05 + dy, el, rotation=90, fontsize=8.45, color='#555', ha='center', va='bottom')
axA.set_xlim(1100, 0); axA.set_ylim(-0.05, ytop + 0.45); axA.set_yticks([])
axA.set_xlabel('Binding Energy (eV)'); axA.set_ylabel('Intensity (norm., offset)')
for s in ('top', 'right', 'left'): axA.spines[s].set_visible(False)

# ---------------- hi-res fits: rows = regions, cols = samples; residual strip under each fit
REG = [('C1s', 'C 1s', ['C–C / C–H', 'C–O', 'O–C–O / N–C=O']),
       ('O1s', 'O 1s', ['C=O / N–C=O', 'C–O–C / C–OH']),
       ('N1s', 'N 1s', ['N–C=O', 'NH$_3^+$'])]
BEST = {'O1s': 2, 'N1s': 2}
XLIM = {'C1s': (294, 277), 'O1s': (539, 526), 'N1s': (407, 394)}
XTICK = {'C1s': [292, 288, 284, 280], 'O1s': np.arange(528, 539, 2), 'N1s': np.arange(396, 407, 2)}
for r, (reg, title, names) in enumerate(REG):
    for c, (key, sur, hires, name) in enumerate(SAMPLES):
        sub = gs[r + 2, c].subgridspec(2, 1, height_ratios=[5.2, 1], hspace=0)
        a = fig.add_subplot(sub[0]); ar = fig.add_subplot(sub[1], sharex=a)
        if r == 0 and c == 0: a00 = a
        R = F.results[key][reg]
        f = R['fit'] if reg == 'C1s' else R['fits'][BEST[reg]]
        x, y, bg = R['x'], R['y'], R['bg']
        fitc = bg + f['model']; res = y - fitc
        a.plot(x, y, '.', color='k', ms=2.2, mew=0, label='Raw Data', zorder=4)
        a.plot(x, bg, color='#888', lw=0.6, ls='--', label='Background')
        a.plot(x, fitc, color=FITC, lw=0.7, label='Overall Fit', zorder=3)
        for k, t in enumerate(f['comps']):
            if f['frac'][k] < 0.5: continue
            a.fill_between(x, bg, bg + F.gl(x, *t), facecolor=CFILL[k], alpha=0.42, edgecolor=CFILL[k], lw=0.5,
                           label=f'{names[k]}  {t[0]:.1f} eV  {f["frac"][k]:.0f}%')
        top, base = max(fitc.max(), y.max()), min(bg.min(), y.min())
        a.set_xlim(*XLIM[reg]); a.set_ylim(base - 0.06 * (top - base), base + (2.6 if reg == 'C1s' else 2.3) * (top - base))
        a.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)
        a.yaxis.get_offset_text().set_fontsize(5.5)
        a.tick_params(labelsize=7.8, length=2.5, labelbottom=False)
        a.legend(loc='upper left', frameon=False, fontsize=6.5, handlelength=1.6, labelspacing=0.3, borderaxespad=0.4)
        a.text(0.975, 0.955, title, transform=a.transAxes, ha='right', va='top', fontsize=9.75)
        a.text(0.975, 0.855, f'{title[0]} {AT[key][{"C1s": 0, "N1s": 1, "O1s": 2}[reg]]:.1f} at.%', transform=a.transAxes,
               ha='right', va='top', fontsize=7.15, color='#444')
        # residual (data - fit); RSD = residual standard deviation, chi2 = reduced chi-square with Poisson sigma
        sigma = np.sqrt(np.maximum(y, 1) / dwell(hires)[reg]); nfree = 2 * len(f['comps']) + 1
        chi2 = ((res / sigma) ** 2).sum() / (len(x) - nfree)
        ar.axhline(0, color='#bbb', lw=0.4)
        ar.plot(x, res, color=RESC, lw=0.6)
        lim = np.abs(res).max() * 1.25; ar.set_ylim(-lim, lim)
        ar.ticklabel_format(axis='y', style='sci', scilimits=(0, 0), useMathText=True)
        ar.yaxis.get_offset_text().set_visible(False)
        ar.set_yticks([]); ar.tick_params(labelsize=7.8, length=2.5)
        ar.set_xticks(XTICK[reg])
        ar.text(0.985, 0.5, f'RSD {res.std():.0f} CPS,  $\\chi^2_\\nu$ {chi2:.1f}', transform=ar.transAxes, ha='right', va='center',
                fontsize=6.5, color='#444', bbox=dict(fc='white', ec='none', pad=0.6, alpha=0.8))
        ar.set_xlabel('Binding Energy (eV)', fontsize=8.45, labelpad=1.5)
        if c == 0:
            a.set_ylabel('Intensity (CPS)', fontsize=8.45); ar.set_ylabel('Res.', fontsize=7.8)
        a.set_title(name, fontsize=11.05, color=COL[key], pad=11)
        if c == 0:
            a.text(-0.2, 1.13, f'({"bcd"[r]}) {title}', transform=a.transAxes, fontsize=13, va='bottom')

# panel (a) label: same figure-x as the (b)/(c)/(d) row labels (axes x = -0.2 of the first column)
xf = fig.transFigure.inverted().transform(a00.transAxes.transform((-0.2, 0)))[0]
fig.text(xf, axA.get_position().y1 + 0.024, '(a) Survey spectra', fontsize=13, va='bottom', ha='left')
fig.savefig(os.path.join(FIG, 'XPS_figure.png'), dpi=300)
fig.savefig(os.path.join(FIG, 'XPS_figure.pdf'))

# ---------------- summary table
rows = []
for key, sur, hires, name in SAMPLES + EXTRA:
    C, N, O = AT[key]
    R = F.results[key]
    c1 = R['C1s']['fit']['frac']; o1 = R['O1s']['fits'][2]['frac']; n1 = R['N1s']['fits'][2]['frac']
    rows.append(dict(sample=name, instrument_label=key, at_source=AT_SRC[key], C_at=C, O_at=O, N_at=N, O_C=O / C, N_C=N / C,
                     C1s_CC=c1[0], C1s_CO=c1[1], C1s_OCO_amide=c1[2],
                     O1s_CO_amide=o1[0], O1s_CO=o1[1], N1s_amide=n1[0], N1s_NH3=n1[1],
                     charge_shift_eV=R['C1s']['shift']))
with open(os.path.join(OUT, 'XPS_table.csv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader()
    for r in rows: w.writerow({k: (f'{v:.2f}' if isinstance(v, float) else v) for k, v in r.items()})
for r in rows:
    print(' | '.join(f'{v:.2f}' if isinstance(v, float) else v for v in r.values()))
print('saved XPS_figure.png / .pdf / XPS_table.csv')
