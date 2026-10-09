"""Careful C 1s / O 1s / N 1s fitting of the three hi-res PHI datasets in ../raw (MultiPak csv exports where
they exist, otherwise the .spe files; see ../../README.md for provenance).
Shirley background, GL(30) peaks (KherveFitting's GL definition: (1-m)*Gauss + m*Lorentz),
shared FWHM within a region, positions bounded. Model selection by comparing 1/2/3-peak fits."""
import os, sys, json, struct
import numpy as np
from scipy.optimize import least_squares
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
import re
XPS = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
D = os.path.join(XPS, 'raw')   # PHI csv exports

def read_csv_blocks(path):
    blocks, cur, name = {}, [], None
    for line in open(path, encoding='latin-1'):
        s = line.strip()
        if not s:
            continue
        parts = s.split(',')
        if len(parts) == 2 and all(re.fullmatch(r'-?\d+(\.\d+)?', p.strip()) for p in parts):
            cur.append((float(parts[0]), float(parts[1])))
        else:
            if cur and name:
                blocks[name] = np.array(cur)
            cur = []
            if re.fullmatch(r'[A-Z][a-z]?\d[spdf]', s) or s == 'Su1s':
                name = s
    if cur and name:
        blocks[name] = np.array(cur)
    for k, a in blocks.items():
        keep = ~((a[:, 0] == 0) & (a[:, 1] == 0))
        blocks[k] = a[keep]
    return blocks

def read_spe(path):
    """PHI MultiPak .spe: text header up to EOFH, then a 16-byte block header, one 96-byte record per
    region (n points, byte count, byte offset) and float32 counts. Same dict of (BE, counts) blocks as
    read_csv_blocks; checked point-for-point against the 260904 csv exports."""
    raw = open(path, 'rb').read()
    j = raw.index(b'\n', raw.index(b'EOFH')) + 1
    regs = [l.split() for l in raw[:j].decode('latin-1').splitlines() if l.startswith('SpectralRegDef:')]
    body = raw[j:]; blocks = {}
    for k in range(struct.unpack('<4I', body[:16])[1]):
        h = struct.unpack('<24I', body[16 + 96 * k:16 + 96 * (k + 1)])
        y = np.frombuffer(body[h[20]:h[20] + h[19]], '<f4').astype(float)
        name, step, start = regs[k][3], float(regs[k][6]), float(regs[k][7])
        blocks[name] = np.column_stack([start + step * np.arange(h[5]), y])
    return blocks

def load(rel):
    path = os.path.join(XPS, rel)
    return read_spe(path) if rel.endswith('.spe') else read_csv_blocks(path)

OUT = os.path.dirname(os.path.abspath(__file__))
SAMPLES = [('WT untreated', 'raw/WT_untreated_hires.csv'),
           ('WT NaOH treated', 'raw/WT_NaOH_treated_hires.spe'),
           ('BslA NaOH treated', 'raw/BslA_NaOH_treated_hires.csv')]
CO_MAX = {'BslA NaOH treated', 'WT NaOH treated'}        # C 1s maximum is the C-O component, not C-C
CROP_OFF = {'WT NaOH treated': -2.0}                     # 261006 session charged ~2 eV more; slide the O/N windows
M_GL = 0.30   # GL(30)

def gl(x, pos, h, fwhm, m=M_GL):
    s = fwhm / (2 * np.sqrt(2 * np.log(2)))
    g = np.exp(-0.5 * ((x - pos) / s) ** 2)
    l = 1.0 / (1 + ((x - pos) / (fwhm / 2)) ** 2)
    return h * ((1 - m) * g + m * l)

def shirley(x, y, n_end=5, iters=50):
    """Iterative Shirley background over the full window; endpoints = mean of n_end points at each side.
    x descends (high BE first)."""
    yl, yr = y[:n_end].mean(), y[-n_end:].mean()        # left = high BE, right = low BE
    b = np.full_like(y, yr)
    for _ in range(iters):
        net = y - b
        # integral from low-BE end up to each point (x descending => cumulative from the right)
        cum = np.cumsum(net[::-1])[::-1]
        b = yr + (yl - yr) * cum / cum[0]
    return b

def fit_region(x, ynet, peaks, share_fwhm=True, fwhm0=1.5, fwhm_bounds=(0.9, 2.6)):
    """peaks: list of dict(pos, lo, hi, name). Returns params, covariance-based sigmas, model, stats."""
    n = len(peaks)
    scale = ynet.max(); ynet = ynet / scale
    p0, lo, hi = [], [], []
    for p in peaks:
        p0 += [p['pos'], max(ynet[(x > p['lo']) & (x < p['hi'])].max(), 1.0) * 0.8]
        lo += [p['lo'], 0.0]; hi += [p['hi'], ynet.max() * 1.5]
    if share_fwhm:
        p0 += [fwhm0]; lo += [fwhm_bounds[0]]; hi += [fwhm_bounds[1]]
    else:
        p0 += [fwhm0] * n; lo += [fwhm_bounds[0]] * n; hi += [fwhm_bounds[1]] * n
    w = 1 / np.sqrt(np.maximum(ynet + np.median(ynet[ynet > 0]) * 0.1, 1))   # ~Poisson weights

    def unpack(q):
        out = []
        for i in range(n):
            f = q[2 * n] if share_fwhm else q[2 * n + i]
            out.append((q[2 * i], q[2 * i + 1], f))
        return out
    def model(q):
        return sum(gl(x, *t) for t in unpack(q))
    def resid(q):
        return (model(q) - ynet) * w
    r = least_squares(resid, p0, bounds=(lo, hi), x_scale='jac', max_nfev=20000, ftol=1e-12, xtol=1e-12)
    J = r.jac; dof = max(len(x) - len(r.x), 1)
    chi2 = (r.fun ** 2).sum() / dof
    try:
        cov = np.linalg.inv(J.T @ J) * chi2
        sig = np.sqrt(np.diag(cov))
    except np.linalg.LinAlgError:
        sig = np.full_like(r.x, np.nan)
    comps = [(p, h * scale, f) for p, h, f in unpack(r.x)]
    step = abs(np.median(np.diff(x)))
    areas = np.array([gl(x, *t).sum() * step for t in comps])
    rms = np.sqrt(((model(r.x) - ynet) ** 2).mean()) / ynet.max() * 100
    return dict(params=r.x, sig=sig, comps=comps, areas=areas, frac=100 * areas / areas.sum(),
                model=model(r.x) * scale, chi2=chi2, rms=rms, n=n, share=share_fwhm)

def crop(block, lo, hi):
    x, y = block[:, 0], block[:, 1]
    m = (x >= lo) & (x <= hi)
    return x[m], y[m]

results = {}
for label, base in SAMPLES:
    b = load(base); off = CROP_OFF.get(label, 0.0)
    res = {}
    # ---------- C 1s: 3 peaks, C-C / C-O / O-C-O(+amide), fixed offsets like the user's Kherve constraints
    x, y = crop(b['C1s'], 273, 292)
    bg = shirley(x, y); yn = y - bg
    # raw C-C guess = lowest-BE resolvable feature
    xmax = x[np.argmax(yn)]
    guessA = xmax - 1.7 if label in CO_MAX else xmax
    pk = [dict(name='C-C', pos=guessA, lo=guessA - 0.4, hi=guessA + 0.4),
          dict(name='C-O', pos=guessA + 1.7, lo=guessA + 1.3, hi=guessA + 2.1),
          dict(name='O-C-O/N-C=O', pos=guessA + 3.2, lo=guessA + 2.8, hi=guessA + 3.6)]
    fc = fit_region(x, yn, pk, fwhm0=1.5)
    shift = 284.8 - fc['comps'][0][0]
    fc['comps'] = [(p + shift, h, f) for p, h, f in fc['comps']]
    res['C1s'] = dict(x=x + shift, y=y, bg=bg, fit=fc, names=[p['name'] for p in pk], shift=shift)

    # ---------- O 1s: compare 1, 2, 3 components
    x, y = crop(b['O1s'], 524 + off, (538 if label != 'BslA NaOH treated' else 537) + off)
    x = x + shift
    bg = shirley(x, y); yn = y - bg
    fits = {}
    fits[1] = fit_region(x, yn, [dict(name='O 1s', pos=532.6, lo=531.5, hi=533.5)], fwhm0=1.7)
    fits[2] = fit_region(x, yn, [dict(name='C=O / N-C=O', pos=531.5, lo=531.0, hi=532.0),
                                 dict(name='C-O-C / C-OH', pos=532.8, lo=532.4, hi=533.2)], fwhm0=1.6)
    fits[3] = fit_region(x, yn, [dict(name='C=O / N-C=O', pos=531.4, lo=530.9, hi=531.9),
                                 dict(name='C-O-C / C-OH', pos=532.7, lo=532.3, hi=533.1),
                                 dict(name='O-C-O / H2O', pos=533.8, lo=533.4, hi=534.4)], fwhm0=1.5)
    res['O1s'] = dict(x=x, y=y, bg=bg, fits=fits)

    # ---------- N 1s: compare 1, 2 components
    x, y = crop(b['N1s'], 393 + off, 406 + off)
    x = x + shift
    bg = shirley(x, y); yn = y - bg
    fits = {}
    fits[1] = fit_region(x, yn, [dict(name='N-C=O amide', pos=400.0, lo=399.5, hi=400.5)], fwhm0=1.7)
    fits[2] = fit_region(x, yn, [dict(name='N-C=O amide', pos=400.0, lo=399.6, hi=400.4),
                                 dict(name='NH3+ / R-NH2+', pos=401.7, lo=401.2, hi=402.3)], fwhm0=1.6)
    res['N1s'] = dict(x=x, y=y, bg=bg, fits=fits)
    results[label] = res

# ---------- report
def report(label, res):
    print(f'\n===== {label}   (charge shift {res["C1s"]["shift"]:+.2f} eV so that C-C = 284.80)')
    fc = res['C1s']['fit']; sh = res['C1s']['shift']
    print(f'C 1s  FWHM {fc["comps"][0][2]:.2f} eV   RMS resid {fc["rms"]:.2f}% of max')
    for nm, (p, h, f), a, fr in zip(res['C1s']['names'], fc['comps'], fc['areas'], fc['frac']):
        print(f'   {nm:16s} {p:7.2f} eV   {fr:5.1f}%')
    for reg in ('O1s', 'N1s'):
        print(f'{reg[:-2]} {reg[-2:]}')
        for n, f in res[reg]['fits'].items():
            print(f'   {n}-peak model: chi2_red {f["chi2"]:.3f}  RMS {f["rms"]:.2f}%  FWHM {f["comps"][0][2]:.2f}')
            names = {1: ['single'], 2: ['C=O/N-C=O', 'C-O-C/C-OH'], 3: ['C=O/N-C=O', 'C-O-C/C-OH', 'O-C-O/H2O']}[n] if reg == 'O1s' \
                else {1: ['N-C=O amide'], 2: ['N-C=O amide', 'NH3+']}[n]
            for nm, (p, h, fw), fr, sp in zip(names, f['comps'], f['frac'], f['sig'][0::2]):
                print(f'        {nm:14s} {p:7.2f} ± {sp:.2f} eV   {fr:5.1f}%')

for label, res in results.items():
    report(label, res)

# ---------- figure: rows = samples, cols = C1s / O1s(best=2) / N1s(best chosen below)
BEST = {'O1s': 2, 'N1s': 2}
fig, axes = plt.subplots(len(results), 3, figsize=(15, 3.7 * len(results)))
fills = ['#2a78d6', '#1baf7a', '#eb6834', '#eda100']
for i, (label, res) in enumerate(results.items()):
    for j, reg in enumerate(['C1s', 'O1s', 'N1s']):
        a = axes[i, j]; r = res[reg]
        f = r['fit'] if reg == 'C1s' else r['fits'][BEST[reg]]
        x, y, bg = r['x'], r['y'], r['bg']
        a.plot(x, y, 'k.', ms=3, label='data')
        a.plot(x, bg, color='#999', lw=1, ls='--', label='Shirley')
        for k, t in enumerate(f['comps']):
            a.fill_between(x, bg, bg + gl(x, *t), color=fills[k], alpha=0.35, lw=0)
            a.text(t[0], bg[np.argmin(abs(x - t[0]))] + t[1] * 1.03, f'{t[0]:.2f}\n{f["frac"][k]:.0f}%', ha='center', va='bottom', fontsize=8)
        a.plot(x, bg + f['model'], color='#333', lw=1.2, label='fit')
        a.set_xlim(x.max(), x.min()); a.set_ylim(top=(bg + f['model']).max() * 1.28)
        a.set_title(f'{label}  {reg[:-2]} {reg[-2:]}   FWHM {f["comps"][0][2]:.2f}  RMS {f["rms"]:.1f}%', fontsize=10, loc='left')
        a.set_yticks([]); a.set_xlabel('Binding energy (eV)')
        for s in ('top', 'right', 'left'): a.spines[s].set_visible(False)
        if i == 0 and j == 0: a.legend(frameon=False, fontsize=8, loc='upper left')
fig.tight_layout(); fig.savefig(os.path.join(OUT, 'hires_fits.png'), dpi=130)

# dump numbers for the later figure
dump = {}
for label, res in results.items():
    dump[label] = dict(shift=res['C1s']['shift'],
                       C1s=[dict(name=n, pos=float(c[0]), fwhm=float(c[2]), frac=float(fr)) for n, c, fr in zip(res['C1s']['names'], res['C1s']['fit']['comps'], res['C1s']['fit']['frac'])],
                       O1s={n: [dict(pos=float(c[0]), fwhm=float(c[2]), frac=float(fr)) for c, fr in zip(f['comps'], f['frac'])] for n, f in res['O1s']['fits'].items()},
                       N1s={n: [dict(pos=float(c[0]), fwhm=float(c[2]), frac=float(fr)) for c, fr in zip(f['comps'], f['frac'])] for n, f in res['N1s']['fits'].items()})
json.dump(dump, open(os.path.join(OUT, 'hires_fits.json'), 'w'), indent=1)
print('\nsaved hires_fits.png / .json')
