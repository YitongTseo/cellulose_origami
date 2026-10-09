"""Texture descriptors for the four TM3000 micrographs in cellulose_closeups/ (Figure 2 of the manuscript).

For each image: pixel size from the embedded scale bar; slow illumination shading removed with a wide
Gaussian high-pass; then
  * grey-level contrast      = std / mean of the flattened image (per cent)
  * orientation coherence    = median structure-tensor coherence in 32-px windows (0 isotropic, 1 aligned)
  * correlation length       = lag at which the radially averaged autocorrelation of the flattened image
                               first falls to 1/e (micrometres) -- the characteristic feature size
Writes SEM_texture.csv and SEM_autocorrelation.png (radial autocorrelation profiles)."""
import os, csv
import numpy as np
from PIL import Image
from scipy import ndimage as ndi
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
IMG = os.path.join(HERE, '..', '..', '..', 'cellulose_closeups')
FIG = os.path.join(HERE, '..', 'figures')
# (file, display name, magnification, scale-bar label in um). Scale-bar length in px is read from the footer.
IMAGES = [('03.28.25_SEM_Control.jpg', 'Untreated WT BC', '300x', 300),
          ('03.28.25_SEM_BslA_engineered.jpg', 'NaOH-treated BslA-activated BC', '300x', 300),
          ('03.28.25_SEM_Control_2.jpg', 'Untreated WT BC', '1500x', 50),
          ('03.28.25_SEM_BslA_engineered_2.jpg', 'NaOH-treated BslA-activated BC', '1200x', 50)]
COL = {'Untreated WT BC': '#F5B041', 'NaOH-treated BslA-activated BC': '#5DADE2'}
FOOTER = 955          # image rows below this are the instrument footer

def scale_bar_px(im):
    bar = im[962:975, :] < 60
    cols = np.where(bar.mean(0) > 0.5)[0]
    return cols.max() - cols.min() + 1

def radial_autocorr(f):
    f = f - f.mean()
    F = np.fft.fft2(f)
    ac = np.fft.fftshift(np.fft.ifft2(np.abs(F) ** 2).real); ac /= ac.max()
    cy, cx = np.array(ac.shape) // 2
    yy, xx = np.indices(ac.shape); r = np.hypot(yy - cy, xx - cx).astype(int)
    prof = np.bincount(r.ravel(), ac.ravel()) / np.bincount(r.ravel())
    return prof[:min(ac.shape) // 2]

rows, profiles = [], []
for fn, name, mag, bar_um in IMAGES:
    im = np.asarray(Image.open(os.path.join(IMG, fn)).convert('L')).astype(float)
    um_px = bar_um / scale_bar_px(im)
    img = im[:FOOTER, :]
    g = ndi.gaussian_filter(img, 3)
    flat = g - ndi.gaussian_filter(img, 0.1 * img.shape[1])          # high-pass: remove shading
    contrast = 100 * flat.std() / g.mean()
    gx, gy = ndi.sobel(ndi.gaussian_filter(img, 1.5), 1), ndi.sobel(ndi.gaussian_filter(img, 1.5), 0)
    w = 32
    A, B, C = (ndi.uniform_filter(gx * gx, w), ndi.uniform_filter(gy * gy, w), ndi.uniform_filter(gx * gy, w))
    coherence = np.median(np.sqrt((A - B) ** 2 + 4 * C ** 2) / (A + B + 1e-9))
    prof = radial_autocorr(flat)
    L = np.argmax(prof < 1 / np.e) * um_px
    rows.append(dict(file=fn, sample=name, magnification=mag, um_per_px=round(um_px, 3),
                     contrast_pct=round(contrast, 1), orientation_coherence=round(coherence, 2), correlation_length_um=round(L, 1)))
    profiles.append((name, mag, um_px, prof, L))
    print(f'{name:32s} {mag:6s} {um_px:.3f} um/px | contrast {contrast:4.1f}% | coherence {coherence:.2f} | correlation length {L:4.1f} um')

with open(os.path.join(HERE, 'SEM_texture.csv'), 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.8))
for ax, want in zip(axes, ('300x', 'hi')):
    for name, mag, um_px, prof, L in profiles:
        if (mag == '300x') != (want == '300x'): continue
        lag = np.arange(len(prof)) * um_px
        m = lag <= (60 if want == '300x' else 12)
        ax.plot(lag[m], prof[m], color=COL[name], lw=1.2, label=f'{name} ({mag}), $L$ = {L:.1f} µm')
        ax.plot(L, 1 / np.e, 'o', color=COL[name], ms=4)
    ax.axhline(1 / np.e, color='#999', lw=0.5, ls='--')
    ax.text(ax.get_xlim()[1] * 0.98, 1 / np.e + 0.03, '1/e', ha='right', fontsize=7, color='#666')
    ax.set_xlabel('Lag (µm)'); ax.set_ylabel('Autocorrelation'); ax.set_ylim(-0.25, 1.02)
    ax.legend(frameon=False, fontsize=6.5, loc='upper right')
    ax.set_title('A   300×' if want == '300x' else 'B   1200× / 1500×', loc='left', fontsize=9, fontweight='bold')
    for s in ('top', 'right'): ax.spines[s].set_visible(False)
fig.tight_layout()
fig.savefig(os.path.join(FIG, 'SEM_autocorrelation.png'), dpi=300); fig.savefig(os.path.join(FIG, 'SEM_autocorrelation.pdf'))
print('saved SEM_texture.csv, SEM_autocorrelation.png/.pdf')
