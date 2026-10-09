# XPS text drafts (caption, methods, results) — working draft

Sample naming below follows the instrument labels (Day 0, Control, 48 hrs). Bracketed text marks things to confirm before submission.

---

## Figure caption

**Figure X. Surface chemistry of bacterial cellulose (BC) films by X-ray photoelectron spectroscopy.**
(A) Survey spectra of untreated wild-type BC (Day 0, Control) and NaOH-treated BslA-expressing BC (48 hrs), normalised to the C 1s peak height and offset for clarity. Dotted lines mark the principal core levels; sodium and phosphorus are detected only on the BslA sample. (B–D) High-resolution C 1s, O 1s and N 1s regions with Shirley backgrounds (dashed) and Gaussian–Lorentzian component fits. C 1s is resolved into C–C/C–H (284.8 eV), C–O (286.2–286.6 eV) and O–C–O/N–C=O (288.0 eV); O 1s into C=O/amide (531.3–531.5 eV) and C–O–C/C–OH (532.5–532.8 eV); N 1s into amide N–C=O (399.8–400.0 eV) and protonated amine NH₃⁺ (401.2–402.1 eV). Component areas are given as percentages of the respective element. Binding energies are charge-referenced to the C–C component at 284.8 eV. The BslA sample shows a 2.3-fold higher N/C ratio (0.12 vs 0.05) and a more cellulose-like O/C ratio (0.59 vs 0.31) than the wild-type samples, consistent with a surface-localised protein layer on purified cellulose.

---

## Methods paragraph

**X-ray photoelectron spectroscopy.** XPS was performed on a PHI VersaProbe II (Physical Electronics) at the MIT.nano Surface Analysis Facility using a monochromatic Al Kα source (1486.6 eV, 50 W, 200 µm spot) at a 45° take-off angle, with dual-beam (electron and Ar⁺) charge neutralisation. Dried BC films were mounted on glass slides [confirm: mounting/adhesive] and analysed without further preparation. Survey spectra were acquired at a pass energy of 187.85 eV with a 0.8 eV step (3 sweeps); high-resolution C 1s, O 1s and N 1s spectra at 23.5 eV pass energy with a 0.1 eV step. Surface atomic concentrations were determined from the survey spectra in MultiPak (PHI) using instrument-specific relative sensitivity factors. High-resolution spectra were processed in [KherveFitting v1.80 / custom Python routine]: a Shirley background was subtracted from each region, and components were fitted with Gaussian–Lorentzian (GL(30)) line shapes with a common full width at half maximum within each region. Binding energies were charge-referenced by setting the aliphatic C–C/C–H component of the C 1s envelope to 284.8 eV (shifts of +1.8 to +2.4 eV). C 1s components were constrained to C–O at +1.3 to +2.1 eV and O–C–O/N–C=O at +2.8 to +3.6 eV relative to C–C. The number of O 1s and N 1s components was chosen by comparing residuals of one-, two- and three-component models; additional components were retained only where they reduced the root-mean-square residual. One survey and one set of high-resolution regions were acquired per sample [n = 1 spot per sample; state explicitly].

---

## Results text (working draft)

### Surface composition
XPS surveys detected carbon, oxygen and nitrogen on all films, with trace silicon, phosphorus and sodium (Table X, Fig. XA). The two wild-type samples (Day 0 and Control) were indistinguishable in composition (C 73 %, O 23 %, N 4 %; O/C = 0.31–0.32, N/C = 0.05). Their O/C ratio is far below that of pure cellulose (0.83), and their C 1s envelopes are dominated by aliphatic C–C/C–H carbon (61–71 % of C), indicating that the outermost ~10 nm of the untreated pellicles is covered by a hydrocarbon-rich overlayer, most plausibly residual cell envelope, lipid and culture-medium components that a purification step would remove. The BslA-expressing sample (48 hrs) differed markedly: oxygen rose to 34 % and nitrogen to 7.2 %, giving O/C = 0.59 and N/C = 0.12. Its C 1s envelope is dominated by C–O carbon (53 %), the signature of the anhydroglucose backbone, with a reduced C–C/C–H component (28 %) and an enlarged O–C–O/N–C=O component (19 %). Because cellulose contains no nitrogen, the N 1s signal is a direct measure of surface protein; the 2.3-fold increase in N/C on the BslA sample relative to wild type is consistent with a BslA layer at the film surface.

### Chemical state of nitrogen and oxygen
The N 1s peak was centred at 399.8–400.0 eV in all samples, characteristic of amide nitrogen (peptide bonds). On the BslA sample the envelope was narrower (FWHM 1.44 vs 1.8–2.0 eV) and carried a distinct high-binding-energy shoulder at 401.2 eV (13 % of N), assigned to protonated amine groups (lysine/arginine side chains); this shoulder was weak (5 %) or absent on the wild-type films. The O 1s envelopes of all three samples were well described by two components, a major C–O–C/C–OH component at 532.5–532.8 eV (69–73 % of O) and a minor C=O/amide component at 531.3–531.5 eV (27–31 %); a third component (acetal or adsorbed water, ~533.5 eV) did not improve the fit. The O 1s line shape was thus essentially identical across samples; the samples differ in the *amount* of oxygen-bearing surface (O/C), not in the type of oxygen functionality.

### Caveats to state
- Only one analysis spot per sample was measured; the differences reported (N/C 0.05 vs 0.12; O/C 0.31 vs 0.59) are large relative to typical spot-to-spot variability on such films (~10–20 % relative), but replicate spots would strengthen the claim.
- Silicon (3.5 at.%, Si 2p ≈ 101.6 eV, consistent with siloxane) was detected on the Day 0 sample and is attributed to handling contamination.
- Sodium (0.5 at.%) and phosphorus (1.4 at.%) were detected only on the BslA sample, consistent with residual NaOH/buffer from its treatment.
- [Prep history: the intended "NaOH-washed wild-type" arm could not be identified with confidence. The data support treating Day 0 and Control as two untreated wild-type replicates. If the "strips" surveys (O/C ≈ 0.5, N ≈ 0.5–1 %) correspond to NaOH-washed wild-type material, they would fill that arm; note that no high-resolution regions were acquired for them.]
- The reviewer's Cl 2p control for antibiotic carryover: Cl 2p was not detected above noise in any survey (≤0.05 at.%).

### One-line summary for the response to reviewers
XPS confirms a surface-localised, amide-rich nitrogen signal (N/C = 0.12) on the BslA-expressing cellulose that is 2.3× that of wild-type films (N/C = 0.05), together with a more cellulose-like O/C ratio, supporting deposition of BslA on the purified BC surface.
