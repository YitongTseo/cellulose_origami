# FTIR tables (working draft)

## Table 1. ATR-FTIR band assignments used in Figure X

| Wavenumber (cm⁻¹) | Assignment | Origin | Ref. |
|---|---|---|---|
| 3340 | ν(O–H), hydrogen-bonded hydroxyl | cellulose | Oh 2005; Kačuráková 2002 |
| 2900 | ν(C–H), CH and CH₂ | cellulose | Oh 2005 |
| 1740 | ν(C=O), ester carbonyl | lipid / medium residue | Naumann 1991 |
| 1650 | amide I, ν(C=O) of peptide bond (overlaps δ(H–O–H) of adsorbed water) | protein | Barth 2007 |
| 1540 | amide II, δ(N–H) + ν(C–N) | protein | Barth 2007 |
| 1430 | δ(CH₂) scissoring, "crystallinity band" (cellulose I); shifts to ~1420 in cellulose II / amorphous | cellulose | Nelson & O'Connor 1964; Oh 2005 |
| 1372 | δ(C–H) bending | cellulose | Nelson & O'Connor 1964 |
| 1315 | ω(CH₂) wagging | cellulose | Oh 2005 |
| 1160 | ν_as(C–O–C), β-1,4-glycosidic link | cellulose | Kačuráková 2002 |
| 1108 | ν(C–O), ring, antisymmetric | cellulose | Kačuráková 2002 |
| 1055 / 1030 | ν(C–O) at C3 / C6 | cellulose | Kačuráková 2002 |
| 898 | C–O–C deformation at β-glycosidic linkage, "amorphous band" | cellulose | Nelson & O'Connor 1964 |
| 750 / 710 | Iα / Iβ allomorph marker bands (weak; not resolved here) | cellulose | Sugiyama 1991; Yamamoto 1996 |

Suggested references (check against your bibliography):
- Nelson, M. L.; O'Connor, R. T. *J. Appl. Polym. Sci.* **1964**, 8, 1311–1324 (TCI, LOI; 1430/898).
- Oh, S. Y. et al. *Carbohydr. Res.* **2005**, 340, 2376–2391 (cellulose I/II assignments).
- Kačuráková, M. et al. *Carbohydr. Res.* **2002**, 337, 1145–1153 (cellulose C–O region).
- Barth, A. *Biochim. Biophys. Acta* **2007**, 1767, 1073–1101 (amide bands).
- Naumann, D. et al. *Nature* **1991**, 351, 81–82 (bacterial FTIR: ester C=O, amide).
- Castro, C. et al. *Carbohydr. Polym.* **2011**, 84, 96–102 (bacterial cellulose FTIR, indices).
- Sugiyama, J. et al. *Macromolecules* **1991**, 24, 2461–2466; Yamamoto, H. et al. *Cellulose* **1996**, 3, 229–242 (Iα/Iβ).

## Table 2. FTIR-derived indices (spectra normalised to the 1030 cm⁻¹ band; peak heights above local linear baselines)

Spectra: `raw/WT_untreated`, `raw/WT_NaOH_treated`, `raw/BslA_NaOH_treated` (`.SPA` instrument files and `.csv` exports).

| Sample | δ(CH₂) position (cm⁻¹) | TCI (A₁₃₇₂/A₂₉₀₀) | LOI (A₁₄₃₀/A₈₉₈) | HBI (A₃₄₀₀/A₁₃₂₀) | Amide I / A₁₀₃₀ | Amide II / A₁₀₃₀ | Amide I position (cm⁻¹) |
|---|---|---|---|---|---|---|---|
| WT untreated | 1428 | 0.22 | 2.76 | 1.97 | 0.27 | 0.11 | 1650 |
| WT NaOH treated | 1428 | 0.29 | 3.16 | 2.28 | 0.16 | 0.06 | 1650 |
| BslA NaOH treated | **1422** | 0.07* | n.d.* | 1.30* | **0.39** | **0.14** | **1631** |

\* Cellulose bands in the BslA NaOH treated spectrum are broadened by the protein overlayer; the 898 and 1372 bands cannot be baselined reliably, so these indices are not comparable and should be reported as "not determined" or with this caveat.

TCI: total crystallinity index; LOI: lateral order index; HBI: hydrogen-bond intensity (Nelson & O'Connor 1964; Oh 2005). Note: indices from ATR spectra are relative, not absolute crystallinities; compare within this dataset only.

### How to read Table 2
- **Amide I / A₁₀₃₀** is the protein-loading number: WT NaOH treated (0.16) < WT untreated (0.27) < BslA NaOH treated (0.39). Amide II, which is free of the adsorbed-water overlap at ~1640 cm⁻¹, gives the same order (0.06 < 0.11 < 0.14). The order matches the XPS N/C ordering, but the NaOH wash removes only ~40% of the amide signal here versus ~6× less surface nitrogen in XPS (ATR samples ~1 µm deep, XPS ~10 nm).
- **TCI** rises from 0.22 to 0.29 on NaOH treatment of WT, consistent with removal of non-cellulosic material; LOI and HBI rise slightly as well.
- **δ(CH₂) position**: 1428 cm⁻¹ for both WT samples, i.e. unchanged by NaOH. BslA NaOH treated sits at 1422; a 6 cm⁻¹ downshift with no new O–H shoulder at ~3440–3480 cm⁻¹ is more consistent with band broadening from the protein overlayer / reduced order than with a cellulose I → II transition. This directly addresses the reviewer's question about NaOH-induced allomorph change: no cellulose II signature is observed in the NaOH-treated WT sample.
- **Amide I position** (maximum after removing a linear baseline between 1800 and 1480 cm⁻¹; marked in Figure panel C): 1650 cm⁻¹ on both WT samples but 1631 cm⁻¹ on BslA NaOH treated. A ~1630 cm⁻¹ component is visible as a shoulder on the WT bands too (the WT NaOH treated band is nearly flat between 1650 and 1630); on the BslA sample it becomes the dominant component, with ~1650 left as the shoulder. So the protein on that sample has a different (more β-sheet-like) band shape, not just more of it.
