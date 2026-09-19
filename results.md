# AVIRIS Cuprite: Spectral Unmixing & Mineral Identification

**Method:** VCA endmember extraction → FCLS abundance mapping → windowed Spectral Angle Distance (SAD) mineral ID against the USGS spectral library.

---

## Summary

This project unmixes the AVIRIS Cuprite benchmark scene into 12–14 spectral endmembers and identifies the physical minerals each one likely represents. A corrupted data band was found and removed early on (cutting reconstruction error by 2.5x), the endmember count was chosen by triangulating three independent methods, and the final mineral identification uses a windowed SAD approach that gives real confidence signals — not just a single angle — for each match. Several endmembers land on minerals well-documented at Cuprite (alunite, kaolinite, montmorillonite, buddingtonite, hematite), while a few remain genuinely low-confidence and are reported as such rather than forced into a match.

---

## 1. Data and Preprocessing

- **Scene:** AVIRIS Cuprite, NV (flown 1997) — the standard hyperspectral unmixing benchmark, covering a well-studied hydrothermal alteration site.
- **Shape:** 250 × 190 pixels, reshaped from a `(188, 47500)` `.mat` file.
- **Value range:** 707–8929 (mean 3511) — raw, not normalized reflectance. All error is reported as *relative* RMSE for this reason.
- **Band-0 fix:** one corrupted band was distorting an endmember and inflating error. Trimming it (188 → 186 bands) dropped reconstruction RMSE from **171.6 → 67.0**.
- **Instrument artifact:** synchronized dips appear in every endmember at ~bands 100 and 135 — this is the known seam between AVIRIS's internal spectrometer modules, a real instrument characteristic, not a defect.
- **Band bookkeeping:** the `.mat` file's `SlectBands` field records exactly which of the original 224 AVIRIS channels survive in the 188/186-band cube — this later gave an exact (not approximate) alignment to the USGS library.

## 2. Endmember Count and Unmixing Quality

Three methods triangulated to **p = 14**: an ambiguous SVD elbow, literature convention for this scene (12–16), and an HFC virtual-dimensionality estimate (13–16 depending on significance threshold).

FCLS abundance maps passed the sum-to-one check almost exactly (0.9999998–1.0000002) at both p=14 and p=12.

| | p = 14 | p = 12 |
|---|---|---|
| Reconstruction RMSE | 66.96 | 67.98 |
| **Relative RMSE** | **1.91%** | **1.94%** |

Dropping 2 endmembers costs only ~0.03 percentage points of accuracy — evidence that some of the p=14 solution is redundant (most likely shade/sunlit duplicates of the same material, a known VCA behavior). **The final mineral ID below uses p=12**, since the more parsimonious set is nearly as accurate and easier to interpret one-to-one.

## 3. Mineral Identification Method

Endmembers were compared to the USGS splib07 spectral library (AVIRIS-1997 convolution, matching the scene's flight year), searching `ChapterM_Minerals` and `ChapterS_SoilsAndMixtures`.

**Why windowed, not whole-spectrum:** an initial pass comparing full 186-band spectra produced suspiciously tight, low-variance angles and a couple of reference files winning implausibly often across unrelated endmembers — a sign that noisy/flat regions of the spectrum (and library files with many bad channels) were distorting the comparison. The final method instead restricts SAD to two physically diagnostic windows and adds three sanity checks:

- **SWIR window (2.00–2.45 μm, 45 bands)** — clay/alteration mineral absorption features (Al-OH, etc.)
- **VNIR window (0.40–1.30 μm, 97 bands)** — iron-oxide absorption features
- **Max band depth** — how strong the endmember's own absorption feature is in that window; near-zero means the spectrum is essentially flat there, and any "top match" is not meaningful (explicitly flagged as such below).
- **Library median angle** — the typical angle across the *whole* library, giving a baseline for "distance expected by chance" in that window.
- **Gap to #2** — how much better the top match is than the runner-up; a near-zero gap means the ID is ambiguous between two candidates, not wrong.
- **`[Cuprite-reported]`** tag — cross-references whether that specific mineral has been independently reported at Cuprite in the literature.

## 4. Results (p = 12)

| EM | SWIR top match (°) | VNIR top match (°) | Read |
|---|---|---|---|
| 1 | Chalcedony — 14.1° [C] | *featureless, unreliable* | Likely silica (chalcedony), SWIR only |
| 2 | Olivine — 32.9° (weak) | Siderite/Goethite — 20.7–21.4° [C] | Iron-bearing phase; VNIR more diagnostic |
| 3 | Kaolinite/Halloysite — 5.3° [C] | Hematite — 9.0° [C] | **Confident:** hematite-stained kaolinite clay |
| 4 | Alunite (K) — 6.6° [C] | weak, no tag | **Confident:** alunite (SWIR-diagnostic mineral) |
| 5 | Montmorillonite — 17.7° [C] | Montmorillonite — 11.3° [C] | **Confident, cross-validated** both windows |
| 6 | Buddingtonite — 13.6° [C] | weak, ambiguous | **Confident:** buddingtonite (SWIR-diagnostic) |
| 7 | Kaolinite/Alunite mix — 8.6° [C] | weak but same family — 24.5° [C] | **Confident:** kaolinite, cross-validated |
| 8 | Chalcedony — 14.4° [C] | Opalized Tuff — 22.6° [C] | Silica phase (opal/chalcedony), both windows agree |
| 9 | Montmorillonite — 16.0° [C] | Montmorillonite — 5.6° [C] | **Confident, cross-validated**, strong in VNIR |
| 10 | Alunite — 16.1° [C] | Goethite+Jarosite+Quartz mix — 11.8° [C] | Likely mixed alunite + iron-oxide alteration zone |
| 11 | Pyrite — 38.8° (weak) | Gibbsite — 25.1° (weak) | **Low confidence** — no strong feature either window |
| 12 | Carbon Black — 39.9° (featureless) | Gibbsite — 19.1° (weak) | **Low confidence** — possibly non-diagnostic/mixed endmember |

*[C] = independently reported at Cuprite in the literature — a strong cross-check that these aren't spurious library matches.*

**Takeaways:**
- 7 of 12 endmembers (3, 4, 5, 6, 7, 9, and to a lesser extent 8, 10) land on well-documented Cuprite alteration minerals, several confirmed in *both* spectral windows independently — alunite, kaolinite, montmorillonite, buddingtonite, hematite, and silica phases.
- Endmember 3 and 10 show *different* minerals dominating each window (e.g. kaolinite in SWIR, hematite in VNIR) — this reads as a genuine mixed-material pixel (hematite-stained clay, or alunite + iron-oxide zone) rather than a contradiction, since these mineral pairs commonly co-occur at Cuprite.
- Endmembers 11 and 12 have low band depth and high angles in both windows — these are honestly reported as unreliable rather than forced to a best-available guess.
- The earlier full-spectrum method's suspicious repeat-matches (one library file dominating many unrelated endmembers) do **not** recur here, supporting the diagnosis that the issue was a whole-spectrum/channel-dropout artifact, now resolved by windowing.

## 5. Remaining Limitations

- VNIR window discarded 449 of 1,485 library spectra for gaps/bad values (vs. 16 in SWIR) — VNIR results are on a smaller, potentially less representative library subset.
- Endmembers 11 and 12 remain unresolved; worth checking whether they correspond to shade/background pixels rather than a distinct material.
- Shade/brightness endmember redundancy (Section 2) has not been corrected at the pixel level (e.g. via L2 shade-normalization before VCA) — only indirectly evidenced via the p=12 ablation.

## 6. Key Numbers

| Item | Value |
|---|---|
| Cube shape (post-trim) | 250 × 190 × 186 |
| Final endmember count | 12 |
| Relative RMSE (p=12) | 1.94% |
| SWIR window | 2.00–2.45 μm, 45 bands |
| VNIR window | 0.40–1.30 μm, 97 bands |
| Confidently identified minerals | Alunite, Kaolinite, Montmorillonite, Buddingtonite, Hematite, Silica (opal/chalcedony) |
| Low-confidence endmembers | 11, 12 |
