"""
mineral_id.py
=============
Identify which USGS reference minerals your 14 VCA endmembers correspond to,
via Spectral Angle Distance (SAD).

Pipeline:
  1. Load your endmember matrix E (bands, p) -- from VCA, on the trimmed
     186-band grid.
  2. Load the USGS splib07b AVIRIS-convolved wavelength file (224 channels).
  3. Load every mineral spectrum from ChapterM_Minerals (and optionally
     ChapterS_SoilsAndMixtures as a fallback).
  4. Align your 186 bands to the USGS 224-channel grid EXACTLY, using the
     'SlectBands' field already present in your original cuprite.mat (the
     1-indexed original AVIRIS band numbers that were kept when that .mat
     file was built) -- no wavelength matching or literature-range
     guessing needed.
  5. Compute SAD between every endmember and every (aligned) mineral
     spectrum, report the top-3 matches per endmember, and save a CSV.

Fill in the CONFIG block below, then run:
    python mineral_id.py
"""

import os
import re
import glob
import numpy as np
import pandas as pd
import scipy.io

# ============================== CONFIG ======================================

# Path to the USGS ASCII library folder you downloaded and unzipped, e.g.
# ASCIIdata_splib07b_cvAVIRISc1997/  (recommended -- closest year to the
# 1997 Cuprite flight) or ASCIIdata_splib07b_cvAVIRISc2014/
USGS_ROOT = "ASCIIdata_splib07b_cvAVIRISc1997"

# Wavelength file lives at the top level of that folder. The exact filename
# varies by download; this glob finds it automatically -- if it matches more
# than one file, set WAVELENGTH_FILE explicitly instead.
WAVELENGTH_FILE = "ASCIIdata_splib07b_cvAVIRISc1997/s07_AV97_AVIRIS_1997_Wavelengths_(microns)_224ch.txt" # e.g. "ASCIIdata_splib07b_cvAVIRISc1997/splib07a_Wavelengths_AVIRIS1997_0.37-2.5_microns.txt"

# Mineral chapter folder(s) to search, in priority order. SAD results will
# note which chapter each match came from.
CHAPTER_DIRS = ["ChapterM_Minerals", "ChapterS_SoilsAndMixtures"]

# Your VCA endmember matrix, shape (186, 14) -- bands x endmembers, on the
# TRIMMED band set (first 2 bands already removed, per your pipeline).
# Point this at wherever you saved E, e.g. np.save('E.npy', E) after VCA.
ENDMEMBER_NPY = "E.npy"

# --- Band alignment ---
# Your original cuprite.mat has a 'SlectBands' field: a (188,1) uint8 array
# of the ORIGINAL 1-indexed AVIRIS band numbers (1-224) that were kept when
# this .mat file was assembled. This gives EXACT band alignment to the USGS
# 224-channel grid -- no wavelength matching or literature-range guessing
# needed. Point MAT_FILE at your original file.
MAT_FILE = "Cuprite/Cuprite_data_R188.mat"
SELECT_BANDS_KEY = "SlectBands"  # key name

# How many bands were trimmed off the FRONT of the 188 by your own pipeline
# (the band-0 corrupted-outlier fix: cube_3d_trimmed = cube_3d[:, :, 2:]).
N_FRONT_TRIMMED = 2

TOP_N_MATCHES = 3
OUTPUT_CSV = "endmember_mineral_matches.csv"

# ============================================================================


def load_usgs_ascii(filepath):
    """Load a single USGS ASCII spectrum or wavelength file.
    First line is a title/metadata line and is skipped. Bad/deleted
    channels are marked -1.23e34 and converted to NaN."""
    with open(filepath, "r") as f:
        lines = f.readlines()
    values = np.array([float(line.strip()) for line in lines[1:] if line.strip()])
    values[values < -1e30] = np.nan
    return values



def load_mineral_library():
    """Load every mineral spectrum from CHAPTER_DIRS.
    Returns dict: name -> (spectrum array, source_chapter)."""
    spectra = {}
    for chapter in CHAPTER_DIRS:
        chapter_path = os.path.join(USGS_ROOT, chapter)
        if not os.path.isdir(chapter_path):
            print(f"  (skipping missing chapter folder: {chapter_path})")
            continue
        txt_files = sorted(glob.glob(os.path.join(chapter_path, "*.txt")))
        for fpath in txt_files:
            name = os.path.splitext(os.path.basename(fpath))[0]
            try:
                spec = load_usgs_ascii(fpath)
            except Exception as e:
                print(f"  WARNING: failed to load {fpath}: {e}")
                continue
            spectra[name] = (spec, chapter)
        print(f"  Loaded {len(txt_files)} spectra from {chapter}")
    return spectra


def dedupe_repeat_measurements(mineral_spectra):
    """Collapse obvious repeat-measurement variants (filenames ending in
    'odd1', 'odd2', etc.) by keeping only the first occurrence per base name,
    so top-N results aren't dominated by the same physical sample measured
    twice. This is a heuristic, not authoritative -- inspect results."""
    seen_bases = {}
    deduped = {}
    for name, (spec, chapter) in mineral_spectra.items():
        base = re.sub(r"_odd\d+$", "", name, flags=re.IGNORECASE)
        if base in seen_bases:
            continue
        seen_bases[base] = name
        deduped[name] = (spec, chapter)
    n_removed = len(mineral_spectra) - len(deduped)
    if n_removed:
        print(f"  Deduped {n_removed} apparent repeat-measurement files "
              f"(kept one per base sample name).")
    return deduped


def align_by_select_bands(n_your_bands, n_usgs_bands=224):
    """Exact alignment using the 'SlectBands' field from your original
    cuprite.mat. That field is a (188,1) array of the 1-indexed original
    AVIRIS band numbers (1-224) that were kept when the .mat file was built.
    Converting to 0-indexed and dropping your own N_FRONT_TRIMMED front
    bands gives the exact position of each of your working bands within
    the USGS library's 224-channel grid -- no matching/guessing required."""
    mat = scipy.io.loadmat(MAT_FILE)
    if SELECT_BANDS_KEY not in mat:
        raise KeyError(
            f"'{SELECT_BANDS_KEY}' not found in {MAT_FILE}. Available keys: "
            f"{[k for k in mat.keys() if not k.startswith('__')]}"
        )
    select_bands_1indexed = mat[SELECT_BANDS_KEY].flatten().astype(int)
    print(f"  SlectBands: {len(select_bands_1indexed)} entries, "
          f"range {select_bands_1indexed.min()}-{select_bands_1indexed.max()} "
          f"(1-indexed, out of {n_usgs_bands})")

    select_bands_0indexed = select_bands_1indexed - 1
    if select_bands_0indexed.max() >= n_usgs_bands:
        raise ValueError(
            f"SlectBands contains an index ({select_bands_0indexed.max() + 1}) "
            f"outside the 1-{n_usgs_bands} range -- double-check n_usgs_bands "
            f"matches your USGS wavelength file's channel count."
        )

    align_idx = select_bands_0indexed[N_FRONT_TRIMMED:]
    print(f"  After dropping your {N_FRONT_TRIMMED} front-trimmed bands: "
          f"{len(align_idx)} indices remain.")

    if len(align_idx) != n_your_bands:
        raise ValueError(
            f"Alignment produced {len(align_idx)} indices but your endmember "
            f"matrix has {n_your_bands} bands. Check N_FRONT_TRIMMED and "
            f"confirm SlectBands corresponds to the same original 188-band "
            f"cube your VCA endmembers were extracted from."
        )
    return align_idx


def sad(a, b):
    """Spectral Angle Distance in degrees, ignoring NaNs in either spectrum."""
    mask = ~(np.isnan(a) | np.isnan(b))
    if mask.sum() < 10:
        return np.nan  # too few valid overlapping channels to trust
    a, b = a[mask], b[mask]
    cos_angle = np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b))
    cos_angle = np.clip(cos_angle, -1, 1)
    return np.degrees(np.arccos(cos_angle))


def interpret(sad_deg):
    if np.isnan(sad_deg):
        return "no valid overlap"
    if sad_deg < 10:
        return "strong match"
    if sad_deg < 20:
        return "plausible, uncertain"
    return "weak match"


def main():
    print("Loading endmembers...")
    E = np.load(ENDMEMBER_NPY)  # (bands, p)
    n_bands, p = E.shape
    print(f"  E shape: {E.shape} ({p} endmembers, {n_bands} bands)")

    print("Loading USGS wavelength grid...")
    wl_path = WAVELENGTH_FILE
    usgs_wavelengths = load_usgs_ascii(wl_path)
    print(f"  USGS grid: {len(usgs_wavelengths)} channels, "
          f"{usgs_wavelengths[0]:.3f}-{usgs_wavelengths[-1]:.3f} microns "
          f"(from {wl_path})")

    print("Loading mineral library...")
    mineral_spectra = load_mineral_library()
    print(f"  Total spectra loaded: {len(mineral_spectra)}")
    mineral_spectra = dedupe_repeat_measurements(mineral_spectra)

    print("Aligning bands (using SlectBands from cuprite.mat -- exact)...")
    align_idx = align_by_select_bands(n_bands, n_usgs_bands=len(usgs_wavelengths))
    alignment_mode = "select_bands-exact"

    # Optional sanity check: print the wavelengths at your two known
    # instrument-seam locations (~band 100, ~135 in your ORIGINAL 188-band
    # indexing, before your 2-band front trim) to confirm the alignment
    # lands somewhere physically sensible (AVIRIS module boundaries).
    for seam_band_188 in (100, 135):
        seam_band_186 = seam_band_188 - N_FRONT_TRIMMED
        if 0 <= seam_band_186 < n_bands:
            wl = usgs_wavelengths[align_idx[seam_band_186]]
            print(f"  Sanity check: your original band {seam_band_188} "
                  f"aligns to {wl:.3f} microns in the USGS grid.")

    print(f"Computing SAD ({alignment_mode} alignment)...")
    results_rows = []
    for i in range(p):
        endmember_spectrum = E[:, i]
        scores = {}
        for name, (spec, chapter) in mineral_spectra.items():
            aligned_ref = spec[align_idx]
            scores[name] = (sad(endmember_spectrum, aligned_ref), chapter)

        ranked = sorted(scores.items(), key=lambda kv: (np.isnan(kv[1][0]), kv[1][0]))
        top_matches = ranked[:TOP_N_MATCHES]

        print(f"\nEndmember {i+1}:")
        for rank, (name, (score, chapter)) in enumerate(top_matches, start=1):
            label = interpret(score)
            print(f"  #{rank}: {name} ({chapter}) -- SAD={score:.2f} deg [{label}]")
            results_rows.append({
                "endmember": i + 1,
                "rank": rank,
                "mineral": name,
                "chapter": chapter,
                "sad_degrees": score,
                "interpretation": label,
                "alignment_mode": alignment_mode,
            })

    df = pd.DataFrame(results_rows)
    df.to_csv(OUTPUT_CSV, index=False)
    print(f"\nSaved full results to {OUTPUT_CSV}")


if __name__ == "__main__":
    main()