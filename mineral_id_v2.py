"""
mineral_id_v2.py
================
Feature-based mineral ID for the 14 VCA endmembers.

Why v2: raw SAD on full-range spectra is dominated by overall brightness/slope,
so almost every bright library spectrum lands within a few degrees of every
endmember (hence "strong match" everywhere, halite/pectolite/etc. showing up).
v2 compares the *absorption features* instead:

  1. Restrict to a diagnostic window (default SWIR 2.00-2.45 um, where the
     Cuprite alteration minerals -- alunite, kaolinite, muscovite, buddingtonite,
     montmorillonite, calcite -- have their diagnostic Al-OH / NH4 / CO3 bands).
  2. Continuum-remove each spectrum (divide by its upper convex hull) and
     convert to band depth, 1 - CR.  Flat spectra -> ~0 everywhere.
  3. SAD between band-depth vectors.  This angle is ~0 only if the absorption
     positions AND relative depths match, so it's a far stricter test.
  4. Require full band coverage in the window (no NaN-masking on a different
     subset per library spectrum -- that made scores non-comparable in v1).
  5. Report, per endmember: top-N (one per mineral name), the gap to the next
     candidate, where #1 sits relative to the library median, and the
     endmember's own feature strength (tiny = featureless -> don't trust any ID).

Run:  python mineral_id_v2.py
"""

import os
import re
import glob
import numpy as np
import pandas as pd
import scipy.io

# ============================== CONFIG ======================================
USGS_ROOT = "ASCIIdata_splib07b_cvAVIRISc1997"
WAVELENGTH_FILE = os.path.join(
    USGS_ROOT, "s07_AV97_AVIRIS_1997_Wavelengths_(microns)_224ch.txt")
CHAPTER_DIRS = ["ChapterM_Minerals", "ChapterS_SoilsAndMixtures"]

ENDMEMBER_NPY = "E.npy"                       # (186, 14)
MAT_FILE = "Cuprite/Cuprite_data_R188.mat"
SELECT_BANDS_KEY = "SlectBands"
N_FRONT_TRIMMED = 2

# Diagnostic windows to evaluate. Each is run separately and saved to its own
# CSV. SWIR is the main one for Cuprite; VNIR (Fe3+ features) helps separate
# jarosite / hematite / goethite.
WINDOWS = {
    "SWIR": (2.00, 2.45),
    "VNIR": (0.40, 1.30),
}

TOP_N_MATCHES = 3
ONE_PER_MINERAL = True    # collapse several samples of the same mineral
MIN_FEATURE_DEPTH = 0.02  # endmember max band depth below this -> "featureless"

# Minerals commonly mapped at Cuprite (Swayze/Clark Tetracorder maps etc.).
# Used only as a soft flag in the output, never to filter.
CUPRITE_REPORTED = [
    "alunite", "kaolin", "dickite", "halloysite", "muscov", "illite",
    "montmoril", "buddingtonite", "chalcedony", "opal", "calcite",
    "jarosite", "hematite", "goethite", "nontronite", "pyrophyl",
    "sphene", "andradite",
]
# ============================================================================


def load_usgs_ascii(filepath):
    with open(filepath, "r") as f:
        lines = f.readlines()
    v = np.array([float(l.strip()) for l in lines[1:] if l.strip()])
    v[v < -1e30] = np.nan
    return v


def load_library():
    spectra = {}
    for chapter in CHAPTER_DIRS:
        path = os.path.join(USGS_ROOT, chapter)
        if not os.path.isdir(path):
            print(f"  (skipping missing {path})")
            continue
        files = sorted(glob.glob(os.path.join(path, "*.txt")))
        for fp in files:
            name = os.path.splitext(os.path.basename(fp))[0]
            try:
                spectra[name] = (load_usgs_ascii(fp), chapter)
            except Exception as e:
                print(f"  WARNING: {fp}: {e}")
        print(f"  Loaded {len(files)} spectra from {chapter}")
    return spectra


def align_idx_from_mat(n_bands, n_usgs):
    mat = scipy.io.loadmat(MAT_FILE)
    sb = mat[SELECT_BANDS_KEY].flatten().astype(int) - 1
    if sb.max() >= n_usgs:
        raise ValueError("SlectBands index outside USGS grid")
    idx = sb[N_FRONT_TRIMMED:]
    if len(idx) != n_bands:
        raise ValueError(f"alignment gives {len(idx)} bands, E has {n_bands}")
    return idx


def upper_hull(x, y):
    """Indices of the upper convex hull of points (x, y), x increasing."""
    hull = []
    for i in range(len(x)):
        while len(hull) >= 2:
            i1, i2 = hull[-2], hull[-1]
            cross = (x[i2] - x[i1]) * (y[i] - y[i1]) - (y[i2] - y[i1]) * (x[i] - x[i1])
            if cross >= 0:        # i2 lies on/below the chord -> drop it
                hull.pop()
            else:
                break
        hull.append(i)
    return np.array(hull)


def band_depth(wl, r):
    """1 - continuum-removed reflectance. Returns None if r has NaNs/<=0."""
    if np.any(~np.isfinite(r)) or np.any(r <= 0):
        return None
    h = upper_hull(wl, r)
    cont = np.interp(wl, wl[h], r[h])
    return 1.0 - r / cont


def angle_deg(a, b):
    na, nb = np.linalg.norm(a), np.linalg.norm(b)
    if na == 0 or nb == 0:
        return np.nan
    return np.degrees(np.arccos(np.clip(a @ b / (na * nb), -1, 1)))


def mineral_key(name):
    # s07_AV97_<Mineral...>_<sample>_... -> first token after prefix, lowercased,
    # stripped of trailing punctuation/fractions (e.g. "Alunite.5+Musc" -> "alunite")
    tok = name.split("_")[2] if name.count("_") >= 2 else name
    return re.split(r"[.+\-(]", tok)[0].lower()


def is_cuprite_reported(name):
    n = name.lower()
    return any(k in n for k in CUPRITE_REPORTED)


def run_window(label, lo, hi, E, lib, wl_aligned, align_idx):
    in_win = (wl_aligned >= lo) & (wl_aligned <= hi)
    # drop bands where the endmembers themselves are NaN
    in_win &= np.all(np.isfinite(E), axis=1)
    wl = wl_aligned[in_win]
    print(f"\n=== {label} window {lo:.2f}-{hi:.2f} um: {in_win.sum()} bands ===")
    if in_win.sum() < 8:
        print("  too few bands, skipping")
        return []

    # library band-depth vectors (only spectra with FULL coverage of the window)
    lib_bd, skipped = {}, 0
    for name, (spec, chap) in lib.items():
        bd = band_depth(wl, spec[align_idx][in_win])
        if bd is None:
            skipped += 1
            continue
        lib_bd[name] = (bd, chap)
    print(f"  library spectra usable: {len(lib_bd)} (skipped {skipped} "
          f"with gaps/non-positive values in window)")

    rows = []
    for i in range(E.shape[1]):
        bd_e = band_depth(wl, E[in_win, i])
        strength = np.nanmax(bd_e) if bd_e is not None else np.nan
        scores = sorted(
            ((angle_deg(bd_e, bd), name, chap) for name, (bd, chap) in lib_bd.items()),
            key=lambda t: (np.isnan(t[0]), t[0]))
        all_angles = np.array([s[0] for s in scores if np.isfinite(s[0])])
        median = np.median(all_angles)

        picked, seen = [], set()
        for s in scores:
            k = mineral_key(s[1])
            if ONE_PER_MINERAL and k in seen:
                continue
            seen.add(k)
            picked.append(s)
            if len(picked) == TOP_N_MATCHES + 1:
                break
        gap = picked[1][0] - picked[0][0] if len(picked) > 1 else np.nan

        flag = "FEATURELESS - ID unreliable" if strength < MIN_FEATURE_DEPTH else ""
        print(f"\nEndmember {i+1}: max band depth={strength:.3f}  "
              f"library median angle={median:.1f} deg  gap #1->#2={gap:.1f} deg  {flag}")
        for r, (a, name, chap) in enumerate(picked[:TOP_N_MATCHES], 1):
            cup = " [Cuprite-reported]" if is_cuprite_reported(name) else ""
            print(f"  #{r}: {a:5.1f} deg  {name} ({chap}){cup}")
            rows.append(dict(window=label, endmember=i + 1, rank=r, mineral=name,
                             chapter=chap, angle_deg=a, gap_to_next=gap if r == 1 else np.nan,
                             library_median_deg=median, endmember_max_depth=strength,
                             cuprite_reported=bool(cup), flag=flag))
    return rows


def main():
    E = np.load(ENDMEMBER_NPY).astype(float)
    n_bands, p = E.shape
    print(f"E: {E.shape}")
    wl_usgs = load_usgs_ascii(WAVELENGTH_FILE)
    align_idx = align_idx_from_mat(n_bands, len(wl_usgs))
    wl_aligned = wl_usgs[align_idx]
    print(f"Aligned wavelengths: {wl_aligned[0]:.3f}-{wl_aligned[-1]:.3f} um")
    lib = load_library()

    all_rows = []
    for label, (lo, hi) in WINDOWS.items():
        all_rows += run_window(label, lo, hi, E, lib, wl_aligned, align_idx)
    df = pd.DataFrame(all_rows)
    df.to_csv("endmember_mineral_matches_v2.csv", index=False)
    print("\nSaved endmember_mineral_matches_v2.csv")


if __name__ == "__main__":
    main()