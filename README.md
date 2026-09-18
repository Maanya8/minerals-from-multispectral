# Cuprite Hyperspectral Mineral Mapping

A CPU-friendly hyperspectral unmixing workflow for the AVIRIS Cuprite scene. The project extracts spectral endmembers with Vertex Component Analysis (VCA), estimates per-pixel material abundances with Fully Constrained Least Squares (FCLS), and compares extracted spectra with mineral spectra from the USGS Spectral Library.

## What the pipeline does

1. Loads the Cuprite AVIRIS hyperspectral cube.
2. Estimates the number of spectral endmembers with Harsanyi-Farrand-Chang (HFC) virtual dimensionality.
3. Extracts representative endmember spectra with VCA.
4. Estimates abundance maps under non-negativity and sum-to-one constraints using FCLS.
5. Measures reconstruction error for the linear mixing model.
6. Matches endmembers to USGS mineral and soil spectra using continuum-removed absorption features in the SWIR and VNIR regions.

The linear mixing model is:

```text
Y = E A + N
```

where `Y` is the observed hyperspectral data, `E` contains endmember spectra, `A` contains abundance fractions, and `N` is residual noise.

## Project layout

```text
Cuprite/                                  Cuprite AVIRIS .mat data and metadata
ASCIIdata_splib07b_cvAVIRISc1997/         USGS spectral library ASCII spectra
file_open.py                              Cuprite data loader
vca_fcls.py                               HFC, VCA, FCLS, and reconstruction pipeline
mineral_id.py                             Initial spectral-library matching workflow
mineral_id_v2.py                          Feature-based mineral matching workflow
mineral_graphs.py                         Abundance and dominant-mineral visualizations
E.npy                                     Extracted endmember spectra
abundance_maps_p12.npy                    FCLS abundance maps
visualizations/                           Generated plots
results.md                                Current experiment results
steps.md                                  Project background and planned extensions
```

## Requirements

- Python 3.10 or newer
- NumPy
- SciPy
- pandas
- Matplotlib
- `pysptools`

Install the dependencies in a virtual environment:

```powershell
py -m venv .venv
.\.venv\Scripts\Activate.ps1
python -m pip install numpy scipy pandas matplotlib pysptools
```

The dataset folders are expected to remain at the paths shown in the project layout. The scripts use relative paths, so run commands from the repository root.

## Running the pipeline

Extract endmembers, calculate abundance maps, save `E.npy`, and generate the VCA plot:

```powershell
python vca_fcls.py
```

Run feature-based mineral identification against the USGS library:

```powershell
python mineral_id_v2.py
```

The matching script evaluates SWIR and VNIR diagnostic windows and writes ranked matches to `endmember_mineral_matches_v2_p12.csv`.

Generate spatial abundance and dominant-mineral visualizations:

```powershell
python mineral_graphs.py
```

## Outputs

Typical outputs include:

- `E.npy`: extracted endmember spectra.
- `abundance_maps_p12.npy` or `abundance_maps_p14.npy`: abundance fractions for each endmember at each image pixel.
- `endmember_mineral_matches_v2_p12.csv`: ranked USGS library matches and feature diagnostics.
- `visualizations/vca_extracted_endmembers_p12.png`: extracted spectra.
- `abundance_maps_grid_p12.png`: one abundance map per endmember.
- `dominant_mineral_map_p12.png`: the highest-abundance endmember at each pixel.

Abundance fractions should be non-negative and sum to approximately 1 for each pixel. The reconstruction RMSE printed by `vca_fcls.py` is a useful baseline for how well the linear model explains the scene.

## Current results

The current experiment estimated 14 endmembers with HFC and evaluated a 12-endmember VCA/FCLS result. The saved 12-endmember abundance maps have shape `(250, 190, 12)`, with abundance sums approximately equal to 1. The reported reconstruction RMSE is approximately `67.98` for the 12-endmember run. See [`results.md`](results.md) for the captured output and comparison with the 14-endmember run.

## Limitations and next steps

- VCA assumes that sufficiently pure pixels exist in the scene; heavily mixed pixels may require a minimum-volume method such as MVSA or MVES.
- The linear mixing model does not represent intimate, nonlinear mixtures well.
- Mineral names are spectral-library matches, not definitive geological labels. Featureless endmembers and close spectral matches should be interpreted cautiously.
- Quantitative abundance accuracy requires benchmark or synthetic data with known ground truth.
- A useful extension is to compare the linear baseline with a bilinear model or an autoencoder-based unmixing method.

## Data attribution

The project uses the Cuprite AVIRIS benchmark scene and spectra from the USGS Spectral Library. Keep the included dataset documentation and original attribution with the data when redistributing the project.
