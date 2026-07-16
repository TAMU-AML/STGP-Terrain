# STGP-Terrain-Aware-Power-Curve

Reproducibility report for the paper "Spatio-Temporal Gaussian Process for Building Terrain-Incorporating Wind Power Curves" by Ahmadreza Chokhachian, V. Roshan Joseph, and Yu Ding.

This repository is also available via GitHub.

## Datasets

The experiments use wind turbine SCADA datasets collected from a utility-scale wind farm, containing measurements from **66 turbines** recorded at **10-minute intervals**.

### Turbine Dataset

SCADA measurements from **2017 and 2018** for all turbines. Each CSV file corresponds to a single turbine and includes measurements such as:

- wind_speed
- temperature
- wind_direction
- turbulence_intensity
- power

Each turbine file contains approximately **40,000–50,000 observations**.

### Terrain Data

A separate CSV file contains terrain features for all **66 turbines**. The terrain variables include:

- slope
- rix
- ridge

These terrain features are used as spatial inputs in the proposed **spatio-temporal Gaussian process (STGP) model**.

### data/processed_data/

Contains the outputs of Stage 1 and Stage 2 of the STGP pipeline (support points and thinned TwinGP outputs). The `support` package used in Stage 1 has been archived from CRAN and normally needs to be installed locally. To avoid installation issues during reproduction, we provide the outputs of these stages directly in this folder.

## Repository Structure

```
STGP-Terrain-Aware-Power-Curve/
│
├──requirements.txt/
│
├──requirements_r.R/
│
├── data/
│   ├── Turbine_i_2017.csv
│   ├── Turbine_i_2018.csv
│   ├── terrain_features.csv
│   ├── turbine_locations.csv
│   └── processed_data/
│
├── run_table23.py          # cross-platform runner (Tables 2–3)
├── run_table4.py           # cross-platform runner (Table 4)
├── run_table23.sh          # Linux / macOS convenience wrapper
├── run_table4.sh           # Linux / macOS convenience wrapper
├── run_table23.bat         # Windows convenience wrapper
├── run_table4.bat          # Windows convenience wrapper
│
├── code/
│   ├── supplement
│   ├── Figures
│      ├── fig_4.R
│      ├── fig_5.R
│      ├── fig_7.R
│   ├── update_final_results.R
│   ├── update_final_results.py
│   ├── Table2-Table3(STGP).R
│   ├── Table2-Table3(twinGP+Binning).R
│   ├── Table2-Table3(NN).ipynb
│   ├── Table2-Table3(BNN).ipynb
│   ├── Table2-Table3(XGBoost).ipynb
│   ├── Table4(STGP).R
│   ├── Table4(twinGP+Binning).R
│   ├── Table4(NN).ipynb
│   ├── Table4(BNN).ipynb
│   ├── Table4(Binning-hetGP).R
│   └── Table4(XGBoost).ipynb
│
├── results/
│   ├── final_results.csv
│   └── intermediate/
│
├── LICENSE
└── README.md
```

**results/final_results.csv**
Stores the aggregated results of tables presented in the main paper (Tables 2, 3, 4). Running each method updates the corresponding row in this file.

**results/intermediate/**
Runtime logs, turbine-level prediction errors, and other non-aggregated outputs from the main-paper methods.

**results/figures/**
Figures produced by `code/Figures/` (`fig_4.R`, `fig_5.R`, `fig_7.R`).

**results/supplement/**
Tables produced by `code/supplement/` (S1, S4, S5) for the Supplemental Material.

## Code

The implementation of the proposed method and the benchmark methods is written in **R** and **Python**.

### R implementations
- STGP (proposed method)
- Binning
- hetGP
- TwinGP

### Python implementations
- Multi-layer neural network
- Bayesian neural network
- XGBoost

## Dependencies

Required R packages: `dplyr`, `data.table`, `twingp`, `readr`, `tidyr`, `hetGP`, `ALEPlot`

Required Python packages: `numpy`, `pandas`, `torch`, `torchbnn`, `tensorflow`, `lightgbm`, `xgboost`, `scikit-learn`, `jupyter`, `nbconvert`

The following software must be installed and available on your system `PATH` before running the reproduction scripts.

| Software | Version tested | Notes |
|----------|-----------------|-------|
| R | ≥ 4.0 | `Rscript` must be on PATH |
| Python | ≥ 3.9 | `python` or `python3` must be on PATH |
| Jupyter | any recent | install via `pip install notebook nbconvert` |

## Note on Runtime Performance (OpenBLAS)

We observed that the runtime of several methods, including **STGP**, is significantly faster when **OpenBLAS** is used for linear algebra operations. If OpenBLAS is not configured in your R installation, we recommend installing it following the instructions at:

https://github.com/david-cortes/R-openblas-in-windows

Without OpenBLAS, the runtime of some methods may be **up to two times slower**.

## Reproducibility Workflow

1. **Download the package from GitHub.**

   Download the repository as a ZIP file and extract it to any directory, or clone it with Git:

```bash
   git clone https://github.com/TAMU-AML/STGP-Terrain.git
```

2. **Install dependencies.**

```bash
   # R
   Rscript requirements_r.R

   # Python
   pip install -r requirements.txt
```

3. **Open a terminal in the project folder.**

   - Windows: open the folder in File Explorer, click the address bar, type `cmd`, press Enter (or right-click inside the folder and select "Open in Terminal").
   - macOS / Linux: open Terminal and `cd` into the extracted folder.

4. **Reproduce Tables 2 and 3.**

   - Windows: double-click `run_table23.bat`, or run `python run_table23.py`
   - macOS / Linux: run `bash run_table23.sh`, or `python3 run_table23.py`

5. **Reproduce Table 4.**

   - Windows: double-click `run_table4.bat`, or run `python run_table4.py`
   - macOS / Linux: run `bash run_table4.sh`, or `python3 run_table4.py`

6. **Reproduce Figures 4, 5, and 7.**

```bash
   Rscript code/Figures/fig_4.R
   Rscript code/Figures/fig_5.R
   Rscript code/Figures/fig_7.R
```

7. **Reproduce Supplemental Table S1.**

```bash
   python code/supplement/S1.py
   Rscript code/supplement/S1.R
```

8. **Reproduce Supplemental Table S4.**

```bash
   Rscript code/supplement/S4.R
```

9. **Reproduce Supplemental Table S5.**

```bash
   Rscript code/supplement/S5.R
```

## Reproduction Reference Table

| Which results to reproduce | Code file(s) | Output | Run time |
|---|---|---|---|
| Table 2 | `code/Table2-Table3(XGBoost).ipynb`, `code/Table2-Table3(NN).ipynb`, `code/Table2-Table3(BNN).ipynb`, `code/Table2-Table3(twinGP+Binning).R`, `code/Table2-Table3(STGP).R` (run together via `run_table23.py` / `.sh` / `.bat`) | `results/final_results.csv` (Table 2 rows), `results/intermediate/` | ~18 hours total, measured (produced together with Table 3) |
| Table 3 | Same scripts as Table 2 above (2017 and 2018 predictions are produced in the same run) | `results/final_results.csv` (Table 3 rows), `results/intermediate/` | Included in the runtime for Table 2 (no separate run required) |
| Table 4 | `code/Table4(XGBoost).ipynb`, `code/Table4(NN).ipynb`, `code/Table4(BNN).ipynb`, `code/Table4(twinGP+Binning).R`, `code/Table4(Binning-hetGP).R`, `code/Table4(STGP).R` (run together via `run_table4.py` / `.sh` / `.bat`) | `results/final_results.csv` (Table 4 rows), `results/intermediate/` | ~3 hours |
| Figure 4 | `code/Figures/fig_4.R` | `results/figures/fig4_lk_ok_rmse_boxplot.pdf`, `results/figures/fig4_turbine60_lk_ok.pdf` | ~5 minutes |
| Figure 5 | `code/Figures/fig_5.R` | `results/figures/fig5_temporal_ale.pdf`, `results/figures/fig5_spatial_ale.pdf`, `results/figures/fig5_ale_ranges.csv` | ~1 hour |
| Figure 7 | `code/Figures/fig_7.R` | `results/figures/fig7_underperformance.pdf`, `results/figures/fig7_underperformance_data.csv` | ~5 hours |
| Table S1 | `code/supplement/S1.py`, `code/supplement/S1.R` | `results/supplement/table_s1.csv` | ~41 hours |
| Table S4 | `code/supplement/S4.R` | `results/supplement/table_s4.csv` | ~28 hours |
| Table S5 | `code/supplement/S5.R` | `results/supplement/table_s5.csv` | ~5 hours |

All outputs land in `results/final_results.csv` (Tables 2–4), `results/figures/` (Figures 4, 5, 7), and `results/supplement/` (Tables S1, S4, S5).

### Note on BHM Implementation

We do not provide code for the **BHM (Bayesian Hierarchical Model)** method in this repository. The same leave-one-out experiment is already implemented in:

https://github.com/TAMU-AML/BHM-Terrain-Paper

In addition, the BHM method is significantly slower than the other methods considered in this paper. Including it in the reproduction scripts would substantially increase the time required to reproduce the results.
