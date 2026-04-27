[![Paper DOI](https://img.shields.io/badge/Paper-DOI-blue)](https://doi.org/10.XXXX/XXXXX)
[![Data DOI](https://img.shields.io/badge/Data-ESS--DIVE-green)](https://doi.org/10.XXXX/XXXXX)
# Imnavait-45-years

**Data and model code for:**

> Mukherjee, N., Gao, B., Coon, E.T., Shuai, P., Hill, D., Neilson, B.T., Kling, G.W., Chen, J., and Cardenas, M.B. — *Seasonal to decadal changes in supra-permafrost soil and aquifer thermal hydrology: broader, deeper, and longer activity.* Submitted to *AGU Advances*.

A field-calibrated, coupled surface–subsurface thermal-hydrological model (ATS v1.6-dev) driven by 44 years of climate forcing (1981–2024) along a hillslope-to-riparian transect at Imnavait Creek, North Slope of Alaska.

---

## Repository Structure

```
.
├── data/
│   ├── forcing_daac_ABoVE_snowmodel_1981_2020.h5        # NASA ABoVE SnowModel met forcing, 1981–2020
│   ├── forcing_daac_ABoVE_snowmodel_2019_2024.h5        # Extended met forcing, 2019–2024
│   ├── forcing_smoothed4spinUp.h5                        # Smoothed forcing for steady-state spin-up
│   ├── column_data_freezeup.h5                           # 1D column freeze-up output
│   ├── column_data_spinup_smoothedMet.h5                 # 1D column spin-up output
│   ├── mesh_column_imnavait.exo                          # 1D column mesh (Exodus II)
│   ├── mesh_transect_imnavait.exo                        # 2D transect mesh (Exodus II)
│   ├── IMP_dtw_dti.csv                                   # Piezometer depth-to-water and depth-to-ice table
│   ├── IMP-DTW-time-series.csv                           # Water table time series (IMP-01 to IMP-19)
│   ├── SM_Data_daily_Imnavait.csv                        # SNOTEL daily snow/met data (Site ID: 968)
│   ├── thaw_depth_vis_1981_2024.csv                      # Modeled thaw depth, 1981–2024
│   ├── max_summer_enthalpy_1981_2024.csv                 # Annual maximum summer subsurface enthalpy
│   └── TRods_2023-2024/                                  # Thermistor array data (T.rod.X®)
│       ├── IMP-05_2023.csv / IMP-05_2024.csv
│       ├── IMP-09_2023.csv / IMP-09_2024.csv
│       ├── IMP-14_2023.csv / IMP-14_2024.csv
│       └── IMP-18_2023.csv / IMP-18_2024.csv
│
├── inputs/                                               # ATS XML input files
│   ├── column_freezeup_imnavait.xml                      # Step 1: 1D freeze-up
│   ├── column_spinup_imnavait_smoothedMet.xml            # Step 2: 1D steady-state spin-up
│   ├── column_spinup_imnavait_realMet_daac_ABoVE_snowmodel_1981_2020.xml   # Step 3: calibration, 1981–2020
│   ├── column_spinup_imnavait_realMet_daac_ABoVE_snowmodel_2019_2024.xml   # Step 3: calibration, 2019–2024
│   ├── transect_run_imnavait_realMet_daac_ABoVE_snowmodel_1981_2020.xml    # Step 4: 2D run, 1981–2020
│   └── transect_run_imnavait_realMet_daac_ABoVE_snowmodel_2019_2024.xml    # Step 4: 2D run, 2019–2024
│
├── runs/                                                 # Checkpoints and observation output
│   ├── run01_column_freezeup/
│   ├── run02_column_spinup_smoothedmet/
│   ├── run03_column_transient_calibration_SnowExp_0p24_swe_1p0_2011-2020/
│   ├── run03_column_transient_calibration_SnowExp_0p26_swe_1p0_2011-2020/  # ← selected calibration
│   ├── run03_column_transient_calibration_SnowExp_0p26_swe*0p2_2019-2024/
│   ├── run03_column_transient_calibration_SnowExp_0p28_swe_1p2_2011-2020/
│   ├── run03_column_transient_calibration_SnowExp_0p30_swe_1p5_2011-2020/
│   ├── run04_transect_transient_1981-2020/
│   └── run_04_transect_transient_2019-2024/
│
├── scripts/
│   ├── IC45y_Figures.ipynb        # Reproduces all manuscript figures
│   ├── gen_column_imnavait.py     # Generates 1D column mesh
│   ├── gen_imnavait_transect.py   # Generates 2D transect mesh
│   ├── column_data.py             # Processes 1D column output
│   ├── smooth_met.py              # Generates smoothed met forcing for spin-up
│   └── utils/
│       ├── analysis.py            # Theil-Sen slope, Mann-Kendall, ZCP diagnostics
│       ├── data_io.py             # HDF5/CSV I/O helpers
│       ├── mesh_utils.py          # Mesh generation utilities
│       └── plotting.py            # Figure styling helpers
│
└── flmd.csv                       # File-level metadata
```

---

## Model Workflow

The simulation follows four sequential steps:

| Step | Input XML | Description |
|------|-----------|-------------|
| 1 | `column_freezeup_imnavait.xml` | 1D column: establish frozen initial state |
| 2 | `column_spinup_imnavait_smoothedMet.xml` | 1D column: 100-year steady-state spin-up with smoothed forcing |
| 3 | `column_spinup_imnavait_realMet_*.xml` | 1D column: transient calibration against SNOTEL (snow compaction exponent sweep) |
| 4 | `transect_run_imnavait_realMet_*.xml` | Full 2D hillslope-riparian transient run, 1981–2024 |

The selected calibration uses **snow compaction exponent = 0.26, SWE scale = 1.0** (`run03_column_transient_calibration_SnowExp_0p26_swe_1p0_2011-2020`). The calibrated 1D pressure and temperature profiles are mapped uniformly to all 109 columns of the 2D transect mesh as initial conditions before Step 4.

---

## Requirements

### ATS
- **ATS v1.6-dev**, commit `414f81aa`
- Source: [https://github.com/amanzi/ats](https://github.com/amanzi/ats) (BSD License)

### Python
```
python >= 3.9
numpy, scipy, pandas, xarray, h5py, matplotlib, pymannkendall
```

---

## Running the Model

```bash
# Step 1 — freeze-up
ats inputs/column_freezeup_imnavait.xml >out.log 2>out.err &

# Step 2 — steady-state spin-up
ats inputs/column_spinup_imnavait_smoothedMet.xml >out.log 2>out.err &

# Step 3 — 1D transient calibration (1981–2020) or (2019-2024)
ats inputs/column_spinup_imnavait_realMet_daac_ABoVE_snowmodel_<xxxx>_<xxxx>.xml >out.log 2>out.err &

# Step 4 — 2D transient run (1981–2020)
mpirun -np 8 ats inputs/transect_run_imnavait_realMet_daac_ABoVE_snowmodel_1981_2020.xml >out.log 2>out.err &

# Step 4 — 2D transient run (2019–2024 extension)
mpirun -np 8 ats inputs/transect_run_imnavait_realMet_daac_ABoVE_snowmodel_2019_2024.xml >out.log 2>out.err &
```

Outputs (checkpoints and observation `.dat` files) are written to the corresponding `runs/` subdirectory.

---

## Reproducing Figures

```bash
cd scripts
jupyter notebook IC45y_Figures.ipynb
```

Data are read from `../data/` and model outputs from `../runs/`.

---

## Field Site

**Imnavait Creek**, North Slope of Alaska (~68°36′39″N, 149°18′56″W), ~10 km east of Toolik Field Station.

- 109 m hillslope-to-riparian transect instrumented with 18 piezometers (IMP-01 to IMP-19)
- Thermistor arrays (T.rod.X®) at IMP-05, 09, 14, 18
- SNOTEL station (Site ID: 968) for snow depth, density, and surface temperature

---

## Climate Forcing

**NASA DAAC ABoVE SnowModel** (Liston et al., 2023) — bias-corrected, 3 km resolution, daily meteorological reconstruction.

> Liston, G.E., Reinking, A.K., & Boleman, N.T. (2023). *Daily SnowModel Outputs Covering the ABoVE Core Domain, 3-km Resolution, 1980–2020*. ORNL DAAC. https://doi.org/10.3334/ORNLDAAC/2105

---

## Data Availability

All data and model outputs will be archived in **ESS-DIVE** upon manuscript acceptance. Permanent DOIs will be provided in the published version.

---

## Citation

```bibtex
@article{mukherjee_inreview,
  author  = {Mukherjee, Neelarun and Gao, Bo and Coon, Ethan T. and Shuai, Pin and
             Hill, Devon and Neilson, Bethany T. and Kling, George W. and
             Chen, Jingyi and Cardenas, M. Bayani},
  title   = {Seasonal to decadal changes in supra-permafrost soil and aquifer
             thermal hydrology: broader, deeper, and longer activity},
  journal = {AGU Advances},
  year    = {in review}
}
```

---

## Funding

- U.S. DOE Office of Science, BER Environmental System Science Program (DE-SC0024091)
- NASA Terrestrial Hydrology Program (80NSSC18K0983)
- NASA FINESST
- U.S. NSF Arctic grants: ARC-1204220, DEB-1026843, 1637459, 0639805, 2224743, PLR-1504006, OPP-1107593, 1936759

---

## Contact

**Neelarun Mukherjee** — [neelarun@utexas.edu](mailto:neelarun@utexas.edu)  
Department of Earth and Planetary Sciences, The University of Texas at Austin
