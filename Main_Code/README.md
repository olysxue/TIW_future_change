# Main Figure Scripts (NCL v6.6.2)

This folder contains NCL scripts used to generate **Main Figures 1–4** in the manuscript:

**"Consistent projected weakening of Tropical Instability Waves under greenhouse warming"**, submitted to *Nature Communications*.

All scripts were developed and tested using **NCL version 6.6.2**. Each script corresponds to one main figure and relies on CESM-UHR simulation outputs.

---

## Script Overview and Data Sources

| Script | Figure | Description | Required Data |
|--------|--------|-------------|---------------|
| `Figure1_CESM_TIW_intensity.ncl` | Fig. 1 | CESM-UHR-simulated TIW intensity and structure compared to reanalysis | four TIW metrics in CESM-UHR: EKE, SST, SSH, and VVEL fields|
| `Figure2_plot_CESM_PD_2xCO2_4xCO2_zonal_ccurrent_TIW_intensity.ncl` | Fig. 2 |Projected changes of tropical mean background state modulating in CESM-UHR | CESM-UHR monthly mean outputs under PD, 2xCO2, and 4xCO2: wind stress, temp, zonal currents |
| `Figure3_CESM_seasonal_zonal_current_TIW_intensity_response.ncl` | Fig. 3 | Seasonal evolution of background conditions and TIW intensity under different CO2 simulations in CESM-UHR | four TIW metrics in CESM-UHR: EKE, SST, SSH, and VVEL fields|
| `Figure4_TIW_EKE_budget_analysis.ncl` | Fig. 4 | TIW-EKE budget analysis under different CO2 simulations in CESM-UHR | Barotrpic energy conversion rate and Baroclinic energy conversion rate in CESM-UHR |

---
Due to file size constraints, only selected `.nc` are included in this GitHub repo.

- **Full processed datasets of CESM-UHR**  are available at Zenodo:[doi:10.5281/zenodo.15454130](https://doi.org/10.5281/zenodo.15454130)
- **E3SM datasets** are available from [https://aims2.llnl.gov/] 
- **CMIP HighresMIP datasets** are available from [https://esgf-node.llnl.gov/projects/cmip6/]
- **GLORYS12** is available at [https://data.marine.copernicus.eu/product/GLOBAL_MULTIYEAR_PHY_001_030/description]

---

## How to Run

Each script can be executed independently in a terminal with NCL:

ncl Figure*********.ncl
