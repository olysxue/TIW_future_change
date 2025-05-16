# Extended Data Figure Scripts (NCL, NCL v6.6.2)

This folder contains NCL scripts used to generate **Extended Data Figures 1–10** in the manuscript:

**"Consistent projected weakening of Tropical Instability Waves under greenhouse warming"**, submitted to *Nature Communications*.

All scripts were developed and tested using **NCL version 6.6.2**. Each script corresponds directly to one Extended Data Figure and is designed to reproduce the results presented in the manuscript.

## File Overview

| Script | Corresponding Figure | Description |
|--------|----------------------|-------------|
| `Extended_Figure1_plot_CMEMS_TIW_intensity.ncl` | Extended Data Fig. 1 | Compares TIW intensity (EKE, VVEL, SST, SSH) in GLORYS12, CESM-UHR, and E3SM models |
| `Extended_Figure2_VVEL_PD_2xCO2_4xCO2.ncl` | Extended Data Fig. 2 | Projects TIW metrics (EKE, SST, SSH) in E3SMv1 and v2 under piControl and 4xCO2 |
| `Extended_Figure3_plot_E3SM_SST_zonal_current_mean_seasonal_cycle.ncl` | Extended Data Fig. 3 | Compares SST and zonal currents between GLORYS12 and climate models |
| `Extended_Figure4_plot_CMIP6_HighResMIP_SST_variance_future_change.ncl` | Extended Data Fig. 4 | SST-based TIW intensity in CMIP6 HighResMIP models under CO₂ forcing |
| `Extended_Figure5_plot_E3SM_mean_state_change.ncl` | Extended Data Fig. 5 | Mean-state changes (SST, currents, shear) in E3SMv1 and v2 under 4xCO2 |
| `Extended_Figure6_plot_CESM_UHR_background_change_4season.ncl` | Extended Data Fig. 6 | Seasonal SST, wind, zonal currents, and their changes in CESM-UHR |
| `Extended_Figure7_plot_E3SM_seasonal_cycle_TIWEKE_SST_SSH.ncl` | Extended Data Fig. 7 | Seasonal evolution of TIW metrics in E3SMv1 and v2 (MAM vs SON) |
| `Extended_Figure8_plot_energy_conversion_seasonal_cycle.ncl` | Extended Data Fig. 8 | Seasonal evolution of barotropic/baroclinic energy conversions (BTR, BCR) |
| `Extended_Figure9_plot_CESM_TIW_heat_momentum_flux_feedback.ncl` | Extended Data Fig. 9 | TIW-induced heat flux convergence (HFC) and momentum flux convergence (MFC) |
| `Extended_Figure10_different_EKE_calculation_method_comparison_CESM_E3SMv1_v2.ncl` | Extended Data Fig. 10 | Comparison of direct vs indirect TIW-EKE computation methods in CESM and E3SM |

## How to Run

Each script can be executed independently using the command line:

ncl Extended_Figure1_plot_CMEMS_TIW_intensity.ncl
