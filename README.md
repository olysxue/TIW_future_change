# Code for "Consistent Projected changes of Tropical Instability Waves under greenhouse warming".

This repository contains the code supporting the manuscript:

## 🗂️ Repository Structure

- `Main_code/`  
  Scripts to produce the main figures1-4 in the main text.
  
- `Extended_data_code/`  
  Scripts for Extended Data Figures 1–10.

- `data_postprocess/CESM-UHR/`  
  CESM-UHR  data postprocessing scripts, including daily climatology removal, 10–50-day filtering, TIW intensity calculation, eddy kinetic energy decomposition, barotropic/baroclinic energy conversion, and TIW flux diagnostics.

- `data_postprocess/E3SM/`  
  data postprocessing Scripts for E3SMv1 and E3SMv2 outputs (e.g., filtering, TIW metrics).

- `data_postprocess/HighresMIP/`  
  postprocessing of CMIP6 HighResMIP models (SST-based TIW intensity).

- `data_postprocess/GLORYS12/`  
 TIW metric calculations from the GLORYS12 ocean reanalysis product.

---
## 📦 Data Availability

The CESM-UHR simulations and full processed datasets used to generate the figures in this manuscript are archived on Zenodo:

▶ **Zenodo Archive (TIW-EKE, BTR/BCR, filtered fields, and flux diagnostics)**  
[https://doi.org/10.5281/zenodo.15454130](https://doi.org/10.5281/zenodo.15454130)


- 🔹 **CESM_2xCO2.zip**, **CESM_4xCO2.zip**, **CESM_PD.zip**  
  Raw CESM-UHR simulation outputs used for the diagnostics (daily SST, SSH, wind stress, 3D temp and currents, etc.).
  
- 🔹 **CESM_TIW_intensity.zip**  
  TIW intensity metrics including:
  - Eddy kinetic energy (TIW-EKE)
  - Std of SST, SSH, and vertical velocity (VVEL) band-passed at 10–50 days

- 🔹 **CESM_energy.zip**  
  Barotropic (BTR) and baroclinic (BCR) energy conversion rates computed in CESM-UHR, supporting TIW energetics analysis.

- 🔹 **CESM_feedback.zip**  
  Monthly TIW-induced feedbacks:
  - Heat flux convergence
  - Momentum flux convergence

- 🔹 **E3SM_TIW_intensity.zip**  
   - Eddy kinetic energy (TIW-EKE)
   - Std of SST and SSH band-passed at 10–50 days

- 🔹 **GLORYS12_TIW_metrics.zip**  
   TIW intensity metrics including:
  - Eddy kinetic energy (TIW-EKE)
  - Std of SST, SSH, and vertical velocity (VVEL) band-passed at 10–50 days

All variables processed at daily resolution using FFT bandpass (10–50 days) and converted to monthly std or variance maps.

▶ **Raw model output sources:**
- CESM-UHR: [https://climatedata.ibs.re.kr](https://climatedata.ibs.re.kr)  
- E3SMv1/v2: [https://aims2.llnl.gov/](https://aims2.llnl.gov/)  
- CMIP6 HighResMIP: [https://esgf-node.llnl.gov/projects/cmip6/](https://esgf-node.llnl.gov/projects/cmip6/)  
- GLORYS12: [https://marine.copernicus.eu](https://marine.copernicus.eu)
---

## 🛠️ Dependencies

- [**NCL**](https://www.ncl.ucar.edu/) (tested with version **6.6.2**)  
- [**CDO**](https://code.mpimet.mpg.de/projects/cdo) (Climate Data Operators)
- [**NCO**](https://nco.sourceforge.net/) (NetCDF Operators; used for variable merging, dimension trimming, etc.)

---

## 📦 Citation

If you use this code, please cite:

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15462389.svg)](https://doi.org/10.5281/zenodo.15462389)

Xue, A. (2025). *Consistent Projected changes of Tropical Instability Waves under Greenhouse Warming*.  
Zenodo. https://doi.org/10.5281/zenodo.15462389

---
## 📧 Contact
For questions or assistance with the code and data, please contact:

Aoyun Xue
Postdoctoral Researcher
📩 aoyunxue@ucsb.edu
University of California, Santa Barbara
