# Code for "Consistent Projected weakening of Tropical Instability Waves under greenhouse warming" submitted to *Nature Communications*.

This repository contains the code supporting the manuscript:

## 🗂️ Repository Structure

- `Main_code/`  
  Scripts to produce the main figures in the main text.
  
- `Extended_data_code/`  
  Scripts for Extended Data Figures 1–10.

- `data_postprocess/CESM-UHR/`  
  CESM-UHR postprocessing scripts, including daily climatology removal, 10–50-day filtering, TIW intensity calculation, eddy kinetic energy decomposition, barotropic/baroclinic energy conversion, and TIW flux diagnostics.

- `data_postprocess/E3SM/`  
  Scripts for E3SMv1 and E3SMv2 outputs (e.g., filtering, TIW metrics, BTR/BCR).

- `data_postprocess/HighresMIP/`  
  Processing of CMIP6 HighResMIP models (e.g., zonal averaging, TIW indices).

- `data_postprocess/GLORYS12/`  
  Filtering and metric extraction from the GLORYS12 ocean reanalysis product.

---
## 📦 Data Availability

Due to size limitations, only selected sample NetCDF files are hosted in this repository.  
The full processed datasets used to generate the figures in this manuscript are archived on Zenodo:

▶ **Zenodo Archive (TIW-EKE, BTR/BCR, filtered fields, and flux diagnostics)**  
[https://doi.org/10.5281/zenodo.XXXXXXXX](https://doi.org/10.5281/zenodo.XXXXXXXX)

Includes:

- CESM-UHR: Four TIW metrics, barotropic/baroclinic energy conversion rates (BTR, BCR), heat & momentum flux convergence
- E3SMv1/v2: TIW-EKE, SST and SSH-based TIW metrics
- GLORYS12: TIW intensity diagnostics from EKE, SST, SSH, and velocity reanalysis fields

▶ **Raw model output sources:**
- CESM-UHR: [https://climatedata.ibs.re.kr](https://climatedata.ibs.re.kr)  
- E3SMv1/v2: [https://aims2.llnl.gov/](https://aims2.llnl.gov/)  
- CMIP6 HighResMIP: [https://esgf-node.llnl.gov/projects/cmip6/](https://esgf-node.llnl.gov/projects/cmip6/)  
- GLORYS12: [https://marine.copernicus.eu](https://marine.copernicus.eu)
---

## 🛠️ Dependencies

- [**NCL**](https://www.ncl.ucar.edu/) (tested with version **6.6.2**)  
- [**CDO**](https://code.mpimet.mpg.de/projects/cdo) (Climate Data Operators)

---

## 📧 Contact
For questions or assistance with the code and data, please contact:

Aoyun Xue
Postdoctoral Researcher
📩 aoyunxue@ucsb.edu
University of California, Santa Barbara
