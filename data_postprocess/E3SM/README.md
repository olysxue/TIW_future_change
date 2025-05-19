# E3SMv1/v2 TIW Diagnostics: SST, SSH, and Eddy Kinetic Energy

This folder contains scripts and workflows for analyzing **Tropical Instability Waves (TIWs)** in daily outputs from **E3SMv1** and **E3SMv2** simulations.
- Sea Surface Temperature (SST)
- Sea Surface Height (SSH)
- TIW-induced Eddy Kinetic Energy (TIW-EKE) via kinetic energy decomposition

The analysis includes two idealized CMIP6 scenarios:
- `piControl` (pre-industrial control)
- `abrupt-4xCO2` (abrupt quadrupled CO2 forcing)

---

## 🔁 Processing Workflow

## 🔁 Workflow Overview

### 🧩 Step 1: Variable Extraction
- Extract daily SST, SSH, UVEL, VVEL from raw E3SM outputs
- Concatenate multi-year files into unified NetCDFs
- Script: `Step1_extract_E3SM_concat_var.sh`

---

### 🧮 Step 2: Remove Mean Seasonal Cycle
- Compute daily climatology (multi-year mean for each calendar day)
- Subtract to obtain daily SST and SSH anomalies
- Scripts:
  - `Step2_calculate_E3SM_daily_anom_picontrol_no_detrend.sh`
  - `Step2_calculate_E3SM_daily_anom_abrupt_4xCO2_moving_20year_detrend.sh`

---

### 🔍 Step 3: Band-Pass Filtering (10–50 Days)
- Apply 10–50 day Fourier band-pass filter using NCL (`ezfftf`)
- Isolate TIW-band variability in SST and SSH
- Scripts:
  - `Step3_calculate_E3SM_picontrol_10-50fft.ncl`
  - `Step3_calculate_E3SM_abrupt_4xCO2_10-50fft.ncl`

---

### 📈 Step 4: Estimate TIW Intensity & EKE

#### 4.1: Standard Deviation of Filtered SST/SSH (TIW Intensity)
- Compute 7-pentad (~35 day) running std dev → monthly TIW intensity
- Scripts:
  - `Step4_1_calculate_TIW_SST_SSH_intensity_pentad.sh`
  - `Step4_2_E3SM_TIW_intensity_pentad_to_monthly.ncl`

#### 4.2: Compute TIW-EKE via KE Decomposition
- Use pentad UVEL/VVEL to calculate:
  - Total KE (TKE) from daily velocities
  - Background KE (LSKE) from monthly mean velocities
  - TIW-EKE = TKE − LSKE
- Script: `Step4_3_calculate_E3SM_TIW_EKE.ncl`

---


## 📂 Script Summary

| File | Description | Software |
|------|-------------|----------|
| `Step1_extract_E3SM_concat_var.sh` | Combine UVEL/VVEL/SST/SSH from raw E3SM output | cdo nco |
| `Step2_*` | Compute SST/SSH anomalies (with/without detrending) | cdo |
| `Step3_*` | Apply 10–50 day bandpass filter | ncl |
| `Step4_1_*` | Calculate TIW SST/SSH intensity (stddev) | cdo |
| `Step4_2_*` | Convert intensity to monthly scale | ncl |
| `Step4_3_*` | Compute TIW-EKE from KE decomposition | ncl |

---

## 📌 Notes

- This workflow **does not compute KE from u/v directly**, and **does not use daily velocities**.
- No energy conversion or flux divergence terms are included due to limited variable availability.
- Data are regridded to 1° and masked to the equatorial Pacific for comparison across models and resolutions.

---

## 📥 Data Availability

All processed outputs are available at:

▶ **Zenodo Archive**  
[https://doi.org/10.5281/zenodo.15454130](https://doi.org/10.5281/zenodo.15454130)

Includes:
- SST/SSH anomalies and TIW intensity
- Monthly TIW-EKE fields based on KE decomposition

Raw model outputs are available via:  
▶ [https://aims2.llnl.gov/](https://aims2.llnl.gov/)

---

## 📧 Contact

**Aoyun Xue**  
Postdoctoral Researcher  
University of California, Santa Barbara  
📧 aoyunxue@ucsb.edu
