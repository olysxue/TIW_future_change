# E3SMv1/v2 TIW Diagnostics: SST, SSH, and Eddy Kinetic Energy

This folder contains scripts and workflows for analyzing **Tropical Instability Waves (TIWs)** in daily outputs from **E3SMv1** and **E3SMv2** simulations, focusing on SST, SSH, and eddy kinetic energy (EKE) derived from kinetic energy decomposition.

The analysis includes two idealized CMIP6 scenarios:
- `piControl` (pre-industrial control)
- `abrupt-4xCO2` (quadrupled CO2 forcing)

---

## 🔁 Processing Workflow

### Step 1: Variable Extraction
- From the raw E3SM outputs, extract the following fields:
  - Daily: **SST**, **SSH**, and **total kinetic energy (TKE)**  
    *(Note: TKE is precomputed, not derived in this workflow)*  
- No u or v extraction from daily files, we cannot calculate TIW-EKE directly based on 10-50-day bandpass filtered velocity fields.

---

### Step 2: Remove Mean Seasonal Cycle
- Compute daily climatology of SST and SSH using a multi-year mean.
- Subtract the climatology to obtain daily anomalies
---

### Step 3: Band-Pass Filtering (10–50 Days)
- Apply a 10-50-day Fourier band-pass filter to the SST and SSH anomalies to isolate TIW-band variability.
- Output: Fourier Band-passed SST and SSH fields

---

### Step 4: Estimate TIW Intensity
- Compute the temporal standard deviation (RMS) of filtered SST and SSH:
- Used as a proxy for TIW strength.

---

### Step 5: Compute TIW-EKE from Precomputed KE Fields
- Decompose total kinetic energy into large-scale and eddy components:
Goal: Quantify TIW-EKE
Method:
TKE = 0.5ρ(u² + v²)
LSKE = 0.5ρ(ū² + v̄²)
TIW-EKE = TKE − LSKE

---

## 📂 Script Overview

| Script | Description |
|--------|-------------|
| `Step2_E3SM_SST_filter_10-50.ncl` | Band-pass filter for SST anomalies |
| `Step2_E3SM_SSH_filter_10-50.ncl` | Band-pass filter for SSH anomalies |
| `calc_TIW_intensity.ncl` | Compute TIW intensity from filtered SST and SSH |
| `calc_TIW_EKE_E3SM.sh` | Subtract background KE from monthly TKE to get TIW-EKE |

---

## 📌 Notes

- This workflow **does not compute KE from u/v directly**, and **does not use daily velocities**.
- The analysis is performed for both `piControl` and `abrupt-4xCO₂` in E3SMv1 and v2.
- No energy conversion or flux divergence terms are included due to limited variable availability.

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
