# CESM-UHR Tropical Instability Wave (TIW) Data Processing Workflow

This folder contains NCL and cdo scripts used to process CESM-UHR simulation output for TIW-related diagnostics. The full workflow includes anomaly extraction, band-pass filtering, TIW intensity estimation, energy conversion calculations, and flux diagnostics.

All scripts were tested with **NCL6.6.2 & cdo 2.0**.

---

## 🔁 Processing Workflow Overview

### **Step 1: Remove Daily Climatology**
- **Goal**: Remove mean seasonal cycle to isolate inter-daily variability
- **Method**: Compute daily climatology from multiyear mean and subtract it
- **Input**: CESM-UHR daily outputs (SST, SSH, UVEL, VVEL)
- **Output**: Daily anomalies (`*_anom.nc`)

### **Step 2: Apply Fourier 10–50 Day Band-Pass Filter**
- **Goal**: Isolate the TIW signal in the 10–50 day band
- **Method**: Lanczos band-pass filtering on daily anomalies
- **Scripts**:
  - `CESM_UHR_SST_filter_10–50.ncl`
  - `CESM_UHR_SSH_filter_10–50.ncl`
  - `CESM_UHR_UVEL_filter_10–50.ncl`
  - `CESM_UHR_VVEL_filter_10–50.ncl`

### **Step 3: Compute TIW Intensity Metrics**
- **Goal**: Quantify TIW amplitude using filtered signal variance
- **Variables**: SST, SSH, UVEL, VVEL (band-pass filtered)
- **Method**: Calculate local standard deviation (RMS)

### **Step 4: Calculate Energy Conversion Rates**
- **Goal**: Diagnose TIW energy sources
- **Terms computed**:
  - Barotropic energy conversion (BTR): \[-⟨u'v'⟩ ∂⟨u⟩/∂y\]
  - Baroclinic energy conversion (BCR): \[-g/ρ₀ ⟨w'T'⟩/⟨∂T/∂z⟩\]
- **Input**: Filtered and mean-state fields of velocity and temperature

### **Step 5: Estimate TIW-Induced Fluxes**
- **Goal**: Assess the dynamical impact of TIWs on large-scale momentum and heat budgets
- **Terms**:
  - TIW heat flux (u'T', v'T' w'T') and its convergence
  - TIW momentum flux (u'u', u'v', v'v') and its convergence

---

## 📂 Script List

| Script | Description |
|--------|-------------|
| `CESM_UHR_SST_filter_10–50.ncl` | Band-pass filter SST anomalies |
| `CESM_UHR_SSH_filter_10–50.ncl` | Band-pass filter SSH anomalies |
| `CESM_UHR_UVEL_filter_10–50.ncl` | Band-pass filter zonal velocity anomalies |
| `CESM_UHR_VVEL_filter_10–50.ncl` | Band-pass filter meridional velocity anomalies |
| _(Other scripts not shown here may handle climatology removal, flux calculation, etc.)_ |

---

## 📌 Notes

- All diagnostics are calculated for the **equatorial Pacific** (approx. 5°S–5°N, 160°–90°W).
- Filtering and variance are applied on a **daily timescale**, after removing climatology.
- You may need to modify file paths and variable names to match your local setup.

---

## 📥 Data Access

- Processed outputs are archived at:  
  ▶ [https://doi.org/10.5281/zenodo.15454130](https://doi.org/10.5281/zenodo.15454130)  

---

## 📧 Contact

For questions or data access assistance, please contact:

**Aoyun Xue**  
Postdoctoral Researcher  
University of California, Santa Barbara  
📧 aoyunxue@ucsb.edu

