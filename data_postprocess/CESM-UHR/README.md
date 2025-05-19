# CESM-UHR Tropical Instability Wave (TIW) Data Processing Workflow

This folder contains NCL and CDO scripts used to process CESM-UHR simulation output for TIW-related diagnostics. The full workflow includes anomaly extraction, band-pass filtering, TIW intensity estimation, energy conversion calculations, eddy flux diagnostics, and kinetic energy decomposition.

All scripts were tested with **NCL 6.6.2** and **CDO 2.0**.

---

## ♻ Processing Workflow Overview

### **Step 1: Remove Daily Climatology**
- **Goal**: Remove the mean seasonal cycle to isolate high-frequency variability
- **Method**: Compute daily climatology from multiyear mean and subtract it
- **Input**: CESM-UHR daily outputs (SST, SSH, UVEL, VVEL, WVEL, TEMP, density)
- **Output**: Daily anomalies (`*_anom.nc`)
- **Script**: `Step1_Calculate_daily_anom.sh`

---

### **Step 2: Apply 10–50-Day Fourier Band-Pass Filter**
- **Goal**: Isolate TIW signals
- **Method**: Lanczos filter on daily anomalies
- **Scripts**:
  - `Step2_CESM_UHR_SST_filter_10–50.ncl`
  - `Step2_CESM_UHR_SSH_filter_10–50.ncl`
  - `Step2_CESM_UHR_UVEL_filter_10–50.ncl`
  - `Step2_CESM_UHR_VVEL_filter_10–50.ncl`
  - `Step2_CESM_UHR_WVEL_filter_10–50.ncl`
  - `Step2_CESM_UHR_density_10–50_filter.ncl`

---

### **Step 3: Compute TIW Intensity Metrics**
- **Goal**: Quantify TIW intensity
- **Method**: Compute TIW-EKE, standard deviation of filtered SST, SSH, VVEL
- **Script**: `Step3_calculate_TIW_metric_CESM-UHR.ncl`

---

### **Step 4: Calculate Energy Conversion Rates**
- **Goal**: Identify TIW energy sources
- **Terms**:
  - Barotropic conversion (BTR): \[-⟨u'v'⟩ ∂⟨u⟩/∂y\]
  - Baroclinic conversion (BCR): \[-⟨ρ'w'\⟩/g\]
- **Scripts**:
  - `Step4_calculate_APE.ncl`
  - `Step4_calculate_CKE.ncl`

---

### **Step 5: Estimate Eddy Heat and Momentum Flux Convergence**
- **Goal**: Evaluate TIW impact on heat and momentum transport
- **Diagnostics**:
  - Eddy heat flux: (u'T', v'T', w'T') and their convergence
  - Eddy momentum flux: (u'u', u'v', v'v') and their convergence
- **Scripts**:
  - `Step5_calculate_heat_flux_momentum_flux_convergence.ncl`
  - `Step5_calculate_eddy_heat_momentum_flux.sh`

---

### **Step 6: Decompose Kinetic Energy (KE)**
- **Goal**: Quantify TIW-EKE
- **Method**:
  - TKE = 0.5ρ(u² + v²)
  - LSKE = 0.5ρ(ū² + v̄²)
  - TIW-EKE = TKE − LSKE
- **Script**: `Step6_calculate_TIW_EKE_decomposition.sh`

---

## 📂 Script List

| Script | Description |
|--------|-------------|
| `Step1_Calculate_daily_anom.sh` | Remove daily climatology (CDO) |
| `Step2_CESM_UHR_SST_filter_10–50.ncl` | 10-50-day Band-pass filter SST anomalies (NCL) |
| `Step2_CESM_UHR_SSH_filter_10–50.ncl` | 10-50-day Band-pass filter SSH anomalies (NCL) |
| `Step2_CESM_UHR_UVEL_filter_10–50.ncl` | 10-50-day Band-pass filter zonal velocity anomalies (NCL)|
| `Step2_CESM_UHR_VVEL_filter_10–50.ncl` | 10-50-day Band-pass filter meridional velocity anomalies (NCL)|
| `Step2_CESM_UHR_WVEL_filter_10–50.ncl` | 10-50-day Band-pass filter vertical velocity anomalies(NCL) |
| `Step2_CESM_UHR_density_10–50_filter.ncl` | 10-50-day Band-pass filter density anomalies (NCL)|
| `Step3_calculate_TIW_metric_CESM-UHR.ncl` | Compute TIW-EKE and SST/SSH/VVEL-based metrics (NCL)|
| `Step4_calculate_APE.ncl` | Compute baroclinic energy conversion rate (NCL) |
| `Step4_calculate_CKE.ncl` | Compute barotropic energy conversion rate (NCL) |
| `Step5_calculate_heat_flux_momentum_flux_convergence.ncl` | Calculate eddy flux and convergence (NCL) |
| `Step5_calculate_eddy_heat_momentum_flux.sh` | Calculate eddy flux (CDO) |
| `Step6_calculate_TIW_EKE_decomposition.sh` | Decompose total KE to extract TIW-EKE (CDO) |

---

## 📌 Notes

- Diagnostics are evaluated over the **equatorial eastern Pacific (10°S–10°N, 180°–80°W)**
- Filtering and diagnostics are applied on daily anomalies
- Update variable names and file paths as needed for your system

---

## 📥 Data Access

- CESM-HR simulation output and Processed TIW data is available at:  
  ▶ [https://doi.org/10.5281/zenodo.15454130](https://doi.org/10.5281/zenodo.15454130)

**Download file**:  
- 🔹 **CESM_TIW_intensity.zip**
   Four TIW metrics 
- 🔹 **CESM_energy.zip**
   Barotropic (BTR) and baroclinic (BCR) energy conversion rates
- 🔹 **CESM_feedback.zip**  Monthly TIW-induced feedbacks
- 🔹 **CESM_2xCO2.zip**, **CESM_4xCO2.zip**, **CESM_PD.zip**  
   Raw CESM-UHR simulation outputs used for the diagnostics above.

---

## 📧 Contact

For questions, please contact:

**Aoyun Xue**  
Postdoctoral Researcher  
University of California, Santa Barbara  
📧 aoyunxue@ucsb.edu
