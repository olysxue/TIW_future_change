
# E3SMv1 and E3SMv2 Tropical Instability Wave (TIW) Diagnostics

This folder contains NCL and CDO scripts used to process **E3SMv1** and **E3SMv2** simulation outputs to calculate TIW intensity. The focus is on SST and SSH variability and eddy kinetic energy (EKE) derived from KE decomposition.

All scripts were tested with **NCL 6.6.2** and **CDO 2.0**.

---

## 🔁 Processing Workflow Overview

### **Step 1: Remove Daily Climatology**
- **Goal**: Remove the mean seasonal cycle to isolate inter-daily variability
- **Input**: Daily SST and SSH from E3SM simulations
- **Output**: Daily SST and SSH anomalies (`*_anom.nc`)

---

### **Step 2: Apply 10–50 Day Band-Pass Filter**
- **Goal**: Isolate the TIW-band signals
- **Method**: fourier 10-50-day bandpass filter applied to daily SST and SSH anomalies
- **Output**: Band-passed SST and SSH fields (`*_fft.nc`)

---

### **Step 3: Compute TIW Intensity Metrics**
- **Goal**: Quantify the amplitude of TIWs
- **Method**: Compute the temporal standard deviation (RMS) of band-passed SST and SSH
- **Output**: TIW intensity for SST and SSH

---

### **Step 4: Compute TIW Eddy Kinetic Energy (TIW-EKE)**
- **Goal**: Estimate mesoscale kinetic energy maintained by TIWs via decomposition
- **Method**:
  - Compute **total KE** from daily horizontal velocity:  
    **TKE = 0.5 × ρ × (u² + v²)**
  - Compute **large scale KE** from monthly mean velocity:  
    **LSKE = 0.5 × ρ × (ū² + v̄²)**
  - Subtract to get TIW-EKE:  
    **TIW-EKE = TKE − LSKE**
- **Assumption**: Constant reference density ρ = **1024 kg/m³**
- **Units**: Joules per cubic meter (J·m⁻³)
- **Tools**: CDO + shell script (`calc_TIW_EKE_E3SM.sh`)
- **Output**: Monthly TIW-EKE

---

## 📂 Script List

| Script | Description |
|--------|-------------|
| `Step2_E3SM_SST_filter_10-50.ncl` | Band-pass filter SST anomalies |
| `Step2_E3SM_SSH_filter_10-50.ncl` | Band-pass filter SSH anomalies |
| `calc_TIW_intensity.ncl` | Compute stddev of filtered SST and SSH |
| `calc_TIW_EKE_E3SM.sh` | Compute TIW-EKE via KE decomposition using CDO, with ρ₀ = 1024 kg/m³ |

---

## 📌 Notes

- This folder only includes **SST**, **SSH**, and **TIW-EKE** diagnostics.
- Energy conversion terms (BTR, BCR) and flux-based diagnostics were not computed for E3SM.
- The region analyzed is the **equatorial Pacific** (approx. 10°S–10°N, 160°E–80°W).
- Scenarios include **PD**, **2xCO₂**, and **4xCO₂**, consistent with CESM-UHR.

---

## 📥 Data Access

All processed outputs are archived at:

▶ **Zenodo Archive**  
[https://doi.org/10.5281/zenodo.15454130](https://doi.org/10.5281/zenodo.15454130)

Includes:
- TIW intensity metrics (SST, SSH)
- TIW-EKE (via KE decomposition, using ρ₀ = 1024 kg/m³) for E3SMv1 and E3SMv2

Raw model outputs are available at:  
▶ [https://aims2.llnl.gov/](https://aims2.llnl.gov/)

---

## 📧 Contact

**Aoyun Xue**  
Postdoctoral Researcher  
University of California, Santa Barbara  
📧 aoyunxue@ucsb.edu
