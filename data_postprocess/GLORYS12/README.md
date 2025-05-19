# GLORYS12 Tropical Instability Wave (TIW) Diagnostics

This folder contains NCL scripts used to extract and quantify **Tropical Instability Wave (TIW)** signals from the **GLORYS12 (CMEMS)** ocean reanalysis. The analysis targets equatorial Pacific variability in daily **SST**, **SSH**, horizontal velocity fields, and eddy kinetic energy (EKE).

---

## 🌊 Dataset

- **Source**: CMEMS GLORYS12V1 global reanalysis (1/12° horizontal resolution)
- **Temporal resolution**: Daily
- **Variables analyzed**: `thetao` (temperature), `uo`, `vo` (zonal/meridional velocity), `zos` (SSH)
- **Spatial domain**: 10°S–10°N, 180°–80°W

---

## 🔁 Workflow Summary

### ✅ Step 1: Daily Anomaly Calculation
- Remove daily climatology to obtain sub-seasonal anomalies
- Applied to: SST, SSH, `uo`, `vo`, and `thetao`

**Script**:  
- `step1_calculate_GLORYS12_daily_anom.sh`

---

### ✅ Step 2: Band-Pass Filtering (10–50 Days)
- Apply FFT-based 10–50 day band-pass filter to isolate TIW-band variability
- Separate filtering for:
  - Surface fields: SST, SSH, `uo`, `vo`
  - Subsurface fields: `thetao`, `uo`, `vo`

**Scripts**:  
- `step2_GLORYS12_10-50_fft_sst_ssh.ncl`  
- `step2_GLORYS12_10-50_fft_thetao_depth_new.ncl`  
- `step2_GLORYS12_10-50_fft_u_depth_new.ncl`  
- `step2_GLORYS12_10-50_fft_v_depth_new.ncl`

---

### ✅ Step 3: TIW Metrics
- **TIW intensity**: standard deviation of filtered SST, SSH, or velocity anomalies  
- **TIW EKE**:  
  - Surface: from filtered `u′`, `v′`  
  - Full-field decomposition: from monthly KE difference  
    \[
    \text{TIW-EKE} = \frac{1}{2} \overline{u^2 + v^2} - \frac{1}{2} (\bar{u}^2 + \bar{v}^2)
    \]

**Scripts**:  
- `step3_calculate_TIW_metric_GLORYS12_surface.ncl`  
- `Step4_calculate_TIW_EKE_GLORYS12_decomposition.ncl`

---

## 📂 Output

- TIW metrics:
  - SST-based metric
  - SSH-based metric
  - TIW-EKE direct calculation
  - TIW-EKE indirect calculation

---

## 📌 Notes

- Band-pass filtering is performed using **discrete Fourier transform (FFT)** with a 10–50 day window.
- SSH anomalies (`zos_anom`) are included alongside SST for surface diagnostics.
- These diagnostics are intended for comparison with model-based TIW simulations (e.g., CESM, E3SM).

---

## 📥 Data Access

All processed outputs are archived at:

▶ **Zenodo Archive**  
[https://doi.org/10.5281/zenodo.15454130](https://doi.org/10.5281/zenodo.15454130)

---

## 📧 Contact

**Aoyun Xue**  
Postdoctoral Researcher  
University of California, Santa Barbara  
📧 aoyunxue@ucsb.edu
