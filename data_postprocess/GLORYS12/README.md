# GLORYS12 Tropical Instability Wave (TIW) Diagnostics

This folder contains NCL scripts used to extract and quantify **Tropical Instability Wave (TIW)** signals from the **GLORYS12 (CMEMS)** ocean reanalysis. The analysis targets equatorial Pacific variability in daily **SST**, **SSH**, horizontal velocity fields, and eddy kinetic energy (EKE).

---

## 🌊 Dataset

- **Source**: CMEMS GLORYS12V1 global reanalysis (1/12° horizontal resolution)
- **Temporal resolution**: Daily
- **Variables analyzed**: `thetao` (temperature), `uo`, `vo` (zonal/meridional velocity), `zos` (SSH)
- **Spatial domain**: 10°S–10°N, 120°E–80°W

---

## 🔁 Processing Workflow

### Step 1: Compute Daily Anomalies
- Remove the daily climatology to isolate sub-seasonal variability
- Applied to: SST, SSH, zonal and meridional velocity (uo, vo), and subsurface temperature

Scripts:
- Embedded in: `CMEMS_10-50_fft_*.ncl`

---

### Step 2: Apply 10–50 Day Band-Pass Filter (Fourier)
- Use FFT-based band-pass filter to isolate TIW-band (10–50 day) signals
- Applied separately for:
  - Surface: SST, SSH, u, v
  - Subsurface: temperature (`thetao`), u, v at selected depths

Scripts:
- `CMEMS_10-50_fft_u_surface.ncl`
- `CMEMS_10-50_fft_u_surface_new.ncl`
- `CMEMS_10-50_fft_v_depth_new.ncl`
- `CMEMS_10-50_fft_u_depth_new.ncl`
- `CMEMS_10-50_fft_t_surface.ncl`
- `CMEMS_10-50_fft_thetao_depth_new.ncl`

---

### Step 3: Compute TIW Intensity and Eddy Kinetic Energy
- **TIW intensity**: temporal standard deviation (RMS) of band-passed anomalies in SST, SSH, and velocities  
  \[
  \sigma(x,y) = \sqrt{ \frac{1}{N} \sum_t \left( X'_{\mathrm{bp}}(t) \right)^2 }
  \]

- **Eddy kinetic energy (EKE)**: computed from filtered u', v' as:  
  \[
  \mathrm{EKE} = \frac{1}{2} \, \rho_0 \left( u'^2 + v'^2 \right), \quad \rho_0 = 1024 \, \mathrm{kg/m^3}
  \]

Scripts:
- `calculate_TIW_std_variance_CMEMS_surface.ncl`
- `calculate_TIW_std_variance_CMEMS_depth.ncl`
- `calculate_TIW_variance_stddev_CMEMS_fourier.ncl`

---

## 📂 Output

- TIW-band **variance** and **standard deviation** for:
  - SST and SSH
  - Horizontal velocity (uo, vo)
  - Subsurface temperature
- TIW-band **EKE fields** at surface and selected subsurface depths  
- Dimensions: 2D fields (lon × lat) and vertical sections (lon/lat × depth)

---

## 📌 Notes

- Band-pass filtering is performed using **discrete Fourier transform (FFT)** with a 10–50 day window.
- EKE is computed using filtered velocity anomalies only; no mean-flow KE included.
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
