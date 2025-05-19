
# 🌊 HighResMIP Tropical Instability Wave (TIW) Diagnostics

This folder contains scripts to compute **Tropical Instability Wave (TIW)** metrics from multiple CMIP6 **HighResMIP** models. The focus is on **sea surface temperature (SST)** variability at intraseasonal timescales (10–50 days) along the equatorial Pacific.

The metrics derived include:
- **Band-passed SST anomalies (10–50 day)**
- **Monthly standard deviation (TIW intensity)**

All outputs are derived using [CDO](https://code.mpimet.mpg.de/projects/cdo) tools and standardized filenames.

---

## 📅 Model Time Periods

| Model              | Scenario       | Years        |
|--------------------|----------------|--------------|
| AWI-CM-1-1-MR      | piControl      | 2851–2900    |
| AWI-CM-1-1-MR      | abrupt-4xCO2   | 1951–2000    |
| CNRM-CM6-1-HR      | piControl      | 2056–2105    |
| CNRM-CM6-1-HR      | abrupt-4xCO2   | 1950–1999    |
| GFDL-CM4           | piControl      | 0601–0650    |
| GFDL-CM4           | abrupt-4xCO2   | 0101–0150    |
| HadGEM3-GC31-MM    | piControl      | 2300–2349    |
| HadGEM3-GC31-MM    | abrupt-4xCO2   | 1950–1999    |

Each model's `tos_Oday_*.nc` file should be extracted and subsetted over the equatorial Pacific domain (e.g., 10°S–10°N, 120°E–80°W).

---

## 🔁 Processing Workflow

The following steps are applied to each model's daily SST data:

✅ Step 1: Select time range
```bash
cdo selyear,<start>/<end> input.nc temp1.nc
---
✅ Step 2: Set missing values to zero
cdo setmisstoc,0 temp1.nc temp2.nc
---
✅ Step 3: Apply 10–50 day band-pass filter
cdo bandpass,7.3,36.5 -del29feb temp2.nc output_10-50fft.nc
---
✅ Step 4: Compute monthly standard deviation (TIW intensity)
cdo monstd output_10-50fft.nc output_10-50fft_std.nc
---
📌 Notes
Only tos (SST) is used; velocity and SSH fields are not available for HighResMIP daily outputs.

Leap days (Feb 29) are removed to ensure a regular calendar before FFT.

No detrending is applied — the bandpass filter removes low-frequency signals.

All processing is done using CDO.

---

## 📦 Data Sources

### ✅ Raw Data

- **Source**: [CMIP6 HighResMIP archive](https://esgf-node.llnl.gov/projects/cmip6/)
- **Access**: ESGF search + `esgf-pyclient` or Wget script
- **Variable used**: `tos` (sea surface temperature)
- **Frequency**: Daily (`Oday` table)
- **Spatial domain**: Eastern Equatorial Pacific (typically extracted as 10°S–10°N, 180°–80°W)
---
### 📥 Processed Outputs

All processed outputs (bandpassed SST and TIW intensity fields) are publicly available at:

▶ **Zenodo archive**  
[https://doi.org/10.5281/zenodo.15454130](https://doi.org/10.5281/zenodo.15454130)
---
👤 Contact
Aoyun Xue
Postdoctoral Researcher
University of California, Santa Barbara
📧 aoyunxue@ucsb.edu
