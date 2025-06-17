CESM SSP585 Tropical Instability Wave (TIW) Data Processing Workflow

This directory contains CDO scripts and documentation for processing both historical (20th century) and SSP585 (21st century) CESM outputs, focusing on Tropical Instability Wave (TIW) intensity and feedback diagnostics.

Requirements: CDO 2.0+, Bash Shell

♻ Processing Workflow Overview

Step 1: Merge Historical & SSP585 Experiments

Goal: Combine historical (1850–2005) and SSP585 (2006–2100) anomaly fields into a continuous time series
Method: Merge time dimensions after preparing each dataset
Input: Multi-member historical and SSP585 outputs (e.g., SST, SSH, UVEL, VVEL)
Output: combined.nc
Script: Step1_merge_hist_ssp585.sh

cdo mergetime raw/hist/*.nc raw/ssp585/*.nc intermediate/combined.nc

Step 2: Region Selection & Variable Extraction

Goal: Extract the tropical Pacific domain (180°–300°E, 10°S–10°N) and select relevant variables
Method: Use CDO to subset spatial domain and variables
Script: Step2_extract_region_variables.sh

cdo sellonlatbox,180,300,-10,10 combined.nc intermediate/region.nc
cdo selvar,SST,SSH,UVEL,VVEL intermediate/region.nc intermediate/tiw_vars.nc

Step 3: Remove Climatology & Detrend

Goal: Remove seasonal cycle and linear trend to obtain pure interannual anomalies
Method: Compute seasonal climatology, subtract, then detrend
Script: Step3_climatology_detrend.sh

cdo seasonmean intermediate/tiw_vars.nc intermediate/season_clim.nc
cdo sub intermediate/tiw_vars.nc intermediate/season_clim.nc intermediate/anom.nc
cdo detrend intermediate/anom.nc intermediate/anom_detrend.nc

Step 4: Band-Pass Filtering (10–50 Days)

Goal: Isolate TIW signals in the 10–50 day band
Method: Apply CDO's bandpass filter
Script: Step4_bandpass_filter.sh

cdo setmisstonn intermediate/anom_detrend.nc intermediate/anom_filled.nc  # fill missing
cdo bandpass,10,50 intermediate/anom_filled.nc intermediate/anom_bp.nc

Step 5: Compute TIW Intensity

Goal: Quantify TIW activity via monthly standard deviation
Method: Use CDO's monstd operator
Script: Step5_compute_TIW_intensity.sh

cdo monstd intermediate/anom_bp.nc output/tiw_intensity.nc

Step 6: Trend & Ensemble Mean Analysis

Goal: Calculate ensemble mean across members and apply long-term smoothing (20-year running mean)
Method: CDO's ensmean and runmean operators
Script: Step6_ensemble_trend.sh

cdo ensmean output/tiw_intensity_*.nc output/tiw_intensity_ensmean.nc
cdo runmean,240 output/tiw_intensity_ensmean.nc output/tiw_intensity_20yr.nc

📂 Script List

Script

Description

Step1_merge_hist_ssp585.sh

Merge historical and SSP585 time series (CDO)

Step2_extract_region_variables.sh

Subset region and extract variables (CDO)

Step3_climatology_detrend.sh

Remove seasonal cycle and detrend (CDO)

Step4_bandpass_filter.sh

10–50 day bandpass filter (CDO)

Step5_compute_TIW_intensity.sh

Monthly STD for TIW intensity (CDO)

Step6_ensemble_trend.sh

Ensemble mean and 20-year smoothing (CDO)

📌 Notes

Region: Equatorial Pacific (10°S–10°N, 180°–300°E)

Input: Raw grid data with missing values filled or no missing values

Modify scripts to include additional variables or different time bands as needed

📥 Data Access

Access CESM SSP585 and historical raw data via institutional server or DOI link.

📧 Contact

Aoyun XuePostdoctoral Researcher, UCSBaoyunxue@ucsb.edu


