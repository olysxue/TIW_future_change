#!/bin/bash
# Author: Aoyun Xue
# Purpose: Compute daily anomalies and 20-yr moving detrended versions
# Format: consistent with sst_abrupt_4xCO2_0001-0150_remap_equator_moving_20yr_detrend.nc

# === abrupt-4xCO2 only ===
SCEN="abrupt-4xCO2"
STARTYR="0001"
ENDYR="0150"
REGION="remap_equator"
SUFFIX="moving_20yr_detrend"

VARS=("sst" "ssh")

for VAR in "${VARS[@]}"; do

  echo "Processing $VAR for $SCEN..."

  # === Input ===
  INFILE="/glade/derecho/scratch/aoyunxue/E3SMv2/${SCEN}/daily/${VAR}_${SCEN}_${STARTYR}-${ENDYR}_${REGION}.nc"

  # === Temporary files ===
  CLIM="${VAR}_${SCEN}_${STARTYR}-${ENDYR}_${REGION}_clim.nc"
  ANOM="${VAR}_${SCEN}_${STARTYR}-${ENDYR}_${REGION}_anom.nc"
  TREND="${VAR}_${SCEN}_${STARTYR}-${ENDYR}_${REGION}_20yrtrend.nc"

  # === Final output ===
  OUTFILE="${VAR}_${SCEN}_${STARTYR}-${ENDYR}_${REGION}_${SUFFIX}.nc"

  # Step 1: Daily climatology
  echo "  ➤ Computing daily climatology..."
  cdo -ydaymean $INFILE $CLIM

  # Step 2: Subtract to get anomaly
  echo "  ➤ Computing daily anomalies..."
  cdo -ydaysub $INFILE $CLIM $ANOM

  # Step 3: Run 20-yr sliding window (~7300 days)
  echo "  ➤ Running 20-year mean filter..."
  cdo runmean,7300 $ANOM $TREND

  # Step 4: Subtract sliding trend
  echo "  ➤ Subtracting to get final detrended anomaly..."
  cdo sub $ANOM $TREND $OUTFILE

  echo "Output written: $OUTFILE"

  # (Optional) cleanup
  rm -f $CLIM $ANOM $TREND

done

echo " Finished all variables for $SCEN."
