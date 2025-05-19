#!/bin/bash
# Author: Aoyun Xue
# Purpose: Compute daily anomalies only (no detrending)
# For: picontrol

SCEN="picontrol"
STARTYR="0001"
ENDYR="0500"
REGION="remap_equator"

VARS=("sst" "ssh")

for VAR in "${VARS[@]}"; do

  echo "Processing $VAR for $SCEN..."

  # === Input ===
  INFILE="/glade/derecho/scratch/aoyunxue/E3SMv2/${SCEN}/daily/${VAR}_${SCEN}_${STARTYR}-${ENDYR}_${REGION}.nc"

  # === Output ===
  CLIM="${VAR}_${SCEN}_${STARTYR}-${ENDYR}_${REGION}_clim.nc"
  ANOM="${VAR}_${SCEN}_${STARTYR}-${ENDYR}_${REGION}_anom.nc"

  # Step 1: Daily climatology
  echo "  ➤ Computing daily climatology..."
  cdo -ydaymean $INFILE $CLIM

  # Step 2: Subtract to get daily anomaly
  echo "  ➤ Computing daily anomaly..."
  cdo -ydaysub $INFILE $CLIM $ANOM

  echo " Output written: $ANOM"
  rm -f $CLIM

done

echo "Done: daily anomalies (no detrend) for $SCEN"
