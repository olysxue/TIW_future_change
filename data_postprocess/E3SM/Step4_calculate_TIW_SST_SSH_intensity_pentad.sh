#!/bin/bash
# Author: Aoyun Xue
# Purpose: Calculate TIW intensity (runstd) and climatology, and interpolate to 1° grid
# Inputs follow your naming standard

# === Setup ===
SCENS=("picontrol" "abrupt_4xCO2")
VARS=("sst" "ssh")
GRID="r360x180"   # 1-degree grid
WIN=7             # sliding window size (7 pentads ≈ 35 days)

for SCEN in "${SCENS[@]}"; do

  # === Year range and filename suffix ===
  if [ "$SCEN" == "picontrol" ]; then
    YEARS="0001-0500"
    DETREND=""  # no suffix
    TAG="10-50fft"
  else
    YEARS="0001-0150"
    DETREND="_moving_20yr_detrend"
    TAG="10-50fft"
  fi

  for VAR in "${VARS[@]}"; do
    echo " Processing: $VAR - $SCEN"

    # === Filenames ===
    INFILE="${VAR}_${SCEN}_${YEARS}_remap_equator${DETREND}_${TAG}.nc"
    OUT_STD="${VAR}_${SCEN}_${YEARS}_remap_equator${DETREND}_${TAG}_1mon_runstd.nc"
    OUT_STD_MEAN="${VAR}_mon_${SCEN}_${YEARS}_remap_equator${DETREND}_${TAG}_1mon_runstd_mean.nc"
    OUT_STD_1DEG="${VAR}_${SCEN}_${YEARS}_remap_equator${DETREND}_${TAG}_1mon_runstd_1degree.nc"
    OUT_MEAN_1DEG="${VAR}_mon_${SCEN}_${YEARS}_remap_equator${DETREND}_${TAG}_1mon_runstd_mean_1degree.nc"

    # === Step 1: sliding stddev
    echo "  ➤ runstd..."
    cdo runstd,${WIN} $INFILE $OUT_STD

    # === Step 2: timmean
    echo "  ➤ timmean..."
    cdo timmean $OUT_STD $OUT_STD_MEAN

    # === Step 3: interpolate to 1°
    echo "  ➤ remap to 1 degree..."
    cdo remapbil,$GRID $OUT_STD $OUT_STD_1DEG
    cdo remapbil,$GRID $OUT_STD_MEAN $OUT_MEAN_1DEG

    echo " Finished: $OUT_STD_1DEG and $OUT_MEAN_1DEG"
  done
done

