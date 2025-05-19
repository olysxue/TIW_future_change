#!/bin/bash
# Author: Aoyun Xue
# Purpose: Compute TIW-EKE from CESM/E3SM daily velocity fields
# Scenarios: PD (0111-0140), 2xCO2 / 4xCO2 (0141-0170)

LATLON="10S_10N.180_80W"
LEVEL="levidx0-50m"
RES="025degree"

for SCEN in PD 2xCO2 4xCO2; do
  if [ "$SCEN" == "PD" ]; then
    YR="0111-0140"
  else
    YR="0141-0170"
  fi

  echo "Processing scenario: $SCEN | Period: $YR"

  # === Input files ===
  UFILE="UVEL.${SCEN}.day.${YR}.${LATLON}.${LEVEL}_${RES}.nc"
  VFILE="VVEL.${SCEN}.day.${YR}.${LATLON}.${LEVEL}_${RES}.nc"

  # === Output files ===
  TKE_OUT="TKE.${SCEN}.mon.${YR}.${LATLON}.${LEVEL}_${RES}_new.nc"
  KE_BG="KE.${SCEN}.mon.${YR}.${LATLON}.${LEVEL}_new_method_${RES}.nc"
  EKE_OUT="TIW_EKE.${SCEN}.mon.${YR}.${LATLON}.${LEVEL}_${RES}.nc"

  # === Step 1: Compute daily total KE = 0.5 * (u² + v²)
  echo "  ➤ Computing daily total KE for $SCEN..."
  cdo -expr,"TKE=0.5*(UVEL*UVEL + VVEL*VVEL)" -merge $UFILE $VFILE tke_daily_${SCEN}.nc

  # === Step 2: Monthly mean of daily total KE
  echo "  ➤ Averaging total KE monthly..."
  cdo monmean tke_daily_${SCEN}.nc $TKE_OUT

  # === Step 3: Monthly mean of UVEL & VVEL
  echo "  ➤ Calculating monthly mean UVEL/VVEL..."
  cdo monmean $UFILE UVEL_mon_${SCEN}.nc
  cdo monmean $VFILE VVEL_mon_${SCEN}.nc

  # === Step 4: Background KE from monthly mean flow
  echo "  ➤ Computing background KE from monthly mean flow..."
  cdo -expr,"KE=0.5*(UVEL*UVEL + VVEL*VVEL)" -merge UVEL_mon_${SCEN}.nc VVEL_mon_${SCEN}.nc $KE_BG

  # === Step 5: Subtract background KE from total KE
  echo "  ➤ Deriving TIW-EKE..."
  cdo sub $TKE_OUT $KE_BG $EKE_OUT

  # Clean up
  rm -f tke_daily_${SCEN}.nc UVEL_mon_${SCEN}.nc VVEL_mon_${SCEN}.nc

  echo "✅ Finished $SCEN → Output: $EKE_OUT"
done
