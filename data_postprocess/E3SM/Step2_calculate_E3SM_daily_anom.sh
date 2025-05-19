#!/bin/bash
# Author: Aoyun Xue
# Purpose: Compute SST and SSH daily anomalies using cdo
# Scenarios: picontrol (0111–0140), abrupt-4xCO2 (0001–0150)

SCENARIOS=("picontrol" "abrupt-4xCO2")

for SCEN in "${SCENARIOS[@]}"; do
  echo "Processing scenario: $SCEN"

  # === Input file name by scenario ===
  if [ "$SCEN" == "picontrol" ]; then
    T_RANGE="00010101-05001231"
  else
    T_RANGE="00010101-01501231"
  fi

  # === Input paths ===
  SST_IN="/glade/derecho/scratch/aoyunxue/CO2/E3SM/${SCEN}/tos/tos_day_E3SM-1-0_${SCEN}_gr_${T_RANGE}.nc"
  SSH_IN="/glade/derecho/scratch/aoyunxue/CO2/E3SM/${SCEN}/zos/zos_day_E3SM-1-0_${SCEN}_gr_${T_RANGE}.nc"

  # === Output file names ===
  SST_CLIM="tos_day_clim_${SCEN}.nc"
  SSH_CLIM="zos_day_clim_${SCEN}.nc"
  SST_ANOM="tos_day_anom_${SCEN}.nc"
  SSH_ANOM="zos_day_anom_${SCEN}.nc"

  # === SST ===
  echo "  ➤ Calculating SST climatology..."
  cdo -ydaymean $SST_IN $SST_CLIM
  echo "  ➤ Calculating SST anomaly..."
  cdo -ydaysub  $SST_IN $SST_CLIM $SST_ANOM

  # === SSH ===
  echo "  ➤ Calculating SSH climatology..."
  cdo -ydaymean $SSH_IN $SSH_CLIM
  echo "  ➤ Calculating SSH anomaly..."
  cdo -ydaysub  $SSH_IN $SSH_CLIM $SSH_ANOM

  echo " Done: $SCEN → $SST_ANOM and $SSH_ANOM"
done

echo "All scenarios processed."
