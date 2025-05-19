#!/bin/bash

# === Configuration ===
BASEDIR="/glade/derecho/scratch/aoyunxue/FROM_CHEYENNE/TIW_data/CESM_HR"
OUTPUT_DIR="$BASEDIR/flux/"
mkdir -p $OUTPUT_DIR

SCENS=("PD" "2xCO2" "4xCO2")

for scen in "${SCENS[@]}"; do
  echo "Processing $scen..."

  UVEL="$BASEDIR/$scen/UVEL/CESM.${scen}.day.UVEL_0_50m.anom.10-50.fft.nc"
  VVEL="$BASEDIR/$scen/VVEL/CESM.${scen}.day.VVEL_0_50m.anom.10-50.fft.nc"
  TEMP="$BASEDIR/$scen/TEMP/CESM.${scen}.day.TEMP_0_50m.anom.10-50.fft.nc"

  # --- Heat flux components (without ρ * cp) ---
  cdo mulc,-1.0 mul $UVEL $TEMP $OUTPUT_DIR/CESM2_HR_daily_heat_flux_u_${scen}.nc
  cdo mulc,-1.0 mul $VVEL $TEMP $OUTPUT_DIR/CESM2_HR_daily_heat_flux_v_${scen}.nc

  # --- Momentum flux components (without density) ---
  cdo mulc,-1.0 mul $UVEL $VVEL $OUTPUT_DIR/CESM2_HR_daily_momentum_flux_uv_${scen}.nc
  cdo mulc,-1.0 mul $UVEL $UVEL $OUTPUT_DIR/CESM2_HR_daily_momentum_flux_uu_${scen}.nc
  cdo mulc,-1.0 mul $VVEL $VVEL $OUTPUT_DIR/CESM2_HR_daily_momentum_flux_vv_${scen}.nc
  cdo
