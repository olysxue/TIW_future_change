#!/bin/bash
# step1_process_4models_highresMIP_TIW.sh
# Author: Aoyun Xue
# Description: Extract TIW-band SST variability (10–50 day bandpass + monthly std) from HighResMIP models

# Target models
MODELS=("AWI-CM-1-1-MR" "CNRM-CM6-1-HR" "GFDL-CM4" "HadGEM3-GC31-MM")
SCENARIOS=("piControl" "abrupt-4xCO2")

# Manually specified time ranges
declare -A START
declare -A END

START["AWI-CM-1-1-MR_piControl"]=2851
END["AWI-CM-1-1-MR_piControl"]=2900
START["AWI-CM-1-1-MR_abrupt-4xCO2"]=1951
END["AWI-CM-1-1-MR_abrupt-4xCO2"]=2000

START["CNRM-CM6-1-HR_piControl"]=2056
END["CNRM-CM6-1-HR_piControl"]=2105
START["CNRM-CM6-1-HR_abrupt-4xCO2"]=1950
END["CNRM-CM6-1-HR_abrupt-4xCO2"]=1999

START["GFDL-CM4_piControl"]=0601
END["GFDL-CM4_piControl"]=0650
START["GFDL-CM4_abrupt-4xCO2"]=0101
END["GFDL-CM4_abrupt-4xCO2"]=0150

START["HadGEM3-GC31-MM_piControl"]=2300
END["HadGEM3-GC31-MM_piControl"]=2349
START["HadGEM3-GC31-MM_abrupt-4xCO2"]=1950
END["HadGEM3-GC31-MM_abrupt-4xCO2"]=1999

# Output folder
OUTDIR="TIW_SST_STD_outputs"
mkdir -p $OUTDIR

# Loop over all model/scenario pairs
for MODEL in "${MODELS[@]}"; do
  for SCEN in "${SCENARIOS[@]}"; do
    KEY="${MODEL}_${SCEN}"
    Y1=${START[$KEY]}
    Y2=${END[$KEY]}

    INPUT="tos_Oday_${MODEL}_${SCEN}_r1i1p1f1_gn_${Y1}-${Y2}_equator.nc"
    if [[ ! -f $INPUT ]]; then
      echo " Input file not found: $INPUT"
      continue
    fi

    echo "Processing: $MODEL | $SCEN | $Y1–$Y2"

    # Step 1: Select years
    cdo selyear,$Y1/$Y2 $INPUT temp1.nc

    # Step 2: Replace missing values with 0
    cdo setmisstoc,0 temp1.nc temp2.nc
    rm -f temp1.nc

    # Step 3: Apply 10–50 day bandpass filter
    FFT_OUT="$OUTDIR/tos_Oday_${MODEL}_${SCEN}_r1i1p1f1_gn_${Y1}-${Y2}_equator_10-50fft.nc"
    cdo bandpass,7.3,36.5 -del29feb temp2.nc $FFT_OUT
    rm -f temp2.nc

    # Step 4: Compute monthly std dev (TIW intensity)
    STD_OUT="$OUTDIR/tos_Oday_${MODEL}_${SCEN}_r1i1p1f1_gn_${Y1}-${Y2}_equator_10-50fft_std.nc"
    cdo monstd $FFT_OUT $STD_OUT

    echo  Output saved: $STD_OUT"
    echo "-----------------------------------"
  done
done

echo "All 4 HighResMIP models processed!"
