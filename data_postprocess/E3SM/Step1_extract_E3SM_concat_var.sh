#!/usr/bin/env bash

# locations of files
basedir="/glade/derecho/scratch/aoyunxue/CO2/E3SM/picontrol/"
# Script to post-process E3SM output

This repository contains some scripts to do manual CMORization of output from the E3SMv1 and v2 simulations. 

Contents of the repository:
- extract_concat_var.sh
  Shell script to extract daily ocean variables from annual output files, concatenate into a single file, rename according to CMOR conventions


# locations of output files
OUTPUT_DIR=/glade/derecho/scratch/aoyunxue/CO2/E3SM/picontrol/

# files to split
list=`ls ${basedir}/mpaso.hist.am.highFrequencyOutput*nc`

cd ${basedir}/

for hf in ${list}
do

   # Get portion of filename associated with the year
   froot=${hf%.nc}
   fnam=${froot##*Output.}
   echo ${fnam}

   # Extract SSH, SST information
   # (writing raw casename into the temporary filename so it will be retained in the netCDF metadata)
   echo ncks -v temperatureAtSurface ${hf} sst_picontrol_${fnam}.nc
   ncks -v temperatureAtSurface ${hf} sst_picontrol_${fnam}.nc

   echo ncks -v ssh ${hf} ssh_picontrol_${fnam}.nc
   ncks -v ssh ${hf} ssh_picontrol_${fnam}.nc

   echo ncks -v kineticEnergyAtSurface ${hf} KE_picontrol_${fnam}.nc
   ncks -v kineticEnergyAtSurface ${hf} KE_picontrol_${fnam}.nc

done

# Concatenate year files into a single file
ncrcat sst*.nc tos_day_E3SM-1-0_picontrol_gr_00010101-05001231.nc
ncrcat ssh*.nc zos_day_E3SM-1-0_picontrol_gr_00010101-05001231.nc
ncrcat KE*.nc KE_day_E3SM-1-0_picontrol_gr_00010101-05001231.nc

# Rename concatenated variables to CMORized names
ncrename -v temperatureAtSurface,tos tos_day_E3SM-1-0_picontrol_gr_00010101-05001231.nc
ncrename -v ssh,zos zos_day_E3SM-1-0_picontrol_gr_00010101-05001231.nc
ncrename -v kineticEnergyAtSurface,KE KE_day_E3SM-1-0_picontrol_gr_00010101-05001231.nc
# Make new directory if it doesn't already exist
mkdir -p /glade/derecho/scratch/aoyunxue/CO2/E3SM/picontrol/tos
mkdir -p /glade/derecho/scratch/aoyunxue/CO2/E3SM/picontrol/zos
mkdir -p /glade/derecho/scratch/aoyunxue/CO2/E3SM/picontrol/KE

# Copy the output into the new location
cp tos_day_E3SM-1-0_picontrol_gr_00010101-05001231.nc /glade/derecho/scratch/aoyunxue/CO2/E3SM/picontrol/tos
cp zos_day_E3SM-1-0_picontrol_gr_00010101-05001231.nc /glade/derecho/scratch/aoyunxue/CO2/E3SM/picontrol/zos
cp KE_day_E3SM-1-0_picontrol_gr_00010101-05001231.nc  /glade/derecho/scratch/aoyunxue/CO2/E3SM/picontrol/KE
# Remove temporary files
rm sst*.nc
rm ssh*.nc
rm KE*.nc
