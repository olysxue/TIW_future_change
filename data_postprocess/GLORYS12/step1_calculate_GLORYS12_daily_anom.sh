#!/bin/bash
# Author: Aoyun Xue
# Purpose: Compute daily anomalies for GLORYS12 variables (uo, vo, thetao, zos)
# Domain: 10°S–10°N, 180°W–80°W | Period: 1993–2023 | Resolution: 0.25°


# === uo ===
echo "  ➤ Processing uo..."
cdo -ydaymean cmems_uo_180.00W-80.00W_10.00S-10.00N_surface_1993-2023_025degree.nc \
               cmems_uo_180.00W-80.00W_10.00S-10.00N_surface_1993-2023_025degree_clm.nc
cdo -ydaysub  cmems_uo_180.00W-80.00W_10.00S-10.00N_surface_1993-2023_025degree.nc \
              cmems_uo_180.00W-80.00W_10.00S-10.00N_surface_1993-2023_025degree_clm.nc \
              cmems_uo_180.00W-80.00W_10.00S-10.00N_surface_1993-2023_025degree_anom.nc

# === vo ===
echo "  ➤ Processing vo..."
cdo -ydaymean cmems_vo_180.00W-80.00W_10.00S-10.00N_surface_1993-2023_025degree.nc \
               cmems_vo_180.00W-80.00W_10.00S-10.00N_surface_1993-2023_025degree_clm.nc
cdo -ydaysub  cmems_vo_180.00W-80.00W_10.00S-10.00N_surface_1993-2023_025degree.nc \
              cmems_vo_180.00W-80.00W_10.00S-10.00N_surface_1993-2023_025degree_clm.nc \
              cmems_vo_180.00W-80.00W_10.00S-10.00N_surface_1993-2023_025degree_anom.nc

# === thetao (SST) ===
echo "  ➤ Processing thetao (SST)..."
cdo -ydaymean cmems_thetao_180.00W-80.00W_10.00S-10.00N_surface_1993-2023_025degree.nc \
               cmems_thetao_180.00W-80.00W_10.00S-10.00N_surface_1993-2023_025degree_clm.nc
cdo -ydaysub  cmems_thetao_180.00W-80.00W_10.00S-10.00N_surface_1993-2023_025degree.nc \
              cmems_thetao_180.00W-80.00W_10.00S-10.00N_surface_1993-2023_025degree_clm.nc \
              cmems_thetao_180.00W-80.00W_10.00S-10.00N_surface_1993-2023_025degree_anom.nc

# === zos (SSH) ===
echo "  ➤ Processing zos (SSH)..."
cdo -ydaymean cmems_zos_180.00W-80.00W_10.00S-10.00N_1993-2023_025degree.nc \
               cmems_zos_180.00W-80.00W_10.00S-10.00N_1993-2023_025degree_clm.nc
cdo -ydaysub  cmems_zos_180.00W-80.00W_10.00S-10.00N_1993-2023_025degree.nc \
              cmems_zos_180.00W-80.00W_10.00S-10.00N_1993-2023_025degree_clm.nc \
              cmems_zos_180.00W-80.00W_10.00S-10.00_1993-2023_025degree_anom.nc
