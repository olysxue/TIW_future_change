  # UVEL
  cd /glade/scratch/aoyunxue/TIW_data/CESM2_HR/PD/VVEL/
  cdo -ydaymean UVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc UVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_clm.nc    
  cdo -ydaysub UVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc UVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_clm.nc UVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_anom.nc
  # VVEL
  cd /glade/scratch/aoyunxue/TIW_data/CESM2_HR/PD/VVEL/
  cdo -ydaymean VVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc VVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_clm.nc
  cdo -ydaysub VVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc VVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_clm.nc VVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_anom.nc
    # WVEL
  cd /glade/scratch/aoyunxue/TIW_data/CESM2_HR/PD/WVEL/
  cdo -ydaymean WVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc WVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_clm.nc
  cdo -ydaysub WVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc WVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_clm.nc WVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_anom.nc
  #temp
  cd /glade/scratch/aoyunxue/TIW_data/CESM2_HR/PD/temp/
  cdo -ydaymean TEMP.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc TEMP.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_clm.nc
  cdo -ydaysub TEMP.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc TEMP.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_clm.nc TEMP.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_anom.nc
  #density
  cd /glade/scratch/aoyunxue/TIW_data/CESM2_HR/PD/density/
  cdo -ydaymean density.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc density.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_clm.nc
  cdo -ydaysub density.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc density.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_clm.nc density.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_anom.nc
   #SST
  cd /glade/scratch/aoyunxue/TIW_data/CESM2_HR/PD/SST/
  cdo -ydaymean SST.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc SST.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_clm.nc
  cdo -ydaysub SST.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc SST.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_clm.nc SST.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_anom.nc
  #SSH
  cd /glade/scratch/aoyunxue/TIW_data/CESM2_HR/PD/SSH/
  cdo -ydaymean SSH.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc SSH.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_clm.nc
  cdo -ydaysub SSH.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc SSH.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_clm.nc SSH.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_anom.nc