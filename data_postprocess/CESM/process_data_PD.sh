cd /glade/scratch/aoyunxue/TIW_data/CESM2_HR/PD/UVEL/
cdo -sellevidx,1/6 UVEL.day.0111-0140.10S_10N.180_80W.levidx01_21.nc UVEL.day.0111-0140.10S_10N.180_80W.levidx01_06.nc
cdo remapbil,025degree_120-300.txt UVEL.day.0111-0140.10S_10N.180_80W.levidx01_06.nc UVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc
rm -rf UVEL.day.0111-0140.10S_10N.180_80W.levidx01_06.nc
# calculate anom
cdo -ydaymean UVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc UVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_clm.nc
cdo -ydaysub UVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc UVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_clm.nc UVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_anom.nc
# VVEL
cd /glade/scratch/aoyunxue/TIW_data/CESM2_HR/PD/VVEL/                                                         
cdo -sellevidx,1/6 VVEL.day.0111-0140.10S_10N.180_80W.levidx01_21.nc VVEL.day.0111-0140.10S_10N.180_80W.levidx01_06.nc
cdo remapbil,025degree_120-300.txt VVEL.day.0111-0140.10S_10N.180_80W.levidx01_06.nc VVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc
rm -rf VVEL.day.0111-0140.10S_10N.180_80W.levidx01_06.nc
cdo -ydaymean VVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc VVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_clm.nc
cdo -ydaysub VVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc VVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_clm.nc VVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_anom.nc
#temp
cd /glade/scratch/aoyunxue/TIW_data/CESM2_HR/PD/TEMP/
cdo -sellevidx,1/6 TEMP.day.0111-0140.10S_10N.180_80W.levidx01_21.nc TEMP.day.0111-0140.10S_10N.180_80W.levidx01_06.nc
cdo remapbil,025degree_120-300.txt TEMP.day.0111-0140.10S_10N.180_80W.levidx01_06.nc TEMP.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc
rm -rf TEMP.day.0111-0140.10S_10N.180_80W.levidx01_06.nc
cdo -ydaymean TEMP.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc TEMP.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_clm.nc
cdo -ydaysub TEMP.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc TEMP.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_clm.nc TEMP.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_anom.nc
