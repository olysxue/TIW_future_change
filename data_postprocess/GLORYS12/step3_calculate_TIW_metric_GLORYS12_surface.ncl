load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/glade/u/home/aoyunxue/libfun/my_utilities.ncl"
begin
;read  CMEMS 
;REFERRENCE LAT AND LON
dens = 1025 
datadir   = "/glade/derecho/scratch/aoyunxue/reanalysis/CMEMS/"
fil       = systemfunc("ls " + datadir + "GLOBAL-REANALYSIS-PHY-001-030-DAILY_vo_sample.nc")
f         = addfile(fil, "r")
lat       = f->latitude({-10:10})
lon       = f->longitude({-180:-60})
;==============================================================================================================
datadir2  =  "/glade/derecho/scratch/aoyunxue/reanalysis/CMEMS/anom/"
fil2      = systemfunc("ls " + datadir2 + "GLOBAL-REANALYSIS-PHY-001-030-DAILY_thetao_surface_anom_10-50_fft.nc")
 f2       = addfile(fil2, "r")
 ta       = f2->ta
 printVarSummary(ta)
 ta!0 = "time"
 ta&time=ispan(0,9860,1)
 ta&time@units = "days since 1993-01-01 00:00:00"
 ta!1 = "lat"
 ta!2 = "lon"
 ta&lat = lat
 ta&lon = lon 
 printVarSummary(ta)
;#######################################################
 datadir02  =  "/glade/derecho/scratch/aoyunxue/reanalysis/CMEMS/anom/"
 fil02      = systemfunc("ls " + datadir02 + "GLOBAL-REANALYSIS-PHY-001-030-DAILY_uo_surface_anom_10-50_fft.nc")
 f02       = addfile(fil02, "r")
 ua        = f02->ua
 printVarSummary(ua)
 ua!0 = "time"
 ua&time=ispan(0,9860,1)
 ua&time@units = "days since 1993-01-01 00:00:00"
 ua!1 = "lat"
 ua!2 = "lon"
 ua&lat = lat
 ua&lon = lon 
 printVarSummary(ua)
;######################################################
 datadir002  =  "/glade/derecho/scratch/aoyunxue/reanalysis/CMEMS/anom/"
 fil002      = systemfunc("ls " + datadir002 + "GLOBAL-REANALYSIS-PHY-001-030-DAILY_vo_surface_anom_10-50_fft.nc")
 f002       = addfile(fil002, "r")
 va        = f002->va
 printVarSummary(va)
 va!0 = "time"
 va&time=ispan(0,9860,1)
 va&time@units = "days since 1993-01-01 00:00:00"
 va!1 = "lat"
 va!2 = "lon"
 va&lat = lat
 va&lon = lon 
 printVarSummary(va)
;######################################################

  ua_std_CMEMS = calculate_monthly_values(ua, "std", 0, False)
  printVarSummary(ua_std_CMEMS)
  ua_var_CMEMS = calculate_monthly_values(ua, "var", 0, False)
  printVarSummary(ua_var_CMEMS)

  va_std_CMEMS = calculate_monthly_values(va, "std", 0, False)
  printVarSummary(va_std_CMEMS)
  va_var_CMEMS = calculate_monthly_values(va, "var", 0, False)
  printVarSummary(va_var_CMEMS)

  ta_std_CMEMS = calculate_monthly_values(ta, "std", 0, False)
  printVarSummary(ta_std_CMEMS)
  ta_var_CMEMS = calculate_monthly_values(ta, "var", 0, False)
  printVarSummary(ta_var_CMEMS)
   
  EKE_CMEMS = kinetic_energy(ua,va,dens)
  EKE_CMEMS_mon = calculate_monthly_values(EKE_CMEMS, "avg", 0, False)
;#####################################################
delete([/ua, va, ta/])
;--------------------------output-------------------------------
   outdir     = "/glade/derecho/scratch/aoyunxue/reanalysis/CMEMS/"

    out_path02 = outdir + "CMEMS_tiw_rms_var_ua_10-50_fft_surface.nc" 
      if(isfilepresent(out_path02)) then
         system("rm -f " + out_path02)
      end if
    fout02 = addfile(out_path02, "c")  
    fout02->tiw_std  = ua_std_CMEMS
    fout02->tiw_var  = ua_var_CMEMS

    out_path002 = outdir + "CMEMS_tiw_rms_var_va_10-50_fft_surface.nc"
      if(isfilepresent(out_path002)) then
         system("rm -f " + out_path002)
      end if
    fout002 = addfile(out_path002, "c")
    fout002->tiw_std  = va_std_CMEMS
    fout002->tiw_var  = va_var_CMEMS

    out_path0002 = outdir + "CMEMS_tiw_rms_var_ta_10-50_fft_surface.nc"
      if(isfilepresent(out_path0002)) then
         system("rm -f " + out_path0002)
      end if
    fout0002 = addfile(out_path0002, "c")
    fout0002->tiw_std  = ta_std_CMEMS
    fout0002->tiw_var  = ta_var_CMEMS

    out_path00002 = outdir + "CMEMS_EKE_10-50_fft_surface_mon.nc"
      if(isfilepresent(out_path00002)) then
         system("rm -f " + out_path00002)
      end if
    fout00002 = addfile(out_path00002, "c")
    fout00002->EKE = EKE_CMEMS_mon

delete([/ta_std_CMEMS, ta_var_CMEMS, ua_std_CMEMS, ua_var_CMEMS, va_std_CMEMS, va_var_CMEMS, EKE_CMEMS/])
end 

