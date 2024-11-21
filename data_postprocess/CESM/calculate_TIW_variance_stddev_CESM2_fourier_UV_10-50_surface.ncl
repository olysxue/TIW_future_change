load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "/glade/u/home/aoyunxue/libfun/my_utilities.ncl"

begin
dens = 1025
;==============================read data====================================
 datadir1  = "/glade/scratch/aoyunxue/TIW_data/CESM2_HR/PD/UVEL/"
 fil1      = systemfunc("ls " + datadir1 + "CESM2.PD.day.UVEL_surface.anom.10-50.fft.nc")
 f1        = addfile(fil1, "r")
 ua        = f1->ua
 printVarSummary(ua)

 datadir2  = "/glade/scratch/aoyunxue/TIW_data/CESM2_HR/PD/VVEL/"
 fil2      = systemfunc("ls " + datadir2 + "CESM2.PD.day.VVEL_surface.anom.10-50.fft.nc")
 f2        = addfile(fil2, "r")
 va        = f2->va
 printVarSummary(va)

;=====================calculate tiw stddev==================
  ua_std_CESM2 = calculate_monthly_values(ua, "std", 0, False)
  printVarSummary(ua_std_CESM2)
  ua_var_CESM2 = calculate_monthly_values(ua, "var", 0, False)
  printVarSummary(ua_var_CESM2)

  va_std_CESM2 = calculate_monthly_values(va, "std", 0, False)
  printVarSummary(va_std_CESM2)
  va_var_CESM2 = calculate_monthly_values(va, "var", 0, False)
  printVarSummary(va_var_CESM2)

   EKE_CESM2 = kinetic_energy(ua,va,dens)
   EKE_CESM2_mon = calculate_monthly_values(EKE_CESM2, "avg", 0, False)
;============================================================
  outdir     = "/glade/scratch/aoyunxue/TIW_data/CESM2_HR/PD/"

    out_path1 = outdir + "CESM2_HR_TIW_rms_var_ua_10-50_fft_surface_PD.nc" 
      if(isfilepresent(out_path1)) then
         system("rm -f " + out_path1)
      end if
    fout1 = addfile(out_path1, "c")  
    fout1->tiw_std  = ua_std_CESM2
    fout1->tiw_var  = ua_var_CESM2

    out_path2 = outdir + "CESM2_HR_TIW_rms_var_va_10-50_fft_surface_PD.nc"
      if(isfilepresent(out_path2)) then
         system("rm -f " + out_path2)
      end if
    fout2 = addfile(out_path2, "c")
    fout2->tiw_std  = va_std_CESM2
    fout2->tiw_var  = va_var_CESM2

    out_path3 = outdir + "CESM2_TIW_EKE_10-50_fft_surface_PD.nc"
      if(isfilepresent(out_path3)) then
         system("rm -f " + out_path3)
      end if
    fout3 = addfile(out_path3, "c")
    fout3->EKE = EKE_CESM2_mon

delete([/ ua_std_CESM2, ua_var_CESM2, va_std_CESM2, va_var_CESM2, EKE_CESM2/])
end 


