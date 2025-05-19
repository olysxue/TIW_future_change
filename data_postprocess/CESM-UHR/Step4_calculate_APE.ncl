load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
; read referrence lat lon 
datadir0 = "/glade/derecho/scratch/aoyunxue/FROM_CHEYENNE/TIW_data/CESM_HR/PD/TEMP/"
fil0 = systemfunc("ls " + datadir0 + "TEMP.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree.nc")
f0   = addfile(fil0,"r")
TEMP = f0->TEMP
;*******************************************
;  parameter set
;varname = "EKE"
 dens= 1025   ; kg/m3 reference density of seawater
 nave = 31
 g = 9.8 ; unit m/s2
;*******************************************
; PD 
; read TIW u & v 

 datadir1 = "/glade/derecho/scratch/aoyunxue/FROM_CHEYENNE/TIW_data/CESM_HR/PD/WVEL/"
 datadir2 = "/glade/derecho/scratch/aoyunxue/FROM_CHEYENNE/TIW_data/CESM_HR/2xCO2/WVEL/"
 datadir3 = "/glade/derecho/scratch/aoyunxue/FROM_CHEYENNE/TIW_data/CESM_HR/4xCO2/WVEL/"

   fil1  = systemfunc("ls "+datadir1 + "CESM.PD.day.WVEL.01-22.anom.10-50.fft.nc")
   fil2  = systemfunc("ls "+datadir2 + "CESM.2xCO2.day.WVEL.01-22.anom.10-50.fft.nc")
   fil3  = systemfunc("ls "+datadir3 + "CESM.4xCO2.day.WVEL.01-22.anom.10-50.fft.nc")
 
   f1    = addfile(fil1, "r")
   f2    = addfile(fil2, "r")
   f3    = addfile(fil3, "r")

   wa0_PD     = f1->wa(:,0:6,:,:) 
   wa0_2xCO2  = f2->wa(:,0:6,:,:)
   wa0_4xCO2  = f3->wa(:,0:6,:,:)
   dsizes = dimsizes(wa0_PD)
   ntim   = dsizes(0)
   nlev   = dsizes(1)
   nlat   = dsizes(2) 
   nlon   = dsizes(3)
    wa_PD = new((/ntim,nlev-1, nlat, nlon/),typeof(wa0_PD))
    wa_2xCO2 = new((/ntim,nlev-1, nlat, nlon/),typeof(wa0_2xCO2))
    wa_4xCO2 = new((/ntim,nlev-1, nlat, nlon/),typeof(wa0_4xCO2))
 
do i = 0, nlev-2
     wa_PD(:,i,:,:) = (wa0_PD(:,i,:,:) + wa0_PD(:,i+1,:,:))/2.
     wa_2xCO2(:,i,:,:) = (wa0_2xCO2(:,i,:,:) + wa0_2xCO2(:,i+1,:,:))/2.
     wa_4xCO2(:,i,:,:) = (wa0_4xCO2(:,i,:,:) + wa0_4xCO2(:,i+1,:,:))/2.
end do 
copy_VarMeta(TEMP, wa_PD)
copy_VarMeta(TEMP, wa_2xCO2)
copy_VarMeta(TEMP, wa_4xCO2)
;------------------------------------------
; read density anom
  ;=======================================
    datadir01 = "/glade/derecho/scratch/aoyunxue/FROM_CHEYENNE/TIW_data/CESM_HR/PD/density/"
    datadir02 = "/glade/derecho/scratch/aoyunxue/FROM_CHEYENNE/TIW_data/CESM_HR/2xCO2/density/"
    datadir03 = "/glade/derecho/scratch/aoyunxue/FROM_CHEYENNE/TIW_data/CESM_HR/4xCO2/density/"
   fil01  = systemfunc("ls "+datadir01 + "CESM.PD.day.density.anom.10-50.fft.nc")
   fil02  = systemfunc("ls "+datadir02 + "CESM.2xCO2.day.density.anom.10-50.fft.nc")
   fil03  = systemfunc("ls "+datadir03 + "CESM.4xCO2.day.density.anom.10-50.fft.nc")
 
   f01    = addfile(fil01, "r")
   f02    = addfile(fil02, "r")
   f03    = addfile(fil03, "r")

   ta_PD     = f01->density
   ta_2xCO2  = f02->density
   ta_4xCO2  = f03->density
 copy_VarMeta(TEMP, ta_PD)
 copy_VarMeta(TEMP, ta_2xCO2)
 copy_VarMeta(TEMP, ta_4xCO2)
   printVarSummary(ta_PD)   


;=====================================read ave======================================
; parameters set 
  g = 9.8 ; m/s2
  BCR0_PD = -1. * g  * ta_PD * wa_PD
  BCR0_2xCO2 =-1.* g  * ta_2xCO2 * wa_2xCO2
  BCR0_4xCO2 = -1.* g  * ta_4xCO2 * wa_4xCO2
  BCR_PD = runave_n_Wrap(BCR0_PD, 61,0,0)
  BCR_2xCO2 = runave_n_Wrap(BCR0_2xCO2, 61,0,0)
  BCR_4xCO2 = runave_n_Wrap(BCR0_4xCO2, 61,0,0)
  copy_VarMeta(ta_PD, BCR_PD)
  copy_VarMeta(ta_2xCO2, BCR_2xCO2)
  copy_VarMeta(ta_4xCO2, BCR_4xCO2)
;=====================================================================================
; output APE

 ;  opdir1    = "/Users/xueaoyun/Documents/work/TIW_future/data/energy/"
    opdir1    = "/glade/derecho/scratch/aoyunxue/FROM_CHEYENNE/TIW_data/CESM_HR/energy/"

    fname1    = "TIW_APE_PD.nc"
    system("/bin/rm -f "+ opdir1 + fname1)
    f1      = addfile(opdir1 + fname1,"c")    ; write netCDF file
    f1->APE = BCR_PD 
  
    fname2    = "TIW_APE_2xCO2.nc"
    system("/bin/rm -f "+ opdir1 + fname2)
    f2      = addfile(opdir1 + fname2,"c")    ; write netCDF file
    f2->APE = BCR_2xCO2 

    fname3    = "TIW_APE_4xCO2.nc"
    system("/bin/rm -f "+ opdir1 + fname3)
    f3      = addfile(opdir1 + fname3,"c")    ; write netCDF file
    f3->APE = BCR_4xCO2 

end 
