load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 

begin
;========================================================================================================
 datadir1  = "/glade/derecho/scratch/aoyunxue/reanalysis/CMEMS/raw/vo/depth/"
 fils1     = systemfunc("ls " + datadir1+ "cmems_vo_180.00W-80.00W_10.00S-10.00N_0-180m_1993-2023_025degree_anom.nc")
 f1        = addfile(fils1, "r")
 vo        = f1->vo(:,:,:,:)
 printVarSummary(vo)
    vo!0   = "time"
    vo!1   = "depth"
    vo!2   = "lat"
    vo!3   = "lon"
    vo&time= f1->time
    vo&depth = f1->depth
    vo&lat = f1->lat
    vo&lon = f1->lon
 ;=======================================================================================================
  va_time = vo(depth|:,lat|:,lon|:,time|:)
    fca   = 1./10.
    fcb   = 1./50.  
    cf1   = ezfftf(va_time)                            ; sub period (2, ntim/2)
   fcf1   = fspan(0, 0.5, dimsizes(va_time&time)/2)    ; 
   ifcf1  = ind(fcf1.lt.fcb.or.fcf1.gt.fca)            ; badpass filtering 
   cf1(:,:,:,:,ifcf1) = 0.0                              ; set coef to zero
   vanom  = ezfftb(cf1, cf1@xbar)
   copy_VarMeta(va_time,vanom)
   delete(va_time)
   delete(vo)
   va  = vanom(time|:,depth|:,lat|:,lon|:)
   printVarSummary(va)
   delete(vanom)
;========================================================================================================
  outdir1   = "/glade/derecho/scratch/aoyunxue/reanalysis/CMEMS/raw/vo/depth/"
  out_path1 = outdir1 + "cmems_vo_180.00W-80.00W_10.00S-10.00N_0-180m_1993-2023_025degree_anom_10-50_fft.nc" 
      if(isfilepresent(out_path1)) then
         system("rm -f " + out_path1)
      end if
    fout1 = addfile(out_path1, "c")  
    fout1->va  = va
    fout1->depth = vo&depth
    fout1->lat = vo&lat
    fout1->lon = vo&lon
;========================================================================================================
end 

