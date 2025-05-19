load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 

begin

;========================================================================================================
 datadir1  = "/glade/derecho/scratch/aoyunxue/reanalysis/CMEMS/anom/"
 fils1     = systemfunc("ls " + datadir1+ "GLOBAL-REANALYSIS-PHY-001-030-DAILY_thetao_surface_anom_10S-10N.nc")
 f1        = addfile(fils1, "r")
 thetao        = f1->thetao(:,:,:)
 printVarSummary(thetao)
    thetao!0   = "time"
    thetao!1   = "lat"
    thetao!2   = "lon"
    thetao&time= f1->time
    thetao&lat = f1->latitude
    thetao&lon = f1->longitude
 ;=======================================================================================================
  ta_time = thetao(lat|:,lon|:,time|:)
    fca   = 1./10.
    fcb   = 1./50.  
    cf1   = ezfftf(ta_time)                            ; sub period (2, ntim/2)
   fcf1   = fspan(0, 0.5, dimsizes(ta_time&time)/2)    ; 
   ifcf1  = ind(fcf1.lt.fcb.or.fcf1.gt.fca)            ; badpass filtering 
   cf1(:,:,:,ifcf1) = 0.0                              ; set coef to zero
   tanom  = ezfftb(cf1, cf1@xbar)
   copy_VarMeta(ta_time,tanom)
   delete(ta_time)
   delete(thetao)
   ta  = tanom(time|:,lat|:,lon|:)
   printVarSummary(ta)
   delete(tanom)
;========================================================================================================
  outdir1   = "/glade/derecho/scratch/aoyunxue/reanalysis/CMEMS/anom/"
  out_path1 = outdir1 + "GLOBAL-REANALYSIS-PHY-001-030-DAILY_thetao_surface_anom_10-50_fft.nc" 
      if(isfilepresent(out_path1)) then
         system("rm -f " + out_path1)
      end if
    fout1 = addfile(out_path1, "c")  
    fout1->ta  = ta 
;========================================================================================================
end 

