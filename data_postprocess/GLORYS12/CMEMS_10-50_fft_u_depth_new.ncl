load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 

begin
;========================================================================================================
 datadir1  = "/glade/derecho/scratch/aoyunxue/reanalysis/CMEMS/raw/uo/depth/"
 fils1     = systemfunc("ls " + datadir1+ "cmems_uo_180.00W-80.00W_10.00S-10.00N_1993-2023_025degree_anom.nc")
 f1        = addfile(fils1, "r")
 uo        = f1->uo(:,:,:,:)
 printVarSummary(uo)
    uo!0   = "time"
    uo!1   = "depth"
    uo!2   = "lat"
    uo!3   = "lon"
    uo&time= f1->time
    uo&depth = f1->depth
    uo&lat = f1->lat
    uo&lon = f1->lon
 ;=======================================================================================================
  ua_time = uo(depth|:,lat|:,lon|:,time|:)
    fca   = 1./10.
    fcb   = 1./50.  
    cf1   = ezfftf(ua_time)                            ; sub period (2, ntim/2)
   fcf1   = fspan(0, 0.5, dimsizes(ua_time&time)/2)    ; 
   ifcf1  = ind(fcf1.lt.fcb.or.fcf1.gt.fca)            ; badpass filtering 
   cf1(:,:,:,:,ifcf1) = 0.0                              ; set coef to zero
   uanom  = ezfftb(cf1, cf1@xbar)
   copy_VarMeta(ua_time,uanom)
   delete(ua_time)
   delete(uo)
   ua  = uanom(time|:,depth|:,lat|:,lon|:)
   printVarSummary(ua)
   delete(uanom)
;========================================================================================================
  outdir1   = "/glade/derecho/scratch/aoyunxue/reanalysis/CMEMS/raw/uo/depth/"
  out_path1 = outdir1 + "cmems_uo_180.00W-80.00W_10.00S-10.00N_1993-2023_025degree_anom_10-50_fft.nc" 
      if(isfilepresent(out_path1)) then
         system("rm -f " + out_path1)
      end if
    fout1 = addfile(out_path1, "c")  
    fout1->ua  = ua 
    fout1->depth = uo&depth
    fout1->lat = uo&lat
    fout1->lon = uo&lon
;========================================================================================================
end 

