load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 

begin
;========================================================================================================
 datadir1  = "/glade/derecho/scratch/aoyunxue/reanalysis/CMEMS/raw/thetao/025degree/"
 fils1     = systemfunc("ls " + datadir1+ "cmems_thetao_180.00W-80.00W_10.00S-10.00N_0-180m_1993-2024_025degree_anom.nc")
 f1        = addfile(fils1, "r")
 thetao        = f1->thetao(:,:,:,:)
 printVarSummary(thetao)
    thetao!0   = "time"
    thetao!1   = "depth"
    thetao!2   = "lat"
    thetao!3   = "lon"
    thetao&time= f1->time
    thetao&depth = f1->depth
    thetao&lat = f1->lat
    thetao&lon = f1->lon
 ;=======================================================================================================
   ta_time = thetao(depth|:,lat|:,lon|:,time|:)
    fca   = 1./10.
    fcb   = 1./50.  
    cf1   = ezfftf(ta_time)                            ; sub period (2, ntim/2)
   fcf1   = fspan(0, 0.5, dimsizes(ta_time&time)/2)    ; 
   ifcf1  = ind(fcf1.lt.fcb.or.fcf1.gt.fca)            ; badpass filtering 
   cf1(:,:,:,:,ifcf1) = 0.0                              ; set coef to zero
   tanom  = ezfftb(cf1, cf1@xbar)
   copy_VarMeta(ta_time,tanom)
   delete(ta_time)
   delete(thetao)
   ta  = tanom(time|:,depth|:,lat|:,lon|:)
   printVarSummary(ta)
   delete(tanom)
;========================================================================================================
  outdir1   = "/glade/derecho/scratch/aoyunxue/reanalysis/CMEMS/raw/thetao/025degree/"
  out_path1 = outdir1 + "cmems_thetao_180.00W-80.00W_10.00S-10.00N_1993-2024_025degree_anom_10-50_fft.nc" 
      if(isfilepresent(out_path1)) then
         system("rm -f " + out_path1)
      end if
    fout1 = addfile(out_path1, "c")  
    fout1->thetao  = ta
    fout1->depth = thetao&depth
    fout1->lat = thetao&lat
    fout1->lon = thetao&lon
;========================================================================================================
end 

