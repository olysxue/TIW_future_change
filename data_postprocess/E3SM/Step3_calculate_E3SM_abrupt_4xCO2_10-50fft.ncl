 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
 
 begin
 ; read tos
   dir1 = "/glade/derecho/scratch/aoyunxue/E3SMv2/abrupt4xCO2/daily/"
  fil0  = systemfunc("ls " + dir1+"sst_abrupt_4xCO2_0001-0150_remap_equator_new.nc")
   f0   = addfile(fil0,"r")
    sst = f0->sst
    
  dir1 = "/glade/derecho/scratch/aoyunxue/E3SMv2/abrupt4xCO2/daily/"
  fil1 = systemfunc("ls " + dir1+"sst_abrupt_4xCO2_0001-0150_remap_equator_moving_20yr_detrend.nc")
   f1   = addfile(fil1,"r")
    ta = f1->sst
    copy_VarMeta(sst,ta) 
    printVarSummary(ta)
    ta_time=ta(lat|:,lon|:,time|:)
     ;===================================================================
      fca     = 1./2.
      fcb     = 1./10.                                              ; end freq
      cf      = ezfftf (ta_time)                      ; sub period (2, ntim/2)
      fcf     = fspan(0, 0.5, dimsizes(ta_time&time)/2)
      ifcf    = ind(fcf.lt.fcb.or.fcf.gt.fca)
      cf(:,:,:,ifcf) = 0.0                                              ; set coef to zero
      tanom1   = ezfftb (cf, cf@xbar)
      copy_VarMeta(ta_time,tanom1)
      ta1=tanom1(time|:,lat|:,lon|:)                                                                                                                                                        
      printVarSummary(ta1)
      delete(cf)
      delete(ta_time)
      delete(tanom1)
     ;===================================================================
     ; read zos
   dir2 = "/glade/derecho/scratch/aoyunxue/E3SMv2/abrupt4xCO2/daily/"
   fil2 = systemfunc("ls " + dir2 +"ssh_abrupt_4xCO2_0001-0150_remap_equator.nc")
   f2   = addfile(fil2,"r")
   ssh = f2->ssh  

   dir3 = "/glade/derecho/scratch/aoyunxue/E3SMv2/abrupt4xCO2/daily/"
   fil3 = systemfunc("ls " + dir3 +"ssh_abrupt_4xCO2_0001-0150_remap_equator_moving_20yr_detrend.nc")
   f3   = addfile(fil3,"r")
    za  = f3->ssh  
   copy_VarMeta(ssh,za)
   printVarSummary(za)
   za_time=za(lat|:,lon|:,time|:)

     fca     = 1./2.
     fcb     = 1./10.                                              ; end freq
     cf      = ezfftf(za_time)                      ; sub period (2, ntim/2)
     fcf     = fspan(0, 0.5, dimsizes(za_time&time)/2)
     ifcf    = ind(fcf.lt.fcb.or.fcf.gt.fca)
     cf(:,:,:,ifcf) = 0.0                                              ; set coef to zero
     zanom1   = ezfftb (cf, cf@xbar)
     copy_VarMeta(za_time,zanom1)
     za1 = zanom1(time|:,lat|:,lon|:)
     printVarSummary(za1)
     delete(cf)
     delete(za_time)
     delete(zanom1)
;====================================================================

fout=addfile("/glade/derecho/scratch/aoyunxue/E3SMv2/abrupt4xCO2/daily/sst_abrupt_4xCO2_0001-0150_remap_equator_moving_20yr_detrend_10-50fft.nc","c")
fout->sst = ta1
  
fout=addfile("/glade/derecho/scratch/aoyunxue/E3SMv2/abrupt4xCO2/daily/ssh_abrupt_4xCO2_0001-0150_remap_equator_moving_20yr_detrend_10-50fft.nc","c")
fout->ssh = za1
end 
