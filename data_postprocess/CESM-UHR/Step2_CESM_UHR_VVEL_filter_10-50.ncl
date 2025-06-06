load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 

begin
 dir1 = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/PD/VVEL/"
 dir2 = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/2xCO2/VVEL/"
 dir3 = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/4xCO2/VVEL/"

 fil1 = systemfunc("ls " + dir1+"VVEL.day.0111-0140.10S_10N.180_80W.levidx01_06_025degree_anom.nc")
 fil1 = systemfunc("ls " + dir2+"VVEL.day.0141-0170.10S_10N.180_80W.levidx01_06_025degree_anom.nc") ;2xCO2
 fil1 = systemfunc("ls " + dir3+"VVEL.day.0141-0170.10S_10N.180_80W.levidx01_06_025degree_anom.nc") ;4xCO2

 f1   = addfile(fil1,"r")
 ta = f1->VVEL
   printVarSummary(ta)
   ;dims = dimsizes(ta)
   ;ntim = dims(0)
   ;ta!0="time"
   ;ta&time=ispan(0,ntim-1,1)
   ta_time=ta(z_t|:,lat|:,lon|:,time|:)
  ;===================================================================
   fca     = 1./10.
   fcb     = 1./50.                                              ; end freq
   cf      = ezfftf (ta_time)                      ; sub period (2, ntim/2)
   fcf     = fspan(0, 0.5, dimsizes(ta_time&time)/2) 
   ifcf    = ind(fcf.lt.fcb.or.fcf.gt.fca)
   cf(:,:,:,:,ifcf) = 0.0                                              ; set coef to zero
   tanom1   = ezfftb (cf, cf@xbar)
   copy_VarMeta(ta_time,tanom1)
   ta1=tanom1(time|:,z_t|:,lat|:,lon|:)
   printVarSummary(ta1)
   delete(cf)
  ;===================================================================
     fca1     = 1./10.
     fcb1     = 1./60.                                              ; end freq
     cf      = ezfftf (ta_time)                      ; sub period (2, ntim/2)
     fcf     = fspan(0, 0.5, dimsizes(ta_time&time)/2) 
     ifcf1   = ind(fcf.lt.fcb1.or.fcf.gt.fca1)
     cf(:,:,:,:,ifcf1) = 0.0                                              ; set coef to zero
     tanom2   = ezfftb (cf, cf@xbar)
     copy_VarMeta(ta_time,tanom2)
     ta2=tanom2(time|:,z_t|:,lat|:,lon|:)
     printVarSummary(ta2)
     delete(cf)
;=====================================================================
   delete(ta_time)
   delete(ta)
   delete([/tanom1, tanom2/])
;   fca1     = 1./20.
;   fcb1     = 1./40.                                              ; end freq
;   cf1      = ezfftf (ta_time)                      ; sub period (2, ntim/2)
;   fcf1     = fspan(0, 0.5, dimsizes(ta_time&time)/2) 
;   ifcf1    = ind(fcf1.lt.fcb1.or.fcf1.gt.fca1)
;   cf1(:,:,:,ifcf1) = 0.0                                              ; set coef to zero
;   tanom1   = ezfftb (cf1, cf1@xbar)
;   copy_VarMeta(ta_time,tanom1)
;   ta2=tanom1(time|:,lat|:,lon|:)
;  
     fout=addfile("/glade/scratch/aoyunxue/TIW_data/CESM_HR/PD/VVEL/CESM.PD.day.VVEL.anom.10-50.fft.nc","c")
     fout->va=ta1
     fout=addfile("/glade/scratch/aoyunxue/TIW_data/CESM_HR/2xCO2/VVEL/CESM.2xCO2.day.VVEL.anom.10-50.fft.nc","c")
     fout->va=ta1
     fout=addfile("/glade/scratch/aoyunxue/TIW_data/CESM_HR/4xCO2/VVEL/CESM.4xCO2.day.VVEL.anom.10-50.fft.nc","c")
     fout->va=ta1

 end
