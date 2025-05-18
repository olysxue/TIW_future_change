load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 

begin
 dir1 = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/PD/SST/"
; dir2 = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/2xCO2/SST/"
 ;dir3 = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/4xCO2/SST/"

 fil1 = systemfunc("ls " + dir1+"SST.day.0111-0140_equator_025degree_10S_10N_anom.nc")
; fil1 = systemfunc("ls " + dir2+"SST.day.0141-0170_equator_025degree_10S_10N_anom.nc")
; fil1 = systemfunc("ls " + dir3+"SST.day.0141-0170_equator_025degree_10S_10N_anom.nc")

 f1   = addfile(fil1,"r")
 ta = f1->SST
   printVarSummary(ta)
   dims = dimsizes(ta)
   ;ntim = dims(0)
   ;ta!0="time"
   ;ta&time=ispan(0,ntim-1,1)
   ta_time=ta(lat|:,lon|:,time|:)
  ;===================================================================
   fca     = 1./10.
   fcb     = 1./50.                                              ; end freq
   cf      = ezfftf (ta_time)                      ; sub period (2, ntim/2)
   fcf     = fspan(0, 0.5, dimsizes(ta_time&time)/2) 
   ifcf    = ind(fcf.lt.fcb.or.fcf.gt.fca)
   cf(:,:,:,ifcf) = 0.0                                              ; set coef to zero
   tanom1   = ezfftb (cf, cf@xbar)
   copy_VarMeta(ta_time,tanom1)
   ta1=tanom1(time|:,lat|:,lon|:)
   printVarSummary(ta1)
   delete(cf)
  ;===================================================================
     fca1     = 1./10.
     fcb1     = 1./60.                                              ; end freq
     cf      = ezfftf (ta_time)                      ; sub period (2, ntim/2)
     fcf     = fspan(0, 0.5, dimsizes(ta_time&time)/2) 
     ifcf1   = ind(fcf.lt.fcb1.or.fcf.gt.fca1)
     cf(:,:,:,ifcf1) = 0.0                                              ; set coef to zero
     tanom2   = ezfftb (cf, cf@xbar)
     copy_VarMeta(ta_time,tanom2)
     ta2=tanom2(time|:,lat|:,lon|:)
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
     fout=addfile("/glade/scratch/aoyunxue/TIW_data/CESM_HR/PD/SST/CESM.PD.day.10S_10N_SST.anom.10-50.fft.nc","c")
     fout->ta=ta1
 ;    fout=addfile("/glade/scratch/aoyunxue/TIW_data/CESM_HR/2xCO2/SST/CESM.2xCO2.day.10S_10N_SST.anom.10-50.fft.nc","c")
 ;    fout->ta=ta1  
 ;    fout=addfile("/glade/scratch/aoyunxue/TIW_data/CESM_HR/4xCO2/SST/CESM.4xCO2.day.10S_10N_SST.anom.10-50.fft.nc","c")
 ;    fout->ta=ta1

 end
