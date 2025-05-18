load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"  
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl" 
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl" 

begin
  dir1 = "/glade/scratch/aoyunxue/TIW_data/CESM2_HR/PD/density/"
; dir2 = "/glade/scratch/aoyunxue/TIW_data/CESM2_HR/2xCO2/density/"
 ;dir3 = "/glade/scratch/aoyunxue/TIW_data/CESM2_HR/4xCO2/density/"

 fil1 = systemfunc("ls " + dir1+"density.day.0111-0140_equator_025degree_10S_10N_anom.nc")
; fil1 = systemfunc("ls " + dir2+"density.day.0141-0170_equator_025degree_10S_10N_anom.nc")
; fil1 = systemfunc("ls " + dir3+"density.day.0141-0170_equator_025degree_10S_10N_anom.nc")

 f1   = addfile(fil1,"r")
 rhoa = f1->density
   printVarSummary(rhoa)
   dims = dimsizes(rhoa)
   ;ntim = dims(0)
   ;rhoa!0="time"
   ;rhoa&time=ispan(0,ntim-1,1)
   rhoa_time=rhoa(lat|:,lon|:,time|:)
  ;===================================================================
   fca     = 1./10.
   fcb     = 1./50.                                              ; end freq
   cf      = ezfftf (rhoa_time)                      ; sub period (2, ntim/2)
   fcf     = fspan(0, 0.5, dimsizes(rhoa_time&time)/2) 
   ifcf    = ind(fcf.lt.fcb.or.fcf.gt.fca)
   cf(:,:,:,ifcf) = 0.0                                              ; set coef to zero
   rhoanom1   = ezfftb (cf, cf@xbar)
   copy_VarMeta(rhoa_time,rhoanom1)
   rhoa1=rhoanom1(time|:,lat|:,lon|:)
   printVarSummary(rhoa1)
   delete(cf)
  ;===================================================================
     fca1     = 1./10.
     fcb1     = 1./60.                                              ; end freq
     cf      = ezfftf (rhoa_time)                      ; sub period (2, ntim/2)
     fcf     = fspan(0, 0.5, dimsizes(rhoa_time&time)/2) 
     ifcf1   = ind(fcf.lt.fcb1.or.fcf.gt.fca1)
     cf(:,:,:,ifcf1) = 0.0                                              ; set coef to zero
     rhoanom2   = ezfftb (cf, cf@xbar)
     copy_VarMeta(rhoa_time,rhoanom2)
     rhoa2=rhoanom2(time|:,lat|:,lon|:)
     printVarSummary(rhoa2)
     delete(cf)
;=====================================================================
   delete(rhoa_time)
   delete(rhoa)
   delete([/rhoanom1, rhoanom2/])
;   fca1     = 1./20.
;   fcb1     = 1./40.                                              ; end freq
;   cf1      = ezfftf (rhoa_time)                      ; sub period (2, ntim/2)
;   fcf1     = fspan(0, 0.5, dimsizes(rhoa_time&time)/2) 
;   ifcf1    = ind(fcf1.lt.fcb1.or.fcf1.gt.fca1)
;   cf1(:,:,:,ifcf1) = 0.0                                              ; set coef to zero
;   rhoanom1   = ezfftb (cf1, cf1@xbar)
;   copy_VarMeta(rhoa_time,rhoanom1)
;   rhoa2=rhoanom1(time|:,lat|:,lon|:)
;  
     fout=addfile("/glade/scratch/aoyunxue/TIW_data/CESM2_HR/PD/density/CESM2.PD.day.10S_10N_density.anom.10-50.fft.nc","c")
     fout->rhoa=rhoa1
 ;    fout=addfile("/glade/scratch/aoyunxue/TIW_data/CESM2_HR/2xCO2/density/CESM2.2xCO2.day.10S_10N_density.anom.10-50.fft.nc","c")
 ;    fout->rhoa=rhoa1  
 ;    fout=addfile("/glade/scratch/aoyunxue/TIW_data/CESM2_HR/4xCO2/density/CESM2.4xCO2.day.10S_10N_density.anom.10-50.fft.nc","c")
 ;    fout->rhoa=rhoa1

 end
