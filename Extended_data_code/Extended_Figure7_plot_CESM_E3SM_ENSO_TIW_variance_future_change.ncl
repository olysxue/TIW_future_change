load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

function get_title_position(plot)
    local vpy, vpx, vph
    begin
      getvalues plot
          "vpYF"       : vpy
          "vpXF"       : vpx
          "vpHeightF"  : vph
      end getvalues
      pos = (/vpx+0.01,vpy+0.01/)
      return(pos)
    end

    function get_title_position1(plot)
    local vpy, vpx, vph
    begin
      getvalues plot
          "vpYF"       : vpy
          "vpXF"       : vpx
          "vpHeightF"  : vph
      end getvalues
      pos = (/vpx+0.19,vpy+0.04/)
      return(pos)
    end

undef("readdata1")
function readdata1(ncdir, filelist, varname, scale)
local dims,temp,f
begin
   f= addfile(ncdir + filelist(0), "r")
   temp=f->$varname$(:,:,:)
   dims=getfilevardimsizes(f, varname)
   out=new((/2,dims(0),dims(1),dims(2)/), typeof(temp))
   
   do i=0,1
   f:= addfile(ncdir + filelist(i), "r")
   out(i,:,:,:)=f->$varname$(:,:,:)
   end do 
   out=out*scale
   copy_VarMeta(temp,out(0,:,:,:))
   return(out)
end

begin 
  ; read nino34 index 
  ncdir1  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/monthly/PD/"
  ncdir2  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/monthly/2xCO2/"
  ncdir3  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/monthly/4xCO2/"
  
  fil1    = systemfunc("ls " + ncdir1 + "TEMP.mon.0071-0140.20S_20N.120E_80W.levidx01_21_1degree_detrend_anom.nc")
  fil2    = systemfunc("ls " + ncdir2 + "TEMP.mon.0071-0170.20S_20N.120E_80W.levidx01_21_1degree_detrend_anom.nc")
  fil3    = systemfunc("ls " + ncdir3 + "TEMP.mon.0071-0170.20S_20N.120E_80W.levidx01_21_1degree_detrend_anom.nc")
    f1    = addfile(fil1, "r")
    f2    = addfile(fil2, "r")
    f3    = addfile(fil3, "r") 
    CESM_TEMP_PD     = f1->TEMP(:,0,:,:)
    CESM_TEMP_2xCO2  = f2->TEMP(:,0,:,:)
    CESM_TEMP_4xCO2  = f3->TEMP(:,0,:,:)
    printVarSummary(CESM_TEMP_PD)
    printVarSummary(CESM_TEMP_2xCO2)

    lat1   = f1->lat({-5:5})
    wgtx1  = latRegWgt(lat1, "float", 0)
    CESM_nino3_PD    = dim_standardize_n_Wrap(wgt_areaave_Wrap(CESM_TEMP_PD(480:839,{-5:5},{210:270}),wgtx1,1.0,0),0,0)  
    CESM_nino3_2xCO2 = dim_standardize_n_Wrap(wgt_areaave_Wrap(CESM_TEMP_2xCO2(840:1199,{-5:5},{210:270}),wgtx1,1.0,0),0,0)  
    CESM_nino3_4xCO2 = dim_standardize_n_Wrap(wgt_areaave_Wrap(CESM_TEMP_4xCO2(840:1199,{-5:5},{210:270}),wgtx1,1.0,0),0,0) 
    
    CESM_nino34_PD    = dim_standardize_n_Wrap(wgt_areaave_Wrap(CESM_TEMP_PD(480:839,{-5:5},{190:240}),wgtx1,1.0,0),0,0)  
    CESM_nino34_2xCO2 = dim_standardize_n_Wrap(wgt_areaave_Wrap(CESM_TEMP_2xCO2(840:1199,{-5:5},{190:240}),wgtx1,1.0,0),0,0)  
    CESM_nino34_4xCO2 = dim_standardize_n_Wrap(wgt_areaave_Wrap(CESM_TEMP_4xCO2(840:1199,{-5:5},{190:240}),wgtx1,1.0,0),0,0) 

    ntime = dimsizes(CESM_nino3_2xCO2)

    CESM_TEMP_PD_stat    = dim_stat4_n(CESM_TEMP_PD,0)
    CESM_TEMP_2xCO2_stat = dim_stat4_n(CESM_TEMP_2xCO2,0)
    CESM_TEMP_4xCO2_stat = dim_stat4_n(CESM_TEMP_4xCO2,0)
    CESM_TEMP_PD_stat!0  = "var"
    CESM_TEMP_PD_stat!1  = "lat"
    CESM_TEMP_PD_stat!2  = "lon"
    CESM_TEMP_PD_stat&var  = ispan(0,3,1)
    CESM_TEMP_PD_stat&lat  = CESM_TEMP_PD&lat
    CESM_TEMP_PD_stat&lon  = CESM_TEMP_PD&lon
    copy_VarMeta(CESM_TEMP_PD_stat,CESM_TEMP_2xCO2_stat)
    copy_VarMeta(CESM_TEMP_PD_stat,CESM_TEMP_4xCO2_stat)
    printVarSummary(CESM_TEMP_PD_stat)
    printVarSummary(CESM_TEMP_2xCO2_stat)
    printVarSummary(CESM_TEMP_4xCO2_stat)
 
   ; CESM_TEMP_PD_var    = dim_stddev_n_Wrap(CESM_TEMP_PD, 0)
   ; CESM_TEMP_2xCO2_var = dim_stddev_n_Wrap(CESM_TEMP_2xCO2, 0)
   ; CESM_TEMP_4xCO2_var = dim_stddev_n_Wrap(CESM_TEMP_4xCO2, 0)

    CESM_TEMP_PD_var    = CESM_TEMP_PD_stat(1,:,:)
    CESM_TEMP_2xCO2_var = CESM_TEMP_2xCO2_stat(1,:,:)
    CESM_TEMP_4xCO2_var = CESM_TEMP_4xCO2_stat(1,:,:)

    CESM_TEMP_PD_skew    = CESM_TEMP_PD_stat(2,:,:)
    CESM_TEMP_2xCO2_skew = CESM_TEMP_2xCO2_stat(2,:,:)
    CESM_TEMP_4xCO2_skew = CESM_TEMP_4xCO2_stat(2,:,:)
    
    CESM_pattern1 = smth9_Wrap(CESM_TEMP_PD_var, 0.50, 0.25, False)
    CESM_pattern2 = smth9_Wrap(CESM_TEMP_2xCO2_var, 0.50, 0.25, False)
    CESM_pattern3 = smth9_Wrap(CESM_TEMP_4xCO2_var, 0.50, 0.25, False)
 
    CESM_pattern01 = smth9_Wrap(CESM_TEMP_PD_skew, 0.50, 0.25, False)
    CESM_pattern02 = smth9_Wrap(CESM_TEMP_2xCO2_skew, 0.50, 0.25, False)
    CESM_pattern03 = smth9_Wrap(CESM_TEMP_4xCO2_skew, 0.50, 0.25, False)
    delete([/fil1, fil2, fil3/])
    delete([/f1, f2, f3/])
;=======================================================================
ncdir  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/mean/"
; read SST var
fil1    = systemfunc("ls " + ncdir + "CESM2_HR_SST_std_var_10-50_fft_PD.nc")
fil2    = systemfunc("ls " + ncdir + "CESM2_HR_SST_std_var_10-50_fft_2xCO2.nc")
fil3    = systemfunc("ls " + ncdir + "CESM2_HR_SST_std_var_10-50_fft_4xCO2.nc")
  f1    = addfile(fil1, "r")
  f2    = addfile(fil2, "r")
  f3    = addfile(fil3, "r") 
  tiw_SST_PD     = f1->tiw_std(:,:,:)
  tiw_SST_2xCO2  = f2->tiw_std(:,:,:)
  tiw_SST_4xCO2  = f3->tiw_std(:,:,:)
  printVarSummary(tiw_SST_PD)

  tiw_SST_PD_clm    = clmMonTLL(tiw_SST_PD)
  tiw_SST_2xCO2_clm = clmMonTLL(tiw_SST_2xCO2)
  tiw_SST_4xCO2_clm = clmMonTLL(tiw_SST_4xCO2)

  tiw_SST_PD_anom     = calcMonAnomTLL(tiw_SST_PD, tiw_SST_PD_clm)
  tiw_SST_2xCO2_anom  = calcMonAnomTLL(tiw_SST_2xCO2, tiw_SST_2xCO2_clm)
  tiw_SST_4xCO2_anom  = calcMonAnomTLL(tiw_SST_4xCO2, tiw_SST_4xCO2_clm)

  tiw_SST_PD_var    = dim_variance_n_Wrap(tiw_SST_PD_anom, 0)
  tiw_SST_2xCO2_var = dim_variance_n_Wrap(tiw_SST_2xCO2_anom, 0)
  tiw_SST_4xCO2_var = dim_variance_n_Wrap(tiw_SST_4xCO2_anom, 0)

  ; regression pattern
  tiw_SST_PD_variance  = new((/dimsizes(tiw_SST_PD_anom&lat),dimsizes(tiw_SST_PD_anom&lon)/),"float")   
  tiw_SST_PD_variance!0="lat"
  tiw_SST_PD_variance!1="lon"
  tiw_SST_PD_variance&lat = tiw_SST_PD_anom&lat
  tiw_SST_PD_variance&lon = tiw_SST_PD_anom&lon
  tiw_SST_PD_variance  = (/regCoef(CESM_nino3_PD,tiw_SST_PD_anom(lat|:,lon|:,time|:))/)

  tiw_SST_2xCO2_variance  = new((/dimsizes(tiw_SST_PD_anom&lat),dimsizes(tiw_SST_PD_anom&lon)/),"float")   
  tiw_SST_2xCO2_variance!0="lat"
  tiw_SST_2xCO2_variance!1="lon"
  tiw_SST_2xCO2_variance&lat = tiw_SST_2xCO2_anom&lat
  tiw_SST_2xCO2_variance&lon = tiw_SST_2xCO2_anom&lon
  tiw_SST_2xCO2_variance  = (/regCoef(CESM_nino3_2xCO2,tiw_SST_2xCO2_anom(lat|:,lon|:,time|:))/)

  tiw_SST_4xCO2_variance  = new((/dimsizes(tiw_SST_PD_anom&lat),dimsizes(tiw_SST_PD_anom&lon)/),"float")   
  tiw_SST_4xCO2_variance!0="lat"
  tiw_SST_4xCO2_variance!1="lon"
  tiw_SST_4xCO2_variance&lat = tiw_SST_4xCO2_anom&lat
  tiw_SST_4xCO2_variance&lon = tiw_SST_4xCO2_anom&lon
  tiw_SST_4xCO2_variance  = (/regCoef(CESM_nino3_4xCO2,tiw_SST_4xCO2_anom(lat|:,lon|:,time|:))/)

 ; tiw_SST_PD_variance     = dim_variance_n_Wrap(tiw_SST_PD_anom, 0)*100.
 ; tiw_SST_2xCO2_variance  = dim_variance_n_Wrap(tiw_SST_2xCO2_anom, 0)*100.
 ; tiw_SST_4xCO2_variance  = dim_variance_n_Wrap(tiw_SST_4xCO2_anom, 0)*100.

  copy_VarMeta(tiw_SST_PD(0,:,:), tiw_SST_PD_variance)
  copy_VarMeta(tiw_SST_2xCO2(0,:,:), tiw_SST_2xCO2_variance)
  copy_VarMeta(tiw_SST_4xCO2(0,:,:), tiw_SST_4xCO2_variance)

; calculate index 
CESM_SST_PD    = wgt_areaave_Wrap(tiw_SST_PD_anom(:,{0:5},{210:250}),1.0,1.0,0)  
CESM_SST_2xCO2 = wgt_areaave_Wrap(tiw_SST_2xCO2_anom(:,{0:5},{210:250}),1.0,1.0,0)  
CESM_SST_4xCO2 = wgt_areaave_Wrap(tiw_SST_4xCO2_anom(:,{0:5},{210:250}),1.0,1.0,0)  
;-----------------------------------------------------------------------------------
; read SSH var
fil1    = systemfunc("ls " + ncdir + "CESM2_HR_SSH_std_var_10-50_fft_PD.nc")
fil2    = systemfunc("ls " + ncdir + "CESM2_HR_SSH_std_var_10-50_fft_2xCO2.nc")
fil3    = systemfunc("ls " + ncdir + "CESM2_HR_SSH_std_var_10-50_fft_4xCO2.nc")
  f1    = addfile(fil1, "r")
  f2    = addfile(fil2, "r")
  f3    = addfile(fil3, "r") 
  tiw_SSH_PD     = f1->tiw_std(:,:,:)
  tiw_SSH_2xCO2  = f2->tiw_std(:,:,:)
  tiw_SSH_4xCO2  = f3->tiw_std(:,:,:)

  printVarSummary(tiw_SSH_PD)

  tiw_SSH_PD_clm    = clmMonTLL(tiw_SSH_PD)
  tiw_SSH_2xCO2_clm = clmMonTLL(tiw_SSH_2xCO2)
  tiw_SSH_4xCO2_clm = clmMonTLL(tiw_SSH_4xCO2)

  tiw_SSH_PD_anom     = calcMonAnomTLL(tiw_SSH_PD, tiw_SSH_PD_clm)
  tiw_SSH_2xCO2_anom  = calcMonAnomTLL(tiw_SSH_2xCO2, tiw_SSH_2xCO2_clm)
  tiw_SSH_4xCO2_anom  = calcMonAnomTLL(tiw_SSH_4xCO2, tiw_SSH_4xCO2_clm)
      
  CESM_SSH_PD    = wgt_areaave_Wrap(tiw_SSH_PD_anom(:,{3:8},{210:250}),1.0,1.0,0)  
  CESM_SSH_2xCO2 = wgt_areaave_Wrap(tiw_SSH_2xCO2_anom(:,{3:8},{210:250}),1.0,1.0,0)  
  CESM_SSH_4xCO2 = wgt_areaave_Wrap(tiw_SSH_4xCO2_anom(:,{3:8},{210:250}),1.0,1.0,0)  
  print(CESM_SSH_4xCO2)

  tiw_SSH_PD_var    = dim_variance_n_Wrap(tiw_SSH_PD_anom, 0)
  tiw_SSH_2xCO2_var = dim_variance_n_Wrap(tiw_SSH_2xCO2_anom, 0)
  tiw_SSH_4xCO2_var = dim_variance_n_Wrap(tiw_SSH_4xCO2_anom, 0)

; regression pattern
  tiw_SSH_PD_variance  = new((/dimsizes(tiw_SSH_PD_anom&lat),dimsizes(tiw_SSH_PD_anom&lon)/),"float")   
  tiw_SSH_PD_variance!0="lat"
  tiw_SSH_PD_variance!1="lon"
  tiw_SSH_PD_variance&lat = tiw_SSH_PD_anom&lat
  tiw_SSH_PD_variance&lon = tiw_SSH_PD_anom&lon
  tiw_SSH_PD_variance  = (/regCoef(CESM_nino3_PD,tiw_SSH_PD_anom(lat|:,lon|:,time|:))/)

  tiw_SSH_2xCO2_variance  = new((/dimsizes(tiw_SSH_PD_anom&lat),dimsizes(tiw_SSH_PD_anom&lon)/),"float")   
  tiw_SSH_2xCO2_variance!0="lat"
  tiw_SSH_2xCO2_variance!1="lon"
  tiw_SSH_2xCO2_variance&lat = tiw_SSH_PD_anom&lat
  tiw_SSH_2xCO2_variance&lon = tiw_SSH_PD_anom&lon
  tiw_SSH_2xCO2_variance  = (/regCoef(CESM_nino3_2xCO2,tiw_SSH_2xCO2_anom(lat|:,lon|:,time|:))/)

  tiw_SSH_4xCO2_variance  = new((/dimsizes(tiw_SSH_PD_anom&lat),dimsizes(tiw_SSH_PD_anom&lon)/),"float")   
  tiw_SSH_4xCO2_variance!0="lat"
  tiw_SSH_4xCO2_variance!1="lon"
  tiw_SSH_4xCO2_variance&lat = tiw_SSH_PD_anom&lat
  tiw_SSH_4xCO2_variance&lon = tiw_SSH_PD_anom&lon
  tiw_SSH_4xCO2_variance  = (/regCoef(CESM_nino3_4xCO2,tiw_SSH_4xCO2_anom(lat|:,lon|:,time|:))/)


  CESM_pattern001 = smth9_Wrap(tiw_SST_PD_variance, 0.50, 0.25, False)
  CESM_pattern002 = smth9_Wrap(tiw_SST_2xCO2_variance, 0.50, 0.25, False)
  CESM_pattern003 = smth9_Wrap(tiw_SST_4xCO2_variance, 0.50, 0.25, False)

  CESM_pattern0001 = smth9_Wrap(tiw_SSH_PD_variance, 0.50, 0.25, False)
  CESM_pattern0002 = smth9_Wrap(tiw_SSH_2xCO2_variance, 0.50, 0.25, False)
  CESM_pattern0003 = smth9_Wrap(tiw_SSH_4xCO2_variance, 0.50, 0.25, False)
  printMinMax(CESM_pattern001, 0)
  printMinMax(CESM_pattern003, 0)

  printMinMax(CESM_pattern0001, 0)
  printMinMax(CESM_pattern0003, 0)

  CESM_SSH_PD    = wgt_areaave_Wrap(tiw_SSH_PD_anom(:,{3:8},{210:250}),1.0,1.0,0)  
  CESM_SSH_2xCO2 = wgt_areaave_Wrap(tiw_SSH_2xCO2_anom(:,{3:8},{210:250}),1.0,1.0,0)  
  CESM_SSH_4xCO2 = wgt_areaave_Wrap(tiw_SSH_4xCO2_anom(:,{3:8},{210:250}),1.0,1.0,0)  
  printVarSummary(CESM_SSH_4xCO2)
  do i = 0,dimsizes(CESM_SSH_PD)-1 
    if ((CESM_nino34_4xCO2(i) .le. 0.6) .and. (CESM_nino34_4xCO2(i) .ge. -0.6)) then
            if (abs(CESM_SSH_4xCO2(i)) .ge. 1.2) then
                    CESM_SSH_4xCO2(i) = CESM_SSH_4xCO2@_FillValue
            end if 
    end if 
 end do 
;=======================================================================
; read E3SM nino34 index 
 ncdir01  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/E3SM/monthly/"
 ncdir02  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/E3SM/monthly/"
 
 fil01    = systemfunc("ls " + ncdir01 + "thetao_Omon_E3SM-1-0_piControl_r1i1p1f1_gr_0001-0080_equator_0-200m.nc")
 fil02    = systemfunc("ls " + ncdir02 + "thetao_Omon_E3SM-1-0_abrupt-4xCO2_r1i1p1f1_gr_0001-0080_equator_0-200m.nc")
   f01    = addfile(fil01, "r")
   f02    = addfile(fil02, "r")
   E3SM_TEMP_PD0     = f01->thetao(50*12:79*12+11,0,:,:)
   E3SM_TEMP_4xCO20  = f02->thetao(50*12:79*12+11,0,:,:)
   E3SM_TEMP_PD         = dtrend_msg_n(E3SM_TEMP_PD0&time, E3SM_TEMP_PD0, True, False, 0)
   E3SM_TEMP_4xCO2      = dtrend_msg_n(E3SM_TEMP_4xCO20&time, E3SM_TEMP_4xCO20, True, False, 0)
   copy_VarCoords(E3SM_TEMP_PD0, E3SM_TEMP_PD)
   copy_VarCoords(E3SM_TEMP_4xCO20, E3SM_TEMP_4xCO2)
   E3SM_TEMP_PD_clm    = clmMonTLL(E3SM_TEMP_PD)
   E3SM_TEMP_4xCO2_clm = clmMonTLL(E3SM_TEMP_4xCO2)
   E3SM_TEMP_PD_anom    = calcMonAnomTLL(E3SM_TEMP_PD, E3SM_TEMP_PD_clm)
   E3SM_TEMP_4xCO2_anom = calcMonAnomTLL(E3SM_TEMP_4xCO2,E3SM_TEMP_4xCO2_clm)
   lat01   = f01->lat({-5:5})
   wgtx01  = latRegWgt(lat01, "float", 0)
   E3SM_nino3_PD    = dim_standardize_n_Wrap(wgt_areaave_Wrap(E3SM_TEMP_PD_anom(:,{-5:5},{210:270}),wgtx01,1.0,0),0,0)  
   E3SM_nino3_4xCO2 = dim_standardize_n_Wrap(wgt_areaave_Wrap(E3SM_TEMP_4xCO2_anom(:,{-5:5},{210:270}),wgtx01,1.0,0),0,0)  
  
   E3SM_nino34_PD    = dim_standardize_n_Wrap(wgt_areaave_Wrap(E3SM_TEMP_PD_anom(:,{-5:5},{190:240}),wgtx01,1.0,0),0,0)  
   E3SM_nino34_4xCO2 = dim_standardize_n_Wrap(wgt_areaave_Wrap(E3SM_TEMP_4xCO2_anom(:,{-5:5},{190:240}),wgtx01,1.0,0),0,0)  

   E3SM_TEMP_PD_stat    = dim_stat4_n(E3SM_TEMP_PD_anom,0)
   E3SM_TEMP_4xCO2_stat = dim_stat4_n(E3SM_TEMP_4xCO2_anom,0)
   E3SM_TEMP_PD_stat!0  = "var"
   E3SM_TEMP_PD_stat!1  = "lat"
   E3SM_TEMP_PD_stat!2  = "lon"
   E3SM_TEMP_PD_stat&var  = ispan(0,3,1)
   E3SM_TEMP_PD_stat&lat  = E3SM_TEMP_PD&lat
   E3SM_TEMP_PD_stat&lon  = E3SM_TEMP_PD&lon
   copy_VarMeta(E3SM_TEMP_PD_stat,E3SM_TEMP_4xCO2_stat)
   printVarSummary(E3SM_TEMP_PD_stat)
   printVarSummary(E3SM_TEMP_4xCO2_stat)

 ; E3SM_TEMP_PD_var    = dim_stddev_n_Wrap(E3SM_TEMP_PD_anom,0)
 ; E3SM_TEMP_4xCO2_var = dim_stddev_n_Wrap(E3SM_TEMP_4xCO2_anom,0)

   E3SM_TEMP_PD_var    = E3SM_TEMP_PD_stat(1,:,:)
   E3SM_TEMP_4xCO2_var = E3SM_TEMP_4xCO2_stat(1,:,:)

   printMinMax(E3SM_TEMP_PD_var, 0)
   E3SM_TEMP_PD_skew    = E3SM_TEMP_PD_stat(2,:,:)
   E3SM_TEMP_4xCO2_skew = E3SM_TEMP_4xCO2_stat(2,:,:)

   E3SM_pattern1 = smth9_Wrap(E3SM_TEMP_PD_var, 0.50, 0.25, False)
   E3SM_pattern2 = smth9_Wrap(E3SM_TEMP_4xCO2_var, 0.50, 0.25, False)

   E3SM_pattern01 = smth9_Wrap(E3SM_TEMP_PD_skew, 0.50, 0.25, False)
   E3SM_pattern02 = smth9_Wrap(E3SM_TEMP_4xCO2_skew, 0.50, 0.25, False)
;----------------------------------------------------------------------------------------------------------------------------
;read TIW SST var 
ncdir  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/E3SM/mean/"
fil1    = systemfunc("ls " + ncdir + "tos_day_E3SM-1-0_picontrol_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mon.nc")
fil2    = systemfunc("ls " + ncdir + "tos_day_E3SM-1-0_4xCO2_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mon.nc")

  f1    = addfile(fil1, "r")
  f2    = addfile(fil2, "r")
  SST0_PD     = f1->tos(50*12:79*12+11,:,:)
  SST0_4xCO2  = f2->tos(50*12:79*12+11,:,:)
  SST0_PD!0 = "time"
  SST0_4xCO2!0 = "time"
  SST0_PD&time = E3SM_TEMP_PD0&time
  SST0_4xCO2&time = E3SM_TEMP_4xCO20&time
 ;SST00_PD         = dtrend_msg_n(SST0_PD&time,    SST0_PD, True, False, 0)
 ;SST00_4xCO2      = dtrend_msg_n(SST0_4xCO2&time, SST0_4xCO2, True, False, 0)
  SST_PD      = sqrt(SST0_PD)
  SST_4xCO2   = sqrt(SST0_4xCO2)
  copy_VarCoords(SST0_PD, SST_PD)
  copy_VarCoords(SST0_4xCO2, SST_4xCO2)
  SST_PD_clm    = clmMonTLL(SST_PD)
  SST_4xCO2_clm = clmMonTLL(SST_4xCO2)
  delete([/SST0_PD, SST0_4xCO2/])
  SST_PD_anom     = calcMonAnomTLL(SST_PD, SST_PD_clm)
  SST_4xCO2_anom  = calcMonAnomTLL(SST_4xCO2, SST_4xCO2_clm)

  tiw_SST_PD_var0    = dim_variance_n_Wrap(SST_PD_anom, 0)
  tiw_SST_4xCO2_var0 = dim_variance_n_Wrap(SST_4xCO2_anom, 0)

 ncdir  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/E3SM/mean/"
 fil1    = systemfunc("ls " + ncdir + "zos_day_E3SM-1-0_picontrol_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mon.nc")
 fil2    = systemfunc("ls " + ncdir + "zos_day_E3SM-1-0_4xCO2_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mon.nc")

   f1    = addfile(fil1, "r")
   f2    = addfile(fil2, "r")
   SSH0_PD     = f1->zos(50*12:79*12+11,:,:)
   SSH0_4xCO2  = f2->zos(50*12:79*12+11,:,:)  
   SSH0_PD!0 = "time"
   SSH0_4xCO2!0 = "time"
   SSH0_PD&time    = E3SM_TEMP_PD0&time
   SSH0_4xCO2&time = E3SM_TEMP_4xCO20&time
 ; SSH0_PD         = dtrend_msg_n(SSH0_PD&time,    SSH0_PD, True, False, 0)
 ; SSH0_4xCO2      = dtrend_msg_n(SSH0_4xCO2&time, SSH0_4xCO2, True, False, 0)

   SSH_PD      = sqrt(SSH0_PD)*100
   SSH_4xCO2   = sqrt(SSH0_4xCO2)*100
   copy_VarCoords(SSH0_PD, SSH_PD)
   copy_VarCoords(SSH0_4xCO2, SSH_4xCO2)
   delete([/SSH0_PD, SSH0_4xCO2/])
   SSH_PD_clm    = clmMonTLL(SSH_PD)
   SSH_4xCO2_clm = clmMonTLL(SSH_4xCO2)
 
   SSH_PD_anom     = calcMonAnomTLL(SSH_PD, SSH_PD_clm)
   SSH_4xCO2_anom  = calcMonAnomTLL(SSH_4xCO2, SSH_4xCO2_clm)

   tiw_SSH_PD_var0    = dim_variance_n_Wrap(SSH_PD_anom, 0)
   tiw_SSH_4xCO2_var0 = dim_variance_n_Wrap(SSH_4xCO2_anom, 0)

   ;regression patterns
SST_PD_variance  = new((/dimsizes(SST_PD_anom&lat),dimsizes(SST_PD_anom&lon)/),"double")   
SST_PD_variance!0="lat"
SST_PD_variance!1="lon"
SST_PD_variance&lat = SST_PD_anom&lat
SST_PD_variance&lon = SST_PD_anom&lon
SST_PD_variance  = (/regCoef(E3SM_nino34_PD,SST_PD_anom(lat|:,lon|:,time|:))/)

SST_4xCO2_variance  = new((/dimsizes(SST_PD_anom&lat),dimsizes(SST_PD_anom&lon)/),"double")   
SST_4xCO2_variance!0="lat"
SST_4xCO2_variance!1="lon"
SST_4xCO2_variance&lat = SST_4xCO2_anom&lat
SST_4xCO2_variance&lon = SST_4xCO2_anom&lon
SST_4xCO2_variance  = (/regCoef(E3SM_nino34_4xCO2,SST_4xCO2_anom(lat|:,lon|:,time|:))/)

SSH_PD_variance  = new((/dimsizes(SSH_PD_anom&lat),dimsizes(SSH_PD_anom&lon)/),"double")   
SSH_PD_variance!0="lat"
SSH_PD_variance!1="lon"
SSH_PD_variance&lat = SSH_PD_anom&lat
SSH_PD_variance&lon = SSH_PD_anom&lon
SSH_PD_variance  = (/regCoef(E3SM_nino34_PD,SSH_PD_anom(lat|:,lon|:,time|:))/)

SSH_4xCO2_variance  = new((/dimsizes(SSH_PD_anom&lat),dimsizes(SSH_PD_anom&lon)/),"double")   
SSH_4xCO2_variance!0="lat"
SSH_4xCO2_variance!1="lon"
SSH_4xCO2_variance&lat = SSH_4xCO2_anom&lat
SSH_4xCO2_variance&lon = SSH_4xCO2_anom&lon
SSH_4xCO2_variance  = (/regCoef(E3SM_nino34_4xCO2,SSH_4xCO2_anom(lat|:,lon|:,time|:))/)

;SST_PD_variance     = dim_variance_n_Wrap(SST_PD_anom, 0)
;SST_4xCO2_variance  = dim_variance_n_Wrap(SST_4xCO2_anom, 0)
;printVarSummary(SST_PD_variance)
copy_VarMeta(SST_PD(0,:,:), SST_PD_variance)
copy_VarMeta(SST_4xCO2(0,:,:), SST_4xCO2_variance)

copy_VarMeta(SSH_PD(0,:,:), SSH_PD_variance)
copy_VarMeta(SSH_4xCO2(0,:,:), SSH_4xCO2_variance)

;SSH_PD_variance     = dim_variance_n_Wrap(SSH_PD_anom, 0)
;SSH_4xCO2_variance  = dim_variance_n_Wrap(SSH_4xCO2_anom, 0)

E3SM_pattern001 = smth9_Wrap(SST_PD_variance, 0.50, 0.25, False)
E3SM_pattern002 = smth9_Wrap(SST_4xCO2_variance, 0.50, 0.25, False)

E3SM_pattern0001 = smth9_Wrap(SSH_PD_variance, 0.50, 0.25, False)
E3SM_pattern0002 = smth9_Wrap(SSH_4xCO2_variance, 0.50, 0.25, False)

printMinMax(E3SM_pattern001, 0)
printMinMax(E3SM_pattern002, 0)

printMinMax(E3SM_pattern0001, 0)
printMinMax(E3SM_pattern0002, 0)


E3SM_SST_PD    = wgt_areaave_Wrap(SST_PD_anom(:,{0:5},{200:240}),1.0,1.0,0)  
E3SM_SST_4xCO2 = wgt_areaave_Wrap(SST_4xCO2_anom(:,{0:5},{210:250}),1.0,1.0,0)   

E3SM_SSH_PD    = wgt_areaave_Wrap(SSH_PD_anom(:,{3:8},{210:250}),1.0,1.0,0)  
E3SM_SSH_4xCO2 = wgt_areaave_Wrap(SSH_4xCO2_anom(:,{3:8},{210:240}),1.0,1.0,0)    

;do i = 0,dimsizes(E3SM_SSH_PD)-1 
;if (E3SM_nino34_4xCO2(i) .le. (-0.5)) then
;        if (E3SM_SSH_4xCO2(i) .le. 0.0) then
;                E3SM_SSH_4xCO2(i) = E3SM_SSH_4xCO2@_FillValue
;        end if 
;end if 
;end do 
print(E3SM_SSH_4xCO2)
printMinMax(E3SM_SSH_4xCO2, 0)
E3SM_SSH_PD = dtrend(E3SM_SSH_PD, False)
E3SM_SSH_4xCO2 = dtrend(E3SM_SSH_4xCO2, False)
;================================================================
; CESM scatter
val1 = new ((/3,360/),"float")
val1@_FillValue =-999
val1(0,:) = where(CESM_nino3_PD .lt. -0.50, CESM_nino3_PD, val1@_FillValue)
val1(1,:) = where(CESM_nino3_PD .ge. -0.50 .and. CESM_nino3_PD .lt. 0.50, CESM_nino3_PD, val1@_FillValue)
val1(2,:) = where(CESM_nino3_PD .ge.  0.50, CESM_nino3_PD, val1@_FillValue)


val2 = new ((/3,360/),"float")
val2@_FillValue =-999
val2(0,:) = where(CESM_nino3_4xCO2 .lt. -0.50, CESM_nino3_4xCO2, val2@_FillValue)
val2(1,:) = where(CESM_nino3_4xCO2 .ge. -0.50 .and. CESM_nino3_4xCO2 .lt. 0.50, CESM_nino3_4xCO2, val2@_FillValue)
val2(2,:) = where(CESM_nino3_4xCO2 .ge.  0.50, CESM_nino3_4xCO2, val2@_FillValue)

;trend 
val01 = new ((/2,360/),"float")
val01@_FillValue =-999
val01(0,:) = where(CESM_nino3_PD .lt. -0.5, CESM_nino3_PD, val01@_FillValue)
val01(1,:) = where(CESM_nino3_PD .ge. 0.5,  CESM_nino3_PD, val01@_FillValue)
rc_SST1     =  regline(val01(0,:),CESM_SST_PD)
rc_SST01    =  regline(val01(1,:),CESM_SST_PD)
trend_SST1  =  rc_SST1*(val01(0,:) -rc_SST1@xave)  + rc_SST1@yave
trend_SST01 =  rc_SST01*(val01(1,:)-rc_SST01@xave) + rc_SST01@yave

val02 = new ((/2,360/),"float")
val02@_FillValue =-999
val02(0,:) = where(CESM_nino3_4xCO2 .lt. -0.5, CESM_nino3_4xCO2, val02@_FillValue)
val02(1,:) = where(CESM_nino3_4xCO2 .ge. 0.5,  CESM_nino3_4xCO2, val02@_FillValue)
rc_SST2     =  regline(val02(0,:),CESM_SST_4xCO2)
rc_SST02    =  regline(val02(1,:),CESM_SST_4xCO2)
trend_SST2  =  rc_SST2*(val02(0,:) -rc_SST2@xave)  + rc_SST2@yave
trend_SST02 =  rc_SST02*(val02(1,:)-rc_SST02@xave) + rc_SST02@yave


rc_SSH1     =  regline(val01(0,:),CESM_SSH_PD)
rc_SSH01    =  regline(val01(1,:),CESM_SSH_PD)
trend_SSH1  =  rc_SSH1*(val01(0,:) -rc_SSH1@xave)  + rc_SSH1@yave
trend_SSH01 =  rc_SSH01*(val01(1,:)-rc_SSH01@xave) + rc_SSH01@yave

rc_SSH2     =  regline(val02(0,:),CESM_SSH_4xCO2)
rc_SSH02    =  regline(val02(1,:),CESM_SSH_4xCO2)
trend_SSH2  =  rc_SSH2*(val02(0,:) -rc_SSH2@xave)  + rc_SSH2@yave
trend_SSH02 =  rc_SSH02*(val02(1,:)-rc_SSH02@xave) + rc_SSH02@yave


  CESM_SST_PD_stat    = dim_stat4_n(CESM_SST_PD,0)
  CESM_SST_4xCO2_stat = dim_stat4_n(CESM_SST_4xCO2,0)
  CESM_SST_PD_var     = dim_variance_n_Wrap(CESM_SST_PD, 0)
  CESM_SST_4xCO2_var  = dim_variance_n_Wrap(CESM_SST_4xCO2, 0)
  CESM_SST_PD_skew    = CESM_SST_PD_stat(2)
  CESM_SST_4xCO2_skew = CESM_SST_4xCO2_stat(2)

  CESM_SSH_PD_stat    = dim_stat4_n(CESM_SSH_PD,0)
  CESM_SSH_4xCO2_stat = dim_stat4_n(CESM_SSH_4xCO2,0)
  CESM_SSH_PD_var     = dim_variance_n_Wrap(CESM_SSH_PD, 0)
  CESM_SSH_4xCO2_var  = dim_variance_n_Wrap(CESM_SSH_4xCO2, 0)
  CESM_SSH_PD_skew    = CESM_SSH_PD_stat(2)
  CESM_SSH_4xCO2_skew = CESM_SSH_4xCO2_stat(2)

;------------------------------------------------------------------------------------
; E3SMv1 
val001 = new ((/3,360/),"float")
val001@_FillValue =-999
val001(0,:) = where(E3SM_nino3_PD .lt. -0.50, E3SM_nino3_PD, val001@_FillValue)
val001(1,:) = where(E3SM_nino3_PD .ge. -0.50 .and. E3SM_nino3_PD .lt. 0.50, E3SM_nino3_PD, val001@_FillValue)
val001(2,:) = where(E3SM_nino3_PD .ge.  0.50, E3SM_nino3_PD, val001@_FillValue)

val002 = new ((/3,360/),"float")
val002@_FillValue =-999
val002(0,:) = where(E3SM_nino3_4xCO2 .lt. -0.50, E3SM_nino3_4xCO2, val002@_FillValue)
val002(1,:) = where(E3SM_nino3_4xCO2 .ge. -0.50 .and. E3SM_nino3_4xCO2 .lt. 0.50, E3SM_nino3_4xCO2, val002@_FillValue)
val002(2,:) = where(E3SM_nino3_4xCO2 .ge.  0.50, E3SM_nino3_4xCO2, val002@_FillValue)

;trend 
val0001 = new ((/2,360/),"float")
val0001@_FillValue =-999
val0001(0,:) = where(E3SM_nino3_PD .lt. -0.5, E3SM_nino3_PD, val0001@_FillValue)
val0001(1,:) = where(E3SM_nino3_PD .ge. 0.5, E3SM_nino3_PD, val0001@_FillValue)
rc1_SST1     =  regline(val0001(0,:),E3SM_SST_PD)
rc1_SST01    =  regline(val0001(1,:),E3SM_SST_PD)
trend1_SST1  =  rc1_SST1*(val0001(0,:) -rc1_SST1@xave)  + rc1_SST1@yave
trend1_SST01 =  rc1_SST01*(val0001(1,:)-rc1_SST01@xave) + rc1_SST01@yave

val0002 = new ((/2,360/),"float")
val0002@_FillValue =-999
val0002(0,:) = where(E3SM_nino3_4xCO2 .lt. -0.5, E3SM_nino3_4xCO2, val0002@_FillValue)
val0002(1,:) = where(E3SM_nino3_4xCO2 .ge. 0.5, E3SM_nino3_4xCO2,  val0002@_FillValue)
rc2_SST2     =  regline(val0002(0,:),E3SM_SST_4xCO2)
rc2_SST02    =  regline(val0002(1,:),E3SM_SST_4xCO2)
trend2_SST2  =   rc2_SST2*(val0002(0,:) -rc2_SST2@xave)  + rc2_SST2@yave
trend2_SST02 =   rc2_SST02*(val0002(1,:)-rc2_SST02@xave) + rc2_SST02@yave

rc1_SSH1     =  regline(val0001(0,:),E3SM_SSH_PD)
rc1_SSH01    =  regline(val0001(1,:),E3SM_SSH_PD)
trend1_SSH1  =  rc_SSH1*(val0001(0,:) -rc1_SSH1@xave)  + rc1_SSH1@yave
trend1_SSH01 =  rc_SSH01*(val0001(1,:)-rc1_SSH01@xave) + rc1_SSH01@yave

rc2_SSH2     =  regline(val0002(0,:),E3SM_SSH_4xCO2)
rc2_SSH02    =  regline(val0002(1,:),E3SM_SSH_4xCO2)
trend2_SSH2  =  rc2_SSH2*(val0002(0,:) -rc2_SSH2@xave)  + rc2_SSH2@yave
trend2_SSH02 =  rc2_SSH02*(val0002(1,:)-rc2_SSH02@xave) + rc2_SSH02@yave


E3SM_SST_PD_stat    = dim_stat4_n(E3SM_SST_PD,0)
E3SM_SST_4xCO2_stat = dim_stat4_n(E3SM_SST_4xCO2,0)
E3SM_SST_PD_var     = dim_variance_n_Wrap(E3SM_SST_PD, 0)
E3SM_SST_4xCO2_var  = dim_variance_n_Wrap(E3SM_SST_4xCO2, 0)
E3SM_SST_PD_skew    = E3SM_SST_PD_stat(2)
E3SM_SST_4xCO2_skew = E3SM_SST_4xCO2_stat(2)

E3SM_SSH_PD_stat    = dim_stat4_n(E3SM_SSH_PD,0)
E3SM_SSH_4xCO2_stat = dim_stat4_n(E3SM_SSH_4xCO2,0)
E3SM_SSH_PD_var     = dim_variance_n_Wrap(E3SM_SSH_PD, 0)
E3SM_SSH_4xCO2_var  = dim_variance_n_Wrap(E3SM_SSH_4xCO2, 0)
E3SM_SSH_PD_skew    = E3SM_SSH_PD_stat(2)
E3SM_SSH_4xCO2_skew = E3SM_SSH_4xCO2_stat(2)
print(E3SM_SSH_PD_skew)
print(E3SM_SSH_4xCO2_skew)
;=============================================================
  wks = gsn_open_wks("pdf","/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/revise_figures/plot_PD_4xCO2_CESM_E3SM_TIW_ENSO_future_variability")
  wks_type = "eps"
  wks_type@wkWidth  = 0.8
  wks_type@wkHeight = 0.2    
;=============================================================================================
  gsn_define_colormap(wks,"BlueDarkRed18")
 ;  gsn_define_colormap(wks,"BlueWhiteOrangeRed")
  res                      = True         
  res@gsnDraw              = False             ; don't draw yet
  res@gsnFrame             = False             ; don't advance frame yet
  res@gsnAddCyclic         = False    
  res@gsnMaximize          = False              ; don't maximize

  res@cnFillOn             = True               ; turn on color for contours
  res@cnLinesOn            = False               ; turn off contour lines
  res@cnLineLabelsOn       = False              ; turn off contour line labels
  res@cnMonoLineColor      = True
  res@cnSmoothingOn        = True
  res@cnSmoothingTensionF  = 5
  res@gsnContourNegLineDashPattern = 1
  res@gsnContourZeroLineThicknessF =2
  res@cnInfoLabelOn        = False                ; turn off info 
 ; res@gsnSpreadColorStart   = 130
 ; res@gsnSpreadColorEnd     = 230
  res@cnLevelSelectionMode = "ExplicitLevels"  
 ; res@cnLevels             =fspan(0.1,1.0,10)
 ; res@cnFillColors          = (/130,140,150,160,170,180,190,200,210,220,230/)
  res@cnLevels             =fspan(0.4,1.0,7)
  res@cnFillColors          = (/0,11,12,13,14,15,16,17/)
  res@mpFillOn             = False
  res@mpLandFillColor      = "mediumpurple1"             ; color of land
  res@mpGeophysicalLineColor = "gray30"
  res@mpGeophysicalLineThicknessF = 1.2
  res@mpMaxLatF            = 10
  res@mpMinLatF            = -10
  res@mpMaxLonF            = 280
  res@mpMinLonF            = 120
  res@mpCenterLonF         = 230
  res@mpShapeMode          = "FreeAspect"
  res@lbLabelBarOn         = False
  res@lbOrientation         = "Vertical"
  res@lbBoxEndCapStyle         = "TriangleBothEnds"
  res@tmXBOn   = True
  res@tmYLMode = "Explicit"
 ; res@tmYLValues = (/-10,-8,-6,-4,-2,0,2,4,6,8,10/)
 ; res@tmYLLabels = (/"10~S~o~N~S","8~S~o~N~S","6~S~o~N~S","4~S~o~N~S","2~S~o~N~S","EQ","2~S~o~N~N","4~S~o~N~N","6~S~o~N~N","8~S~o~N~N","10~S~o~N~N"/)
 ; res@tmYLMinorValues         = ispan(-10,10,2) 
  res@tmYLValues = (/-10,-5,0,5,10/)
  res@tmYLLabels = (/"10~S~o~N~S","5~S~o~N~S","EQ","5~S~o~N~N","10~S~o~N~N"/)
  res@tmYLMinorValues         = ispan(-10,10,5) 
  res@tmXBMode = "Explicit"
  res@tmXBValues = (/120,150,180,210,240,270/)
  res@tmXBLabels = (/"","150~S~o~N~E","180","150~S~o~N~W","120~S~o~N~W","90~S~o~N~W"/)  
  res@tmXBLabelFontHeightF  = 0.008
  res@tmYLLabelFontHeightF  = 0.008
  res@tmXBLabelFontThicknessF  = 2.0
  res@tmYLLabelFontThicknessF  = 2.0 
  res@tmBorderThicknessF       = 2.0
  res@tmXBMajorLengthF        = -0.006
  res@tmYLMajorLengthF        = -0.006
  res@tmXBMinorLengthF        = -0.003
  res@tmYLMinorLengthF        = -0.003 
  res@tmXBMajorThicknessF     = 1.5
  res@tmYLMajorThicknessF     = 1.5
  res@tmXBMinorThicknessF     = 1.0
  res@tmYLMinorThicknessF     = 1.0 
  res@tmXBLabelDeltaF         = -0.4
  res@tmYLLabelDeltaF         = -0.4
  res@tmXTOn   = False
  res@tmYROn   = False
  res@tmXMajorGrid                = True
  res@tmXMajorGridLineDashPattern =  2
  res@tmXMajorGridLineColor       = "grey30"
  res@tmXMajorGridThicknessF      = 0.8
  res@tmYMajorGrid                = True
  res@tmYMajorGridLineDashPattern =  2
  res@tmYMajorGridLineColor       = "grey30"
  res@tmYMajorGridThicknessF      = 0.8 
  res@tmYMinorGrid                = False 
  res@tmYMinorGridLineDashPattern =  2
  res@tmYMinorGridLineColor       = "grey30"
  res@tmYMinorGridThicknessF      = 0.8 
  res@tmXMinorGrid                = False 
  res@tmXMinorGridLineDashPattern =  2
  res@tmXMinorGridLineColor       = "grey30"
  res@tmXMinorGridThicknessF      = 0.8 
  res@gsnRightString        = ""
  res@gsnRightStringFontHeightF  = 0.008
  res@vpHeightF            = 0.08    ; Changes the aspect ratio
  res@vpWidthF             = 0.18
  res@vpYF    = 0.95
  res@vpXF    = 0.05+0.03      ; Make plots wider
  res@pmLabelBarWidthF     = 0.03               ; default is shorter
  res@pmLabelBarHeightF    = 0.08
  ;res@pmLabelBarOrthogonalPosF = 0.10
 ; res@pmLabelBarParallelPosF   = -0.00
  res@lbLabelFontHeightF       =0.008
  res@gsnLeftString        = ""
  res@gsnLeftStringFontHeightF  = 0.008
  res@gsnCenterString = "PD  ENSO variance"
  res@gsnCenterStringFontHeightF  = 0.010
  plot1  = gsn_csm_contour_map(wks,CESM_pattern1,res)
  res@tmYLLabels = (/"","","","",""/)
  res@gsnCenterString = ""
  res@gsnLeftString        = ""
  res@vpYF    = 0.95
  res@vpXF    = 0.25+0.03      ; Make plots wider
  res@lbLabelBarOn         = True
  res@gsnCenterString = "4xCO2  ENSO variance"
  plot2  = gsn_csm_contour_map(wks,CESM_pattern3,res)
  res@tmYLLabels = (/"10~S~o~N~S","5~S~o~N~S","EQ","5~S~o~N~N","10~S~o~N~N"/)
  res@mpMaxLatF            = 10
  res@mpMinLatF            = -5
  res@mpMaxLonF            = 280
  res@mpMinLonF            = 180
 ;res@lbTitleString    = ""
 ;res@lbTitlePosition  = "Bottom"
 ;res@lbTitleDirection = "Across" 
 ;res@lbTitleFontHeightF= 0.008
 ;res@lbTitleOffsetF    = 0.05
  delete(res@cnLevels)
  delete(res@cnFillColors) 
 ;res@cnLevels             =fspan(0.005,0.035,7)
 ;res@cnFillColors          = (/0,11,12,13,14,15,16,17/)
  res@cnLevels             =fspan(-0.1,-0.02,9)
  res@cnFillColors          = (/2,3,4,5,6,7,8,9,10,0/)
  res@vpYF    = 0.83
  res@vpXF    = 0.05+0.03      ; Make plots wider
  res@gsnCenterString = "PD  SST IAV"
  res@gsnCenterStringFontHeightF  = 0.010
  res@lbLabelBarOn         = False 
  plot3  = gsn_csm_contour_map(wks,CESM_pattern001,res)
  res@tmYLLabels = (/"","","","",""/)
  res@vpYF    = 0.83
  res@vpXF    = 0.25+0.03      ; Make plots wider
  res@lbLabelBarOn         = True
  res@gsnCenterString = "4xCO2  SST IAV"
  plot4  = gsn_csm_contour_map(wks,CESM_pattern003,res)
  delete(res@cnLevels)
  delete(res@cnFillColors)  
  res@tmYLLabels = (/"10~S~o~N~S","5~S~o~N~S","EQ","5~S~o~N~N","10~S~o~N~N"/)
  res@cnLevels             =fspan(-0.7,-0.1,7)
  res@cnFillColors          = (/3,4,5,6,7,8,9,0/)
  ;res@cnLevels             =fspan(0.5,3.5,7)
  ;res@cnFillColors          = (/0,11,12,13,14,15,16,17/)
  res@vpYF    = 0.71
  res@vpXF    = 0.05+0.03      ; Make plots wider
  res@gsnCenterString = "PD  SSH IAV"
  res@gsnCenterStringFontHeightF  = 0.010
  res@lbLabelBarOn         = False 
  plot5  = gsn_csm_contour_map(wks,CESM_pattern0001,res)
  res@tmYLLabels = (/"","","","",""/)
  res@vpYF    = 0.71
  res@vpXF    = 0.25+0.03      ; Make plots wider
  res@lbLabelBarOn         = True
  res@gsnCenterString = "4xCO2  SSH IAV"
  plot6  = gsn_csm_contour_map(wks,CESM_pattern0003,res)
;================================================================
  res@cnLevels              = fspan(0.4,2.8,7)
  res@cnFillColors          = (/0,11,12,13,14,15,16,17/)
  res@tmYLLabels = (/"10~S~o~N~S","5~S~o~N~S","EQ","5~S~o~N~N","10~S~o~N~N"/)
  res@vpYF    = 0.95
  res@vpXF    = 0.52+0.03      ; Make plots wider
  res@mpMaxLatF            = 10
  res@mpMinLatF            = -10
  res@mpMaxLonF            = 280
  res@mpMinLonF            = 120
  res@gsnCenterString = "PD  ENSO variance"
  res@gsnCenterStringFontHeightF  = 0.010
  res@lbLabelBarOn         = False 
  plot01  = gsn_csm_contour_map(wks,E3SM_pattern1,res)
  res@tmYLLabels = (/"","","","",""/)
  res@vpYF    = 0.95
  res@vpXF    = 0.72+0.03     ; Make plots wider
  res@lbLabelBarOn         = True
  res@gsnCenterString = "4xCO2  ENSO variance"
  plot02  = gsn_csm_contour_map(wks,E3SM_pattern2,res)
  delete(res@cnLevels)
  delete(res@cnFillColors) 
  res@tmYLLabels = (/"10~S~o~N~S","5~S~o~N~S","EQ","5~S~o~N~N","10~S~o~N~N"/)
  res@cnLevels              = fspan(-0.20,-0.04,9)
  res@cnFillColors          = (/2,3,4,5,6,7,8,9,10,0/)
  ;res@cnLevels             =fspan(0.02,0.14,7)
  ;res@cnFillColors          = (/0,11,12,13,14,15,16,17/)
  res@mpMaxLatF            = 10
  res@mpMinLatF            = -5
  res@mpMaxLonF            = 280
  res@mpMinLonF            = 180
 ; res@lbTitleString    = ""
 ; res@lbTitlePosition  = "Bottom"
 ; res@lbTitleDirection = "Across" 
 ; res@lbTitleFontHeightF= 0.008
 ; res@lbTitleOffsetF    = 0.05

  res@vpYF    = 0.83
  res@vpXF    = 0.52+0.03      ; Make plots wider
  res@gsnCenterString = "PD  SST IAV"
  res@gsnCenterStringFontHeightF  = 0.010
  res@lbLabelBarOn         = False 
  plot03  = gsn_csm_contour_map(wks,E3SM_pattern001,res)
  res@tmYLLabels = (/"","","","",""/)
  res@vpYF    = 0.83
  res@vpXF    = 0.72+0.03      ; Make plots wider
  res@lbLabelBarOn         = True
  res@gsnCenterString = "4xCO2  SST IAV"
  plot04  = gsn_csm_contour_map(wks,E3SM_pattern002,res)
  delete(res@cnLevels)
  delete(res@cnFillColors) 
  res@tmYLLabels = (/"10~S~o~N~S","5~S~o~N~S","EQ","5~S~o~N~N","10~S~o~N~N"/)
  res@cnLevels             =fspan(-0.9,-0.1,9)
  res@cnFillColors          = (/2,3,4,5,6,7,8,9,10,0/)
 ;res@cnLevels             =fspan(0.5,3.5,7)
 ;res@cnFillColors          = (/0,11,12,13,14,15,16,17/)
  res@vpYF    = 0.71
  res@vpXF    = 0.52+0.03      ; Make plots wider
  res@gsnCenterString = "PD  SSH IAV"
  res@gsnCenterStringFontHeightF  = 0.010
  res@lbLabelBarOn         = False 
  plot05  = gsn_csm_contour_map(wks,E3SM_pattern0001,res)
  res@tmYLLabels = (/"","","","",""/)
  res@vpYF    = 0.71
  res@vpXF    = 0.72+0.03      ; Make plots wider
  res@lbLabelBarOn         = True
  res@gsnCenterString = "4xCO2  SSH IAV"
  plot06  = gsn_csm_contour_map(wks,E3SM_pattern0002,res)


plres=True
plres@gsLineThicknessF  = 2.5
; plres@gsFillColor     = -1
plres@gsLineDashPattern= 2
plres@gsLineColor="black"
plres@gsLineOpacityF=100
w1=360-160
e1=360-100

w2=360-160
e2=360-110

s1=0
n1=5

s2=2
n2=6
dum = new(4, graphic)

dum(0)=gsn_add_polyline(wks, plot3, (/w2,e2,e2,w2,w1/), (/s1,s1,n1,n1,s1/), plres)
dum(1)=gsn_add_polyline(wks, plot4, (/w2,e2,e2,w2,w2/), (/s1,s1,n1,n1,s1/), plres)
dum(2)=gsn_add_polyline(wks, plot03, (/w2,e2,e2,w2,w2/), (/s1,s1,n1,n1,s1/), plres)
dum(3)=gsn_add_polyline(wks, plot04, (/w2,e2,e2,w2,w2/), (/s1,s1,n1,n1,s1/), plres)

plres=True
plres@gsLineThicknessF=2.5
; plres@gsFillColor     = -1
plres@gsLineDashPattern= 2
plres@gsLineColor="black"
plres@gsLineOpacityF=100
w1=360-160
e1=360-110

s1=3
n1=8

dum1=new(4, graphic)

dum1(0)=gsn_add_polyline(wks, plot5, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum1(1)=gsn_add_polyline(wks, plot6, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum1(2) =gsn_add_polyline(wks,plot05, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum1(3) =gsn_add_polyline(wks,plot06, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)

;=============================plot skew=======================================
  draw(plot1)
  draw(plot2)
  draw(plot3)
  draw(plot4)
  draw(plot5)
  draw(plot6)

  draw(plot01)
  draw(plot02)
  draw(plot03)
  draw(plot04)
  draw(plot05)
  draw(plot06)
;=========================================================================
 ;pres                 = True
 ;pres@gsnPanelRowSpec = True                   ; tell panel what order to plot
 ;pres@gsnFrame         = False               ; save panel until both ready
 ;pres@gsnPanelLabelBar = True		      ; commom label bar
 ;pres@lbLabelStride    = 1	      ; every 4th label
 ;gsn_panel(wks,plot,(/4,4/),pres)

    pos1 = get_title_position(plot1)         ; 1st row leftmost plot
    pos2 = get_title_position(plot2)         ; 2nd row leftmost plot
    pos3 = get_title_position(plot3)         ; 1st row leftmost plot
    pos4 = get_title_position(plot4)         ; 2nd row leftmost plot
    pos5 = get_title_position(plot5)         ; 1st row leftmost plot
    pos6 = get_title_position(plot6)         ; 2nd row leftmost plot

    pos01 = get_title_position(plot01)         ; 1st row leftmost plot
    pos02 = get_title_position(plot02)         ; 2nd row leftmost plot
    pos03 = get_title_position(plot03)         ; 1st row leftmost plot
    pos04 = get_title_position(plot04)         ; 2nd row leftmost plot
    pos05 = get_title_position(plot05)         ; 1st row leftmost plot
    pos06 = get_title_position(plot06)         ; 2nd row leftmost plot

    pos001 = get_title_position1(plot1)         ; 1st row leftmost plot
    pos002 = get_title_position1(plot01)         ; 2nd row leftmost plot
    txres               = True
    txres@txFontHeightF = 0.015
    txres@txAngleF      = 0.
    txres@txFontThicknessF = 6
    gsn_text_ndc(wks, "~F22~a",  pos1(0), pos1(1), txres)
   ; gsn_text_ndc(wks, "b)",  pos2(0), pos2(1), txres)
    gsn_text_ndc(wks, "~F22~b",  pos3(0), pos3(1), txres)
   ; gsn_text_ndc(wks, "d)",  pos4(0), pos4(1), txres)
    gsn_text_ndc(wks, "~F22~c",  pos5(0), pos5(1), txres)
   ; gsn_text_ndc(wks, "f)",  pos6(0), pos6(1), txres)

    gsn_text_ndc(wks, "~F22~f",  pos01(0), pos01(1), txres)
   ; gsn_text_ndc(wks, "h)",  pos02(0), pos02(1), txres)
    gsn_text_ndc(wks, "~F22~g",  pos03(0), pos03(1), txres)
   ; gsn_text_ndc(wks, "j)",  pos04(0), pos04(1), txres)
    gsn_text_ndc(wks, "~F22~h",  pos05(0), pos05(1), txres)
  ;  gsn_text_ndc(wks, "l)",  pos06(0), pos06(1), txres)
   ;gsn_text_ndc(wks, "g)",  pos03(0),pos03(1), txres)
   ;gsn_text_ndc(wks, "h)",  pos04(0),pos04(1), txres)

    gsn_text_ndc(wks, "CESM-UHR", pos001(0),pos001(1), txres)
    gsn_text_ndc(wks, "E3SMv1",  pos002(0),pos002(1), txres)
;-----------------------------------------------------------------------------
sres             = True
sres@gsnDraw     = False      ; Will draw later, after overlaying
sres@gsnFrame    = False      ; all plots

sres@vpWidthF    = 0.18        ; Make plots wider
sres@vpHeightF   = 0.18


sres@trXMinF     = -4.0
sres@trXMaxF     = 4.0
sres@trYMinF     =-0.4
sres@trYMaxF     = 0.4
sres@gsnYRefLine           = 0.0           ; reference line
sres@gsnYRefLineColor      ="black"
sres@gsnYRefLineDashPattern =1
sres@gsnYRefLineThicknessF = 2.0
sres@gsnXRefLine           = 0.0              ; reference line
sres@gsnXRefLineColor      ="black"
sres@gsnXRefLineDashPattern =1
sres@gsnXRefLineThicknessF=2.0
; sres@tmYLMinorPerMajor = 0
; sres@tmXBMinorPerMajor = 18
;sres@tmBorderThicknessF = 6

sres@tmYLMajorThicknessF =3
sres@tmYLMinorLengthF = 0.003
sres@tmYLMinorLengthF = 0.003
;sres@tmXBMajorThicknessF = 6
sres@tmYLMajorLengthF = 0.01

; res@tmXBLabelAngleF = 45
sres@tmXBLabelFontHeightF = 0.008
sres@tmYLLabelFontHeightF = 0.008
    
;sres@tiYAxisOffsetXF = 0.01
;sres@tiYAxisFontHeightF =0.01

sres@tmYLMode   = "Explicit"          
sres@tmYLValues      = fspan(-1.0,1.0,11)
sres@tmYLLabels      = (/"-1.0","-0.8","-0.6","-0.4","-0.2","0","0.2","0.4","0.6","0.8","1.0"/)
;---Turn off some tickmarks
sres@tmXTOn                 = False    ; bottom off
sres@tmYROn                 = False    ; right off

sres@tiXAxisString = "  Ni"+"n~H-13V2F35~D~FV-2H3~"+"o 34(~S~o~N~C)"
sres@tiXAxisFontHeightF =0.010
sres@tiYAxisString = "SST(~S~o~N~C)"
sres@tiYAxisFontHeightF =0.010

sres@xyMarkLineModes     = (/"Markers", "Markers","Markers"/) ; choose which have markers
sres@xyMarkers           = (/16,16,16/)                    ; choose type of marker 
sres@xyMarkerColors       = (/"blue","gray40","red"/)                  ; Marker color
sres@xyMarkerSizeF       = 0.003                 ; Marker size (default 0.01)
; res@xyDashPatterns      = 1                      ; solid line 
;  res@tmYLFormat          = "f"                    ; not necessary but nicer labels 

; res@tiMainString        = "Output from regline_stats"  ; title
sres@tfPolyDrawOrder   ="PostDraw"
sres@gsnLeftString     =""
sres@gsnLeftStringFontHeightF =0.01
sres@gsnCenterString = "PD Ni"+"n~H-13V2F35~D~FV-2H3~"+"o34 & SST"
sres@gsnCenterStringFontHeightF  = 0.010
sres@vpYF    = 0.57
sres@vpXF    = 0.05+0.03            ; Make plots wider
plot7 = gsn_csm_xy(wks,val1,CESM_SST_PD,sres) ; create plot

  xres = True
  xres@gsnDraw=False
  xres@gsnFrame=False
  xres@xyDashPattern=2
  xres@xyLineThicknessF=2
  xres@xyLineColor      ="green"   
    plot07 = gsn_csm_xy(wks,val01(0,:),trend_SST1,xres)
    plot007 = gsn_csm_xy(wks,val01(1,:),trend_SST01,xres)
    overlay(plot7,plot07)
    overlay(plot7,plot007)
sres@tmYLOn        = False
sres@vpYF    = 0.57
sres@vpXF    = 0.25+0.03            ; Make plots wider
sres@tiYAxisString = ""
sres@gsnCenterString = "4xCO2 Ni"+"n~H-13V2F35~D~FV-2H3~"+"o34 & SST"
sres@gsnCenterStringOrthogonalPosF = 0.020
plot8 = gsn_csm_xy (wks,val2,CESM_SST_4xCO2,sres) ; create plot
plot08 = gsn_csm_xy(wks,val02(0,:),trend_SST2,xres)
plot008 = gsn_csm_xy(wks,val02(1,:),trend_SST02,xres)
overlay(plot8,plot08)
overlay(plot8,plot008)
;===================================================
sres@tmYLOn       = True
sres@trYMinF     =-0.8
sres@trYMaxF     = 0.8
sres@vpYF    = 0.57
sres@vpXF    = 0.52+0.03            ; Make plots wider
sres@tiYAxisString = "SST(~S~o~N~C)"
sres@gsnCenterString = "PD Ni"+"n~H-13V2F35~D~FV-2H3~"+"o34 & SST"
sres@gsnCenterStringOrthogonalPosF = 0.0
plot9 = gsn_csm_xy (wks,val001,E3SM_SST_PD,sres) ; create plot
plot09 = gsn_csm_xy(wks,val0001(0,:),trend1_SST1,xres)
plot009 = gsn_csm_xy(wks,val0001(1,:),trend1_SST01,xres)
overlay(plot9,plot09)
overlay(plot9,plot009)
sres@tmYLOn       = False
sres@vpYF    = 0.57
sres@vpXF    = 0.72+0.03            ; Make plots wider
sres@tiYAxisString = ""
sres@gsnCenterString = "4xCO2 Ni"+"n~H-13V2F35~D~FV-2H3~"+"o34 & SST"
sres@gsnCenterStringOrthogonalPosF = 0.020
plot10 = gsn_csm_xy (wks,val002,E3SM_SST_4xCO2,sres) ; create plot
plot010 = gsn_csm_xy(wks,val0002(0,:),trend2_SST2,xres)
plot0010 = gsn_csm_xy(wks,val0002(1,:),trend2_SST02,xres)
overlay(plot10,plot010)
overlay(plot10,plot0010)
;----------------------------------------------------------------------
sres@tmYLOn       = True
sres@trYMinF     =-3.0
sres@trYMaxF     = 3.0
sres@vpYF    = 0.32
sres@vpXF    = 0.05+0.03            ; Make plots wider
sres@tiYAxisString = "SSH(~S~o~N~C)"
sres@gsnCenterString = "PD Ni"+"n~H-13V2F35~D~FV-2H3~"+"o34 & SSH"
sres@gsnCenterStringOrthogonalPosF = 0.0
sres@tmYLValues      = fspan(-3.0,3.0,11)
sres@tmYLLabels      = (/"-3.0","-2.4","-1.8","-1.2","-0.6","0","0.6","1.2","1.8","2.4","3.0"/)
plot11 = gsn_csm_xy (wks,val1,CESM_SSH_PD,sres) ; create plot
plot011 = gsn_csm_xy(wks,val01(0,:),trend_SSH1,xres)
plot0011 = gsn_csm_xy(wks,val01(1,:),trend_SSH01,xres)
overlay(plot11,plot011)
overlay(plot11,plot0011)

sres@tmYLOn       = False
sres@vpYF    = 0.32
sres@vpXF    = 0.25+0.03            ; Make plots wider
sres@tiYAxisString = ""
sres@gsnCenterString = "4xCO2 Ni"+"n~H-13V2F35~D~FV-2H3~"+"o34 & SSH"
sres@gsnCenterStringOrthogonalPosF = 0.020
plot12 = gsn_csm_xy (wks,val2,CESM_SSH_4xCO2,sres) ; create plot
plot012 = gsn_csm_xy(wks,val02(0,:),trend_SSH2,xres)
plot0012 = gsn_csm_xy(wks,val02(1,:),trend_SSH02,xres)
overlay(plot12,plot012)
overlay(plot12,plot0012)
;----------------------------------------------------------
sres@tmYLValues      = fspan(-4.0,4.0,11)
sres@tmYLLabels      = (/"-4.0","-3.2","-2.4","-1.6","-0.8","0","0.8","1.6","2.4","3.2","4.0"/)
sres@tmYLOn       = True
sres@trYMinF     =-4.0
sres@trYMaxF     = 4.0
sres@vpYF    = 0.32
sres@vpXF    = 0.52+0.03            ; Make plots wider
sres@tiYAxisString = "SSH(~S~o~N~C)"
sres@gsnCenterString = "PD Ni"+"n~H-13V2F35~D~FV-2H3~"+"o34 & SSH"
sres@gsnCenterStringOrthogonalPosF = 0.0
plot13 = gsn_csm_xy (wks,val001,E3SM_SSH_PD,sres) ; create plot
plot013 = gsn_csm_xy(wks,val0001(0,:),trend1_SSH1,xres)
plot0013 = gsn_csm_xy(wks,val0001(1,:),trend1_SSH01,xres)
overlay(plot13,plot013)
overlay(plot13,plot0013)

sres@tmYLOn       = False
sres@vpYF    = 0.32
sres@vpXF    = 0.72+0.03            ; Make plots wider
sres@tiYAxisString = ""
sres@gsnCenterString = "4xCO2 Ni"+"n~H-13V2F35~D~FV-2H3~"+"o34 & SSH"
sres@gsnCenterStringOrthogonalPosF = 0.020
plot14   = gsn_csm_xy (wks,val002,E3SM_SSH_4xCO2,sres) ; create plot
plot014  = gsn_csm_xy(wks,val0002(0,:),trend2_SSH2,xres)
plot0014 = gsn_csm_xy(wks,val0002(1,:),trend2_SSH02,xres)
overlay(plot14,plot014)
overlay(plot14,plot0014)
;---------------------------------------------------------------------
txres                       = True
txres@txPerimOn             = False  
txres@txBackgroundFillColor = "transparent"
txres@txFontHeightF         = 0.010
txres@txFontThicknessF      = 6
txres@txFont              = 21
txres@txFontColor         = "blue"

txres@txFontColor         = "blue"
txid1     = gsn_create_text(wks, "Var =" + " ("+sprintf("%5.2f", CESM_SST_PD_var) +")", txres)
txres@txFontColor         = "red"
txid2     = gsn_create_text(wks, "Skew =" + " ("+sprintf("%5.2f",CESM_SST_PD_skew) +")", txres)

txres@txFontColor         = "blue"
txid3     = gsn_create_text(wks, "Var =" + " ("+sprintf("%5.2f", CESM_SST_4xCO2_var) +")", txres)
txres@txFontColor         = "red"
txid4     = gsn_create_text(wks, "Skew =" + " ("+sprintf("%5.2f",CESM_SST_4xCO2_skew) +")", txres)

txres@txFontColor         = "blue"
txid5     = gsn_create_text(wks, "Var =" + " ("+sprintf("%5.2f", E3SM_SST_PD_var) +")", txres)
txres@txFontColor         = "red"
txid6     = gsn_create_text(wks, "Skew =" + " ("+sprintf("%5.2f",E3SM_SST_PD_skew) +")", txres)

txres@txFontColor         = "blue"
txid7     = gsn_create_text(wks, "Var =" + " ("+sprintf("%5.2f", E3SM_SST_4xCO2_var) +")", txres)
txres@txFontColor         = "red"
txid8     = gsn_create_text(wks, "Skew =" + " ("+sprintf("%5.2f",E3SM_SST_4xCO2_skew) +")", txres)


amres                  = True
amres@amJust           = "TopCenter"
amres@amParallelPosF   =  0.24     ; 0.5 This is the right edge of the plot.
amres@amOrthogonalPosF = -0.49     ; 0.5 This is the bottom edge of the plot.
annoid1 = gsn_add_annotation(plot7, txid1, amres)
amres@amOrthogonalPosF = -0.40     ; 0.5 This is the bottom edge of the plot.
annoid2 = gsn_add_annotation(plot7, txid2, amres)

amres@amParallelPosF   =  0.24     ; 0.5 This is the right edge of the plot.
amres@amOrthogonalPosF = -0.49     ; 0.5 This is the bottom edge of the plot.
annoid3 = gsn_add_annotation(plot8, txid3, amres)
amres@amOrthogonalPosF = -0.40     ; 0.5 This is the bottom edge of the plot.
annoid4 = gsn_add_annotation(plot8, txid4, amres)


amres@amParallelPosF   =  0.24     ; 0.5 This is the right edge of the plot.
amres@amOrthogonalPosF = -0.49     ; 0.5 This is the bottom edge of the plot.
annoid5 = gsn_add_annotation(plot9, txid5, amres)
amres@amOrthogonalPosF = -0.40     ; 0.5 This is the bottom edge of the plot.
annoid6 = gsn_add_annotation(plot9, txid6, amres)

amres@amParallelPosF   =  0.24      ; 0.5 This is the right edge of the plot.
amres@amOrthogonalPosF = -0.49     ; 0.5 This is the bottom edge of the plot.
annoid7 = gsn_add_annotation(plot10, txid7, amres)
amres@amOrthogonalPosF = -0.40     ; 0.5 This is the bottom edge of the plot.
annoid8 = gsn_add_annotation(plot10, txid8, amres)
;============================================================================
txres@txFontColor         = "blue"
txid01     = gsn_create_text(wks, "Var =" + " ("+sprintf("%5.2f", CESM_SSH_PD_var) +")", txres)
txres@txFontColor         = "red"
txid02     = gsn_create_text(wks, "Skew =" + " ("+sprintf("%5.2f",CESM_SSH_PD_skew) +")", txres)

txres@txFontColor         = "blue"
txid03     = gsn_create_text(wks, "Var =" + " ("+sprintf("%5.2f", CESM_SSH_4xCO2_var) +")", txres)
txres@txFontColor         = "red"
txid04     = gsn_create_text(wks, "Skew =" + " ("+sprintf("%5.2f",CESM_SSH_4xCO2_skew) +")", txres)

txres@txFontColor         = "blue"
txid05     = gsn_create_text(wks, "Var =" + " ("+sprintf("%5.2f", E3SM_SSH_PD_var) +")", txres)
txres@txFontColor         = "red"
txid06     = gsn_create_text(wks, "Skew =" + " ("+sprintf("%5.2f",E3SM_SSH_PD_skew) +")", txres)

txres@txFontColor         = "blue"
txid07     = gsn_create_text(wks, "Var =" + " ("+sprintf("%5.2f", E3SM_SSH_4xCO2_var) +")", txres)
txres@txFontColor         = "red"
txid08     = gsn_create_text(wks, "Skew =" + " ("+sprintf("%5.2f",E3SM_SSH_4xCO2_skew) +")", txres)


amres                  = True
amres@amJust           = "TopCenter"
amres@amParallelPosF   =  0.24     ; 0.5 This is the right edge of the plot.
amres@amOrthogonalPosF = -0.49     ; 0.5 This is the bottom edge of the plot.
annoid1 = gsn_add_annotation(plot11, txid01, amres)
amres@amOrthogonalPosF = -0.40     ; 0.5 This is the bottom edge of the plot.
annoid2 = gsn_add_annotation(plot11, txid02, amres)

amres@amParallelPosF   =  0.24     ; 0.5 This is the right edge of the plot.
amres@amOrthogonalPosF = -0.49     ; 0.5 This is the bottom edge of the plot.
annoid3 = gsn_add_annotation(plot12, txid03, amres)
amres@amOrthogonalPosF = -0.40     ; 0.5 This is the bottom edge of the plot.
annoid4 = gsn_add_annotation(plot12, txid04, amres)


amres@amParallelPosF   =  0.24     ; 0.5 This is the right edge of the plot.
amres@amOrthogonalPosF = -0.49     ; 0.5 This is the bottom edge of the plot.
annoid5 = gsn_add_annotation(plot13, txid05, amres)
amres@amOrthogonalPosF = -0.40     ; 0.5 This is the bottom edge of the plot.
annoid6 = gsn_add_annotation(plot13, txid06, amres)

amres@amParallelPosF   =  0.24      ; 0.5 This is the right edge of the plot.
amres@amOrthogonalPosF = -0.49     ; 0.5 This is the bottom edge of the plot.
annoid7 = gsn_add_annotation(plot14, txid07, amres)
amres@amOrthogonalPosF = -0.40     ; 0.5 This is the bottom edge of the plot.
annoid8 = gsn_add_annotation(plot14, txid08, amres)

  getvalues plot7
    "trYMinF": ymin1
    "trYMaxF": ymax1
  end getvalues

  getvalues plot8
    "trYMinF": ymin2
    "trYMaxF": ymax2
  end getvalues

  getvalues plot9
  "trYMinF": ymin3
  "trYMaxF": ymax3
  end getvalues

 getvalues plot10
    "trYMinF": ymin4
    "trYMaxF": ymax4
  end getvalues

  getvalues plot11
    "trYMinF": ymin5
    "trYMaxF": ymax5
  end getvalues

  getvalues plot12
    "trYMinF": ymin6
    "trYMaxF": ymax6
  end getvalues

  getvalues plot13
  "trYMinF": ymin7
  "trYMaxF": ymax7
  end getvalues

 getvalues plot14
    "trYMinF": ymin8
    "trYMaxF": ymax8
  end getvalues
xmin = (/-0.5/)
xmax = (/0.5/)

xbars = (/(/xmin(0), xmax(0), xmax(0), xmin(0), xmin(0)/)/)

ybars1 = (/ (/ymin1,ymin1,ymax1,ymax1,ymin1/)/)
ybars2 = (/ (/ymin2,ymin2,ymax2,ymax2,ymin2/)/)
ybars3 = (/ (/ymin3,ymin3,ymax3,ymax3,ymin3/)/)
ybars4 = (/ (/ymin4,ymin4,ymax4,ymax4,ymin4/)/)
ybars5 = (/ (/ymin5,ymin5,ymax5,ymax5,ymin5/)/)
ybars6 = (/ (/ymin6,ymin6,ymax6,ymax6,ymin6/)/)
ybars7 = (/ (/ymin7,ymin7,ymax7,ymax7,ymin7/)/)
ybars8 = (/ (/ymin7,ymin8,ymax8,ymax8,ymin8/)/)

;;;;;; Create the bars and attach to appropriate plot
gres             = True
gres@gsFillColor = "pink"
gres@tfPolyDrawOrder = "Draw" 
gres@gsFillOpacityF  = 0.60
; gres@gsEdgesOn   = True       ; draw border around polygons
; gres@gsEdgeColor = "lightgray"  
; gres@gsEdgeDashPattern =2
poly1 = gsn_add_polygon(wks,plot7,xbars,ybars1,gres)
poly2 = gsn_add_polygon(wks,plot8,xbars,ybars2,gres)
poly3 = gsn_add_polygon(wks,plot9,xbars,ybars3,gres)
poly4 = gsn_add_polygon(wks,plot10,xbars,ybars4,gres)
poly5 = gsn_add_polygon(wks,plot11,xbars,ybars5,gres)
poly6 = gsn_add_polygon(wks,plot12,xbars,ybars6,gres)
poly7 = gsn_add_polygon(wks,plot13,xbars,ybars7,gres)
poly8 = gsn_add_polygon(wks,plot14,xbars,ybars8,gres)

;txres1 = True ; text mods desired
;txres1@txJust = "CenterLeft" ; Default is "CenterCenter" ; note that a colon is the default function code, but since I have personally 
;; set the default code to be a tilde (~) in my .hluresfile, I manually 
;; reset it here, since the colon makes for a cleaner example. 
;txres1@txFuncCode = "~"
;;**************************** 
;; method 1: place each component manually on the page. I found this 
;; method to be a bit faster to create and easier to read 
;;**************************** 
;txres1@txFontHeightF = 0.011
;txres1@txFontThicknessF =1.5
;
;gsn_text_ndc(wks,"PD SST",0.10,0.91,txres1) 
;gsn_text_ndc(wks,"4xCO2 SST",0.30,0.91,txres1) 
;gsn_text_ndc(wks,"PD SSH",0.55,0.91,txres1) 
;gsn_text_ndc(wks,"4xCO2 SSH",0.75,0.91,txres1) 

draw(plot7)
draw(plot8)
draw(plot9)
draw(plot10)
draw(plot11)
draw(plot12)
draw(plot13)
draw(plot14)


pos7  = get_title_position(plot7)         ; 1st row leftmost plot
pos8  = get_title_position(plot8)         ; 2nd row leftmost plot
pos9  = get_title_position(plot9)         ; 1st row leftmost plot
pos10 = get_title_position(plot10)         ; 2nd row leftmost plot
pos11 = get_title_position(plot11)         ; 2nd row leftmost plot
pos12 = get_title_position(plot12)         ; 2nd row leftmost plot
pos13 = get_title_position(plot13)         ; 2nd row leftmost plot
pos14 = get_title_position(plot14)         ; 2nd row leftmost plot

txres               = True
txres@txFontHeightF = 0.015
txres@txAngleF      = 0.
txres@txFontColor         = "black"
gsn_text_ndc(wks, "~F22~d",  pos7(0), pos7(1), txres)
gsn_text_ndc(wks, "~F22~i",  pos9(0), pos9(1), txres)
gsn_text_ndc(wks, "~F22~e",  pos11(0), pos11(1), txres)
gsn_text_ndc(wks, "~F22~j",  pos13(0), pos13(1), txres)

end 