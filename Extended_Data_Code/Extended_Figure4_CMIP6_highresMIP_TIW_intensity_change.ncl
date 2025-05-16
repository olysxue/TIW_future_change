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
  pos = (/vpx*0.65,vpy-(vph/2.)/)
  return(pos)
end

; plot mean TIW EKE and RMS
begin 
 ncdir1  = "/Volumes/oly2/highresMIP/CO2/std/"

 ;=======================================================================================
 ; read CMCC HR4
 fil1    = systemfunc("ls " + ncdir1 + "tos_Oday_AWI-CM-1-1-MR_piControl_r1i1p1f1_gn_2851-2900_equator_10-50fft_std.nc")
 fil2    = systemfunc("ls " + ncdir1 + "tos_Oday_AWI-CM-1-1-MR_abrupt-4xCO2_r1i1p1f1_gn_1951-2000_equator_10-50fft_std.nc")

   f1    = addfile(fil1, "r")
   f2    = addfile(fil2, "r")
   tiw_AWI_hist0    = f1->tos(:,:,:)
   tiw_AWI_future0  = f2->tos(:,:,:)

 fil1    = systemfunc("ls " + ncdir1 + "tos_Oday_CNRM-CM6-1-HR_piControl_r1i1p1f2_gn_2056-2105_equator_10-50fft_std.nc")
 fil2    = systemfunc("ls " + ncdir1 + "tos_Oday_CNRM-CM6-1-HR_abrupt-4xCO2_r1i1p1f2_gn_1950-1999_equator_10-50fft_std.nc")
  
   f1    = addfile(fil1, "r")
   f2    = addfile(fil2, "r")
   tiw_CNRM_hist0    = f1->tos(:,:,:)
   tiw_CNRM_future0  = f2->tos(:,:,:)
;==========================================================================================
; read GFDL
 fil01    = systemfunc("ls " + ncdir1 + "tos_Oday_GFDL-CM4_piControl_r1i1p1f1_gn_0601-0650_equator_10-50fft_std.nc")
 fil02    = systemfunc("ls " + ncdir1 + "tos_Oday_GFDL-CM4_abrupt-4xCO2_r1i1p1f1_gn_0101-0150_equator_10-50fft_std.nc")

   f01    = addfile(fil01, "r")
   f02    = addfile(fil02, "r")

   tiw_GFDL_hist0   = f01->tos(:,:,:)
   tiw_GFDL_future0   = f02->tos(:,:,:)


   fil01    = systemfunc("ls " + ncdir1 + "tos_Oday_HadGEM3-GC31-MM_piControl_r1i1p1f1_gn_2300-2349_equator_10-50fft_std.nc")
   fil02    = systemfunc("ls " + ncdir1 + "tos_Oday_HadGEM3-GC31-MM_abrupt-4xCO2_r1i1p1f3_gn_1950-1999_equator_10-50fft_std.nc")
     f01    = addfile(fil01, "r")
     f02    = addfile(fil02, "r")
     tiw_Had_hist0   = f01->tos(:,:,:)
     tiw_Had_future0   = f02->tos(:,:,:)

     tiw_AWI_hist      = dim_avg_n_Wrap(tiw_AWI_hist0, 0)
     tiw_AWI_future    = dim_avg_n_Wrap(tiw_AWI_future0, 0)
     
     tiw_CNRM_hist      = dim_avg_n_Wrap(tiw_CNRM_hist0, 0)
     tiw_CNRM_future    = dim_avg_n_Wrap(tiw_CNRM_future0, 0)

     tiw_GFDL_hist      = dim_avg_n_Wrap(tiw_GFDL_hist0, 0)
     tiw_GFDL_future    = dim_avg_n_Wrap(tiw_GFDL_future0, 0)

     tiw_Had_hist      = dim_avg_n_Wrap(tiw_Had_hist0, 0)
     tiw_Had_future    = dim_avg_n_Wrap(tiw_Had_future0, 0)

     tiw1_AWI_hist    = wgt_areaave_Wrap(tiw_AWI_hist0(:,{0:5},{220:270}), 1.0, 1.0, 0)
     tiw1_AWI_future  = wgt_areaave_Wrap(tiw_AWI_future0(:,{0:5},{220:270}), 1.0, 1.0, 0)

     tiw1_CNRM_hist    = wgt_areaave_Wrap(tiw_CNRM_hist0(:,{0:5},{220:270}), 1.0, 1.0, 0)
     tiw1_CNRM_future  = wgt_areaave_Wrap(tiw_CNRM_future0(:,{0:5},{220:270}), 1.0, 1.0, 0)

     tiw1_GFDL_hist    = wgt_areaave_Wrap(tiw_GFDL_hist0(:,{0:5},{220:270}), 1.0, 1.0, 0)
     tiw1_GFDL_future  = wgt_areaave_Wrap(tiw_GFDL_future0(:,{0:5},{220:270}), 1.0, 1.0, 0)
 
     tiw1_Had_hist    = wgt_areaave_Wrap(tiw_Had_hist0(:,{0:5},{220:270}), 1.0, 1.0, 0)
     tiw1_Had_future  = wgt_areaave_Wrap(tiw_Had_future0(:,{0:5},{220:270}), 1.0, 1.0, 0)

;==================================================================================================
stat  = 0                         ; 0=mean; 1=variance; 2=std. dev; ...
n     = 50                        ; sub-sample size;  default is n=N
nBoot = 10000                      ; bootstrap replications
N     = 50
opt   = False                     ; use defaults
if (isvar("n") .and. (n.lt.N)) then
    opt   = True
    opt@sample_size = n           ; default N
end if

 BootStrap_sst1_PD     = bootstrap_stat(tiw1_AWI_hist, stat, nBoot, 0, opt)
 BootStrap_sst1_4xCO2  = bootstrap_stat(tiw1_AWI_future, stat, nBoot, 0, opt)

 BootStrap_sst2_PD     = bootstrap_stat(tiw1_CNRM_hist, stat, nBoot, 0, opt)
 BootStrap_sst2_4xCO2  = bootstrap_stat(tiw1_CNRM_future, stat, nBoot, 0, opt)

 BootStrap_sst3_PD     = bootstrap_stat(tiw1_GFDL_hist, stat, nBoot, 0, opt)
 BootStrap_sst3_4xCO2  = bootstrap_stat(tiw1_GFDL_future, stat, nBoot, 0, opt)

 BootStrap_sst4_PD     = bootstrap_stat(tiw1_Had_hist, stat, nBoot, 0, opt)
 BootStrap_sst4_4xCO2  = bootstrap_stat(tiw1_Had_future, stat, nBoot, 0, opt)

 xBoot1_sst_PD       = BootStrap_sst1_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot1_sstavg_PD    = BootStrap_sst1_PD[1]       ; Mean of bootstrapped estimates
 xBoot1_sststd_PD    = BootStrap_sst1_PD[2]       ; Std dev of bootstrapped estimates
printVarSummary(xBoot1_sstavg_PD)
 xBoot2_sst_PD       = BootStrap_sst2_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot2_sstavg_PD    = BootStrap_sst2_PD[1]       ; Mean of bootstrapped estimates
 xBoot2_sststd_PD    = BootStrap_sst2_PD[2]       ; Std dev of bootstrapped estimates

 xBoot3_sst_PD       = BootStrap_sst3_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot3_sstavg_PD    = BootStrap_sst3_PD[1]       ; Mean of bootstrapped estimates
 xBoot3_sststd_PD    = BootStrap_sst3_PD[2]       ; Std dev of bootstrapped estimates

 xBoot4_sst_PD       = BootStrap_sst4_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot4_sstavg_PD    = BootStrap_sst4_PD[1]       ; Mean of bootstrapped estimates
 xBoot4_sststd_PD    = BootStrap_sst4_PD[2]       ; Std dev of bootstrapped estimates
;==================================================================================================
   xBoot1_sst_4xCO2    = BootStrap_sst1_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot1_sstavg_4xCO2    = BootStrap_sst1_4xCO2[1]       ; Mean of bootstrapped estimates
xBoot1_sststd_4xCO2    = BootStrap_sst1_4xCO2[2]       ; Std dev of bootstrapped estimates
print(xBoot1_sststd_4xCO2)
   xBoot2_sst_4xCO2    = BootStrap_sst2_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot2_sstavg_4xCO2    = BootStrap_sst2_4xCO2[1]       ; Mean of bootstrapped estimates
xBoot2_sststd_4xCO2    = BootStrap_sst2_4xCO2[2]       ; Std dev of bootstrapped estimates

   xBoot3_sst_4xCO2    = BootStrap_sst3_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot3_sstavg_4xCO2    = BootStrap_sst3_4xCO2[1]       ; Mean of bootstrapped estimates
xBoot3_sststd_4xCO2    = BootStrap_sst3_4xCO2[2]       ; Std dev of bootstrapped estimates

   xBoot4_sst_4xCO2    = BootStrap_sst4_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot4_sstavg_4xCO2    = BootStrap_sst4_4xCO2[1]       ; Mean of bootstrapped estimates
xBoot4_sststd_4xCO2    = BootStrap_sst4_4xCO2[2]       ; Std dev of bootstrapped estimates
;==================================================================================================

xBootLow_sst1_PD    = bootstrap_estimate(xBoot1_sst_PD, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_sst1_PD    = bootstrap_estimate(xBoot1_sst_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_sst1_PD    = bootstrap_estimate(xBoot1_sst_PD, 0.975, False)   ; 97.5% upper confidence bound
print(xBootLow_sst1_PD)
print(xBootHi_sst1_PD)
 xBootLow_sst2_PD   = bootstrap_estimate(xBoot2_sst_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_sst2_PD   = bootstrap_estimate(xBoot2_sst_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_sst2_PD   = bootstrap_estimate(xBoot2_sst_PD, 0.975, False)   ; 97.5% upper confidence bound

  xBootLow_sst3_PD  = bootstrap_estimate(xBoot3_sst_PD, 0.025, False)   ;  2.5% lower confidence bound
  xBootMed_sst3_PD  = bootstrap_estimate(xBoot3_sst_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
   xBootHi_sst3_PD  = bootstrap_estimate(xBoot3_sst_PD, 0.975, False)   ; 97.5% upper confidence bound

   xBootLow_sst4_PD = bootstrap_estimate(xBoot4_sst_PD, 0.025, False)   ;  2.5% lower confidence bound
   xBootMed_sst4_PD = bootstrap_estimate(xBoot4_sst_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
    xBootHi_sst4_PD = bootstrap_estimate(xBoot4_sst_PD, 0.975, False)   ; 97.5% upper confidence bound
;===================================================================================================
xBootLow_sst1_4xCO2    = bootstrap_estimate(xBoot1_sst_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_sst1_4xCO2    = bootstrap_estimate(xBoot1_sst_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_sst1_4xCO2    = bootstrap_estimate(xBoot1_sst_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_sst2_4xCO2    = bootstrap_estimate(xBoot2_sst_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_sst2_4xCO2    = bootstrap_estimate(xBoot2_sst_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_sst2_4xCO2    = bootstrap_estimate(xBoot2_sst_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_sst3_4xCO2    = bootstrap_estimate(xBoot3_sst_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_sst3_4xCO2    = bootstrap_estimate(xBoot3_sst_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_sst3_4xCO2    = bootstrap_estimate(xBoot3_sst_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_sst4_4xCO2    = bootstrap_estimate(xBoot4_sst_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_sst4_4xCO2    = bootstrap_estimate(xBoot4_sst_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_sst4_4xCO2    = bootstrap_estimate(xBoot4_sst_4xCO2, 0.975, False)   ; 97.5% upper confidence bound
;====================================================================================================
stdarr1 = new((/2/),"float") ;  SST
stdarr1(0) =  xBoot1_sstavg_PD
stdarr1(1) =  xBoot1_sstavg_4xCO2

stdarr2 = new((/2/),"float") ; SSH
stdarr2(0) =  xBoot2_sstavg_PD
stdarr2(1) =  xBoot2_sstavg_4xCO2

stdarr3 = new((/2/),"float") ; EKE 
stdarr3(0) =  xBoot3_sstavg_PD
stdarr3(1) =  xBoot3_sstavg_4xCO2 

stdarr4 = new((/2/),"float") ; EKE 
stdarr4(0) =  xBoot4_sstavg_PD
stdarr4(1) =  xBoot4_sstavg_4xCO2 
;=====================================================================================================
  wks = gsn_open_wks("pdf","/Users/xueaoyun/Documents/work/TIW_future/final_figures/plot_TIW_mean_intensity_pic_4xCO2_hist_future")
  wks_type = "pdf"
  wks_type@wkWidth  = 0.8
  wks_type@wkHeight = 0.2    
;=============================================================================================
  gsn_define_colormap(wks,"WhiteYellowOrangeRed")
  plot = new(12,graphic)
  res                      = True         
  res@gsnDraw              = False             ; don't draw yet
  res@gsnFrame             = False             ; don't advance frame yet
  res@gsnAddCyclic         = False    
  res@gsnMaximize          = False              ; don't maximize
  res@cnFillOn             = True               ; turn on color for contours
  res@cnLinesOn            = False              ; turn off contour lines
  res@cnLineLabelsOn       = False              ; turn off contour line labels
  res@cnMonoLineColor      = True
  res@cnSmoothingOn        = True
  res@cnSmoothingTensionF  = 5
  res@gsnContourNegLineDashPattern = 1
  res@gsnContourZeroLineThicknessF =2
  res@cnInfoLabelOn        = False                ; turn off info 
  res@cnLevelSelectionMode = "ExplicitLevels"  
 ; res@cnLevels             =(/0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6/)*2
 ; res@cnFillColors          = (/0,10,11,12,13,14,15,16,17,18,19/)
  res@cnLevels             = fspan(0.2,1.0,81)
  res@mpFillOn             = False
  res@mpLandFillColor      = "mediumpurple1"             ; color of land
  res@mpGeophysicalLineColor = "gray30"
  res@mpGeophysicalLineThicknessF = 1.2
  res@mpMaxLatF            = 10
  res@mpMinLatF            = -10
  res@mpMaxLonF            = 280
  res@mpMinLonF            = 180
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
  res@tmXBValues = (/180,200,220,240,260,280/)
  res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W","80~S~o~N~W"/)  
  res@tmXBLabelFontHeightF  = 0.012
  res@tmYLLabelFontHeightF  = 0.012
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
  res@tmYMinorGrid                = True  
  res@tmYMinorGridLineDashPattern =  2
  res@tmYMinorGridLineColor       = "grey30"
  res@tmYMinorGridThicknessF      = 0.8 
  res@tmXMinorGrid                = True     
  res@tmXMinorGridLineDashPattern =  2
  res@tmXMinorGridLineColor       = "grey30"
  res@tmXMinorGridThicknessF      = 0.8 
  res@gsnRightString        = ""
  res@gsnRightStringFontHeightF  = 0.008
  res@vpHeightF            = 0.15    ; Changes the aspect ratio
  res@vpWidthF             = 0.30
  res@vpYF    = 0.90
  res@vpXF    = 0.10      ; Make plots wider
  res@pmLabelBarWidthF     = 0.03               ; default is shorter
  res@pmLabelBarHeightF    = 0.55
  res@pmLabelBarOrthogonalPosF = 0.07
  res@pmLabelBarParallelPosF   = 0.50
  res@lbLabelFontHeightF       =0.010
  res@gsnLeftString        = "~F22~a"
  res@gsnCenterString      = "AWI-HR picontrol "
  res@gsnLeftStringFontHeightF  = 0.015
  res@gsnCenterStringFontHeightF  = 0.012
  plot(0)  = gsn_csm_contour_map(wks,tiw_AWI_hist,res)
  res@gsnLeftString = "~F22~b"
  res@gsnCenterString = "AWI-HR 4xCO2"
 res@vpYF    = 0.90
 res@vpXF    = 0.45      ; Make plots wider
  res@lbLabelBarOn         = False
  plot(1)  = gsn_csm_contour_map(wks,tiw_AWI_future,res)
  res@gsnLeftString = "~F22~d"
  res@gsnCenterString = "CNRM-HR picontrol"

  res@lbLabelBarOn         = False
  res@vpYF    = 0.70
  res@vpXF    = 0.10      ; Make plots wider
  plot(2) = gsn_csm_contour_map(wks,tiw_CNRM_hist,res)
 
  res@gsnLeftString        = "~F22~e" 
  res@gsnCenterString      = "CNRM-HR 4xCO2"
  res@lbLabelBarOn         = False
  res@pmLabelBarWidthF     = 0.03               ; default is shorter
  res@vpYF    = 0.70
  res@vpXF    = 0.45      ; Make plots wider
  res@cnLevelSelectionMode = "ExplicitLevels"  
  plot(3) = gsn_csm_contour_map(wks,tiw_CNRM_future,res)
  res@vpHeightF            = 0.15    ; Changes the aspect ratio
  res@vpWidthF             = 0.30
  res@pmLabelBarWidthF     = 0.03               ; default is shorter
  res@pmLabelBarOrthogonalPosF = 0.07
  res@pmLabelBarParallelPosF   = 0.50
  res@lbLabelFontHeightF       =0.010
  res@lbLabelBarOn         = False
  res@gsnLeftString        = "~F22~g"
  res@gsnCenterString       = "GFDL-HR picontrol"
  res@vpYF    = 0.50
  res@vpXF    = 0.10      ; Make plots wider
  plot(4)  = gsn_csm_contour_map(wks,tiw_GFDL_hist,res)
  res@gsnLeftString        = "~F22~h"
  res@gsnCenterString      = "GFDL-HR 4xCO2"
  res@lbLabelBarOn         = False
  res@vpYF    = 0.50
  res@vpXF    = 0.45      ; Make plots wider
  plot(5)  = gsn_csm_contour_map(wks,tiw_GFDL_future,res)
  res@gsnCenterString      = "HadGEM3-MM picontrol"
  res@gsnLeftString        = "~F22~j"
  res@lbLabelBarOn         = False
  res@vpYF    = 0.30
  res@vpXF    = 0.10      ; Make plots wider
  plot(6) = gsn_csm_contour_map(wks,tiw_Had_hist,res)
  res@gsnCenterString      = "HadGEM3-MM 4xCO2"
  res@gsnLeftString        = "~F22~k" 
  res@lbLabelBarOn         = False
  res@pmLabelBarWidthF     = 0.03               ; default is shorter
  res@cnLevelSelectionMode = "ExplicitLevels"  
  res@vpYF    = 0.30
  res@vpXF    = 0.45      ; Make plots wider
  plot(7) = gsn_csm_contour_map(wks,tiw_Had_future,res)
;=========================================================================
 xyres                       = True                 ; plot mods desired
 xyres@gsnFrame              = False                ; don't advance frame yet
 xyres@gsnDraw               = False
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.6                ; bring bars down to zero
 xyres@trXMinF               = 0                    ; adds space on either end
 xyres@trXMaxF               = 3                   ; of the 1st and last bars
 xyres@vpWidthF              = 0.15
 xyres@vpHeightF             = 0.15

 xyres@gsnXYBarChart         = True                 ; turn on bar chart
 xyres@gsnXYBarChartBarWidth = 0.24                  ; change bar widths
 xyres@tmXBOn                = True                 ; turn off tickmarks at bot
 xyres@tmYLOn                = True
 xyres@tmYROn                = False                 ; turn off tickmarks at bot
 xyres@tmXTOn                = False                 ; turn off tickmarks at bot
 xyres@tmXBMode          = "Explicit"         ; explicit labels
 xyres@tmXBValues        = (/1,2,3/)
 xyres@tmXBLabels = (/"PI","4xCO2","" /)
 xyres@tmXBLabelFontHeightF = 0.012
 xyres@tmYLMode          = "Explicit"         ; explicit labels
 xyres@tmYLValues        = (/0,0.1,0.2,0.3,0.4,0.5,0.6/)
 xyres@tmYLLabels        = (/"0","0.1","0.2","0.3","0.4","0.5","0.6"/)
 xyres@tmYLLabelFontHeightF = 0.012
 xyres@tiMainString          = " "

 xyres@tiMainFontHeightF    = 0.012
 xyres@tiYAxisFontHeightF   = 0.012
 xyres@tiXAxisFontHeightF   = 0.012
 xyres@tmXBLabelFontThicknessF  = 1.5
 xyres@tmYLLabelFontThicknessF  = 1.5
 xyres@tmBorderThicknessF       = 1.5
 xyres@tmXBMajorLengthF         = 0.008
 xyres@tmYLMajorLengthF         = 0.008
 xyres@tmXBMajorOutwardLengthF  = 0.008
 xyres@tmYLMajorOutwardLengthF  = 0.008
 xyres@tmXBMinorLengthF         = 0.005
 xyres@tmYLMinorLengthF         = 0.005
 xyres@tmXBMinorOutwardLengthF  = 0.005
 xyres@tmYLMinorOutwardLengthF  = 0.005
 xyres@tmXBMajorThicknessF      = 1.5
 xyres@tmYLMajorThicknessF      = 1.5
 xyres@tmXBMinorThicknessF      = 1.0
 xyres@tmYLMinorThicknessF      = 1.0
 xyres@tmXBLabelDeltaF          = -0.4
 xyres@tmYLLabelDeltaF          = -0.4
 xyres@tmXTOn        = False ; turn off top tickmarks
 xyres@tmYROn        = False ; turn off right tickmarks
 xyres@tmXTBorderOn  = False ; turn off top axis
 xyres@tmYRBorderOn  = False ; turn off right axis

 xyres@vpXF    = 0.80
 xyres@vpYF    = 0.90
 xyres@gsnXYBarChartColors = (/"cyan3","orange3"/)
 xyres@gsnCenterString       = "  "
 xyres@gsnLeftString         = "~F22~c"
 xyres@gsnLeftStringFontHeightF  = 0.015
 xyres@gsnCenterStringOrthogonalPosF = -0.10
 xyres@gsnCenterStringParallelPosF = 0.475
 plot(8)  = gsn_csm_xy(wks,ispan(1,2,1),stdarr1(:),xyres)   ; advancing the frame
 xyres@vpXF    = 0.80
 xyres@vpYF    = 0.70
 xyres@gsnCenterString       = "   "
 xyres@gsnCenterStringParallelPosF = 0.45
 xyres@gsnLeftString         = "~F22~f"
 plot(9)  = gsn_csm_xy(wks,ispan(1,2,1),stdarr2(:),xyres)   ; advancing the frame
 xyres@vpXF    = 0.80
 xyres@vpYF    = 0.50
 xyres@gsnLeftString         = "~F22~i "
 xyres@gsnCenterString       = "   "
 xyres@gsnCenterStringParallelPosF = 0.45
 plot(10)  = gsn_csm_xy(wks,ispan(1,2,1),stdarr3(:),xyres)   ; advancing the frame
 xyres@vpXF    = 0.80
 xyres@vpYF    = 0.30
 xyres@gsnLeftString         = "~F22~l"
 xyres@gsnCenterString       = " "
 xyres@gsnCenterStringParallelPosF = 0.45
 plot(11)  = gsn_csm_xy(wks,ispan(1,2,1),stdarr4(:),xyres)   ; advancing the frame

  polyres=True
  polyres@gsLineDashPattern=0
  polyres@gsLineColor="black"
  polyres@gsLineThicknessF=1.0
  polyres@tfPolyDrawOrder = "postDraw"

  aumm1    = gsn_add_polyline(wks,plot(8),(/1,1/),(/xBootLow_sst1_PD,xBootHi_sst1_PD/), polyres)
  aumm001  = gsn_add_polyline(wks,plot(8),(/2,2/),(/xBootLow_sst1_4xCO2,xBootHi_sst1_4xCO2/), polyres)

  bumm1    = gsn_add_polyline(wks,plot(9),(/1,1/),(/xBootLow_sst2_PD,xBootHi_sst2_PD/), polyres)
  bumm001  = gsn_add_polyline(wks,plot(9),(/2,2/),(/xBootLow_sst2_4xCO2,xBootHi_sst2_4xCO2/), polyres)

  cumm1    = gsn_add_polyline(wks,plot(10),(/1,1/),(/xBootLow_sst3_PD,xBootHi_sst3_PD/), polyres)
  cumm001  = gsn_add_polyline(wks,plot(10),(/2,2/),(/xBootLow_sst3_4xCO2,xBootHi_sst3_4xCO2/), polyres)

  dumm1    = gsn_add_polyline(wks,plot(11),(/1,1/),(/xBootLow_sst4_PD,xBootHi_sst4_PD/), polyres)
  dumm001  = gsn_add_polyline(wks,plot(11),(/2,2/),(/xBootLow_sst4_4xCO2,xBootHi_sst4_4xCO2/), polyres)

  wid = 0.24/2
  eumm1   = gsn_add_polyline(wks,plot(8),(/1-wid,1+wid/),(/xBootLow_sst1_PD,xBootLow_sst1_PD/), polyres)
  eumm01  = gsn_add_polyline(wks,plot(8),(/1-wid,1+wid/),(/xBootHi_sst1_PD,xBootHi_sst1_PD/), polyres)

  fumm1   = gsn_add_polyline(wks,plot(8),(/2-wid,2+wid/),(/xBootLow_sst1_4xCO2,xBootLow_sst1_4xCO2/), polyres)
  fumm01  = gsn_add_polyline(wks,plot(8),(/2-wid,2+wid/),(/xBootHi_sst1_4xCO2,xBootHi_sst1_4xCO2/), polyres)
;----------------------------------------------------------------------------------------------------
eumm2   = gsn_add_polyline(wks,plot(9),(/1-wid,1+wid/),(/xBootLow_sst2_PD,xBootLow_sst2_PD/), polyres)
eumm02  = gsn_add_polyline(wks,plot(9),(/1-wid,1+wid/),(/xBootHi_sst2_PD,xBootHi_sst2_PD/), polyres)

fumm2   = gsn_add_polyline(wks,plot(9),(/2-wid,2+wid/),(/xBootLow_sst2_4xCO2,xBootLow_sst2_4xCO2/), polyres)
fumm02  = gsn_add_polyline(wks,plot(9),(/2-wid,2+wid/),(/xBootHi_sst2_4xCO2,xBootHi_sst2_4xCO2/), polyres)
;-------------------------------------------------------------------------------
eumm3   = gsn_add_polyline(wks,plot(10),(/1-wid,1+wid/),(/xBootLow_sst3_PD,xBootLow_sst3_PD/), polyres)
eumm03  = gsn_add_polyline(wks,plot(10),(/1-wid,1+wid/),(/xBootHi_sst3_PD,xBootHi_sst3_PD/), polyres)

fumm3   = gsn_add_polyline(wks,plot(10),(/2-wid,2+wid/),(/xBootLow_sst3_4xCO2,xBootLow_sst3_4xCO2/), polyres)
fumm03  = gsn_add_polyline(wks,plot(10),(/2-wid,2+wid/),(/xBootHi_sst3_4xCO2,xBootHi_sst3_4xCO2/), polyres)
;-----------------------------------------------------------------------------------
eumm4   = gsn_add_polyline(wks,plot(11),(/1-wid,1+wid/),(/xBootLow_sst4_PD,xBootLow_sst4_PD/), polyres)
eumm04  = gsn_add_polyline(wks,plot(11),(/1-wid,1+wid/),(/xBootHi_sst4_PD,xBootHi_sst4_PD/), polyres)

fumm4   = gsn_add_polyline(wks,plot(11),(/2-wid,2+wid/),(/xBootLow_sst4_4xCO2,xBootLow_sst4_4xCO2/), polyres)
fumm04  = gsn_add_polyline(wks,plot(11),(/2-wid,2+wid/),(/xBootHi_sst4_4xCO2,xBootHi_sst4_4xCO2/), polyres)
;------------------------------------------------------------------------------------------------------
;=========================================================================
draw(plot(0))
draw(plot(1))
draw(plot(2))
draw(plot(3))
draw(plot(4))
draw(plot(5))
draw(plot(6))
draw(plot(7))
draw(plot(8))
draw(plot(9))
draw(plot(10))
draw(plot(11))
end 
