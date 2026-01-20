load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/Users/xueaoyun/Documents/work/libfun/my_utilities.ncl"

begin
;=======================================================================
  datadir1 = "/Volumes/oly2/CESM/energy/mon/"
  fil1    = systemfunc("ls " + datadir1 + "TIW_CKE_01-06.PD.mon.nc") 
  f1      = addfile(fil1, "r")
  CKE_PD  = f1->CKE_PD
  CKE_PD_mon = dim_avg_n_Wrap(CKE_PD(:,0:5,:,:), 1)
  datadir2 = "/Volumes/oly2/CESM/energy/mon/"
  fil2    = systemfunc("ls " + datadir2 + "TIW_CKE_01-06.2xCO2.mon.nc") 
  f2      = addfile(fil2, "r")
  CKE_2xCO2  = f2->CKE_2xCO2
  CKE_2xCO2_mon = dim_avg_n_Wrap(CKE_2xCO2(:,0:5,:,:), 1)
  datadir3 = "/Volumes/oly2/CESM/energy/mon/"
  fil3    = systemfunc("ls " + datadir3 + "TIW_CKE_01-06.4xCO2.mon.nc") 
  f3      = addfile(fil3, "r")
  CKE_4xCO2  =f3->CKE_4xCO2
  CKE_4xCO2_mon = dim_avg_n_Wrap(CKE_4xCO2(:,0:5,:,:), 1)

  CKE_PD_Clm    = clmMonTLL(CKE_PD_mon)  
  CKE_2xCO2_Clm = clmMonTLL(CKE_2xCO2_mon)
  CKE_4xCO2_Clm = clmMonTLL(CKE_4xCO2_mon)

  CKE_PD_anom    = calcMonAnomTLL(CKE_PD_mon,  CKE_PD_Clm)
  CKE_2xCO2_anom = calcMonAnomTLL(CKE_2xCO2_mon,  CKE_2xCO2_Clm)
  CKE_4xCO2_anom = calcMonAnomTLL(CKE_4xCO2_mon,  CKE_4xCO2_Clm)

  CKE_PD_mean    = dim_avg_n_Wrap(CKE_PD_Clm,0)
  CKE_2xCO2_mean = dim_avg_n_Wrap(CKE_2xCO2_Clm,0)
  CKE_4xCO2_mean = dim_avg_n_Wrap(CKE_4xCO2_Clm,0)
  printVarSummary(CKE_PD_mean)

  CKE_PD_mon1     = wgt_areaave_Wrap(CKE_PD_mon(:,{0:2},{200:250}), 1.0, 1.0, 0)
  CKE_2xCO2_mon1  = wgt_areaave_Wrap(CKE_2xCO2_mon(:,{0:2},{200:250}), 1.0, 1.0, 0)
  CKE_4xCO2_mon1  = wgt_areaave_Wrap(CKE_4xCO2_mon(:,{0:2},{200:250}), 1.0, 1.0, 0)

  CKE_PD_mon2     = wgt_areaave_Wrap(CKE_PD_mon(:,{3:5},{200:250}), 1.0, 1.0, 0)
  CKE_2xCO2_mon2  = wgt_areaave_Wrap(CKE_2xCO2_mon(:,{3:5},{200:250}), 1.0, 1.0, 0)
  CKE_4xCO2_mon2  = wgt_areaave_Wrap(CKE_4xCO2_mon(:,{3:5},{200:250}), 1.0, 1.0, 0)

  stat  = 0                         ; 0=mean; 1=variance; 2=std. dev; ...
  n     = 30                        ; sub-sample size;  default is n=N
  nBoot = 10000                      ; bootstrap replications
  N     = 30
  opt   = False                     ; use defaults
  if (isvar("n") .and. (n.lt.N)) then
      opt   = True
      opt@sample_size = n           ; default N
  end if

 BootStrap_CKE_PD     = bootstrap_stat(CKE_PD_mon1, stat, nBoot, 0, opt)
 BootStrap_CKE_2xCO2  = bootstrap_stat(CKE_2xCO2_mon1, stat, nBoot, 0, opt)
 BootStrap_CKE_4xCO2  = bootstrap_stat(CKE_4xCO2_mon1, stat, nBoot, 0, opt)

 BootStrap_CKE1_PD     = bootstrap_stat(CKE_PD_mon2, stat, nBoot, 0, opt)
 BootStrap_CKE1_2xCO2  = bootstrap_stat(CKE_2xCO2_mon2, stat, nBoot, 0, opt)
 BootStrap_CKE1_4xCO2  = bootstrap_stat(CKE_4xCO2_mon2, stat, nBoot, 0, opt)

;=====================================================
 xBoot_CKE_PD       = BootStrap_CKE_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_CKEavg_PD    = BootStrap_CKE_PD[1]       ; Mean of bootstrapped estimates
 xBoot_CKEstd_PD    = BootStrap_CKE_PD[2]       ; Std dev of bootstrapped estimates

 xBoot_CKE_2xCO2       = BootStrap_CKE_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_CKEavg_2xCO2    = BootStrap_CKE_2xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_CKEstd_2xCO2    = BootStrap_CKE_2xCO2[2]       ; Std dev of bootstrapped estimates

 xBoot_CKE_4xCO2       = BootStrap_CKE_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_CKEavg_4xCO2    = BootStrap_CKE_4xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_CKEstd_4xCO2    = BootStrap_CKE_4xCO2[2]       ; Std dev of bootstrapped estimates

 xBootLow_CKE_PD    = bootstrap_estimate(xBoot_CKE_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_CKE_PD    = bootstrap_estimate(xBoot_CKE_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_CKE_PD    = bootstrap_estimate(xBoot_CKE_PD, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_CKE_2xCO2 = bootstrap_estimate(xBoot_CKE_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_CKE_2xCO2 = bootstrap_estimate(xBoot_CKE_2xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_CKE_2xCO2 = bootstrap_estimate(xBoot_CKE_2xCO2, 0.975, False)   ; 97.5% upper confidence bound
 
 xBootLow_CKE_4xCO2 = bootstrap_estimate(xBoot_CKE_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_CKE_4xCO2 = bootstrap_estimate(xBoot_CKE_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_CKE_4xCO2 = bootstrap_estimate(xBoot_CKE_4xCO2, 0.975, False)   ; 97.5% upper confidence bound
;-------------------------------------------------------------------------
  xBoot_CKE1_PD       = BootStrap_CKE1_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
  xBoot_CKEavg1_PD    = BootStrap_CKE1_PD[1]       ; Mean of bootstrapped estimates
  xBoot_CKEstd1_PD    = BootStrap_CKE1_PD[2]       ; Std dev of bootstrapped estimates
 
  xBoot_CKE1_2xCO2       = BootStrap_CKE1_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
  xBoot_CKEavg1_2xCO2    = BootStrap_CKE1_2xCO2[1]       ; Mean of bootstrapped estimates
  xBoot_CKEstd1_2xCO2    = BootStrap_CKE1_2xCO2[2]       ; Std dev of bootstrapped estimates
 
  xBoot_CKE1_4xCO2       = BootStrap_CKE1_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
  xBoot_CKEavg1_4xCO2    = BootStrap_CKE1_4xCO2[1]       ; Mean of bootstrapped estimates
  xBoot_CKEstd1_4xCO2    = BootStrap_CKE1_4xCO2[2]       ; Std dev of bootstrapped estimates
 
  xBootLow_CKE1_PD    = bootstrap_estimate(xBoot_CKE1_PD, 0.025, False)   ;  2.5% lower confidence bound
  xBootMed_CKE1_PD    = bootstrap_estimate(xBoot_CKE1_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
   xBootHi_CKE1_PD    = bootstrap_estimate(xBoot_CKE1_PD, 0.975, False)   ; 97.5% upper confidence bound
 
  xBootLow_CKE1_2xCO2 = bootstrap_estimate(xBoot_CKE1_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
  xBootMed_CKE1_2xCO2 = bootstrap_estimate(xBoot_CKE1_2xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
   xBootHi_CKE1_2xCO2 = bootstrap_estimate(xBoot_CKE1_2xCO2, 0.975, False)   ; 97.5% upper confidence bound
  
  xBootLow_CKE1_4xCO2 = bootstrap_estimate(xBoot_CKE1_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
  xBootMed_CKE1_4xCO2 = bootstrap_estimate(xBoot_CKE1_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
   xBootHi_CKE1_4xCO2 = bootstrap_estimate(xBoot_CKE1_4xCO2, 0.975, False)   ; 97.5% upper confidence bound
;=========================================================================
  datadir1 = "/Volumes/oly2/CESM/energy/mon/"
  fil1    = systemfunc("ls " + datadir1 + "TIW_APE_01-06.PD.mon.nc") 
  f1      = addfile(fil1, "r")
  APE_PD0  =f1->APE_PD(:,:,:,:)
  APE_PD  = APE_PD0 * 10000.
  copy_VarMeta(APE_PD0, APE_PD)
  APE_PD_mon = dim_avg_n_Wrap(APE_PD(:,0:5,:,:), 1)

  datadir2 = " /Volumes/oly2/CESM/energy/mon/"
  fil2    = systemfunc("ls " + datadir2 + "TIW_APE_01-06.2xCO2.mon.nc") 
  f2      = addfile(fil2, "r")
  APE_2xCO20 = f2->APE_2xCO2
  APE_2xCO2  = APE_2xCO20 * 10000.
  copy_VarMeta(APE_2xCO20, APE_2xCO2)
  APE_2xCO2_mon = dim_avg_n_Wrap(APE_2xCO2(:,0:5,:,:), 1)

  datadir3 = "/Volumes/oly2/CESM/energy/mon/"
  fil3    = systemfunc("ls " + datadir3 + "TIW_APE_01-06.4xCO2.mon.nc") 
  f3      = addfile(fil3, "r")
  APE_4xCO20  =f3->APE_4xCO2
  APE_4xCO2  = APE_4xCO20 * 10000.
  copy_VarMeta(APE_4xCO20, APE_4xCO2)
  APE_4xCO2_mon = dim_avg_n_Wrap(APE_4xCO2(:,0:5,:,:), 1)

  APE_PD_Clm    = clmMonTLL(APE_PD_mon)  
  APE_2xCO2_Clm = clmMonTLL(APE_2xCO2_mon)
  APE_4xCO2_Clm = clmMonTLL(APE_4xCO2_mon)

  APE_PD_anom    = calcMonAnomTLL(APE_PD_mon,   APE_PD_Clm)
  APE_2xCO2_anom = calcMonAnomTLL(APE_2xCO2_mon,APE_2xCO2_Clm)
  APE_4xCO2_anom = calcMonAnomTLL(APE_4xCO2_mon,APE_4xCO2_Clm)

  APE_PD_mean    = dim_avg_n_Wrap(APE_PD_Clm,0)
  APE_2xCO2_mean = dim_avg_n_Wrap(APE_2xCO2_Clm,0)
  APE_4xCO2_mean = dim_avg_n_Wrap(APE_4xCO2_Clm,0)
  printVarSummary(APE_PD_mean)
printMinMax(APE_PD_mean,0)
printMinMax(APE_2xCO2_mean, 0)
printMinMax(APE_4xCO2_mean, 0)
  APE_PD_mon1     = wgt_areaave_Wrap(APE_PD_mon(:,{0:2},{200:260}), 1.0, 1.0, 0)
  APE_2xCO2_mon1  = wgt_areaave_Wrap(APE_2xCO2_mon(:,{0:2},{200:260}), 1.0, 1.0, 0)
  APE_4xCO2_mon1  = wgt_areaave_Wrap(APE_4xCO2_mon(:,{0:2},{200:260}), 1.0, 1.0, 0)

  APE_PD_mon2     = wgt_areaave_Wrap(APE_PD_mon(:,{3:5},{200:250}), 1.0, 1.0, 0)
  APE_2xCO2_mon2  = wgt_areaave_Wrap(APE_2xCO2_mon(:,{3:5},{200:250}), 1.0, 1.0, 0)
  APE_4xCO2_mon2  = wgt_areaave_Wrap(APE_4xCO2_mon(:,{3:5},{200:250}), 1.0, 1.0, 0)
;---------------------------------------------------------------------
  BootStrap_APE_PD     = bootstrap_stat(APE_PD_mon1, stat, nBoot, 0, opt)
  BootStrap_APE_2xCO2  = bootstrap_stat(APE_2xCO2_mon1, stat, nBoot, 0, opt)
  BootStrap_APE_4xCO2  = bootstrap_stat(APE_4xCO2_mon1, stat, nBoot, 0, opt)

  BootStrap_APE1_PD     = bootstrap_stat(APE_PD_mon2, stat, nBoot, 0, opt)
  BootStrap_APE1_2xCO2  = bootstrap_stat(APE_2xCO2_mon2, stat, nBoot, 0, opt)
  BootStrap_APE1_4xCO2  = bootstrap_stat(APE_4xCO2_mon2, stat, nBoot, 0, opt)

  xBoot_APE_PD       = BootStrap_APE_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
  xBoot_APEavg_PD    = BootStrap_APE_PD[1]       ; Mean of bootstrapped estimates
  xBoot_APEstd_PD    = BootStrap_APE_PD[2]       ; Std dev of bootstrapped estimates
 
  xBoot_APE_2xCO2       = BootStrap_APE_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
  xBoot_APEavg_2xCO2    = BootStrap_APE_2xCO2[1]       ; Mean of bootstrapped estimates
  xBoot_APEstd_2xCO2    = BootStrap_APE_2xCO2[2]       ; Std dev of bootstrapped estimates
 
  xBoot_APE_4xCO2       = BootStrap_APE_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
  xBoot_APEavg_4xCO2    = BootStrap_APE_4xCO2[1]       ; Mean of bootstrapped estimates
  xBoot_APEstd_4xCO2    = BootStrap_APE_4xCO2[2]       ; Std dev of bootstrapped estimates

  xBoot_APE1_PD       = BootStrap_APE1_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
  xBoot_APEavg1_PD    = BootStrap_APE1_PD[1]       ; Mean of bootstrapped estimates
  xBoot_APEstd1_PD    = BootStrap_APE1_PD[2]       ; Std dev of bootstrapped estimates
 
  xBoot_APE1_2xCO2       = BootStrap_APE1_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
  xBoot_APEavg1_2xCO2    = BootStrap_APE1_2xCO2[1]       ; Mean of bootstrapped estimates
  xBoot_APEstd1_2xCO2    = BootStrap_APE1_2xCO2[2]       ; Std dev of bootstrapped estimates
 
  xBoot_APE1_4xCO2       = BootStrap_APE1_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
  xBoot_APEavg1_4xCO2    = BootStrap_APE1_4xCO2[1]       ; Mean of bootstrapped estimates
  xBoot_APEstd1_4xCO2    = BootStrap_APE1_4xCO2[2]       ; Std dev of bootstrapped estimates
;-----------------------------------------------------------------------------------------------------------
  xBootLow_APE_PD    = bootstrap_estimate(xBoot_APE_PD, 0.025, False)   ;  2.5% lower confidence bound
  xBootMed_APE_PD    = bootstrap_estimate(xBoot_APE_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
   xBootHi_APE_PD    = bootstrap_estimate(xBoot_APE_PD, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_APE_2xCO2 = bootstrap_estimate(xBoot_APE_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_APE_2xCO2 = bootstrap_estimate(xBoot_APE_2xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_APE_2xCO2 = bootstrap_estimate(xBoot_APE_2xCO2, 0.975, False)   ; 97.5% upper confidence bound
 
 xBootLow_APE_4xCO2 = bootstrap_estimate(xBoot_APE_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_APE_4xCO2 = bootstrap_estimate(xBoot_APE_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_APE_4xCO2 = bootstrap_estimate(xBoot_APE_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

  xBootLow_APE1_PD    = bootstrap_estimate(xBoot_APE1_PD, 0.025, False)   ;  2.5% lower confidence bound
  xBootMed_APE1_PD    = bootstrap_estimate(xBoot_APE1_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
   xBootHi_APE1_PD    = bootstrap_estimate(xBoot_APE1_PD, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_APE1_2xCO2 = bootstrap_estimate(xBoot_APE1_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_APE1_2xCO2 = bootstrap_estimate(xBoot_APE1_2xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_APE1_2xCO2 = bootstrap_estimate(xBoot_APE1_2xCO2, 0.975, False)   ; 97.5% upper confidence bound
 
 xBootLow_APE1_4xCO2 = bootstrap_estimate(xBoot_APE1_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_APE1_4xCO2 = bootstrap_estimate(xBoot_APE1_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_APE1_4xCO2 = bootstrap_estimate(xBoot_APE1_4xCO2, 0.975, False)   ; 97.5% upper confidence bound
;------------------------------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------------------
stdarr1 = new((/3/),"double") ;  SST
stdarr1(0) =  xBoot_CKEavg_PD
stdarr1(1) =  xBoot_CKEavg_2xCO2 
stdarr1(2) =  xBoot_CKEavg_4xCO2

stdarr2 = new((/3/),"double") ;  SST
stdarr2(0) =  xBoot_CKEavg1_PD
stdarr2(1) =  xBoot_CKEavg1_2xCO2 
stdarr2(2) =  xBoot_CKEavg1_4xCO2

stdarr3 = new((/3/),"float") ; SSH
stdarr3(0) =  xBoot_APEavg_PD
stdarr3(1) =  xBoot_APEavg_2xCO2 
stdarr3(2) =  xBoot_APEavg_4xCO2

stdarr4 = new((/3/),"float") ; SSH
stdarr4(0) =  xBoot_APEavg1_PD
stdarr4(1) =  xBoot_APEavg1_2xCO2 
stdarr4(2) =  xBoot_APEavg1_4xCO2

  ;=====================================================================================

pltDir =  "/Users/xueaoyun/Documents/work/TIW_future/revision1/main_figure/"
   pltName = "PD_2xCO2_4xCO2_mean_CKE_APE"
   pltPath= pltDir+"/"+pltName
   pltTypeLocal = "pdf"
  wks  = gsn_open_wks(pltTypeLocal, pltPath)  
 ;gsn_define_colormap(wks,"sunshine_9lev") 
 gsn_define_colormap(wks,"BlueWhiteOrangeRed") 
res               = True                      ; plot mods desired
res@gsnDraw      =  False                   ; do not draw the plot
res@gsnFrame     =  False 
res@gsnSpreadColors  = True                      ; use full colormap
res@gsnSpreadColorStart = 140
res@gsnSpreadColorEnd = 240
res@gsnAddCyclic         = False        ; data not cyclic
;res@cnMissingValFillColor = "black"
;res@gsnScalarContour        = True ;important!
res@cnLevelSelectionMode  = "ExplicitLevels"
res@cnLevels              = fspan(0.1,0.5,41)*80
;print(res@cnLevels)
;res@cnFillColors         = (/2,11,30,50,70,90,110,130,150,170,190,230/)
;res@cnFillColors          = (/0,2,3,4,5,6,7,8,9,10,11/)
res@cnFillOn         = True                     ; turn on color
res@cnLinesOn        = False 
res@cnLineLabelsOn   = False                     ; no contour line labels
res@cnInfoLabelOn    = False                     ; no info label

res@mpMinLatF     = -5
res@mpMaxLatF     = 10
res@mpMinLonF     = 180
res@mpMaxLonF     = 280
res@mpShapeMode = "FreeAspect"
res@mpCenterLonF  =230.
res@tmXBMinorOn = False
res@tmYLMinorOn = False
res@tmYLTickSpacingF = 2
res@tmXBTickSpacingF = 5
res@tmBorderThicknessF = 2
res@tmXBMajorThicknessF= 2
res@tmYLMajorThicknessF= 2
res@tmXBLabelFontHeightF = 0.01
res@tmYLLabelFontHeightF = 0.01
res@tmXBMajorLengthF = -0.01
res@tmYLMajorLengthF = -0.01
res@tmXBLabelDeltaF = -0.6
res@tmYLLabelDeltaF = -0.8
;res@trXMinF= 142

res@tmYLMode        = "Explicit"  
res@tmYLValues      = (/-10,-8,-6,-4,-2,0,2,4,6,8,10/)
res@tmYLLabels      = (/"10~S~o~N~S","8~S~o~N~S","6~S~o~N~S","4~S~o~N~S","2~S~o~N~S","0","2~S~o~N~N","4~S~o~N~N","6~S~o~N~N","8~S~o~N~N","10~S~o~N~N"/)

res@tmXBMode        = "Explicit"  
res@tmXBValues            = (/180,200,220,240,260,280/)
;res@tmXBLabels            = (/"180","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W","80~S~o~N~W"/)
res@tmXBLabels            = (/"","","","","",""/)
res@gsnLeftString  =""
res@gsnRightString =""
res@tiYAxisFontHeightF = 0.015
res@tiXAxisFontHeightF = 0.015
;res@gsnLeftStringFontHeightF = 0.015

res@lbLabelBarOn =False
res@lbOrientation         = "vertical"
res@pmLabelBarWidthF      = 0.06               ; default is shorter
res@pmLabelBarHeightF     = 0.16              ; default is taller
res@pmLabelBarParallelPosF   = 0.5    ; 0.5 This is the right edge of the plot.
res@pmLabelBarOrthogonalPosF = 0.06    ; 0.5 This is the bottom edge of the plot.
res@lbLabelFontHeightF    = 0.01              ; default is HUGE
;res@lbTitleString         = "W/m~S~3~N"
;res@lbTitleFontHeightF = 0.015
res@lbBoxSeparatorLinesOn= False  
res@vpHeightF            = 0.15      ; Changes the aspect ratio
res@vpWidthF             = 0.25
res@vpXF    = 0.05      ; Make plots wider
res@vpYF    = 0.95
res@tmXBLabelsOn = True
plot0 = gsn_csm_contour_map(wks,CKE_PD_mean,res)
;****************************************************
res@tmYLLabelsOn = False
res@vpHeightF            = 0.15      ; Changes the aspect ratio
res@vpWidthF             = 0.25
res@vpXF    = 0.33      ; Make plots wider
res@vpYF    = 0.95
plot1 = gsn_csm_contour_map(wks,CKE_2xCO2_mean,res)
;*****************************************************
res@vpHeightF            = 0.15      ; Changes the aspect ratio
res@vpWidthF             = 0.25
res@vpXF    = 0.61      ; Make plots wider
res@vpYF    = 0.95
res@lbLabelBarOn =True
plot2 = gsn_csm_contour_map(wks,CKE_4xCO2_mean,res)
delete(res@cnLevels)
res@tmXBLabels            = (/"180","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W","80~S~o~N~W"/)
res@tmYLLabelsOn = True
res@tmYLMode        = "Explicit"  
res@tmYLValues      = (/-10,-8,-6,-4,-2,0,2,4,6,8,10/)
res@tmYLLabels      = (/"10~S~o~N~S","8~S~o~N~S","6~S~o~N~S","4~S~o~N~S","2~S~o~N~S","0","2~S~o~N~N","4~S~o~N~N","6~S~o~N~N","8~S~o~N~N","10~S~o~N~N"/)
res@cnLevels              = fspan(0.1,0.6,51)*30
res@lbLabelStride   = 10
res@vpXF    = 0.05      ; Make plots wider
res@vpYF    = 0.75
res@tmXBLabelsOn = True
res@lbLabelBarOn =False
plot3 = gsn_csm_contour_map(wks,APE_PD_mean,res)
;****************************************************
res@tmYLLabelsOn = False
res@vpHeightF            = 0.15      ; Changes the aspect ratio
res@vpWidthF             = 0.25
res@vpXF    = 0.33      ; Make plots wider
res@vpYF    = 0.75
plot4 = gsn_csm_contour_map(wks,APE_2xCO2_mean,res)
;*****************************************************
res@vpHeightF            = 0.15      ; Changes the aspect ratio
res@vpWidthF             = 0.25
res@vpXF    = 0.61      ; Make plots wider
res@vpYF    = 0.75
res@lbLabelBarOn =True
plot5 = gsn_csm_contour_map(wks,APE_4xCO2_mean,res)
;-------------------------------------------------------------------------   

 
  pres                  = True           ; polyline mods desired
  pres@gsLineThicknessF =2.0            ; line thickness   
  pres@gsLineColor      = "black"         ; color of lines
  pres@gsLineDashPattern =1
  pres@tfPolyDrawOrder       = "PostDraw" 
  ;dum = gsn_add_polyline(wks,plot, (/0, 24/) , (/7.5,7.5/) ,pres)

  
  dum01 = gsn_add_polyline(wks,plot0, (/180,280/) , (/0,0/) ,pres)
  dum02 = gsn_add_polyline(wks,plot1, (/180,280/) , (/0,0/) ,pres)
  dum03 = gsn_add_polyline(wks,plot2, (/180,280/) , (/0,0/) ,pres)
  dum04 = gsn_add_polyline(wks,plot3, (/180,280/) , (/0,0/) ,pres)
  dum05 = gsn_add_polyline(wks,plot4, (/180,280/) , (/0,0/) ,pres)
  dum06 = gsn_add_polyline(wks,plot5, (/180,280/) , (/0,0/) ,pres)


  dum001 = gsn_add_polyline(wks,plot0, (/180,280/) , (/5,5/) ,pres)
  dum002 = gsn_add_polyline(wks,plot1, (/180,280/) , (/5,5/) ,pres)
  dum003 = gsn_add_polyline(wks,plot2, (/180,280/) , (/5,5/) ,pres)
  dum004 = gsn_add_polyline(wks,plot3, (/180,280/) , (/5,5/) ,pres)
  dum005 = gsn_add_polyline(wks,plot4, (/180,280/) , (/5,5/) ,pres)
  dum006 = gsn_add_polyline(wks,plot5, (/180,280/) , (/5,5/) ,pres)


  txres                       = True
  txres@txPerimOn             = False  
  txres@txBackgroundFillColor = "White"
  txres@txFontHeightF         = 0.013
  txres@txFontThicknessF      = 1.5
  txid1     = gsn_create_text(wks, "~F22~a    PD-BTR", txres)
  txid2     = gsn_create_text(wks, "~F22~b    2xCO2-BTR", txres)
  txid3     = gsn_create_text(wks, "~F22~c    4xCO2-BTR", txres)
  txid03    = gsn_create_text(wks, "10~S~-6~N~Jm~S~-3~N~s~S~-1~N", txres)
  txid4     = gsn_create_text(wks, "~F22~d    PD-BCR", txres)
  txid5     = gsn_create_text(wks, "~F22~e    2xCO2-BCR", txres)
  txid6     = gsn_create_text(wks, "~F22~f    4xCO2-BCR", txres)
  txid06    = gsn_create_text(wks,  "10~S~-6~N~Jm~S~-3~N~s~S~-1~N", txres)
 
  amres                  = True
  amres@amParallelPosF   = -0.485     ; 0.5 This is the right edge of the plot.
  amres@amOrthogonalPosF = -0.70     ; 0.5 This is the bottom edge of the plot.
  amres@amJust           = "TopLeft"

  annoid1   = gsn_add_annotation(plot0, txid1, amres)
  annoid2   = gsn_add_annotation(plot1, txid2, amres)
  annoid3   = gsn_add_annotation(plot2, txid3, amres) 
  annoid01  = gsn_add_annotation(plot3, txid4, amres)
  annoid02  = gsn_add_annotation(plot4, txid5, amres)
  annoid03  = gsn_add_annotation(plot5, txid6, amres) 

  txres@txFontHeightF         = 0.015
  amres@amParallelPosF   = 0.51     ; 0.5 This is the right edge of the plot.
  amres@amOrthogonalPosF = -0.72     ; 0.5 This is the bottom edge of the plot.
  txres@txFontHeightF         = 0.008
  annoid03 = gsn_add_annotation(plot2, txid03, amres)  
  annoid03 = gsn_add_annotation(plot5, txid06, amres)  

  plres=True
plres@gsLineThicknessF=2.5
; plres@gsFillColor     = -1
plres@gsLineDashPattern= 2
plres@gsLineColor="blue"
plres@gsLineOpacityF=100
w1=360-160
e1=360-110

w2=360-160
e2=360-110


s1=3
n1=5

s2 = 0
n2 = 2

ddum0=gsn_add_polyline(wks, plot0, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum1=gsn_add_polyline(wks, plot1, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum2=gsn_add_polyline(wks, plot2, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum3=gsn_add_polyline(wks, plot3, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum4=gsn_add_polyline(wks, plot4, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum5=gsn_add_polyline(wks, plot5, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)

ddum6=gsn_add_polyline(wks, plot0, (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum7=gsn_add_polyline(wks, plot1, (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum8=gsn_add_polyline(wks, plot2, (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum9=gsn_add_polyline(wks, plot3, (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum10=gsn_add_polyline(wks,plot4, (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum11=gsn_add_polyline(wks,plot5, (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)

draw(plot0)
draw(plot1)
draw(plot2)
draw(plot3)
draw(plot4)
draw(plot5)


xyres                       = True                 ; plot mods desired
 xyres@gsnFrame              = False                ; don't advance frame yet
 xyres@gsnDraw               = False
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 50                 ; bring bars down to zero
 xyres@trXMinF               = 0                    ; adds space on either end
 xyres@trXMaxF               = 4                   ; of the 1st and last bars
 xyres@vpWidthF              = 0.15
 xyres@vpHeightF             = 0.15

 xyres@gsnXYBarChart         = True                 ; turn on bar chart
 xyres@gsnXYBarChartBarWidth = 0.24                  ; change bar widths
 xyres@tmXBOn                = True                 ; turn off tickmarks at bot
 xyres@tmYLOn                = True
 xyres@tmYROn                = False                 ; turn off tickmarks at bot
 xyres@tmXTOn                = False                 ; turn off tickmarks at bot
 xyres@tmXBMode          = "Explicit"         ; explicit labels
 xyres@tmXBValues        = (/1,2,3,4/)
 xyres@tmXBLabels = (/"PD","2xCO2","4xCO2","" /)
 xyres@tmXBLabelFontHeightF = 0.006
 xyres@tmYLMode          = "Explicit"         ; explicit labels
 xyres@tmYLValues        = (/0,0.1,0.2,0.3,0.4,0.5/)*100
 xyres@tmYLLabels        = (/"0","10","20","30","40","50"/)
 xyres@tmYLLabelFontHeightF = 0.008
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
 xyres@tmXTBorderOn  = True ; turn off top axis
 xyres@tmYRBorderOn  = False ; turn off right axis

 xyres@vpXF    = 0.14
 yspace=0.03
 yf=0.53
 xyres@vpYF    = yf
 xyres@gsnXYBarChartColors = (/"cyan3"/)
 xyres@gsnCenterString       = "    BTR_EQ"
 xyres@gsnLeftString         = ""
 xyres@gsnCenterStringOrthogonalPosF = -0.20+0.02
 xyres@gsnCenterStringParallelPosF = 0.45
 plot6 = gsn_csm_xy(wks,ispan(1,3,1),stdarr1(:),xyres)   ; advancing the frame
 xyres@vpXF = xyres@vpWidthF+xyres@vpXF
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 50                 ; bring bars down to zero
 xyres@gsnXYBarChartColors = (/"cyan3"/)
 xyres@gsnCenterString       = "    BTR_NEQ"
 xyres@gsnCenterStringOrthogonalPosF = -0.18+0.02
 xyres@gsnCenterStringParallelPosF = 0.45
 xyres@tmXTBorderOn  = True  ; turn off top axis
 xyres@tmYRBorderOn  = True ; turn off right axis
 xyres@tmYLBorderOn  = False ; turn off right axis
 xyres@tmYLOn        = False 
 xyres@gsnLeftString         = " "
 plot7 = gsn_csm_xy(wks,ispan(1,3,1),stdarr2(:),xyres)   ; advancing the frame

 xyres@gsnXYBarChartColors = (/"pink"/)
 xyres@vpXF = xyres@vpWidthF+xyres@vpXF+0.05
 print("xyres@vpXF="+xyres@vpXF)
 xyres@tmYLBorderOn  = True ; turn off right axis
 xyres@tmYRBorderOn  = False
 xyres@tmYLOn        = True
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 20                 ; bring bars down to zero
 xyres@gsnXYBarChartColors = (/"pink"/)
 xyres@tmYLValues        = (/0,1,2,3,4,5/)*4
 xyres@tmYLLabels        = (/"0","4","8","12","16","20"/)
 xyres@gsnCenterString       = "    BCR_EQ"
 xyres@gsnCenterStringOrthogonalPosF = -0.20+0.02
 xyres@gsnCenterStringParallelPosF = 0.45
 xyres@gsnLeftString         = ""
 plot8 = gsn_csm_xy(wks,ispan(1,3,1),stdarr3(:),xyres)   ; advancing the frame
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 20                ; bring bars down to zero
 xyres@gsnXYBarChartColors = (/"pink"/)
 xyres@vpYF = yf
 xyres@vpXF = xyres@vpWidthF+xyres@vpXF
 print("xyres@vpXF="+xyres@vpXF)
 xyres@gsnCenterString       = "BCR_NEQ"
 xyres@gsnCenterStringOrthogonalPosF = -0.18+0.02
 xyres@gsnCenterStringParallelPosF = 0.51
 xyres@gsnLeftString         = ""
 xyres@tmXTBorderOn  = True ; turn off top axis
 xyres@tmYRBorderOn  = True ; turn off right axis
 xyres@tmYLBorderOn  = False ; turn off right axis
 xyres@tmYLOn        = False 
 plot9 = gsn_csm_xy(wks,ispan(1,3,1),stdarr4(:),xyres)   ; advancing the frame

 
 txres                       = True
 txres@txPerimOn             = False  
 txres@txBackgroundFillColor = "White"
 txres@txFontHeightF         = 0.013
 txres@txFontThicknessF      = 1.5
 txid7     = gsn_create_text(wks, "~F22~g  BTR", txres)
 txid8     = gsn_create_text(wks, "~F22~h  BCR", txres)

 amres                  = True
 amres@amParallelPosF   = -0.485     ; 0.5 This is the right edge of the plot.
 amres@amOrthogonalPosF = -0.74     ; 0.5 This is the bottom edge of the plot.
 amres@amJust           = "TopLeft"

 annoid7 = gsn_add_annotation(plot6, txid7, amres)
 amres@amOrthogonalPosF = -0.73     ; 0.5 This is the bottom edge of the plot.
 annoid8 = gsn_add_annotation(plot8, txid8, amres)

  polyres=True
  polyres@gsLineDashPattern=0
  polyres@gsLineColor="black"
  polyres@gsLineThicknessF=1.0
  polyres@tfPolyDrawOrder = "postDraw"

  aumm1    = gsn_add_polyline(wks,plot6,(/1,1/),(/xBootLow_CKE_PD,xBootHi_CKE_PD/), polyres)
  aumm01   = gsn_add_polyline(wks,plot6,(/2,2/),(/xBootLow_CKE_2xCO2,xBootHi_CKE_2xCO2/), polyres)
  aumm001  = gsn_add_polyline(wks,plot6,(/3,3/),(/xBootLow_CKE_4xCO2,xBootHi_CKE_4xCO2/), polyres)

  bumm1    = gsn_add_polyline(wks,plot7,(/1,1/),(/xBootLow_CKE1_PD,xBootHi_CKE1_PD/), polyres)
  bumm01   = gsn_add_polyline(wks,plot7,(/2,2/),(/xBootLow_CKE1_2xCO2,xBootHi_CKE1_2xCO2/), polyres)
  bumm001  = gsn_add_polyline(wks,plot7,(/3,3/),(/xBootLow_CKE1_4xCO2,xBootHi_CKE1_4xCO2/), polyres)

 cumm1    = gsn_add_polyline(wks,plot8,(/1,1/),(/xBootLow_APE_PD,xBootHi_APE_PD/), polyres)
 cumm01   = gsn_add_polyline(wks,plot8,(/2,2/),(/xBootLow_APE_2xCO2,xBootHi_APE_2xCO2/), polyres)
 cumm001  = gsn_add_polyline(wks,plot8,(/3,3/),(/xBootLow_APE_4xCO2,xBootHi_APE_4xCO2/), polyres)

 dumm1    = gsn_add_polyline(wks,plot9,(/1,1/),(/xBootLow_APE1_PD,xBootHi_APE1_PD/), polyres)
 dumm01   = gsn_add_polyline(wks,plot9,(/2,2/),(/xBootLow_APE1_2xCO2,xBootHi_APE1_2xCO2/), polyres)
 dumm001  = gsn_add_polyline(wks,plot9,(/3,3/),(/xBootLow_APE1_4xCO2,xBootHi_APE1_4xCO2/), polyres)


  wid = 0.24/2
  eumm1   = gsn_add_polyline(wks,plot6,(/1-wid,1+wid/),(/xBootLow_CKE_PD,xBootLow_CKE_PD/), polyres)
  eumm01  = gsn_add_polyline(wks,plot6,(/1-wid,1+wid/),(/xBootHi_CKE_PD,xBootHi_CKE_PD/), polyres)

  fumm1   = gsn_add_polyline(wks,plot6,(/2-wid,2+wid/),(/xBootLow_CKE_2xCO2,xBootLow_CKE_2xCO2/), polyres)
  fumm01  = gsn_add_polyline(wks,plot6,(/2-wid,2+wid/),(/xBootHi_CKE_2xCO2,xBootHi_CKE_2xCO2/), polyres)

  gumm1   = gsn_add_polyline(wks,plot6,(/3-wid,3+wid/),(/xBootLow_CKE_4xCO2,xBootLow_CKE_4xCO2/), polyres)
  gumm01  = gsn_add_polyline(wks,plot6,(/3-wid,3+wid/),(/xBootHi_CKE_4xCO2,xBootHi_CKE_4xCO2/), polyres)
;----------------------------------------------------------------------------------------------------
eumm2   = gsn_add_polyline(wks,plot7,(/1-wid,1+wid/),(/xBootLow_CKE1_PD,xBootLow_CKE1_PD/), polyres)
eumm02  = gsn_add_polyline(wks,plot7,(/1-wid,1+wid/),(/xBootHi_CKE1_PD,xBootHi_CKE1_PD/), polyres)

fumm2   = gsn_add_polyline(wks,plot7,(/2-wid,2+wid/),(/xBootLow_CKE1_2xCO2,xBootLow_CKE1_2xCO2/), polyres)
fumm02  = gsn_add_polyline(wks,plot7,(/2-wid,2+wid/),(/xBootHi_CKE1_2xCO2,xBootHi_CKE1_2xCO2/), polyres)

gumm2   = gsn_add_polyline(wks,plot7,(/3-wid,3+wid/),(/xBootLow_CKE1_4xCO2,xBootLow_CKE1_4xCO2/), polyres)
gumm02  = gsn_add_polyline(wks,plot7,(/3-wid,3+wid/),(/xBootHi_CKE1_4xCO2,xBootHi_CKE1_4xCO2/), polyres)
;-----------------------------------------------------------------------------------------------------
  eumm3   = gsn_add_polyline(wks,plot8,(/1-wid,1+wid/),(/xBootLow_APE_PD,xBootLow_APE_PD/), polyres)
  eumm03  = gsn_add_polyline(wks,plot8,(/1-wid,1+wid/),(/xBootHi_APE_PD,xBootHi_APE_PD/), polyres)

  fumm3   = gsn_add_polyline(wks,plot8,(/2-wid,2+wid/),(/xBootLow_APE_2xCO2,xBootLow_APE_2xCO2/), polyres)
  fumm03  = gsn_add_polyline(wks,plot8,(/2-wid,2+wid/),(/xBootHi_APE_2xCO2,xBootHi_APE_2xCO2/), polyres)

  gumm3   = gsn_add_polyline(wks,plot8,(/3-wid,3+wid/),(/xBootLow_APE_4xCO2,xBootLow_APE_4xCO2/), polyres)
  gumm03  = gsn_add_polyline(wks,plot8,(/3-wid,3+wid/),(/xBootHi_APE_4xCO2,xBootHi_APE_4xCO2/), polyres)

;------------------------------------------------------------------------------------------------------
eumm4   = gsn_add_polyline(wks,plot9,(/1-wid,1+wid/),(/xBootLow_APE1_PD,xBootLow_APE1_PD/), polyres)
eumm04  = gsn_add_polyline(wks,plot9,(/1-wid,1+wid/),(/xBootHi_APE1_PD,xBootHi_APE1_PD/), polyres)

fumm4   = gsn_add_polyline(wks,plot9,(/2-wid,2+wid/),(/xBootLow_APE1_2xCO2,xBootLow_APE1_2xCO2/), polyres)
fumm04  = gsn_add_polyline(wks,plot9,(/2-wid,2+wid/),(/xBootHi_APE1_2xCO2,xBootHi_APE1_2xCO2/), polyres)

gumm4   = gsn_add_polyline(wks,plot9,(/3-wid,3+wid/),(/xBootLow_APE1_4xCO2,xBootLow_APE1_4xCO2/), polyres)
gumm04  = gsn_add_polyline(wks,plot9,(/3-wid,3+wid/),(/xBootHi_APE1_4xCO2,xBootHi_APE1_4xCO2/), polyres)

;====================================================================================================
  draw(plot6)
  draw(plot7)
  draw(plot8)
  draw(plot9)

end 


 
