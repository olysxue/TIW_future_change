;--------------------------------------------------------------
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
begin 
;=============================================================
; read data
;=============================================================
ncdir1  = "/Volumes/oly2/CESM/energy/"
filename1   = "TKE.PD.mon.0111-0140.20S_20N.120E_80W.levidx01_025degree_new_mean.nc"
filename2   = "KE.PD.mon.0111-0140.20S_20N.120E_80W.levidx01_new_method_025degree_mean.nc"
filename3   = "TIW_EKE.PD.mon.0111-0140.20S_20N.120E_80W.levidx01_025degree_mean.nc"
fil1        = systemfunc("ls " + ncdir1 + filename1)
fil2        = systemfunc("ls " + ncdir1 + filename2)
fil3        = systemfunc("ls " + ncdir1 + filename3)
f1          =  addfile(fil1, "r")
f2          =  addfile(fil2, "r")
f3          =  addfile(fil3, "r")
TKE1        = f1->KE(0,0,{-10:10},:)
LSKE1      = f2->UVEL(0,0,{-10:10},:)
; EKE1        = f3->KE(0,0,{-10:10},:)
EKE1 = TKE1 -LSKE1
copy_VarMeta(TKE1, EKE1)
;----------------------------------------------------------------
ncdir2  = "/Volumes/oly2/E3SM/energy/v1/"
filename01   = "KE_picontrol_mon_0001-0500_remap_equator_mean_1degree.nc"
filename02   = "LSKE_0001-0500_remap_equator_mean.nc"
filename03   = "EKE_picontrol_mon_0001-0500_remap_equator_mean.nc"
fil01        = systemfunc("ls " + ncdir2 + filename01)
fil02        = systemfunc("ls " + ncdir2 + filename02)
fil03        = systemfunc("ls " + ncdir2 + filename03)
f01          =  addfile(fil01, "r")
f02          =  addfile(fil02, "r")
f03          =  addfile(fil03, "r")
TKE2_t         = f01->KE(0,:,:)
TKE2        = TKE2_t * 1024.
copy_VarMeta(TKE2_t, TKE2)
LSKE2       = f02->uo(0,:,:)
EKE2         = f03->KE(0,:,:)
;-----------------------------------------------------------------
ncdir3  = "/Volumes/oly2/E3SM/energy/v2/"
filename001   = "KE_picontrol_mon_0001-0500_remap_equator_mean_1degree.nc"
filename002   = "LSKE_0001-0500_remap_equator_mean.nc"
filename003   = "EKE_picontrol_mon_0001-0500_remap_equator_mean.nc"
fil001        = systemfunc("ls " + ncdir3 + filename001)
fil002        = systemfunc("ls " + ncdir3 + filename002)
fil003        = systemfunc("ls " + ncdir3 + filename003)
f001          =  addfile(fil001, "r")
f002          =  addfile(fil002, "r")
f003          =  addfile(fil003, "r")
TKE3_t          = f001->kineticEnergyAtSurface(0,:,:)
TKE3        = TKE3_t * 1024.
copy_VarMeta(TKE3_t, TKE3)
LSKE3         = f002->uo(0,:,:)
EKE3          = f003->kineticEnergyAtSurface(0,:,:)

printMinMax(TKE1, 0)
printMinMax(LSKE1, 0)
printMinMax(EKE1, 0)

printMinMax(TKE2, 0)
printMinMax(LSKE2, 0)
printMinMax(EKE2, 0)

printMinMax(TKE3, 0)
printMinMax(LSKE3, 0)
printMinMax(EKE3, 0)

delete([/filename1, filename2, filename3,filename01, filename02, filename03, filename001, filename002, filename003/])
delete([/fil1, fil2, fil3, fil01, fil02, fil03, fil001, fil002, fil003/])
delete([/f1, f2, f3, f01, f02, f03, f001, f002, f003/])
;=============================================================
;read monthly data 
ncdir1  = "/Volumes/oly2/CESM/energy/"
filename1   = "TKE.PD.mon.0111-0140.20S_20N.120E_80W.levidx01_025degree_new.nc"
filename2   = "KE.PD.mon.0111-0140.20S_20N.120E_80W.levidx01_new_method_025degree.nc"
filename3   = "TIW_EKE.PD.mon.0111-0140.20S_20N.120E_80W.levidx01_025degree.nc"
fil1        = systemfunc("ls " + ncdir1 + filename1)
fil2        = systemfunc("ls " + ncdir1 + filename2)
fil3        = systemfunc("ls " + ncdir1 + filename3)
f1          =  addfile(fil1, "r")
f2          =  addfile(fil2, "r")
f3          =  addfile(fil3, "r")
TKE01       = f1->KE(:,0,{-10:10},:)
LSKE01      = f2->UVEL(:,0,{-10:10},:)
EKE01       = f3->KE(:,0,{-10:10},:)

ncdir2  = "/Volumes/oly2/E3SM/energy/v1/"
filename01   = "KE_picontrol_mon_0001-0500_remap_equator_1degree.nc"
filename02   = "LSKE_0001-0500_remap_equator.nc"
filename03   = "EKE_picontrol_mon_0001-0500_remap_equator.nc"
fil01        = systemfunc("ls " + ncdir2 + filename01)
fil02        = systemfunc("ls " + ncdir2 + filename02)
fil03        = systemfunc("ls " + ncdir2 + filename03)
f01          =  addfile(fil01, "r")
f02          =  addfile(fil02, "r")
f03          =  addfile(fil03, "r")
TKE02_t        = f01->KE(:,:,:)
TKE02        = TKE02_t * 1024.
copy_VarMeta(TKE02_t, TKE02)
LSKE02       = f02->uo(:,:,:)
EKE02        = f03->KE(:,:,:)

ncdir3  = "/Volumes/oly2/E3SM/energy/v2/"
filename001   = "KE_picontrol_mon_0001-0500_remap_equator_1degree.nc"
filename002   = "LSKE_0001-0500_remap_equator.nc"
filename003   = "EKE_picontrol_mon_0001-0500_remap_equator.nc"
fil001        = systemfunc("ls " + ncdir3 + filename001)
fil002        = systemfunc("ls " + ncdir3 + filename002)
fil003        = systemfunc("ls " + ncdir3 + filename003)
f001          =  addfile(fil001, "r")
f002          =  addfile(fil002, "r")
f003          =  addfile(fil003, "r")
TKE03_t       = f001->kineticEnergyAtSurface(:,:,:)
TKE03        = TKE03_t * 1024.
copy_VarMeta(TKE03_t, TKE03)
LSKE03       = f002->uo(:,:,:)
EKE03        = f003->kineticEnergyAtSurface(:,:,:)
;================================================================

; calculate area average 
TKE_PD1    = wgt_areaave_Wrap(TKE01(:,{0:5},{200:250}), 1.0, 1.0, 0)
TKE_PD2    = doubletofloat(wgt_areaave_Wrap(TKE02(:,{0:5},{200:250}), 1.0, 1.0, 0))
TKE_PD3    = doubletofloat(wgt_areaave_Wrap(TKE03(:,{0:5},{200:250}), 1.0, 1.0, 0))

LSKE_PD1    = wgt_areaave_Wrap(LSKE01(:,{0:5},{200:250}), 1.0, 1.0, 0)
LSKE_PD2    = doubletofloat(wgt_areaave_Wrap(LSKE02(:,{0:5},{200:250}), 1.0, 1.0, 0))
LSKE_PD3    = doubletofloat(wgt_areaave_Wrap(LSKE03(:,{0:5},{200:250}), 1.0, 1.0, 0))

EKE_PD1     = wgt_areaave_Wrap(EKE01(:,{0:5},{200:250}), 1.0, 1.0, 0)
EKE_PD2     =  doubletofloat(wgt_areaave_Wrap(EKE02(:,{0:5},{200:250}), 1.0, 1.0, 0))
EKE_PD3     =  doubletofloat(wgt_areaave_Wrap(EKE03(:,{0:5},{200:250}), 1.0, 1.0, 0))

stat  = 0                         ; 0=mean; 1=variance; 2=std. dev; ...
n     = 1800                        ; sub-sample size;  default is n=N
nBoot = 10000                      ; bootstrap replications
N     = 1800
opt   = False                     ; use defaults
if (isvar("n") .and. (n.lt.N)) then
    opt   = True
    opt@sample_size = n           ; default N
end if

BootStrap_TKE1_PD = bootstrap_stat(TKE_PD1, stat, nBoot, 0, opt)
BootStrap_TKE2_PD = bootstrap_stat(TKE_PD2, stat, nBoot, 0, opt)
BootStrap_TKE3_PD = bootstrap_stat(TKE_PD3, stat, nBoot, 0, opt)

BootStrap_LSKE1_PD = bootstrap_stat(LSKE_PD1, stat, nBoot, 0, opt)
BootStrap_LSKE2_PD = bootstrap_stat(LSKE_PD2, stat, nBoot, 0, opt)
BootStrap_LSKE3_PD = bootstrap_stat(LSKE_PD3, stat, nBoot, 0, opt)

BootStrap_EKE1_PD = bootstrap_stat(EKE_PD1, stat, nBoot, 0, opt)
BootStrap_EKE2_PD = bootstrap_stat(EKE_PD2, stat, nBoot, 0, opt)
BootStrap_EKE3_PD = bootstrap_stat(EKE_PD3, stat, nBoot, 0, opt)

xBoot_TKE1_PD       = BootStrap_TKE1_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_TKEavg1_PD    = BootStrap_TKE1_PD[1]       ; Mean of bootstrapped estimates
xBoot_TKEstd1_PD    = BootStrap_TKE1_PD[2]       ; Std dev of bootstrapped estimates

xBoot_TKE2_PD       = BootStrap_TKE2_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_TKEavg2_PD    = BootStrap_TKE2_PD[1]       ; Mean of bootstrapped estimates
xBoot_TKEstd2_PD    = BootStrap_TKE2_PD[2]       ; Std dev of bootstrapped estimates

xBoot_TKE3_PD       = BootStrap_TKE3_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_TKEavg3_PD    = BootStrap_TKE3_PD[1]       ; Mean of bootstrapped estimates
xBoot_TKEstd3_PD    = BootStrap_TKE3_PD[2]       ; Std dev of bootstrapped estimates

xBoot_LSKE1_PD       = BootStrap_LSKE1_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_LSKEavg1_PD    = BootStrap_LSKE1_PD[1]       ; Mean of bootstrapped estimates
xBoot_LSKEstd1_PD    = BootStrap_LSKE1_PD[2]       ; Std dev of bootstrapped estimates

xBoot_LSKE2_PD       = BootStrap_LSKE2_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_LSKEavg2_PD    = BootStrap_LSKE2_PD[1]       ; Mean of bootstrapped estimates
xBoot_LSKEstd2_PD    = BootStrap_LSKE2_PD[2]       ; Std dev of bootstrapped estimates

xBoot_LSKE3_PD       = BootStrap_LSKE3_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_LSKEavg3_PD    = BootStrap_LSKE3_PD[1]       ; Mean of bootstrapped estimates
xBoot_LSKEstd3_PD    = BootStrap_LSKE3_PD[2]       ; Std dev of bootstrapped estimates

xBoot_EKE1_PD       = BootStrap_EKE1_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_EKEavg1_PD    = BootStrap_EKE1_PD[1]       ; Mean of bootstrapped estimates
xBoot_EKEstd1_PD    = BootStrap_EKE1_PD[2]       ; Std dev of bootstrapped estimates

xBoot_EKE2_PD       = BootStrap_EKE2_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_EKEavg2_PD    = BootStrap_EKE2_PD[1]       ; Mean of bootstrapped estimates
xBoot_EKEstd2_PD    = BootStrap_EKE2_PD[2]       ; Std dev of bootstrapped estimates

xBoot_EKE3_PD       = BootStrap_EKE3_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_EKEavg3_PD    = BootStrap_EKE3_PD[1]       ; Mean of bootstrapped estimates
xBoot_EKEstd3_PD    = BootStrap_EKE3_PD[2]       ; Std dev of bootstrapped estimates
;-------------------------------------------------------------
xBootLow_TKE1_PD    = bootstrap_estimate(xBoot_TKE1_PD, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_TKE1_PD    = bootstrap_estimate(xBoot_TKE1_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_TKE1_PD    = bootstrap_estimate(xBoot_TKE1_PD, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_LSKE1_PD   = bootstrap_estimate(xBoot_LSKE1_PD, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_LSKE1_PD   = bootstrap_estimate(xBoot_LSKE1_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_LSKE1_PD   = bootstrap_estimate(xBoot_LSKE1_PD, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_EKE1_PD    = bootstrap_estimate(xBoot_EKE1_PD, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_EKE1_PD    = bootstrap_estimate(xBoot_EKE1_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_EKE1_PD    = bootstrap_estimate(xBoot_EKE1_PD, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------------
xBootLow_TKE2_PD    = bootstrap_estimate(xBoot_TKE2_PD, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_TKE2_PD    = bootstrap_estimate(xBoot_TKE2_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_TKE2_PD    = bootstrap_estimate(xBoot_TKE2_PD, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_LSKE2_PD   = bootstrap_estimate(xBoot_LSKE2_PD, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_LSKE2_PD   = bootstrap_estimate(xBoot_LSKE2_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_LSKE2_PD   = bootstrap_estimate(xBoot_LSKE2_PD, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_EKE2_PD    = bootstrap_estimate(xBoot_EKE2_PD, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_EKE2_PD    = bootstrap_estimate(xBoot_EKE2_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_EKE2_PD    = bootstrap_estimate(xBoot_EKE2_PD, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------------
xBootLow_TKE3_PD    = bootstrap_estimate(xBoot_TKE3_PD, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_TKE3_PD    = bootstrap_estimate(xBoot_TKE3_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_TKE3_PD    = bootstrap_estimate(xBoot_TKE3_PD, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_LSKE3_PD   = bootstrap_estimate(xBoot_LSKE3_PD, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_LSKE3_PD   = bootstrap_estimate(xBoot_LSKE3_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_LSKE3_PD   = bootstrap_estimate(xBoot_LSKE3_PD, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_EKE3_PD    = bootstrap_estimate(xBoot_EKE3_PD, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_EKE3_PD    = bootstrap_estimate(xBoot_EKE3_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_EKE3_PD    = bootstrap_estimate(xBoot_EKE3_PD, 0.975, False)   ; 97.5% upper confidence bound
;=====================================================================================================
stdarr1 = new((/3/),typeof(TKE1)) ;  TKE
stdarr1(0) =  xBoot_TKEavg1_PD
stdarr1(1) =  xBoot_TKEavg2_PD
stdarr1(2) =  xBoot_TKEavg3_PD

stdarr2 = new((/3/),typeof(LSKE1)) ; SSH
stdarr2(0) =  xBoot_LSKEavg1_PD
stdarr2(1) =  xBoot_LSKEavg2_PD
stdarr2(2) =  xBoot_LSKEavg3_PD

stdarr3 = new((/3/),typeof(EKE1)) ; EKE 
stdarr3(0) =  xBoot_EKEavg1_PD
stdarr3(1) =  xBoot_EKEavg2_PD
stdarr3(2) =  xBoot_EKEavg3_PD
;=============================================================
; plot mean TIW EKE and RMS
;=============================================================
outdir = "/Users/xueaoyun/Documents/work/TIW_future/revise_figures/"
wks = gsn_open_wks("pdf",outdir+"Different_EKE_calculation_method_comparison_CESM_E3SMv1_v2")
wks_type = "pdf"
; wks_type@wkWidth  = 0.8
; wks_type@wkHeight = 0.2  
;--------------------------------------
; gsn_define_colormap(wks,"posneg_1")
gsn_define_colormap(wks,"WhiteYellowOrangeRed")

res                      = True         
res@gsnDraw              = False
res@gsnFrame             = False
res@gsnAddCyclic         = False    
res@gsnMaximize          = False
; map set 
res@mpFillOn             = False
res@mpLandFillColor      = "mediumpurple1"
res@mpGeophysicalLineColor = "gray30"
res@mpGeophysicalLineThicknessF = 1.2
res@mpMinLatF            = -5
res@mpMaxLatF            =  10
res@mpMinLonF            = 180
res@mpMaxLonF            = 280
res@mpCenterLonF         = 230
res@mpShapeMode          = "FreeAspect"
; tmset
res@tmXBOn   = True
res@tmXTOn   = False
res@tmYROn   = False
res@tmYLMode = "Explicit"
res@tmYLValues = (/-10,-5,0,5,10/)
res@tmYLLabels = (/"10~S~o~N~S","5~S~o~N~S","EQ","5~S~o~N~N","10~S~o~N~N"/)
res@tmYLMinorValues         = ispan(-10,10,5) 
res@tmXBMode = "Explicit"
res@tmXBValues = (/180,200,220,240,260,280/)
res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W",""/)  
res@tmXBLabelFontHeightF  = 0.008
res@tmYLLabelFontHeightF  = 0.008
res@tmXBLabelFontThicknessF  = 2
res@tmYLLabelFontThicknessF  = 2 
res@tmBorderThicknessF       = 2
res@tmXBMajorLengthF        = -0.006
res@tmYLMajorLengthF        = -0.006
res@tmXBMinorLengthF        = -0.003
res@tmYLMinorLengthF        = -0.003 
res@tmXBMajorThicknessF     = 1.5
res@tmYLMajorThicknessF     = 1.5
res@tmXBMinorThicknessF     = 1
res@tmYLMinorThicknessF     = 1
; res@tmXBLabelDeltaF         = -0.4
; res@tmYLLabelDeltaF         = -0.4

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

; contour set
res@cnFillOn             = True
res@cnLinesOn            = False
res@cnLineLabelsOn       = False
res@cnMonoLineColor      = True
res@cnSmoothingOn        = True
res@cnSmoothingTensionF  = 5
res@gsnContourNegLineDashPattern = 1
res@cnInfoLabelOn        = False
res@cnLevelSelectionMode = "ExplicitLevels" 
res@cnLevels        =fspan(100,200,51)

;labelbar set
res@lbOrientation        = "Vertical"
res@lbBoxEndCapStyle     = "TriangleBothEnds"
res@pmLabelBarWidthF     = 0.03
res@pmLabelBarHeightF    = 0.10
res@pmLabelBarOrthogonalPosF = 0.05
res@pmLabelBarParallelPosF   = 0.50
res@lbLabelFontHeightF       =0.01
res@lbLabelStride        = 5
res@lbBoxLinesOn=False
res@lbTitlePosition  = "Bottom"
res@lbTitleDirection = "Across" 
res@lbTitleFontHeightF= 0.010
res@lbTitleOffsetF    = 0.05
res@lbTitleString     =  "J/m~S~3~N~"
res@lbLabelBarOn      = True 
; gsnset
res@gsnLeftString        = "~F22~a    CESM_TKE"
res@gsnRightString        = ""
res@gsnLeftStringFontHeightF = 0.010
res@gsnLeftStringOrthogonalPosF = 0.1
res@gsnCenterString = "~F22~CESM-UHR"
res@gsnCenterStringFontHeightF = 0.015
res@gsnCenterStringOrthogonalPosF = 0.35
res@gsnSpreadColors      = True
res@gsnSpreadColorStart  = 10
res@gsnSpreadColorEnd    = 230

res@vpHeightF            = 0.10
res@vpWidthF             = 0.20
res@vpXF                 = 0.05
res@vpYF                 = 0.90

plot1= gsn_csm_contour_map_ce(wks,TKE1,res)
res@gsnCenterString = ""
res@gsnCenterStringFontHeightF = 0.015
res@vpHeightF            = 0.10
res@vpWidthF             = 0.20
res@vpXF                 = 0.05
res@vpYF                 = 0.75
res@gsnLeftStringFontHeightF = 0.010
res@lbLabelBarOn      = True
res@gsnLeftString ="~F22~b    CESM LSKE"
res@cnLevels        =fspan(100,200,51)
plot2= gsn_csm_contour_map_ce(wks,LSKE1,res)
delete(res@cnLevels)
res@lbLabelBarOn      = True
res@gsnCenterString = ""
res@gsnCenterStringFontHeightF = 0.015
res@vpHeightF            = 0.10
res@vpWidthF             = 0.20
res@vpXF                 = 0.05
res@vpYF                 = 0.60
res@gsnLeftStringFontHeightF = 0.010
res@gsnLeftString ="~F22~c    CESM EKE"
res@cnLevels        =fspan(20,50,31)
plot3 = gsn_csm_contour_map_ce(wks,EKE1,res)
delete(res@cnLevels)
; plot E3SM

res@tmYLValues = (/-10,-5,0,5,10/)
res@tmYLLabels = (/"","","","",""/)
res@gsnCenterString = "~F22~E3SMv1"
res@gsnCenterStringFontHeightF = 0.015
res@gsnCenterStringOrthogonalPosF = 0.35
res@vpHeightF            = 0.10
res@vpWidthF             = 0.20
res@vpXF                 = 0.30
res@vpYF                 = 0.90
res@gsnLeftString ="~F22~d    E3SMv1 TKE "
res@cnLevels        =fspan(100,200,51)
plot4= gsn_csm_contour_map_ce(wks,TKE2,res)
res@lbLabelBarOn      = True
res@gsnCenterString = ""
res@gsnCenterStringFontHeightF = 0.015
res@vpHeightF            = 0.10
res@vpWidthF             = 0.20
res@vpXF                 = 0.30
res@vpYF                 = 0.75
res@cnLevels        =fspan(100,200,51)
res@gsnLeftString ="~F22~e    E3SMv1 LSKE "
plot5= gsn_csm_contour_map_ce(wks,LSKE2,res)
delete(res@cnLevels)
res@lbLabelBarOn      = True  
res@gsnCenterString = ""
res@gsnCenterStringFontHeightF = 0.015
res@vpHeightF            = 0.10
res@vpWidthF             = 0.20
res@vpXF                 = 0.30
res@vpYF                 = 0.60
res@cnLevels        =fspan(60,100,41)
res@gsnLeftString ="~F22~f    E3SMv1 EKE "
plot6= gsn_csm_contour_map_ce(wks,EKE2,res)
; E3SMv2
delete(res@cnLevels)
res@gsnCenterString = "~F22~E3SMv2"
res@gsnCenterStringFontHeightF = 0.015
res@gsnCenterStringOrthogonalPosF = 0.35
res@vpHeightF            = 0.10
res@vpWidthF             = 0.20
res@vpXF                 = 0.55
res@vpYF                 = 0.90
res@gsnLeftString ="~F22~g    E3SMv2 TKE "
delete(res@cnLevels)
res@cnLevels        =fspan(100,200,51)
plot7= gsn_csm_contour_map_ce(wks,TKE3,res)
res@lbLabelBarOn      = True
res@gsnCenterString = ""
res@gsnCenterStringFontHeightF = 0.015
res@vpHeightF            = 0.10
res@vpWidthF             = 0.20
res@vpXF                 = 0.55
res@vpYF                 = 0.75
res@cnLevels        =fspan(100,200,51)
res@gsnLeftString ="~F22~h    E3SMv2 LSKE "
plot8= gsn_csm_contour_map_ce(wks,LSKE3,res)
delete(res@cnLevels)
res@lbLabelBarOn      = True  
res@gsnCenterString = ""
res@gsnCenterStringFontHeightF = 0.015
res@vpHeightF            = 0.10
res@vpWidthF             = 0.20
res@vpXF                 = 0.55
res@vpYF                 = 0.60
res@cnLevels        =fspan(40,80,41)
res@gsnLeftString ="~F22~i    E3SMv2 EKE "
plot9= gsn_csm_contour_map_ce(wks,EKE3,res)

;=========================================================================
plres=True
plres@gsLineThicknessF=2.5
; plres@gsFillColor     = -1
plres@gsLineDashPattern= 2
plres@gsLineColor="blue"
plres@gsLineOpacityF=100
w1=360-160
e1=360-100

s1=0
n1=5

dum1=gsn_add_polyline(wks, plot1, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum2=gsn_add_polyline(wks, plot2, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum3=gsn_add_polyline(wks, plot3, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum4=gsn_add_polyline(wks, plot4, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum5=gsn_add_polyline(wks, plot5, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum6=gsn_add_polyline(wks, plot6, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum7=gsn_add_polyline(wks, plot7, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum8=gsn_add_polyline(wks, plot8, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum9=gsn_add_polyline(wks, plot9, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)


draw(plot1)
draw(plot2)
draw(plot3)
draw(plot4)
draw(plot5)
draw(plot6)
draw(plot7)
draw(plot8)
draw(plot9)
;================================================================================
; plot bar
 xyres                       = True                 ; plot mods desired
 xyres@gsnFrame              = False                ; don't advance frame yet
 xyres@gsnDraw               = False
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 200.                ; bring bars down to zero
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 4.                   ; of the 1st and last bars
 xyres@vpWidthF              = 0.10
 xyres@vpHeightF             = 0.10
 xyres@vpXF                  = 0.85
 xyres@vpYF                  = 0.90

 xyres@gsnXYBarChart         = True                 ; turn on bar chart
 xyres@gsnXYBarChartBarWidth = 0.24                  ; change bar widths
 xyres@tmXBOn                = True                 ; turn off tickmarks at bot
 xyres@tmYLOn                = True
 xyres@tmYROn                = False                 ; turn off tickmarks at bot
 xyres@tmXTOn                = False                 ; turn off tickmarks at bot
 xyres@tmXBMode          = "Explicit"         ; explicit labels
 xyres@tmXBValues        = (/1,2,3,4/)
 xyres@tmXBLabels = (/"CESM","E3SMv1","E3SMv2",""/)
 xyres@tmXBLabelFontHeightF = 0.006
 xyres@tmXBLabelAngleF   = 45
 xyres@tmYLMode          = "Explicit"         ; explicit labels
 xyres@tmYLValues        = (/0,50,100,150,200/)
 xyres@tmYLLabels        = (/"0","50","100","150","200"/)
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
 xyres@tmXTBorderOn  = False ; turn off top axis
 xyres@tmYRBorderOn  = False ; turn off right axis

 xyres@gsnXYBarChartColors = (/"cyan3"/)
 xyres@gsnCenterString       = "TKE"
 xyres@gsnCenterStringOrthogonalPosF = -0.10
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString  =  "~F22~j"
 xyres@gsnLeftStringFontHeightF = 0.010
 xyres@gsnLeftStringOrthogonalPosF = 0.
 xyres@gsnLeftStringParallelPosF = -0.05
 plot10  = gsn_csm_xy(wks,ispan(1,3,1),stdarr1(:),xyres)   ; advancing the frame
 xyres@vpWidthF              = 0.10
 xyres@vpHeightF             = 0.10
 xyres@vpXF                  = 0.85
 xyres@vpYF                  = 0.75
 xyres@gsnXYBarChartColors = (/"pink"/)
 xyres@gsnCenterString       = "LSKE"
 xyres@gsnCenterStringParallelPosF = 0.45
 xyres@gsnLeftString  = "~F22~k"
 xyres@gsnLeftStringFontHeightF = 0.010
 xyres@gsnLeftStringOrthogonalPosF = 0.
 xyres@gsnLeftStringParallelPosF = -0.05
 plot11  = gsn_csm_xy(wks,ispan(1,3,1),stdarr2(:),xyres)   ; advancing the frame

 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)

 xyres@vpWidthF              = 0.10
 xyres@vpHeightF             = 0.10
 xyres@vpXF                  = 0.85
 xyres@vpYF                  = 0.60
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 100                 ; bring bars down to zero
                       ; turn off right axis
 xyres@gsnXYBarChartColors = (/"chartreuse2"/)
 xyres@tmYLValues        = (/0,20,40,60,80,100/)
 xyres@tmYLLabels        = (/"0","20","40","60","80","100"/)
 xyres@gsnCenterString       = "EKE"
 xyres@gsnCenterStringOrthogonalPosF = -0.12
 xyres@gsnCenterStringParallelPosF = 0.45
 xyres@gsnLeftString  = "~F22~l"
 xyres@gsnLeftStringFontHeightF = 0.010
 xyres@gsnLeftStringOrthogonalPosF = 0.
 plot12  = gsn_csm_xy(wks,ispan(1,3,1),stdarr3(:),xyres)   ; advancing the frame
; plot stddev
  polyres=True
  polyres@gsLineDashPattern=0
  polyres@gsLineColor="black"
  polyres@gsLineThicknessF=1.0
  polyres@tfPolyDrawOrder = "postDraw"

  aumm1    = gsn_add_polyline(wks,plot10,(/1,1/),(/xBootLow_TKE1_PD,xBootHi_TKE1_PD/), polyres)
  aumm01   = gsn_add_polyline(wks,plot10,(/2,2/),(/xBootLow_TKE2_PD,xBootHi_TKE2_PD/), polyres)
  aumm001  = gsn_add_polyline(wks,plot10,(/3,3/),(/xBootLow_TKE3_PD,xBootHi_TKE3_PD/), polyres)

  bumm1    = gsn_add_polyline(wks,plot11,(/1,1/),(/xBootLow_LSKE1_PD,xBootHi_LSKE1_PD/), polyres)
  bumm01   = gsn_add_polyline(wks,plot11,(/2,2/),(/xBootLow_LSKE2_PD,xBootHi_LSKE2_PD/), polyres)
  bumm001  = gsn_add_polyline(wks,plot11,(/3,3/),(/xBootLow_LSKE3_PD,xBootHi_LSKE3_PD/), polyres)

  cumm1    = gsn_add_polyline(wks,plot12,(/1,1/),(/xBootLow_EKE1_PD,xBootHi_EKE1_PD/), polyres)
  cumm01   = gsn_add_polyline(wks,plot12,(/2,2/),(/xBootLow_EKE2_PD,xBootHi_EKE2_PD/), polyres)
  cumm001  = gsn_add_polyline(wks,plot12,(/3,3/),(/xBootLow_EKE3_PD,xBootHi_EKE3_PD/), polyres)

  wid = 0.24/2
  dumm1   = gsn_add_polyline(wks,plot10,(/1-wid,1+wid/),(/xBootLow_TKE1_PD,xBootLow_TKE1_PD/), polyres)
  dumm01  = gsn_add_polyline(wks,plot10,(/1-wid,1+wid/),(/xBootHi_TKE1_PD,xBootHi_TKE1_PD/), polyres)

  eumm1   = gsn_add_polyline(wks,plot10,(/2-wid,2+wid/),(/xBootLow_TKE2_PD,xBootLow_TKE2_PD/), polyres)
  eumm01  = gsn_add_polyline(wks,plot10,(/2-wid,2+wid/),(/xBootHi_TKE2_PD,xBootHi_TKE2_PD/), polyres)

  fumm1   = gsn_add_polyline(wks,plot10,(/3-wid,3+wid/),(/xBootLow_TKE3_PD,xBootLow_TKE3_PD/), polyres)
  fumm01  = gsn_add_polyline(wks,plot10,(/3-wid,3+wid/),(/xBootHi_TKE3_PD,xBootHi_TKE3_PD/), polyres)
;----------------------------------------------------------------------------------------------------
  gumm1   = gsn_add_polyline(wks,plot11,(/1-wid,1+wid/),(/xBootLow_LSKE1_PD,xBootLow_LSKE1_PD/), polyres)
  gumm01  = gsn_add_polyline(wks,plot11,(/1-wid,1+wid/),(/xBootHi_LSKE1_PD,xBootHi_LSKE1_PD/), polyres)

  humm1   = gsn_add_polyline(wks,plot11,(/2-wid,2+wid/),(/xBootLow_LSKE2_PD,xBootLow_LSKE2_PD/), polyres)
  humm01  = gsn_add_polyline(wks,plot11,(/2-wid,2+wid/),(/xBootHi_LSKE2_PD,xBootHi_LSKE2_PD/), polyres)

  iumm1   = gsn_add_polyline(wks,plot11,(/3-wid,3+wid/),(/xBootLow_LSKE3_PD,xBootLow_LSKE3_PD/), polyres)
  iumm01  = gsn_add_polyline(wks,plot11,(/3-wid,3+wid/),(/xBootHi_LSKE3_PD,xBootHi_LSKE3_PD/), polyres)
  ;----------------------------------------------------------------------------------------------------
  lumm1   = gsn_add_polyline(wks,plot12,(/1-wid,1+wid/),(/xBootLow_EKE1_PD,xBootLow_EKE1_PD/), polyres)
  lumm01  = gsn_add_polyline(wks,plot12,(/1-wid,1+wid/),(/xBootHi_EKE1_PD,xBootHi_EKE1_PD/), polyres)

  mumm1   = gsn_add_polyline(wks,plot12,(/2-wid,2+wid/),(/xBootLow_EKE2_PD,xBootLow_EKE2_PD/), polyres)
  mumm01  = gsn_add_polyline(wks,plot12,(/2-wid,2+wid/),(/xBootHi_EKE2_PD,xBootHi_EKE2_PD/), polyres)

  numm1   = gsn_add_polyline(wks,plot12,(/3-wid,3+wid/),(/xBootLow_EKE3_PD,xBootLow_EKE3_PD/), polyres)
  numm01  = gsn_add_polyline(wks,plot12,(/3-wid,3+wid/),(/xBootHi_EKE3_PD,xBootHi_EKE3_PD/), polyres)
;------------------------------------------------------------------------------------------------------

  draw(plot10)
  draw(plot11)
  draw(plot12)

end 