;==========================functions=============================
;--------------------------------------
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
;-----------------------------------------------------------
; read E3SM 
ncdir1  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/E3SM/E3SMv2/"
;============================================================================
fil1 = systemfunc("ls " + ncdir1+"tos_day_E3SM-2-0_0001-0500_remap_equator_10-50fft_new_1mon_std_mean.nc")
f1   = addfile(fil1,"r")
tos1 = f1->sst(0,:,:)
printVarSummary(tos1)
;---------------------------------------------------------
fil2 = systemfunc("ls " + ncdir1 + "sst_abrupt_4xCO2_0001-0160_remap_equator_moving_20yr_detrend_10-50fft_1mon_runstd_mean.nc")
f2   = addfile(fil2,"r")
tos2 = f2->sst(0,:,:)
printVarSummary(tos2)
;===========================================================================
fil01 = systemfunc("ls " + ncdir1+"zos_day_E3SM-2-0_0001-0500_remap_equator_10-50fft_new_1mon_std_mean.nc")
f01   = addfile(fil01,"r")
zos01 = f01->ssh(0,:,:)
zos1  = zos01 *100.
copy_VarCoords(zos01, zos1)
printVarSummary(zos1)
;---------------------------------------------------------
fil02 = systemfunc("ls " + ncdir1 +"ssh_abrupt_4xCO2_0001-0160_remap_equator_moving_20yr_detrend_10-50fft_1mon_runstd_mean.nc")
f02   = addfile(fil02,"r")
zos02 = f02->ssh(0,:,:)
zos2  = zos02 *100.
copy_VarCoords(zos02, zos2)
printVarSummary(zos2)
;===========================================================================
fil001 = systemfunc("ls " + ncdir1+"EKE_picontrol_0001-0500_remap_equator_runave_mean.nc")
f001   = addfile(fil001,"r")
KE1 = f001->kineticEnergyAtSurface(0,:,:)
printVarSummary(KE1)
;---------------------------------------------------------
fil002 = systemfunc("ls " + ncdir1 +"EKE_abrupt_4xCO2_0001-0160_remap_equator_runave_mean.nc")
f002   = addfile(fil002,"r")
KE2 = f002->KE(0,:,:)
printVarSummary(KE2)
;=========================================================
tos_pattern1 = smth9_Wrap(tos1, 0.50, 0.25, False)
tos_pattern2 = smth9_Wrap(tos2, 0.50, 0.25, False)
zos_pattern1 = smth9_Wrap(zos1, 0.50, 0.25, False)
zos_pattern2 = smth9_Wrap(zos2, 0.50, 0.25, False)
KE_pattern1  = smth9_Wrap(KE1, 0.50, 0.25, False)
KE_pattern2  = smth9_Wrap(KE2, 0.50, 0.25, False)

pattern_tos = new((/2, dimsizes(tos_pattern1&lat), dimsizes(tos_pattern1&lon)/), typeof(tos_pattern1))
pattern_tos(0,:,:) = tos_pattern1 
pattern_tos(1,:,:) = tos_pattern2

pattern_zos = new((/2, dimsizes(zos_pattern1&lat), dimsizes(zos_pattern1&lon)/), typeof(zos_pattern1))
pattern_zos(0,:,:) = zos_pattern1 
pattern_zos(1,:,:) = zos_pattern2

pattern_KE = new((/2, dimsizes(KE_pattern1&lat), dimsizes(KE_pattern1&lon)/), typeof(KE_pattern1))
pattern_KE(0,:,:) = KE_pattern1 
pattern_KE(1,:,:) = KE_pattern2
;===========================================================================
outdir = "/Users/xueaoyun/Documents/work/TIW_future/revise_figures/"
wks = gsn_open_wks("pdf",outdir+"plot_E3SMv2_TIW_mean_intensity")
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

res@tmXBOn   = True
res@tmXTOn   = False
res@tmYROn   = False
res@tmYLMode = "Explicit"
; res@tmYLValues = (/-10,-8,-6,-4,-2,0,2,4,6,8,10/)
; res@tmYLLabels = (/"10~S~o~N~S","8~S~o~N~S","6~S~o~N~S","4~S~o~N~S","2~S~o~N~S","EQ","2~S~o~N~N","4~S~o~N~N","6~S~o~N~N","8~S~o~N~N","10~S~o~N~N"/)
; res@tmYLMinorValues         = ispan(-10,10,2) 
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
res@gsnRightString        = ""
res@gsnRightStringFontHeightF  = 0.008

res@cnFillOn             = True
res@cnLinesOn            = False
res@cnLineLabelsOn       = False
res@cnMonoLineColor      = True
res@cnSmoothingOn        = True
res@cnSmoothingTensionF  = 5
res@gsnContourNegLineDashPattern = 1
; res@gsnContourZeroLineThicknessF = 2
res@cnInfoLabelOn        = False
res@cnLevelSelectionMode = "ExplicitLevels"  
res@lbOrientation        = "Vertical"
res@lbBoxEndCapStyle     = "TriangleBothEnds"
res@pmLabelBarWidthF     = 0.03
res@pmLabelBarHeightF    = 0.12
res@pmLabelBarOrthogonalPosF = 0.05
; res@pmLabelBarParallelPosF   = 0.6
res@lbLabelFontHeightF       =0.008
res@gsnLeftString        = ""
res@gsnCenterStringFontHeightF  = 0.014

res@gsnSpreadColorEnd  =230

res@vpHeightF            = 0.12
res@vpWidthF             = 0.20

res@lbTitlePosition  = "Bottom"
res@lbTitleDirection = "Across" 
res@lbTitleFontHeightF= 0.008
res@lbTitleOffsetF    = 0.00
res@lbLabelStride     = 10

plot=new(9,graphic)
dum=plot

; ----->1 plot eke 
panels = (/"~F22~a", "~F22~b", "~F22~d", "~F22~e", "~F22~g", "~F22~h"/)
text=(/"TIW-EKE","TIW-SST","TIW-SSH"/)
yspace=0.03
yspace1 = 0.03
yf=0.9
xf=0.08
res@tiDeltaF=3
res@gsnCenterStringOrthogonalPosF=0.10
res@cnLevels                = fspan(60,80,41)
res@vpYF = yf
res@lbTitleString    = "J/m~S~3~N~"
do i=0,1
res@vpXF  = xf +i*yspace+i*res@vpWidthF

res@tiYAxisOffsetXF= 0.00
res@tiYAxisFontHeightF=0.014
if i.eq.0 then
res@gsnCenterString = "PD"
res@tiYAxisString=text(i)
res@tiMainString=""
else
res@gsnCenterString = "4xCO2"
res@tiYAxisString= ""
res@tiMainString=""
end if

if i.eq.1 then
res@lbLabelBarOn    = True
res@tmXBLabels = (/"","","","","",""/)  
else
res@lbLabelBarOn    = False
res@tmXBLabels = (/"","","","","",""/)
end if
res@gsnLeftString  = panels(i)
res@gsnLeftStringFontHeightF = 0.015
res@gsnLeftStringOrthogonalPosF = 0.01
res@gsnLeftStringParallelPosF  = 0.03
plot(i)  = gsn_csm_contour_map(wks,pattern_KE(i,:,:),res)
end do


; ----->2 plot sst
res@lbTitleString    = "degC"
delete(res@cnLevels)
res@vpYF    = res@vpYF-res@vpHeightF-yspace
res@cnLevels             = fspan(0.2,0.6,41)
res@cnFillPalette         = "WhiteBlueGreenYellowOrangeRed"
; res@cnFillColors         = (/0,12,13,14,15,17,19/)
res@tmYLLabels = (/"","","","",""/)
do i=0,1
res@vpXF    = xf+i*yspace+i*res@vpWidthF
if i.eq.0 then
res@gsnCenterString = ""
res@tiYAxisString=text(i+1)
res@tmYLLabels = (/"10~S~o~N~S","5~S~o~N~S","EQ","5~S~o~N~N","10~S~o~N~N"/)
else
res@gsnCenterString = ""
res@tiYAxisString= ""
end if 
if i.eq.1 then
res@lbLabelBarOn    = True
res@tmXBLabels = (/"","","","","",""/)  
else
res@lbLabelBarOn    = False
end if
res@gsnLeftString  = panels(i+2)
res@gsnLeftStringFontHeightF = 0.015
res@gsnLeftStringOrthogonalPosF = 0.01
res@gsnLeftStringParallelPosF  = 0.03
plot(i+2)  = gsn_csm_contour_map(wks,pattern_tos(i,:,:),res)
end do
; ----->3 plot ssh 
delete(res@cnLevels)
res@lbLabelStride     = 5
res@vpYF    = res@vpYF - res@vpHeightF -yspace
res@lbTitleString    = "10~S~-2~N~ m"
res@cnLevels              = fspan(2,3,21)
res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W",""/)  
do i=0,1
res@vpXF    = xf+i*yspace+i*res@vpWidthF
if i.eq.0 then
res@gsnCenterString = ""
res@tiYAxisString=text(i+2)
res@tmYLLabels =(/"10~S~o~N~S","5~S~o~N~S","EQ","5~S~o~N~N","10~S~o~N~N"/)
else
res@gsnCenterString = ""
res@tiYAxisString= ""
end if
if i.eq.1 then
res@lbLabelBarOn    = True
else
res@lbLabelBarOn    = False

end if
res@gsnLeftString  = panels(i+4)
res@gsnLeftStringFontHeightF = 0.015
res@gsnLeftStringOrthogonalPosF = 0.01
plot(i+4)  = gsn_csm_contour_map(wks,pattern_zos(i,:,:),res)
end do
;--------------------------------------------------------------------------------------------------------------

;=========================================================================
plres=True
plres@gsLineThicknessF=2.5
; plres@gsFillColor     = -1
plres@gsLineDashPattern= 2
plres@gsLineColor="blue"
plres@gsLineOpacityF=100
w1=360-160
e1=360-100

w2=360-160
e2=360-110

w3=360-160
e3=360-110

s1=0
n1=5

s2 = 3
n2 = 6

dum(0)=gsn_add_polyline(wks, plot(0), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum(1)=gsn_add_polyline(wks, plot(1), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum(2)=gsn_add_polyline(wks, plot(2), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum(3)=gsn_add_polyline(wks, plot(3), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum(4)=gsn_add_polyline(wks, plot(4), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
dum(5)=gsn_add_polyline(wks, plot(5), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)


draw(plot(0))
draw(plot(1))
draw(plot(2))
draw(plot(3))
draw(plot(4))
draw(plot(5))
delete([/tos1, tos2/])
delete([/zos1, zos2/])
delete([/KE1,  KE2/])

;====================================================================================
fil1 = systemfunc("ls " + ncdir1+"tos_mon_E3SM-2-0_0001-0500_remap_equator_10-50fft_new_1mon_std_1degree.nc")
f1   = addfile(fil1,"r")
sst1 = f1->sst(:,:,:)
printVarSummary(sst1)
;---------------------------------------------------------
fil2 = systemfunc("ls " + ncdir1 + "sst_mon_abrupt_4xCO2_0001-0160_remap_equator_moving_20yr_detrend_10-50fft_1mon_runstd_1degree.nc")
f2   = addfile(fil2,"r")
sst2 = f2->sst(:,:,:)
printVarSummary(sst2)
;===========================================================================
fil01 = systemfunc("ls " + ncdir1+"zos_mon_E3SM-2-0_0001-0500_remap_equator_10-50fft_new_1mon_std_1degree.nc")
f01   = addfile(fil01,"r")
ssh01 = f01->ssh(:,:,:)
ssh1  = ssh01 *100.
copy_VarCoords(ssh01, ssh1)
printVarSummary(ssh1)
;---------------------------------------------------------
fil02 = systemfunc("ls " + ncdir1 +"ssh_mon_abrupt_4xCO2_0001-0160_remap_equator_moving_20yr_detrend_10-50fft_1mon_runstd_1degree.nc")
f02   = addfile(fil02,"r")
ssh02 = f02->ssh(:,:,:)
ssh2  = ssh02 *100.
copy_VarCoords(ssh02, ssh2)
printVarSummary(ssh2)
;===========================================================================
fil001 = systemfunc("ls " + ncdir1+"EKE_mon_picontrol_0001-0500_remap_equator_runave.nc")
f001   = addfile(fil001,"r")
KE1 = f001->KE(:,:,:)
printVarSummary(KE1)
;---------------------------------------------------------
fil002 = systemfunc("ls " + ncdir1 +"EKE_mon_abrupt_4xCO2_0001-0160_remap_equator_runave.nc")
f002   = addfile(fil002,"r")
KE2 = f002->KE(:,:,:)
printVarSummary(KE2)

sst1_PD    = wgt_areaave_Wrap(sst1(:,{0:5},{200:260}), 1.0, 1.0, 0)
sst1_4xCO2 = wgt_areaave_Wrap(sst2(:,{0:5},{200:260}), 1.0, 1.0, 0)

ssh1_PD    = wgt_areaave_Wrap(ssh1(:,{3:6},{210:250}), 1.0, 1.0, 0)
ssh1_4xCO2 = wgt_areaave_Wrap(ssh2(:,{3:6},{200:240}), 1.0, 1.0, 0)

eke1_PD    = wgt_areaave_Wrap(KE1(:,{0:5},{210:260}), 1.0, 1.0, 0)
eke1_4xCO2 = wgt_areaave_Wrap(KE2(:,{0:5},{200:250}), 1.0, 1.0, 0)

stat  = 0                         ; 0=mean; 1=variance; 2=std. dev; ...
n     = 400                       ; sub-sample size;  default is n=N
nBoot = 10000                      ; bootstrap replications
N     = 10000
opt   = False                     ; use defaults
if (isvar("n") .and. (n.lt.N)) then
    opt   = True
    opt@sample_size = n           ; default N
end if

 BootStrap_sst1_PD     = bootstrap_stat(sst1_PD, stat, nBoot, 0, opt)
 BootStrap_sst1_4xCO2  = bootstrap_stat(sst1_4xCO2, stat, nBoot, 0, opt)

 BootStrap_ssh1_PD     = bootstrap_stat(ssh1_PD, stat, nBoot, 0, opt)
 BootStrap_ssh1_4xCO2  = bootstrap_stat(ssh1_4xCO2, stat, nBoot, 0, opt)

 BootStrap_eke1_PD     = bootstrap_stat(eke1_PD, stat, nBoot, 0, opt)
 BootStrap_eke1_4xCO2  = bootstrap_stat(eke1_4xCO2, stat, nBoot, 0, opt)

 xBoot_sst1_PD       = BootStrap_sst1_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sstavg1_PD    = BootStrap_sst1_PD[1]       ; Mean of bootstrapped estimates
 xBoot_sststd1_PD    = BootStrap_sst1_PD[2]       ; Std dev of bootstrapped estimates

 xBoot_ssh1_PD       = BootStrap_ssh1_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sshavg1_PD    = BootStrap_ssh1_PD[1]       ; Mean of bootstrapped estimates
 xBoot_sshstd1_PD    = BootStrap_ssh1_PD[2]       ; Std dev of bootstrapped estimates

 xBoot_eke1_PD       = BootStrap_eke1_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_ekeavg1_PD    = BootStrap_eke1_PD[1]       ; Mean of bootstrapped estimates
 xBoot_ekestd1_PD    = BootStrap_eke1_PD[2]       ; Std dev of bootstrapped estimates
;=====================================================
    xBoot_sst1_4xCO2    = BootStrap_sst1_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sstavg1_4xCO2    = BootStrap_sst1_4xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_sststd1_4xCO2    = BootStrap_sst1_4xCO2[2]       ; Std dev of bootstrapped estimates

    xBoot_ssh1_4xCO2    = BootStrap_ssh1_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sshavg1_4xCO2    = BootStrap_ssh1_4xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_sshstd1_4xCO2    = BootStrap_ssh1_4xCO2[2]       ; Std dev of bootstrapped estimates

    xBoot_eke1_4xCO2    = BootStrap_eke1_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_ekeavg1_4xCO2    = BootStrap_eke1_4xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_ekestd1_4xCO2    = BootStrap_eke1_4xCO2[2]       ; Std dev of bootstrapped estimates
;=====================================================
 xBootLow_sst1_PD    = bootstrap_estimate(xBoot_sst1_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_sst1_PD    = bootstrap_estimate(xBoot_sst1_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_sst1_PD    = bootstrap_estimate(xBoot_sst1_PD, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_ssh1_PD    = bootstrap_estimate(xBoot_ssh1_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_ssh1_PD    = bootstrap_estimate(xBoot_ssh1_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_ssh1_PD    = bootstrap_estimate(xBoot_ssh1_PD, 0.975, False)   ; 97.5% upper confidence bound
 
 xBootLow_eke1_PD    = bootstrap_estimate(xBoot_eke1_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_eke1_PD    = bootstrap_estimate(xBoot_eke1_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_eke1_PD    = bootstrap_estimate(xBoot_eke1_PD, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------
xBootLow_sst1_4xCO2    = bootstrap_estimate(xBoot_sst1_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_sst1_4xCO2    = bootstrap_estimate(xBoot_sst1_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_sst1_4xCO2    = bootstrap_estimate(xBoot_sst1_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_ssh1_4xCO2    = bootstrap_estimate(xBoot_ssh1_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_ssh1_4xCO2    = bootstrap_estimate(xBoot_ssh1_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_ssh1_4xCO2    = bootstrap_estimate(xBoot_ssh1_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_eke1_4xCO2    = bootstrap_estimate(xBoot_eke1_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_eke1_4xCO2    = bootstrap_estimate(xBoot_eke1_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_eke1_4xCO2    = bootstrap_estimate(xBoot_eke1_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

 stdarr01 = new((/2/),typeof(sst1_PD)) ;  SST
stdarr01(0) =  xBoot_sstavg1_PD
stdarr01(1) =  xBoot_sstavg1_4xCO2
printMinMax(stdarr01, 0)
stdarr02 = new((/2/),typeof(ssh1_PD)) ; SSH
stdarr02(0) =  xBoot_sshavg1_PD
stdarr02(1) =  xBoot_sshavg1_4xCO2

stdarr03 = new((/2/),typeof(eke1_PD)) ; EKE 
stdarr03(0) =  xBoot_ekeavg1_PD
stdarr03(1) =  xBoot_ekeavg1_4xCO2 

 xyres                       = True                 ; plot mods desired
 xyres@gsnFrame              = False                ; don't advance frame yet
 xyres@gsnDraw               = False
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 100.                 ; bring bars down to zero
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 3.                   ; of the 1st and last bars
 xyres@vpWidthF              = 0.10
 xyres@vpHeightF             = 0.12

 xyres@gsnXYBarChart         = True                 ; turn on bar chart
 xyres@gsnXYBarChartBarWidth = 0.24                  ; change bar widths
 xyres@tmXBOn                = True                 ; turn off tickmarks at bot
 xyres@tmYLOn                = True
 xyres@tmYROn                = False                 ; turn off tickmarks at bot
 xyres@tmXTOn                = False                 ; turn off tickmarks at bot
 xyres@tmXBMode          = "Explicit"         ; explicit labels
 xyres@tmXBValues        = (/1,2,3/)
 xyres@tmXBLabels = (/"PD","4xCO2",""/)
 xyres@tmXBLabelFontHeightF = 0.006
 xyres@tmYLMode          = "Explicit"         ; explicit labels
 xyres@tmYLValues        = (/0,20,40,60,80,100/)*1.0
 xyres@tmYLLabels        = (/"0","20","40","60","80","100"/)
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
 xyres@tmXBLabelsOn  = True 
 xyres@vpXF    = xf+2*yspace+2*res@vpWidthF + 2*yspace
 xyres@vpYF    = yf
 xyres@gsnXYBarChartColors = (/"orange3"/)
 xyres@gsnCenterString       = "E3SMv2_EKE"
 xyres@gsnCenterStringOrthogonalPosF = -0.125
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString  = "~F22~c"
 xyres@gsnLeftStringFontHeightF = 0.015
 xyres@gsnLeftStringOrthogonalPosF = 0.
 xyres@gsnLeftStringParallelPosF = -0.05
 plot(6)  = gsn_csm_xy(wks,ispan(1,2,1),stdarr03(:),xyres)   ; advancing the frame
 xyres@vpYF    = xyres@vpYF + xyres@vpWidthF + yspace
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 3.                   ; of the 1st and last bars
 xyres@tmXBValues        = (/1,2,3/)
 xyres@tmXBLabels = (/"PD","4xCO2","" /)
 xyres@tmYLOn        = True
 xyres@tmYRBorderOn  = False                        ; turn off right axis
 xyres@tmYLBorderOn  = True 
 xyres@vpWidthF              = 0.10
 xyres@vpHeightF             = 0.12
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.5                 ; bring bars down to zero                   
 xyres@tmYLOn        = True 
 xyres@tmYLMode          = "Explicit"         ; explicit labels
 xyres@tmYLValues        = (/0,0.1,0.2,0.3,0.4,0.5/)
 xyres@tmYLLabels        = (/"0","0.1","0.2","0.3","0.4","0.5"/)
 xyres@gsnXYBarChartColors = (/"cyan3"/)
 xyres@gsnCenterString       = "E3SMv2_SST"
 xyres@gsnCenterStringOrthogonalPosF = -0.125
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString  = "~F22~f"
 xyres@gsnLeftStringFontHeightF = 0.015
 xyres@gsnLeftStringOrthogonalPosF = 0.
 xyres@gsnLeftStringParallelPosF = -0.05
 xyres@vpXF    = xf+2*yspace+2*res@vpWidthF + 2*yspace
 xyres@vpYF    = yf-res@vpHeightF-yspace
 plot(7)  = gsn_csm_xy(wks,ispan(1,2,1),stdarr01(:),xyres)   ; advancing the frame

 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 3.                   ; of the 1st and last bars
 xyres@tmXBValues        = (/1,2,3/)
 xyres@tmXBLabels = (/"PD","4xCO2","" /)
 xyres@vpXF    = xf+2*yspace+2*res@vpWidthF + 2*yspace
 xyres@vpYF    = yf-2*res@vpHeightF-2*yspace

 xyres@tmYLOn        = True
 xyres@tmYRBorderOn  = False                        ; turn off right axis
 xyres@tmYLBorderOn  = True 
 xyres@vpWidthF              = 0.10
 xyres@vpHeightF             = 0.12
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 4.                 ; bring bars down to zero                   
 xyres@tmYLOn        = True 
 xyres@tmYLMode          = "Explicit"         ; explicit labels
 xyres@tmYLValues        = (/0,1,2,3,4,5/)*1.0
 xyres@tmYLLabels        = (/"0","1","2","3","4","5"/)
 xyres@gsnXYBarChartColors = (/"pink"/)
 xyres@gsnCenterString       = "E3SMv2_SSH"
 xyres@gsnCenterStringOrthogonalPosF = -0.125
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString  = "~F22~i"
 xyres@gsnLeftStringFontHeightF = 0.015
 xyres@gsnLeftStringOrthogonalPosF = 0.
 xyres@gsnLeftStringParallelPosF = -0.05
 plot(8)  = gsn_csm_xy(wks,ispan(1,2,1),stdarr02(:),xyres)   ; advancing the frame
 

  polyres=True
  polyres@gsLineDashPattern=0
  polyres@gsLineColor="black"
  polyres@gsLineThicknessF=1.0
  polyres@tfPolyDrawOrder = "postDraw"

  aumm1    = gsn_add_polyline(wks,plot(6),(/1,1/),(/xBootLow_eke1_PD,xBootHi_eke1_PD/), polyres)
  aumm001  = gsn_add_polyline(wks,plot(6),(/2,2/),(/xBootLow_eke1_4xCO2,xBootHi_eke1_4xCO2/), polyres)

  bumm1    = gsn_add_polyline(wks,plot(7),(/1,1/),(/xBootLow_sst1_PD,xBootHi_sst1_PD/), polyres)
  bumm001  = gsn_add_polyline(wks,plot(7),(/2,2/),(/xBootLow_sst1_4xCO2,xBootHi_sst1_4xCO2/), polyres)

  cumm1    = gsn_add_polyline(wks,plot(8),(/1,1/),(/xBootLow_ssh1_PD,xBootHi_ssh1_PD/), polyres)
  cumm001  = gsn_add_polyline(wks,plot(8),(/2,2/),(/xBootLow_ssh1_4xCO2,xBootHi_ssh1_4xCO2/), polyres)

  wid = 0.24/2
;----------------------------------------------------------------------------------------------------
  dumm1   = gsn_add_polyline(wks,plot(6),(/1-wid,1+wid/),(/xBootLow_eke1_PD,xBootLow_eke1_PD/), polyres)
  dumm01  = gsn_add_polyline(wks,plot(6),(/1-wid,1+wid/),(/xBootHi_eke1_PD,xBootHi_eke1_PD/), polyres)
  
  eumm1   = gsn_add_polyline(wks,plot(6),(/2-wid,2+wid/),(/xBootLow_eke1_4xCO2,xBootLow_eke1_4xCO2/), polyres)
  eumm01  = gsn_add_polyline(wks,plot(6),(/2-wid,2+wid/),(/xBootHi_eke1_4xCO2,xBootHi_eke1_4xCO2/), polyres)
;------------------------------------------------------------------------------------------------------
  fumm1   = gsn_add_polyline(wks,plot(7),(/1-wid,1+wid/),(/xBootLow_sst1_PD,xBootLow_sst1_PD/), polyres)
  fumm01  = gsn_add_polyline(wks,plot(7),(/1-wid,1+wid/),(/xBootHi_sst1_PD,xBootHi_sst1_PD/), polyres)
  
  gumm1   = gsn_add_polyline(wks,plot(7),(/2-wid,2+wid/),(/xBootLow_sst1_4xCO2,xBootLow_sst1_4xCO2/), polyres)
  gumm01  = gsn_add_polyline(wks,plot(7),(/2-wid,2+wid/),(/xBootHi_sst1_4xCO2,xBootHi_sst1_4xCO2/), polyres)
;----------------------------------------------------------------------------------------------------
  humm1   = gsn_add_polyline(wks,plot(8),(/1-wid,1+wid/),(/xBootLow_ssh1_PD,xBootLow_ssh1_PD/), polyres)
  humm01  = gsn_add_polyline(wks,plot(8),(/1-wid,1+wid/),(/xBootHi_ssh1_PD,xBootHi_ssh1_PD/), polyres)
  
  iumm1   = gsn_add_polyline(wks,plot(8),(/2-wid,2+wid/),(/xBootLow_ssh1_4xCO2,xBootLow_ssh1_4xCO2/), polyres)
  iumm01  = gsn_add_polyline(wks,plot(8),(/2-wid,2+wid/),(/xBootHi_ssh1_4xCO2,xBootHi_ssh1_4xCO2/), polyres)
;----------------------------------------------------------------------------------------------------
  draw(plot(6))
  draw(plot(7))
  draw(plot(8))

end 