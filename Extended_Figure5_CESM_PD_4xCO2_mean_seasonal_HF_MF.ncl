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


;=========================================================================================
; plot SST gradient PD and future
begin 
; parameter constant 
PI  = 3.14159265
R   = 6378388.
rad = PI/180.
dens= 1025 ; reference density of seawater
sec = 3600*24*30.42 ;one month

latS   =  -5.
latN   =   10. 
lonL   =  180
lonR   =  280
;-----------------------------------------------------------------------
ncdir1  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/mean/flux/"

; read mean patterns
fil1    = systemfunc("ls " + ncdir1 + "CESM2_HR_monthly_heat_flux.10-50_mean.nc")
fil2    = systemfunc("ls " + ncdir1 + "CESM2_HR_monthly_momentum_flux.10-50_mean.nc")
  f1    = addfile(fil1, "r")
  f2    = addfile(fil2, "r")
 
     hf0_PD     = f1->hf_mon_PD(0,0:1,:,:)
     hf0_2xCO2  = f1->hf_mon_2xCO2(0,0:1,:,:)
     hf0_4xCO2  = f1->hf_mon_4xCO2(0,0:1,:,:)

     hf_PD      = hf0_PD/100.
     hf_2xCO2   = hf0_2xCO2/100.
     hf_4xCO2   = hf0_4xCO2/100.  
     copy_VarMeta(hf0_PD, hf_PD)
     copy_VarMeta(hf0_2xCO2, hf_2xCO2)
     copy_VarMeta(hf0_4xCO2, hf_4xCO2)
     printVarSummary(hf_PD)
     printVarSummary(hf_2xCO2)
     printVarSummary(hf_4xCO2)
     hf_diff   = hf_4xCO2 - hf_PD
     copy_VarMeta(hf_PD, hf_diff)
     hf_PD_mean = dim_avg_n_Wrap(hf_PD, 0)
     hf_2xCO2_mean = dim_avg_n_Wrap(hf_2xCO2, 0)
     hf_4xCO2_mean = dim_avg_n_Wrap(hf_4xCO2, 0)
     hf_diff_mean = dim_avg_n_Wrap(hf_diff, 0)
;=================================================================
     mf0_PD     = f2->mf_mon_PD(0,0:1,:,:)
     mf0_2xCO2  = f2->mf_mon_2xCO2(0,0:1,:,:)
     mf0_4xCO2  = f2->mf_mon_4xCO2(0,0:1,:,:)

     mf_PD      = mf0_PD/100.
     mf_2xCO2   = mf0_2xCO2/100.
     mf_4xCO2   = mf0_4xCO2/100.
     copy_VarMeta(mf0_PD, mf_PD)
     copy_VarMeta(mf0_2xCO2, mf_2xCO2)
     copy_VarMeta(mf0_4xCO2, mf_4xCO2)
     printVarSummary(mf_PD)
     printVarSummary(mf_2xCO2)
     printVarSummary(mf_4xCO2)
     mf_diff   = mf_4xCO2 - mf_PD
     copy_VarMeta(mf_PD, mf_diff)
     mf_PD_mean = dim_avg_n_Wrap(mf_PD, 0)
     mf_2xCO2_mean = dim_avg_n_Wrap(mf_2xCO2, 0)
     mf_4xCO2_mean = dim_avg_n_Wrap(mf_4xCO2, 0)
     mf_diff_mean = dim_avg_n_Wrap(mf_diff, 0)
;===================================================================
     pattern1 = smth9_Wrap(hf_PD_mean, 0.50, 0.25, False)
     pattern2 = smth9_Wrap(hf_2xCO2_mean, 0.50, 0.25, False)
     pattern3 = smth9_Wrap(hf_4xCO2_mean, 0.50, 0.25, False)
     pattern4 = smth9_Wrap(hf_diff_mean, 0.50, 0.25, False)
     pattern_hf = new((/2, dimsizes(pattern1&lat), dimsizes(pattern1&lon)/), typeof(pattern1))
     pattern_hf(0,:,:) = pattern1 
     pattern_hf(1,:,:) = pattern3

     pattern01 = smth9_Wrap(mf_PD_mean, 0.50, 0.25, False)
     pattern02 = smth9_Wrap(mf_2xCO2_mean, 0.50, 0.25, False)
     pattern03 = smth9_Wrap(mf_4xCO2_mean, 0.50, 0.25, False)
     pattern04 = smth9_Wrap(mf_diff_mean, 0.50, 0.25, False)
     pattern_mf = new((/2, dimsizes(pattern1&lat), dimsizes(pattern1&lon)/), typeof(pattern1))
     pattern_mf(0,:,:) = pattern01 
     pattern_mf(1,:,:) = pattern03

;=====================================================================================================
ncdir  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/mean/flux/mon/"
; read SST var
fil1    = systemfunc("ls " + ncdir + "CESM2_HR_monthly_momentum_flux.10-50.1degree.nc")
  f1    = addfile(fil1, "r")
  tiw_mf_PD0     = f1->mf_mon_PD(:,0:1,:,:)
  tiw_mf_2xCO20  = f1->mf_mon_2xCO2(:,0:1,:,:)
  tiw_mf_4xCO20  = f1->mf_mon_4xCO2(:,0:1,:,:)
  tiw_mf_PD      = dim_avg_n_Wrap(tiw_mf_PD0, 1)/100.
  tiw_mf_2xCO2   = dim_avg_n_Wrap(tiw_mf_2xCO20, 1)/100.
  tiw_mf_4xCO2   = dim_avg_n_Wrap(tiw_mf_4xCO20, 1)/100.
  copy_VarMeta(tiw_mf_PD0(:,0,:,:), tiw_mf_PD)
  copy_VarMeta(tiw_mf_2xCO20(:,0,:,:), tiw_mf_2xCO2)
  copy_VarMeta(tiw_mf_4xCO20(:,0,:,:), tiw_mf_4xCO2)
  printVarSummary(tiw_mf_PD)

  tiw_mf_PD_ac1    = clmMonTLL(tiw_mf_PD)
  tiw_mf_2xCO2_ac1 = clmMonTLL(tiw_mf_2xCO2)
  tiw_mf_4xCO2_ac1 = clmMonTLL(tiw_mf_4xCO2)
  printVarSummary(tiw_mf_PD_ac1)

  tiw_mf_PD_ac      = (tiw_mf_PD_ac1(8,:,:)   +tiw_mf_PD_ac1(9,:,:)   +tiw_mf_PD_ac1(10,:,:)   -tiw_mf_PD_ac1(2,:,:)   -tiw_mf_PD_ac1(3,:,:)   -tiw_mf_PD_ac1(4,:,:))/6.
  tiw_mf_2xCO2_ac   = (tiw_mf_2xCO2_ac1(8,:,:)+tiw_mf_2xCO2_ac1(9,:,:)+tiw_mf_2xCO2_ac1(10,:,:)-tiw_mf_2xCO2_ac1(2,:,:)-tiw_mf_2xCO2_ac1(3,:,:)-tiw_mf_2xCO2_ac1(4,:,:))/6.
  tiw_mf_4xCO2_ac   = (tiw_mf_4xCO2_ac1(8,:,:)+tiw_mf_4xCO2_ac1(9,:,:)+tiw_mf_4xCO2_ac1(10,:,:)-tiw_mf_4xCO2_ac1(2,:,:)-tiw_mf_4xCO2_ac1(3,:,:)-tiw_mf_4xCO2_ac1(4,:,:))/6.
  
  copy_VarMeta(tiw_mf_PD(0,:,:), tiw_mf_PD_ac)
  copy_VarMeta(tiw_mf_2xCO2(0,:,:), tiw_mf_2xCO2_ac)
  copy_VarMeta(tiw_mf_4xCO2(0,:,:), tiw_mf_4xCO2_ac)

  pattern001 = smth9_Wrap(tiw_mf_PD_ac, 0.50, 0.25, False)
  pattern002 = smth9_Wrap(tiw_mf_2xCO2_ac, 0.50, 0.25, False)
  pattern003 = smth9_Wrap(tiw_mf_4xCO2_ac, 0.50, 0.25, False)
printMinMax(pattern001, 0)
  pattern_mf_ac = new((/2, dimsizes(pattern001&lat), dimsizes(pattern001&lon)/), typeof(pattern001))
  pattern_mf_ac(0,:,:) = pattern001
  pattern_mf_ac(1,:,:) = pattern003
;=========================================================================================================
  fil2    = systemfunc("ls " + ncdir + "CESM2_HR_monthly_heat_flux.10-50.1degree.nc")
  f2    = addfile(fil2, "r")
  tiw_hf_PD0     = f2->hf_mon_PD(:,0:1,:,:)
  tiw_hf_2xCO20  = f2->hf_mon_2xCO2(:,0:1,:,:)
  tiw_hf_4xCO20  = f2->hf_mon_4xCO2(:,0:1,:,:)
  tiw_hf_PD      = dim_avg_n_Wrap(tiw_hf_PD0, 1)/100.
  tiw_hf_2xCO2   = dim_avg_n_Wrap(tiw_hf_2xCO20, 1)/100.
  tiw_hf_4xCO2   = dim_avg_n_Wrap(tiw_hf_4xCO20, 1)/100.
  copy_VarCoords(tiw_hf_PD0(:,0,:,:), tiw_hf_PD)
  copy_VarCoords(tiw_hf_2xCO20(:,0,:,:), tiw_hf_2xCO2)
  copy_VarCoords(tiw_hf_4xCO20(:,0,:,:), tiw_hf_4xCO2)
  printVarSummary(tiw_hf_PD)

  tiw_hf_PD_ac1    = clmMonTLL(tiw_hf_PD)
  tiw_hf_2xCO2_ac1 = clmMonTLL(tiw_hf_2xCO2)
  tiw_hf_4xCO2_ac1 = clmMonTLL(tiw_hf_4xCO2)
  printVarSummary(tiw_hf_PD_ac1)

  tiw_hf_PD_ac      = (tiw_hf_PD_ac1(8,:,:)   +tiw_hf_PD_ac1(9,:,:)   +tiw_hf_PD_ac1(10,:,:)   -tiw_hf_PD_ac1(2,:,:)   -tiw_hf_PD_ac1(3,:,:)   -tiw_hf_PD_ac1(4,:,:))/6.
  tiw_hf_2xCO2_ac   = (tiw_hf_2xCO2_ac1(8,:,:)+tiw_hf_2xCO2_ac1(9,:,:)+tiw_hf_2xCO2_ac1(10,:,:)-tiw_hf_2xCO2_ac1(2,:,:)-tiw_hf_2xCO2_ac1(3,:,:)-tiw_hf_2xCO2_ac1(4,:,:))/6.
  tiw_hf_4xCO2_ac   = (tiw_hf_4xCO2_ac1(8,:,:)+tiw_hf_4xCO2_ac1(9,:,:)+tiw_hf_4xCO2_ac1(10,:,:)-tiw_hf_4xCO2_ac1(2,:,:)-tiw_hf_4xCO2_ac1(3,:,:)-tiw_hf_4xCO2_ac1(4,:,:))/6.
  
  copy_VarMeta(tiw_hf_PD(0,:,:), tiw_hf_PD_ac)
  copy_VarMeta(tiw_hf_2xCO2(0,:,:), tiw_hf_2xCO2_ac)
  copy_VarMeta(tiw_hf_4xCO2(0,:,:), tiw_hf_4xCO2_ac)

  pattern0001 = smth9_Wrap(tiw_hf_PD_ac, 0.50, 0.25, False)
  pattern0002 = smth9_Wrap(tiw_hf_2xCO2_ac, 0.50, 0.25, False)
  pattern0003 = smth9_Wrap(tiw_hf_4xCO2_ac, 0.50, 0.25, False)
  printMinMax(pattern0001, 0)
  pattern_hf_ac = new((/2, dimsizes(pattern0001&lat), dimsizes(pattern0001&lon)/), typeof(pattern0001))
  pattern_hf_ac(0,:,:) = pattern0001
  pattern_hf_ac(1,:,:) = pattern0003
;==========================================================================================
; plot
  dj      = (/"~F22~a","~F22~b","~F22~c","~F22~d","~F22~e","~F22~f","~F22~g","~F22~h","~F22~i","~F22~j","~F22~k","~F22~l"/) 
  text=(/"PD","4xCO2"/)
  wks = gsn_open_wks("pdf","/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/revise_figures/CESM_PD_4xCO2_mean_seasonal_hf_mf")
  gsn_define_colormap(wks,"BlueWhiteOrangeRed")  
  plot=new(19,graphic)
  res                      = True         
  res@gsnDraw              = False        ; don't draw yet
  res@gsnFrame             = False        ; don't advance frame yet
  res@gsnAddCyclic         = False      
  res@cnFillOn             = True               ; turn on color for contours
  res@cnLinesOn            = False              ; turn off contour lines
  res@mpFillOn             = False
  res@cnLineLabelsOn       = False              ; turn off contour line labels
  res@cnInfoLabelOn        = False                ; turn off info 
  res@mpLandFillColor      = "gray"             ; color of land
  res@mpGeophysicalLineColor      = "gray30"
  res@mpGeophysicalLineThicknessF = 0.7
  res@mpMaxLatF            = latN
  res@mpMinLatF            = latS
  res@mpMaxLonF            = lonR
  res@mpMinLonF            = lonL  
  res@mpCenterLonF         = 200
  res@tmXBOn               = False 
  res@tmXBMode              = "Explicit"
  res@tmXBValues            = (/180,200,220,240,260,280/)
  res@tmXBLabels            = (/" ","","","","",""/)
  res@tmXBMinorValues       = ispan(120,280,20)
  res@tmXTOn   = False
  res@tmYROn   = False
  res@tmYLMode = "Explicit"
  res@tmYLValues = (/-20,-10,-5,0,5,10,20/)
  res@tmYLLabels = (/"20~S~o~N~S","10~S~o~N~S","5~S~o~N~S","0","5~S~o~N~N","10~S~o~N~N","20~S~o~N~N"/)
  res@tmYLMinorValues =ispan(-20,20,5) 
  res@cnSmoothingOn         = True
  res@cnSmoothingTensionF   = 0.05
  res@tmXBLabelFontHeightF     = 0.008
  res@tmYLLabelFontHeightF     = 0.008
  res@tmXBLabelFontThicknessF  = 1.5
  res@tmYLLabelFontThicknessF  = 1.5 
  res@tmBorderThicknessF       = 2.5
  res@tmXBMajorLengthF        = -0.006
  res@tmYLMajorLengthF        = -0.006
  res@tmXBMinorLengthF        = -0.003
  res@tmYLMinorLengthF        = -0.003 
  res@tmXBMajorThicknessF     = 1.5
  res@tmYLMajorThicknessF     = 1.5
  res@tmXBMinorThicknessF     = 1
  res@tmYLMinorThicknessF     = 1
  res@tmXBLabelDeltaF          = -0.4
  res@tmYLLabelDeltaF          = -0.4
  res@cnFillOn              = True                 ; turn on color fill

  res@tmXBOn   = True
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

  res@lbBoxEndCapStyle     = "TriangleBothEnds"
  res@pmLabelBarWidthF     = 0.20
  ; res@pmLabelBarHeightF    = 0.55
  res@pmLabelBarOrthogonalPosF = 0.19
  ; res@pmLabelBarParallelPosF   = 0.6
  res@lbLabelFontHeightF       =0.008

  res@lbTitlePosition  = "Bottom"
  res@lbTitleDirection = "Across" 
  res@lbTitleFontHeightF= 0.010
  res@lbTitleOffsetF    = 0.05

 ;res@cnFillDrawOrder       = "PreDraw"
  res@gsnRightString          = " "
  res@tiMainString            = " "
  res@mpShapeMode             = "FreeAspect" 
  res@vpHeightF            = 0.10
  res@vpWidthF             = 0.20
  res@vpXF                    = 0.08
  res@vpYF                    = 0.9

 yspace1=0.03
 yspace2 = 0.03
 yf=0.9
res@tiDeltaF=3
; ----->1 plot mf
res@lbTitleString    = "10~S~2~N~ m~S~2~N~/s~S~2~N"
res@vpXF    = 0.08
res@cnLevelSelectionMode  = "ExplicitLevels"
res@cnLevels               := fspan(-1.5,1.5,31)
do i=0,1
res@vpYF    = yf-i*yspace1-i*res@vpHeightF
res@tiYAxisString=text(i)
res@tiYAxisOffsetXF=-0.02
res@tiYAxisFontHeightF=0.014

if i.eq.0 then
res@gsnLeftString  = dj(0)
res@gsnLeftStringFontHeightF = 0.012
res@gsnCenterString = "PD MF_mean"
res@gsnCenterStringFontHeightF = 0.012
res@tiMainString=""
res@lbLabelBarOn    = False
res@tmXBLabels = (/"","","","","",""/)
end if

if(i .eq. 1)then
res@gsnLeftString  = dj(1) 
res@gsnLeftStringFontHeightF = 0.012  
res@gsnCenterString = "  4xCO2 MF_mean"
res@lbLabelBarOn    = True
res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W","80~S~o~N~W"/)  
end if

plot(i)  = gsn_csm_contour_map(wks, pattern_mf(i,:,:),res)
end do

delete([/res@tiYAxisString,res@tiYAxisOffsetXF,res@tiYAxisFontHeightF/])

; ----->2 plot hf
res@vpXF    = res@vpWidthF+res@vpXF+yspace2
res@lbTitleString    = "~S~o~N~C m/s"
delete(res@cnLevels)
res@cnLevels                = fspan(-0.05,0.05,21)
res@tmYLLabels = (/"","","","","","",""/)
do i=0,1
res@vpYF    = yf-i*yspace1-i*res@vpHeightF
if i.eq.0 then
res@gsnLeftString  = dj(2)
res@gsnLeftStringFontHeightF = 0.012
res@gsnCenterString = "PD HF_mean"
res@lbLabelBarOn    = False
res@tiYAxisString   = ""
res@tmXBLabels = (/"","","","","",""/)
end if

if i.eq.1 then
res@gsnLeftString  = dj(3) 
res@gsnLeftStringFontHeightF = 0.012  
res@gsnCenterString =  "4xCO2 HF_mean"
res@lbLabelBarOn    = True
res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W","80~S~o~N~W"/)  
end if

plot(2+i)  = gsn_csm_contour_map(wks,pattern_hf(i,:,:),res)
end do
;------------------------------------------------------------
; plot seasonal cycle 

; ----->3  plot mf seasonal cycle
res@lbTitleString    = "10~S~2~N~ m~S~2~N~/s~S~2~N"
res@vpXF    = res@vpWidthF+res@vpXF+2*yspace2
res@cnLevelSelectionMode  = "ExplicitLevels"
res@cnLevels               := fspan(-1.5,1.5,31)
do i=0,1
res@vpYF    = yf-i*yspace1-i*res@vpHeightF
res@tiYAxisString=""
res@tiYAxisFontHeightF=0.014

if i.eq.0 then
res@gsnLeftString  = dj(4)
res@gsnLeftStringFontHeightF = 0.012
res@gsnCenterString = "PD MF_ac"
res@gsnCenterStringFontHeightF = 0.012
res@tiMainString=""
res@lbLabelBarOn    = False
res@tmXBLabels = (/"","","","","",""/)
end if

if(i .eq. 1)then
res@gsnLeftString  = dj(5) 
res@gsnLeftStringFontHeightF = 0.012  
res@gsnCenterString = "  4xCO2 MF_ac"
res@lbLabelBarOn    = True
res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W","80~S~o~N~W"/)  
end if

plot(4+i)  = gsn_csm_contour_map(wks, pattern_mf_ac(i,:,:),res)
end do
;====================================================================================
delete([/res@tiYAxisString,res@tiYAxisOffsetXF,res@tiYAxisFontHeightF/])
delete(res@cnLevels)
; ----->2 plot hf seasonal cycle 
res@vpXF    = res@vpWidthF+res@vpXF+yspace2
res@lbTitleString    = "~S~o~N~C m/s"
res@cnLevels                = fspan(-0.05,0.05,21)
res@tmYLLabels = (/"","","","","","",""/)
do i=0,1
res@vpYF    = yf-i*yspace1-i*res@vpHeightF
if i.eq.0 then
res@gsnLeftString  = dj(6)
res@gsnLeftStringFontHeightF = 0.012
res@gsnCenterString = "PD HF_ac"
  res@lbLabelBarOn    = False
  res@tmXBLabels = (/"","","","","",""/)
end if

if i.eq.1 then
res@gsnLeftString  = dj(7) 
res@gsnLeftStringFontHeightF = 0.012  
res@gsnCenterString =  "4xCO2 HF_ac"
res@lbLabelBarOn    = True
res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W","80~S~o~N~W"/)  
end if

plot(6+i)  = gsn_csm_contour_map(wks,pattern_hf_ac(i,:,:),res)
end do


plres=True
plres@gsLineThicknessF=2
plres@gsLineDashPattern= 2
plres@gsLineColor="blue"
plres@gsLineOpacityF=100

w1=360-160
e1=360-110

w2=360-160
e2=360-100

w4=360-160
e4=360-120

s1 = 2
n1 = 5 

s2=-1
n2=1

dum1=gsn_add_polyline(wks, plot(0), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum2=gsn_add_polyline(wks, plot(1), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum3=gsn_add_polyline(wks, plot(4), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum4=gsn_add_polyline(wks, plot(5), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)

plres@gsLineColor="red"
eum1=gsn_add_polyline(wks, plot(0), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
eum2=gsn_add_polyline(wks, plot(1), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
eum3=gsn_add_polyline(wks, plot(4), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
eum4=gsn_add_polyline(wks, plot(5), (/w4,e4,e4,w4,w4/), (/s2,s2,n2,n2,s2/), plres)

w3=360-160
e3=360-110

s3= 0.
n3= 5.
plres=True
plres@gsLineThicknessF=2
plres@gsLineDashPattern= 2
plres@gsLineColor="blue"
plres@gsLineOpacityF=100

fum1=gsn_add_polyline(wks, plot(2), (/w3,e3,e3,w3,w3/), (/s3,s3,n3,n3,s3/), plres)
fum2=gsn_add_polyline(wks, plot(3), (/w3,e3,e3,w3,w3/), (/s3,s3,n3,n3,s3/), plres)
fum3=gsn_add_polyline(wks, plot(6), (/w3,e3,e3,w3,w3/), (/s3,s3,n3,n3,s3/), plres)
fum4=gsn_add_polyline(wks, plot(7), (/w3,e3,e3,w3,w3/), (/s3,s3,n3,n3,s3/), plres)

txres                       = True
txres@txPerimOn             = False  
txres@txBackgroundFillColor = "transparent"
txres@txFontHeightF         = 0.02
txres@txFontThicknessF      = 4
txres@txFont                = 22
txid1  = gsn_create_text(wks, "Mean State", txres)
txid2  = gsn_create_text(wks, "Seasonality", txres)

amres = True
amres@amParallelPosF   =  0.6    ; This is the center of the plot.
amres@amOrthogonalPosF = -0.8    ; This is the top edge of the plot.
amres@amJust           = "BottomCenter"
annoid1 = gsn_add_annotation(plot(0), txid1, amres)
annoid2 = gsn_add_annotation(plot(4), txid2, amres)

draw(plot(0))
draw(plot(1))
draw(plot(2))
draw(plot(3))
draw(plot(4))
draw(plot(5))
draw(plot(6))
draw(plot(7))

delete([/hf0_PD, hf0_2xCO2, hf0_4xCO2/])
delete([/hf_PD,  hf_2xCO2,  hf_4xCO2/])

delete([/mf0_PD,mf0_2xCO2,mf0_4xCO2/])
delete([/mf_PD, mf_2xCO2, mf_4xCO2/])
;=============================================================
; read data
;=============================================================

ncdir01  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/mean/flux/mon/"
fil01    = systemfunc("ls " + ncdir01 + "CESM2_HR_monthly_heat_flux.10-50.1degree.nc")
fil02    = systemfunc("ls " + ncdir01 + "CESM2_HR_monthly_momentum_flux.10-50.1degree.nc")
  f01    = addfile(fil01, "r")
  f02    = addfile(fil02, "r")
 
   hf0_PD     = f01->hf_mon_PD(:,0:1,:,:)
   hf0_2xCO2  = f01->hf_mon_2xCO2(:,0:1,:,:)
   hf0_4xCO2  = f01->hf_mon_4xCO2(:,0:1,:,:)
   
   mf0_PD     = f02->mf_mon_PD(:,0:1,:,:)
   mf0_2xCO2  = f02->mf_mon_2xCO2(:,0:1,:,:)
   mf0_4xCO2  = f02->mf_mon_4xCO2(:,0:1,:,:)
            
   hf_PD    = dim_avg_n_Wrap(hf0_PD,1)/100.
   hf_2xCO2 = dim_avg_n_Wrap(hf0_2xCO2, 1)/100.
   hf_4xCO2 = dim_avg_n_Wrap(hf0_4xCO2, 1)/100.
   copy_VarMeta(hf0_PD(:,0,:,:), hf_PD)
   copy_VarMeta(hf0_2xCO2(:,0,:,:), hf_2xCO2)
   copy_VarMeta(hf0_4xCO2(:,0,:,:), hf_4xCO2)

   mf_PD    = dim_avg_n_Wrap(mf0_PD,1)/100.
   mf_2xCO2 = dim_avg_n_Wrap(mf0_2xCO2, 1)/100.
   mf_4xCO2 = dim_avg_n_Wrap(mf0_4xCO2, 1)/100.
   copy_VarMeta(mf0_PD(:,0,:,:), mf_PD)
   copy_VarMeta(mf0_2xCO2(:,0,:,:), mf_2xCO2)
   copy_VarMeta(mf0_4xCO2(:,0,:,:), mf_4xCO2)

   hf_PD_ac0      = (hf_PD(8::12,:,:)   + hf_PD(9::12,:,:)   + hf_PD(10::12,:,:)   -hf_PD(2::12,:,:)   -hf_PD(3::12,:,:)   -hf_PD(4::12,:,:))/6.
   hf_2xCO2_ac0   = (hf_2xCO2(8::12,:,:)+ hf_2xCO2(9::12,:,:)+ hf_2xCO2(10::12,:,:)-hf_2xCO2(2::12,:,:)-hf_2xCO2(3::12,:,:)-hf_2xCO2(4::12,:,:))/6.
   hf_4xCO2_ac0   = (hf_4xCO2(8::12,:,:)+ hf_4xCO2(9::12,:,:)+ hf_4xCO2(10::12,:,:)-hf_4xCO2(2::12,:,:)-hf_4xCO2(3::12,:,:)-hf_4xCO2(4::12,:,:))/6.
   copy_VarMeta(hf_PD(0::12,:,:),    hf_PD_ac0)
   copy_VarMeta(hf_2xCO2(0::12,:,:), hf_2xCO2_ac0)
   copy_VarMeta(hf_4xCO2(0::12,:,:), hf_4xCO2_ac0)

   mf_PD_ac0      = (mf_PD(8::12,:,:)   + mf_PD(9::12,:,:)   + mf_PD(10::12,:,:)   -mf_PD(2::12,:,:)   -mf_PD(3::12,:,:)   -mf_PD(4::12,:,:))/6.
   mf_2xCO2_ac0   = (mf_2xCO2(8::12,:,:)+ mf_2xCO2(9::12,:,:)+ mf_2xCO2(10::12,:,:)-mf_2xCO2(2::12,:,:)-mf_2xCO2(3::12,:,:)-mf_2xCO2(4::12,:,:))/6.
   mf_4xCO2_ac0   = (mf_4xCO2(8::12,:,:)+ mf_4xCO2(9::12,:,:)+ mf_4xCO2(10::12,:,:)-mf_4xCO2(2::12,:,:)-mf_4xCO2(3::12,:,:)-mf_4xCO2(4::12,:,:))/6.
   copy_VarMeta(mf_PD(0::12,:,:),    mf_PD_ac0)
   copy_VarMeta(mf_2xCO2(0::12,:,:), mf_2xCO2_ac0)
   copy_VarMeta(mf_4xCO2(0::12,:,:), mf_4xCO2_ac0)
;-----------------------------------------------------------------
hf1_PD    = wgt_areaave_Wrap(hf_PD(:,{0:5},{210:250}), 1.0, 1.0, 0)
hf1_2xCO2 = wgt_areaave_Wrap(hf_2xCO2(:,{0:5},{210:250}), 1.0, 1.0, 0)
hf1_4xCO2 = wgt_areaave_Wrap(hf_4xCO2(:,{0:5},{210:250}), 1.0, 1.0, 0)
;------------------------------------------------------------------
mf1_PD    = wgt_areaave_Wrap(mf_PD(:,{2:5},{200:250}), 1.0, 1.0, 0)
mf1_2xCO2 = wgt_areaave_Wrap(mf_2xCO2(:,{2:5},{200:250}), 1.0, 1.0, 0)
mf1_4xCO2 = wgt_areaave_Wrap(mf_4xCO2(:,{2:5},{200:250}), 1.0, 1.0, 0)

mf2_PD    = wgt_areaave_Wrap(mf_PD(:,{-1:1},{200:260}), 1.0, 1.0, 0)
mf2_2xCO2 = wgt_areaave_Wrap(mf_2xCO2(:,{-1:1},{200:260}), 1.0, 1.0, 0)
mf2_4xCO2 = wgt_areaave_Wrap(mf_4xCO2(:,{-1:1},{200:260}), 1.0, 1.0, 0)
; seasonal cycle
hf1_PD_ac    = wgt_areaave_Wrap(hf_PD_ac0(:,{0:5},{210:250}), 1.0, 1.0, 0)
hf1_2xCO2_ac = wgt_areaave_Wrap(hf_2xCO2_ac0(:,{0:5},{210:250}), 1.0, 1.0, 0)
hf1_4xCO2_ac = wgt_areaave_Wrap(hf_4xCO2_ac0(:,{0:5},{210:250}), 1.0, 1.0, 0)

mf1_PD_ac    = wgt_areaave_Wrap(mf_PD_ac0(:,{2:5},{200:250}), 1.0, 1.0, 0)
mf1_2xCO2_ac = wgt_areaave_Wrap(mf_2xCO2_ac0(:,{2:5},{200:250}), 1.0, 1.0, 0)
mf1_4xCO2_ac = wgt_areaave_Wrap(mf_4xCO2_ac0(:,{2:5},{200:250}), 1.0, 1.0, 0)

mf2_PD_ac    = wgt_areaave_Wrap(mf_PD_ac0(:,{-1:1},{200:260}), 1.0, 1.0, 0)
mf2_2xCO2_ac = wgt_areaave_Wrap(mf_2xCO2_ac0(:,{-1:1},{200:240}), 1.0, 1.0, 0)
mf2_4xCO2_ac = wgt_areaave_Wrap(mf_4xCO2_ac0(:,{-1:1},{200:240}), 1.0, 1.0, 0)
;=======================================================================

stat  = 0                         ; 0=mean; 1=variance; 2=std. dev; ...
n     = 360                        ; sub-sample size;  default is n=N
nBoot = 10000                      ; bootstrap replications
N     = 360
opt   = False                     ; use defaults
if (isvar("n") .and. (n.lt.N)) then
    opt   = True
    opt@sample_size = n           ; default N
end if

 BootStrap_hf1_PD     = bootstrap_stat(hf1_PD, stat, nBoot, 0, opt)
 BootStrap_hf1_2xCO2  = bootstrap_stat(hf1_2xCO2, stat, nBoot, 0, opt)
 BootStrap_hf1_4xCO2  = bootstrap_stat(hf1_4xCO2, stat, nBoot, 0, opt)

 BootStrap_mf1_PD     = bootstrap_stat(mf1_PD, stat, nBoot, 0, opt)
 BootStrap_mf1_2xCO2  = bootstrap_stat(mf1_2xCO2, stat, nBoot, 0, opt)
 BootStrap_mf1_4xCO2  = bootstrap_stat(mf1_4xCO2, stat, nBoot, 0, opt)

 BootStrap_mf2_PD     = bootstrap_stat(mf2_PD, stat, nBoot, 0, opt)
 BootStrap_mf2_2xCO2  = bootstrap_stat(mf2_2xCO2, stat, nBoot, 0, opt)
 BootStrap_mf2_4xCO2  = bootstrap_stat(mf2_4xCO2, stat, nBoot, 0, opt)

 ; seasonal cycle 
BootStrap_hf1_PD_ac     = bootstrap_stat(hf1_PD_ac, stat, nBoot, 0, opt)
BootStrap_hf1_2xCO2_ac  = bootstrap_stat(hf1_2xCO2_ac, stat, nBoot, 0, opt)
BootStrap_hf1_4xCO2_ac  = bootstrap_stat(hf1_4xCO2_ac, stat, nBoot, 0, opt)

BootStrap_mf1_PD_ac     = bootstrap_stat(mf1_PD_ac, stat, nBoot, 0, opt)
BootStrap_mf1_2xCO2_ac  = bootstrap_stat(mf1_2xCO2_ac, stat, nBoot, 0, opt)
BootStrap_mf1_4xCO2_ac  = bootstrap_stat(mf1_4xCO2_ac, stat, nBoot, 0, opt)

BootStrap_mf2_PD_ac     = bootstrap_stat(mf2_PD_ac, stat, nBoot, 0, opt)
BootStrap_mf2_2xCO2_ac  = bootstrap_stat(mf2_2xCO2_ac, stat, nBoot, 0, opt)
BootStrap_mf2_4xCO2_ac  = bootstrap_stat(mf2_4xCO2_ac, stat, nBoot, 0, opt)

;-----------------------------------------------
 xBoot_hf1_PD       = BootStrap_hf1_PD[1]       
 xBoot_mf1_PD       = BootStrap_mf1_PD[1]       
 xBoot_mf2_PD       = BootStrap_mf2_PD[1]           
;=====================================================
xBoot_hf1_2xCO2       = BootStrap_hf1_2xCO2[1]       
xBoot_mf1_2xCO2       = BootStrap_mf1_2xCO2[1]       
xBoot_mf2_2xCO2       = BootStrap_mf2_2xCO2[1]        
;----------------------------------------------------
xBoot_hf1_4xCO2       = BootStrap_hf1_4xCO2[1]       
xBoot_mf1_4xCO2       = BootStrap_mf1_4xCO2[1]       
xBoot_mf2_4xCO2       = BootStrap_mf2_4xCO2[1]       
;========================================================================
 xBoot_hf01_PD       = BootStrap_hf1_PD[0]       
 xBoot_mf01_PD       = BootStrap_mf1_PD[0]       
 xBoot_mf02_PD       = BootStrap_mf2_PD[0]          
;=====================================================
xBoot_hf01_2xCO2       = BootStrap_hf1_2xCO2[0]       
xBoot_mf01_2xCO2       = BootStrap_mf1_2xCO2[0]       
xBoot_mf02_2xCO2       = BootStrap_mf2_2xCO2[0]          
;-----------------------------------------------------
xBoot_hf01_4xCO2       = BootStrap_hf1_4xCO2[0]     
xBoot_mf01_4xCO2       = BootStrap_mf1_4xCO2[0]     
xBoot_mf02_4xCO2       = BootStrap_mf2_4xCO2[0]     
;---------------------------------------------------
xBoot_hf1_PD_ac    = BootStrap_hf1_PD_ac[1]   ; mean value
xBoot_mf1_PD_ac    = BootStrap_mf1_PD_ac[1]   ; mean value
xBoot_mf2_PD_ac    = BootStrap_mf2_PD_ac[1]   ; mean value  
;----------------------------------------------------
xBoot_hf1_2xCO2_ac    = BootStrap_hf1_2xCO2_ac[1]   ; mean value
xBoot_mf1_2xCO2_ac    = BootStrap_mf1_2xCO2_ac[1]   ; mean value
xBoot_mf2_2xCO2_ac    = BootStrap_mf2_2xCO2_ac[1]   ; mean value  
;----------------------------------------------------   
xBoot_hf1_4xCO2_ac      = BootStrap_hf1_4xCO2_ac[1]   ; mean value 
xBoot_mf1_4xCO2_ac      = BootStrap_mf1_4xCO2_ac[1]   ; mean value 
xBoot_mf2_4xCO2_ac      = BootStrap_mf2_4xCO2_ac[1]   ; mean value 
;------------------------------------------------
xBoot_hf01_PD_ac      = BootStrap_hf1_PD_ac[0]      
xBoot_mf01_PD_ac      = BootStrap_mf1_PD_ac[0]      
xBoot_mf02_PD_ac      = BootStrap_mf2_PD_ac[0]     
;------------------------------------------------
xBoot_hf01_2xCO2_ac   = BootStrap_hf1_2xCO2_ac[0]      
xBoot_mf01_2xCO2_ac   = BootStrap_mf1_2xCO2_ac[0]      
xBoot_mf02_2xCO2_ac   = BootStrap_mf2_2xCO2_ac[0] 
;------------------------------------------------
xBoot_hf01_4xCO2_ac   = BootStrap_hf1_4xCO2_ac[0]   
xBoot_mf01_4xCO2_ac   = BootStrap_mf1_4xCO2_ac[0]   
xBoot_mf02_4xCO2_ac   = BootStrap_mf2_4xCO2_ac[0] 
;=======================================================================
xBootLow_hf1_PD    = bootstrap_estimate(xBoot_hf01_PD, 0.05, False)  ;  2.5% lower confidence bound
xBootMed_hf1_PD    = bootstrap_estimate(xBoot_hf01_PD, 0.500, False)  ; 50.0% median of bootstrapped estimates
 xBootHi_hf1_PD    = bootstrap_estimate(xBoot_hf01_PD, 0.95, False)  ; 97.5% upper confidence bound

xBootLow_mf1_PD    = bootstrap_estimate(xBoot_mf01_PD, 0.05, False)  ;  2.5% lower confidence bound
xBootMed_mf1_PD    = bootstrap_estimate(xBoot_mf01_PD, 0.500, False)  ; 50.0% median of bootstrapped estimates
 xBootHi_mf1_PD    = bootstrap_estimate(xBoot_mf01_PD, 0.95, False)  ; 97.5% upper confidence bound

xBootLow_mf2_PD    = bootstrap_estimate(xBoot_mf02_PD, 0.05, False)  ;  2.5% lower confidence bound
xBootMed_mf2_PD    = bootstrap_estimate(xBoot_mf02_PD, 0.500, False)  ; 50.0% median of bootstrapped estimates
 xBootHi_mf2_PD    = bootstrap_estimate(xBoot_mf02_PD, 0.95, False)  ; 97.5% upper confidence bound
;============================================================================
; seasonal cycle
xBootLow_hf1_PD_ac    = bootstrap_estimate(xBoot_hf01_PD_ac, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_hf1_PD_ac    = bootstrap_estimate(xBoot_hf01_PD_ac, 0.500, False)   ; 50.0% median of bootstrapped estimates  
 xBootHi_hf1_PD_ac    = bootstrap_estimate(xBoot_hf01_PD_ac, 0.95, False)   ; 97.5% upper confidence bound 

xBootLow_mf1_PD_ac    = bootstrap_estimate(xBoot_mf01_PD_ac, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_mf1_PD_ac    = bootstrap_estimate(xBoot_mf01_PD_ac, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_mf1_PD_ac    = bootstrap_estimate(xBoot_mf01_PD_ac, 0.95, False)   ; 97.5% upper confidence bound

xBootLow_mf2_PD_ac    = bootstrap_estimate(xBoot_mf02_PD_ac, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_mf2_PD_ac    = bootstrap_estimate(xBoot_mf02_PD_ac, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_mf2_PD_ac    = bootstrap_estimate(xBoot_mf02_PD_ac, 0.95, False)   ; 97.5% upper confidence bound
;===========================================================================
;-----------------------------------------------------
xBootLow_hf1_2xCO2    = bootstrap_estimate(xBoot_hf01_2xCO2, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_hf1_2xCO2    = bootstrap_estimate(xBoot_hf01_2xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_hf1_2xCO2    = bootstrap_estimate(xBoot_hf01_2xCO2, 0.95, False)   ; 97.5% upper confidence bound

xBootLow_mf1_2xCO2    = bootstrap_estimate(xBoot_mf01_2xCO2, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_mf1_2xCO2    = bootstrap_estimate(xBoot_mf01_2xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_mf1_2xCO2    = bootstrap_estimate(xBoot_mf01_2xCO2, 0.95, False)  ; 97.5% upper confidence bound

xBootLow_mf2_2xCO2    = bootstrap_estimate(xBoot_mf02_2xCO2, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_mf2_2xCO2    = bootstrap_estimate(xBoot_mf02_2xCO2, 0.500, False)  ; 50.0% median of bootstrapped estimates
 xBootHi_mf2_2xCO2    = bootstrap_estimate(xBoot_mf02_2xCO2, 0.95, False)   ; 97.5% upper confidence bound

 ;--------seasonal cycle 
xBootLow_hf1_2xCO2_ac    = bootstrap_estimate(xBoot_hf01_2xCO2_ac, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_hf1_2xCO2_ac    = bootstrap_estimate(xBoot_hf01_2xCO2_ac, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_hf1_2xCO2_ac    = bootstrap_estimate(xBoot_hf01_2xCO2_ac, 0.95, False)   ; 97.5% upper confidence bound


xBootLow_mf1_2xCO2_ac    = bootstrap_estimate(xBoot_mf01_2xCO2_ac, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_mf1_2xCO2_ac    = bootstrap_estimate(xBoot_mf01_2xCO2_ac, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_mf1_2xCO2_ac    = bootstrap_estimate(xBoot_mf01_2xCO2_ac, 0.95, False)   ; 97.5% upper confidence bound

xBootLow_mf2_2xCO2_ac    = bootstrap_estimate(xBoot_mf02_2xCO2_ac, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_mf2_2xCO2_ac    = bootstrap_estimate(xBoot_mf02_2xCO2_ac, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_mf2_2xCO2_ac    = bootstrap_estimate(xBoot_mf02_2xCO2_ac, 0.95, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------
xBootLow_hf1_4xCO2    = bootstrap_estimate(xBoot_hf01_4xCO2, 0.05, False)  ;  2.5% lower confidence bound
xBootMed_hf1_4xCO2    = bootstrap_estimate(xBoot_hf01_4xCO2, 0.500, False)  ; 50.0% median of bootstrapped estimates
 xBootHi_hf1_4xCO2    = bootstrap_estimate(xBoot_hf01_4xCO2, 0.95, False)  ; 97.5% upper confidence bound
   
xBootLow_mf1_4xCO2    = bootstrap_estimate(xBoot_mf01_4xCO2, 0.05, False)  ;  2.5% lower confidence bound
xBootMed_mf1_4xCO2    = bootstrap_estimate(xBoot_mf01_4xCO2, 0.500, False)  ; 50.0% median of bootstrapped estimates
 xBootHi_mf1_4xCO2    = bootstrap_estimate(xBoot_mf01_4xCO2, 0.95, False)  ; 97.5% upper confidence bound

xBootLow_mf2_4xCO2    = bootstrap_estimate(xBoot_mf02_4xCO2, 0.05, False) ;  2.5% lower confidence bound
xBootMed_mf2_4xCO2    = bootstrap_estimate(xBoot_mf02_4xCO2, 0.500, False)  ; 50.0% median of bootstrapped estimates
 xBootHi_mf2_4xCO2    = bootstrap_estimate(xBoot_mf02_4xCO2, 0.95, False)  ; 97.5% upper confidence bound
 ;--------seasonal cycle 
 xBootLow_hf1_4xCO2_ac    = bootstrap_estimate(xBoot_hf01_4xCO2_ac, 0.05, False)   ;  2.5% lower confidence bound
 xBootMed_hf1_4xCO2_ac    = bootstrap_estimate(xBoot_hf01_4xCO2_ac, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_hf1_4xCO2_ac    = bootstrap_estimate(xBoot_hf01_4xCO2_ac, 0.95, False)   ; 97.5% upper confidence bound
 
 xBootLow_mf1_4xCO2_ac    = bootstrap_estimate(xBoot_mf01_4xCO2_ac, 0.05, False)   ;  2.5% lower confidence bound
 xBootMed_mf1_4xCO2_ac    = bootstrap_estimate(xBoot_mf01_4xCO2_ac, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_mf1_4xCO2_ac    = bootstrap_estimate(xBoot_mf01_4xCO2_ac, 0.95, False)   ; 97.5% upper confidence bound
 
 xBootLow_mf2_4xCO2_ac    = bootstrap_estimate(xBoot_mf02_4xCO2_ac, 0.05, False)   ;  2.5% lower confidence bound
 xBootMed_mf2_4xCO2_ac    = bootstrap_estimate(xBoot_mf02_4xCO2_ac, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_mf2_4xCO2_ac    = bootstrap_estimate(xBoot_mf02_4xCO2_ac, 0.95, False)   ; 97.5% upper confidence bound
;==========================================================
stdarr1 = new((/3/),"double") ; mf1
stdarr1(0) =  xBoot_mf1_PD
stdarr1(1) =  xBoot_mf1_2xCO2
stdarr1(2) =  xBoot_mf1_4xCO2

stdarr2 = new((/3/),"double") ; mf2
stdarr2(0) =  -1 * xBoot_mf2_PD
stdarr2(1) =  -1 * xBoot_mf2_2xCO2
stdarr2(2) =  -1 * xBoot_mf2_4xCO2

stdarr3 = new((/3/),"double") ; hf1
stdarr3(0) =  xBoot_hf1_PD
stdarr3(1) =  xBoot_hf1_2xCO2
stdarr3(2) =  xBoot_hf1_4xCO2
;======================================
stdarr01 = new((/3/),"double") ; mf1 ;NEQ
stdarr01(0) =  xBoot_mf1_PD_ac
stdarr01(1) =  xBoot_mf1_2xCO2_ac
stdarr01(2) =  xBoot_mf1_4xCO2_ac

stdarr02 = new((/3/),"double") ; mf2
stdarr02(0) =  -1 * xBoot_mf2_PD_ac
stdarr02(1) =  -1 * xBoot_mf2_2xCO2_ac
stdarr02(2) =  -1 * xBoot_mf2_4xCO2_ac

stdarr03 = new((/3/),"double") ; hf1
stdarr03(0) =  xBoot_hf1_PD_ac
stdarr03(1) =  xBoot_hf1_2xCO2_ac
stdarr03(2) =  xBoot_hf1_4xCO2_ac

;=============================================
 xyres                       = True                 ; plot mods desired
 xyres@gsnFrame              = False                ; don't advance frame yet
 xyres@gsnDraw               = False
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 1.0                 ; bring bars down to zero
 xyres@trXMinF               = 0.5                    ; adds space on either end
 xyres@trXMaxF               = 4.                   ; of the 1st and last bars
 xyres@vpWidthF              = 0.10
 xyres@vpHeightF             = 0.10

 xyres@gsnXYBarChart         = True                 ; turn on bar chart
 xyres@gsnXYBarChartBarWidth = 0.24                  ; change bar widths
 xyres@tmXBOn                = True                 ; turn off tickmarks at bot
 xyres@tmXTOn                = False                  ; turn off tickmarks at bot
 xyres@tmXTBorderOn          = True 
 xyres@tmXBMode          = "Explicit"         ; explicit labels
 xyres@tmXBValues        = (/1,2,3,4/)
 xyres@tmXBLabels = (/"PD","2xCO2","4xCO2","" /)
 xyres@tmXBLabelFontHeightF = 0.006
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

 xyres@vpXF    = 0.08
 yspace=0.06
 yf=0.9
 xyres@tmYLOn         = True
 xyres@tmYROn        = False ; turn off right tickmarks
 xyres@tmYRBorderOn  = False ; turn off right axis
 xyres@tmXTBorderOn  = False 
 xyres@vpYF    = yf-yspace-2*res@vpHeightF-0.05
 xyres@gsnXYBarChartColors = (/"orange3"/)
 xyres@gsnCenterString       = "MF (NEQ)"
 xyres@gsnCenterStringOrthogonalPosF = -0.15
 xyres@gsnCenterStringParallelPosF = 0.45
 xyres@gsnLeftString  = dj(8)
 xyres@gsnLeftStringFontHeightF = 0.012
 xyres@gsnLeftStringOrthogonalPosF = 0.
 xyres@gsnLeftStringParallelPosF = -0.05
 xyres@tmYLMode          = "Explicit"         ; explicit labels
 xyres@tmYLValues        = (/0,0.2,0.4,0.6,0.8,1.0/)
 xyres@tmYLLabels        = (/"0","0.2","0.4","0.6","0.8","1.0"/)
 plot(8)  = gsn_csm_xy(wks,ispan(1,3,1),stdarr1(:),xyres)   ; advancing the frame
 xyres@tmYROn        = False ; turn off right tickmarks
 xyres@tmYLOn        = False
 xyres@tmYRBorderOn  = False                        ; turn off right axis
 xyres@tmYLBorderOn  = False  
 xyres@vpXF = xyres@vpWidthF+xyres@vpXF
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 1.0               ; bring bars down to zero
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 3.5   
 xyres@gsnXYBarChartColors = (/"cyan3"/)
 xyres@gsnCenterString       =  " -1xMF (EQ)"
 xyres@gsnCenterStringParallelPosF = 0.55
 xyres@gsnCenterStringOrthogonalPosF = -0.12
 xyres@gsnLeftString  = " "
 xyres@gsnLeftStringFontHeightF = 0.012
 xyres@gsnLeftStringOrthogonalPosF = 0.
 xyres@gsnLeftStringParallelPosF = -0.05
 plot(9)  = gsn_csm_xy(wks,ispan(1,3,1),stdarr2(:),xyres)   ; advancing the frame
 xyres@vpXF = 1.5*xyres@vpWidthF+xyres@vpXF+0.03
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@tmYLOn         = True 
 xyres@tmYROn        = False ; turn off right tickmarks
 xyres@tmYRBorderOn  = False    ; turn off right axis
 xyres@tmYLBorderOn  = True   ; turn off right axis
 xyres@vpWidthF              = 0.10
 xyres@vpHeightF             = 0.10
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.05                  ; bring bars down to zero
 xyres@trXMaxF               = 4.
                       ; turn off right axis
 xyres@gsnXYBarChartColors = (/"orange3"/)
 xyres@tmYLValues        = (/0,0.01,0.02,0.03,0.04,0.05/)
 xyres@tmYLLabels        = (/"0","0.01","0.02","0.03","0.04","0.05"/)
 xyres@gsnCenterString       = " HF_NEQ"
 xyres@gsnCenterStringOrthogonalPosF = -0.15
 xyres@gsnCenterStringParallelPosF = 0.40
 xyres@trXMinF               = 0.5                    ; adds space on either end
 xyres@trXMaxF               = 4.0                   ; of the 1st and last bars
 xyres@gsnLeftString  = dj(9)
 xyres@gsnLeftStringFontHeightF = 0.012
 xyres@gsnLeftStringOrthogonalPosF = 0.
 plot(10)  = gsn_csm_xy(wks,ispan(1,3,1),stdarr3(:),xyres)   ; advancing the frame
;-----------------------------------
 ; plot annual cycle 
 xyres@vpXF    = 1.5*xyres@vpWidthF+xyres@vpXF+2*yspace2
 yspace=0.06
 yf=0.9
 xyres@tmYLOn         = True
 xyres@tmYLBorderOn   = True 
 xyres@tmYROn        = False ; turn off right tickmarks
 xyres@tmYRBorderOn  = False ; turn off right axis
 xyres@vpYF    = yf-yspace-2*res@vpHeightF-0.05
 xyres@gsnXYBarChartColors = (/"orange3"/)
 xyres@gsnCenterString       = "MF (NEQ)"
 xyres@gsnCenterStringOrthogonalPosF = -0.15
 xyres@gsnCenterStringParallelPosF = 0.52
 xyres@gsnLeftString  = dj(10)
 xyres@gsnLeftStringFontHeightF = 0.012
 xyres@gsnLeftStringOrthogonalPosF = 0.
 xyres@gsnLeftStringParallelPosF = -0.05
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 1.0                  ; bring bars down to zero
 xyres@tmYLValues        = (/0,0.2,0.4,0.6,0.8,1.0/)
 xyres@tmYLLabels        = (/"0","0.2","0.4","0.6","0.8","1.0"/)
 plot(11)  = gsn_csm_xy(wks,ispan(1,3,1),stdarr01(:),xyres)   ; advancing the frame
 xyres@tmYROn        = False ; turn off right tickmarks
 xyres@tmYLOn        = False
 xyres@tmYRBorderOn  = False                         ; turn off right axis
 xyres@tmYLBorderOn  = False  
 xyres@vpXF = xyres@vpWidthF+xyres@vpXF
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 1.0                  ; bring bars down to zero
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 3.5   
 xyres@gsnXYBarChartColors = (/"cyan3"/)
 xyres@gsnCenterString       =  " -1xMF (EQ)"
 xyres@gsnCenterStringParallelPosF = 0.55
 xyres@gsnCenterStringOrthogonalPosF = -0.12
 xyres@gsnLeftString  = " "
 xyres@gsnLeftStringFontHeightF = 0.012
 xyres@gsnLeftStringOrthogonalPosF = 0.
 xyres@gsnLeftStringParallelPosF = -0.05
 plot(12)  = gsn_csm_xy(wks,ispan(1,3,1),stdarr02(:),xyres)   ; advancing the frame
 xyres@vpXF = 1.5*xyres@vpWidthF+xyres@vpXF+0.03
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 4. 
 xyres@tmYLOn         = True 
 xyres@tmYROn        = False ; turn off right tickmarks
 xyres@tmYRBorderOn  = False    ; turn off right axis
 xyres@tmYLBorderOn  = True    ; turn off right axis
 xyres@vpWidthF              = 0.10
 xyres@vpHeightF             = 0.10
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.05                  ; bring bars down to zero
                       ; turn off right axis
 xyres@gsnXYBarChartColors = (/"orange3"/)
 xyres@tmYLValues        = (/0,0.01,0.02,0.03,0.04,0.05/)
 xyres@tmYLLabels        = (/"0","0.01","0.02","0.03","0.04","0.05"/)
 xyres@gsnCenterString       = " -1*HF_NEQ"
 xyres@gsnCenterStringOrthogonalPosF = -0.15
 xyres@gsnCenterStringParallelPosF = 0.40
 xyres@trXMinF               = 0.5                    ; adds space on either end
 xyres@trXMaxF               = 4.0                   ; of the 1st and last bars
 xyres@gsnLeftString  = dj(11)
 xyres@gsnLeftStringFontHeightF = 0.012
 xyres@gsnLeftStringOrthogonalPosF = 0.
 plot(13)  = gsn_csm_xy(wks,ispan(1,3,1),stdarr03(:),xyres)   ; advancing the frame


 polyres=True
 polyres@gsLineDashPattern=0
 polyres@gsLineColor="black"
 polyres@gsLineThicknessF=1.0
 polyres@tfPolyDrawOrder = "postDraw"
;
 aumm1    = gsn_add_polyline(wks,plot(8),(/1,1/),(/xBootLow_mf1_PD,   xBootHi_mf1_PD/), polyres)
 aumm01   = gsn_add_polyline(wks,plot(8),(/2,2/),(/xBootLow_mf1_2xCO2,xBootHi_mf1_2xCO2/), polyres)
 aumm001  = gsn_add_polyline(wks,plot(8),(/3,3/),(/xBootLow_mf1_4xCO2,xBootHi_mf1_4xCO2/), polyres)
;
 bumm1    = gsn_add_polyline(wks,plot(9),(/1,1/),(/-1*xBootHi_mf2_PD,   -1*xBootLow_mf2_PD/), polyres)
 bumm01   = gsn_add_polyline(wks,plot(9),(/2,2/),(/-1*xBootHi_mf2_2xCO2,-1*xBootLow_mf2_2xCO2/), polyres)
 bumm001  = gsn_add_polyline(wks,plot(9),(/3,3/),(/-1*xBootHi_mf2_4xCO2,-1*xBootLow_mf2_4xCO2/), polyres)

 cumm1    = gsn_add_polyline(wks,plot(10),(/1,1/),(/xBootHi_hf1_PD,   xBootLow_hf1_PD/), polyres)
 cumm01   = gsn_add_polyline(wks,plot(10),(/2,2/),(/xBootHi_hf1_2xCO2,xBootLow_hf1_2xCO2/), polyres)
 cumm001  = gsn_add_polyline(wks,plot(10),(/3,3/),(/xBootHi_hf1_4xCO2,xBootLow_hf1_4xCO2/), polyres)

; plot seasonal cycle 
 aumm0001    = gsn_add_polyline(wks,plot(11),(/1,1/),(/xBootLow_mf1_PD_ac,   xBootHi_mf1_PD_ac/), polyres)
 aumm00001   = gsn_add_polyline(wks,plot(11),(/2,2/),(/xBootLow_mf1_2xCO2_ac,xBootHi_mf1_2xCO2_ac/), polyres)
 aumm000001   = gsn_add_polyline(wks,plot(11),(/3,3/),(/xBootLow_mf1_4xCO2_ac,xBootHi_mf1_4xCO2_ac/), polyres)
;
 bumm0001    = gsn_add_polyline(wks,plot(12),(/1,1/),(/-1*xBootHi_mf2_PD_ac,   -1*xBootLow_mf2_PD_ac/), polyres)
 bumm00001   = gsn_add_polyline(wks,plot(12),(/2,2/),(/-1*xBootHi_mf2_2xCO2_ac,-1*xBootLow_mf2_2xCO2_ac/), polyres)
 bumm000001  = gsn_add_polyline(wks,plot(12),(/3,3/),(/-1*xBootHi_mf2_4xCO2_ac,-1*xBootLow_mf2_4xCO2_ac/), polyres)

 cumm0001    = gsn_add_polyline(wks,plot(13),(/1,1/),(/xBootHi_hf1_PD_ac,    xBootLow_hf1_PD_ac/), polyres)
 cumm00001   = gsn_add_polyline(wks,plot(13),(/2,2/),(/xBootHi_hf1_2xCO2_ac, xBootLow_hf1_2xCO2_ac/), polyres)
 cumm000001  = gsn_add_polyline(wks,plot(13),(/3,3/),(/xBootHi_hf1_4xCO2_ac, xBootLow_hf1_4xCO2_ac/), polyres)
;
 wid = 0.24/2
 ;-------------------------------------------------
 cum1   = gsn_add_polyline(wks,plot(8),(/1-wid,1+wid/),(/xBootLow_mf1_PD,xBootLow_mf1_PD/), polyres)
 cum01  = gsn_add_polyline(wks,plot(8),(/1-wid,1+wid/),(/xBootHi_mf1_PD,xBootHi_mf1_PD/), polyres)

 cum2   = gsn_add_polyline(wks,plot(8),(/2-wid,2+wid/),(/xBootLow_mf1_2xCO2,xBootLow_mf1_2xCO2/), polyres)
 cum02  = gsn_add_polyline(wks,plot(8),(/2-wid,2+wid/),(/xBootHi_mf1_2xCO2,xBootHi_mf1_2xCO2/), polyres)

 cum3   = gsn_add_polyline(wks,plot(8),(/3-wid,3+wid/),(/xBootLow_mf1_4xCO2,xBootLow_mf1_4xCO2/), polyres)
 cum03  = gsn_add_polyline(wks,plot(8),(/3-wid,3+wid/),(/xBootHi_mf1_4xCO2,xBootHi_mf1_4xCO2/), polyres)

 dum1   = gsn_add_polyline(wks,plot(9),(/1-wid,1+wid/),(/-1*xBootLow_mf2_PD,-1*xBootLow_mf2_PD/), polyres)
 dum01  = gsn_add_polyline(wks,plot(9),(/1-wid,1+wid/),(/-1*xBootHi_mf2_PD, -1*xBootHi_mf2_PD/), polyres)

 dum2   = gsn_add_polyline(wks,plot(9),(/2-wid,2+wid/),(/-1*xBootLow_mf2_2xCO2,-1*xBootLow_mf2_2xCO2/), polyres)
 dum02  = gsn_add_polyline(wks,plot(9),(/2-wid,2+wid/),(/-1*xBootHi_mf2_2xCO2, -1*xBootHi_mf2_2xCO2/), polyres)

 dum3   = gsn_add_polyline(wks,plot(9),(/3-wid,3+wid/),(/-1*xBootLow_mf2_4xCO2,-1*xBootLow_mf2_4xCO2/), polyres)
 dum03  = gsn_add_polyline(wks,plot(9),(/3-wid,3+wid/),(/-1*xBootHi_mf2_4xCO2, -1*xBootHi_mf2_4xCO2/), polyres)

 fum1   = gsn_add_polyline(wks,plot(10),(/1-wid,1+wid/),(/xBootLow_hf1_PD,xBootLow_hf1_PD/), polyres)
 fum01  = gsn_add_polyline(wks,plot(10),(/1-wid,1+wid/),(/xBootHi_hf1_PD, xBootHi_hf1_PD/), polyres)

 fum2   = gsn_add_polyline(wks,plot(10),(/2-wid,2+wid/),(/xBootLow_hf1_2xCO2,xBootLow_hf1_2xCO2/), polyres)
 fum02  = gsn_add_polyline(wks,plot(10),(/2-wid,2+wid/),(/xBootHi_hf1_2xCO2, xBootHi_hf1_2xCO2/), polyres)

 fum3   = gsn_add_polyline(wks,plot(10),(/3-wid,3+wid/),(/xBootLow_hf1_4xCO2,xBootLow_hf1_4xCO2/), polyres)
 fum03  = gsn_add_polyline(wks,plot(10),(/3-wid,3+wid/),(/xBootHi_hf1_4xCO2, xBootHi_hf1_4xCO2/), polyres)

;-------------------------------------------------
; plot seasonal cycle
aum001   = gsn_add_polyline(wks,plot(11),(/1-wid,1+wid/),(/xBootLow_mf1_PD_ac,xBootLow_mf1_PD_ac/), polyres)
aum0001  = gsn_add_polyline(wks,plot(11),(/1-wid,1+wid/),(/xBootHi_mf1_PD_ac,xBootHi_mf1_PD_ac/), polyres)

bum001   = gsn_add_polyline(wks,plot(11),(/2-wid,2+wid/),(/xBootLow_mf1_2xCO2_ac,xBootLow_mf1_2xCO2_ac/), polyres)
bum0001  = gsn_add_polyline(wks,plot(11),(/2-wid,2+wid/),(/xBootHi_mf1_2xCO2_ac,xBootHi_mf1_2xCO2_ac/), polyres)

cum001   = gsn_add_polyline(wks,plot(11),(/3-wid,3+wid/),(/xBootLow_mf1_4xCO2_ac,xBootLow_mf1_4xCO2_ac/), polyres)
cum0001  = gsn_add_polyline(wks,plot(11),(/3-wid,3+wid/),(/xBootHi_mf1_4xCO2_ac,xBootHi_mf1_4xCO2_ac/), polyres)

dum001   = gsn_add_polyline(wks,plot(12),(/1-wid,1+wid/),(/-1*xBootLow_mf2_PD_ac,-1*xBootLow_mf2_PD_ac/), polyres)
dum0001  = gsn_add_polyline(wks,plot(12),(/1-wid,1+wid/),(/-1*xBootHi_mf2_PD_ac, -1*xBootHi_mf2_PD_ac/), polyres)

eum001   = gsn_add_polyline(wks,plot(12),(/2-wid,2+wid/),(/-1*xBootLow_mf2_2xCO2_ac,-1*xBootLow_mf2_2xCO2_ac/), polyres)
eum0001  = gsn_add_polyline(wks,plot(12),(/2-wid,2+wid/),(/-1*xBootHi_mf2_2xCO2_ac, -1*xBootHi_mf2_2xCO2_ac/), polyres)

fum001   = gsn_add_polyline(wks,plot(12),(/3-wid,3+wid/),(/-1*xBootLow_mf2_4xCO2_ac,-1*xBootLow_mf2_4xCO2_ac/), polyres)
fum0001  = gsn_add_polyline(wks,plot(12),(/3-wid,3+wid/),(/-1*xBootHi_mf2_4xCO2_ac, -1*xBootHi_mf2_4xCO2_ac/), polyres)

gum001   = gsn_add_polyline(wks,plot(13),(/1-wid,1+wid/),(/xBootLow_hf1_PD_ac,xBootLow_hf1_PD_ac/), polyres)
gum0001  = gsn_add_polyline(wks,plot(13),(/1-wid,1+wid/),(/xBootHi_hf1_PD_ac, xBootHi_hf1_PD_ac/), polyres)

hum001   = gsn_add_polyline(wks,plot(13),(/2-wid,2+wid/),(/xBootLow_hf1_2xCO2_ac,xBootLow_hf1_2xCO2_ac/), polyres)
hum0001  = gsn_add_polyline(wks,plot(13),(/2-wid,2+wid/),(/xBootHi_hf1_2xCO2_ac, xBootHi_hf1_2xCO2_ac/), polyres)

ium001   = gsn_add_polyline(wks,plot(13),(/3-wid,3+wid/),(/xBootLow_hf1_4xCO2_ac,xBootLow_hf1_4xCO2_ac/), polyres)
ium0001  = gsn_add_polyline(wks,plot(13),(/3-wid,3+wid/),(/xBootHi_hf1_4xCO2_ac, xBootHi_hf1_4xCO2_ac/), polyres)


 ;----------------------------------------------------------------------------------------------------
;------------------------------------------------------------------------------------------------------

  draw(plot(8))
  draw(plot(9))
  draw(plot(10))
  draw(plot(11))
  draw(plot(12))
  draw(plot(13))
end 
