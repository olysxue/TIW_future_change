; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
; line 65: check units of EKE whether cm2/s2 in file
; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
undef("readdata1")
function readdata1(ncdir, filelist, varname, scale)
local dims,temp,f
begin
   f= addfile(ncdir + filelist(0), "r")
   temp=f->$varname$(0,:,:)
   dims=getfilevardimsizes(f, varname)
   out=new((/5,dims(1),dims(2)/), float)
   
   do i=0,2
   f:= addfile(ncdir + filelist(i), "r")
   out(i,:,:)=f->$varname$(0,:,:)
   end do 
   out(3,:,:)=out(1,:,:)-out(0,:,:) ; 2xCO2 - PD
   out(4,:,:)=out(2,:,:)-out(0,:,:) ; 4xCO2 - PD
   out=out/scale

   do i=0,4
   out(i,:,:)  = smth9_Wrap(out(i,:,:), 0.50, 0.25, False)
   end do
  
   copy_VarMeta(temp,out(0,:,:))

   return(out)
end

undef("readdata")
function readdata(ncdir, filelist, varname, scale)
local dims,temp,f
begin
   f= addfile(ncdir + filelist(0), "r")
   temp=f->$varname$(:,:,:)
   dims=getfilevardimsizes(f, varname)
   out=new((/3,dims(0),dims(1),dims(2)/), float)
   
   do i=0,2
   f:= addfile(ncdir + filelist(i), "r")
   out(i,:,:,:)=f->$varname$(:,:,:)
   end do 
   out=out/scale
  
   copy_VarMeta(temp,out(0,:,:,:))
   return(out)
end

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
varname="tiw_std"
ncdir  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/mean/"
scale=100. ; units cm to m by scale
filelist=(/"CESM2_HR_TIW_rms_var_va_10-50_fft_surface_PD_mean.nc","CESM2_HR_TIW_rms_var_va_10-50_fft_surface_2xCO2_mean.nc","CESM2_HR_TIW_rms_var_va_10-50_fft_surface_4xCO2_mean.nc"/)
va0=readdata1(ncdir, filelist, varname, scale)
va0@units="m/s"

scale=100. ; units cm to m by scale
filelist=(/"CESM2_HR_TIW_rms_var_va_10-50_fft_surface_PD.nc","CESM2_HR_TIW_rms_var_va_10-50_fft_surface_2xCO2.nc","CESM2_HR_TIW_rms_var_va_10-50_fft_surface_4xCO2.nc"/)
va=readdata(ncdir, filelist, varname, scale)
va@units="meter"

va1_PD    = wgt_areaave_Wrap(va(0,:,{3:6},{210:250}), 1.0, 1.0, 0)
va1_2xCO2 = wgt_areaave_Wrap(va(1,:,{3:6},{210:250}), 1.0, 1.0, 0)
va1_4xCO2 = wgt_areaave_Wrap(va(2,:,{3:6},{210:250}), 1.0, 1.0, 0)

va2_PD    = wgt_areaave_Wrap(va(0,:,{-2:1},{220:260}), 1.0, 1.0, 0)
va2_2xCO2 = wgt_areaave_Wrap(va(1,:,{-2:1},{220:250}), 1.0, 1.0, 0)
va2_4xCO2 = wgt_areaave_Wrap(va(2,:,{-2:1},{210:250}), 1.0, 1.0, 0)


stat  = 0                         ; 0=mean; 1=variance; 2=std. dev; ...
n     = 30                        ; sub-sample size;  default is n=N
nBoot = 10000                      ; bootstrap replications
N     = 30
opt   = False                     ; use defaults
if (isvar("n") .and. (n.lt.N)) then
    opt   = True
    opt@sample_size = n           ; default N
end if

 BootStrap_va1_PD     = bootstrap_stat(va1_PD, stat, nBoot, 0, opt)
 BootStrap_va1_2xCO2  = bootstrap_stat(va1_2xCO2, stat, nBoot, 0, opt)
 BootStrap_va1_4xCO2  = bootstrap_stat(va1_4xCO2, stat, nBoot, 0, opt)

 BootStrap_va2_PD     = bootstrap_stat(va2_PD, stat, nBoot, 0, opt)
 BootStrap_va2_2xCO2  = bootstrap_stat(va2_2xCO2, stat, nBoot, 0, opt)
 BootStrap_va2_4xCO2  = bootstrap_stat(va2_4xCO2, stat, nBoot, 0, opt)
;==================================================
  xBoot_va1_PD       = BootStrap_va1_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
  xBoot_vaavg1_PD    = BootStrap_va1_PD[1]       ; Mean of bootstrapped estimates
  xBoot_vastd1_PD    = BootStrap_va1_PD[2]       ; Std dev of bootstrapped estimates
  
  xBoot_va2_PD       = BootStrap_va2_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
  xBoot_vaavg2_PD    = BootStrap_va2_PD[1]       ; Mean of bootstrapped estimates
  xBoot_vastd2_PD    = BootStrap_va2_PD[2]       ; Std dev of bootstrapped estimates
;=====================================================
    xBoot_va1_2xCO2    = BootStrap_va1_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_vaavg1_2xCO2    = BootStrap_va1_2xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_vastd1_2xCO2    = BootStrap_va1_2xCO2[2]       ; Std dev of bootstrapped estimates

    xBoot_va2_2xCO2    = BootStrap_va2_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_vaavg2_2xCO2    = BootStrap_va2_2xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_vastd2_2xCO2    = BootStrap_va2_2xCO2[2]       ; Std dev of bootstrapped estimates
;-----------------------------------------------------
    xBoot_va1_4xCO2    = BootStrap_va1_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_vaavg1_4xCO2    = BootStrap_va1_4xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_vastd1_4xCO2    = BootStrap_va1_4xCO2[2]       ; Std dev of bootstrapped estimates

    xBoot_va2_4xCO2    = BootStrap_va2_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_vaavg2_4xCO2    = BootStrap_va2_4xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_vastd2_4xCO2    = BootStrap_va2_4xCO2[2]       ; Std dev of bootstrapped estimates

;=====================================================
 xBootLow_va1_PD    = bootstrap_estimate(xBoot_va1_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_va1_PD    = bootstrap_estimate(xBoot_va1_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_va1_PD    = bootstrap_estimate(xBoot_va1_PD, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_va2_PD    = bootstrap_estimate(xBoot_va2_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_va2_PD    = bootstrap_estimate(xBoot_va2_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_va2_PD    = bootstrap_estimate(xBoot_va2_PD, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------
 xBootLow_va1_2xCO2    = bootstrap_estimate(xBoot_va1_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_va1_2xCO2    = bootstrap_estimate(xBoot_va1_2xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_va1_2xCO2    = bootstrap_estimate(xBoot_va1_2xCO2, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_va2_2xCO2    = bootstrap_estimate(xBoot_va2_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_va2_2xCO2    = bootstrap_estimate(xBoot_va2_2xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_va2_2xCO2    = bootstrap_estimate(xBoot_va2_2xCO2, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------
 xBootLow_va1_4xCO2    = bootstrap_estimate(xBoot_va1_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_va1_4xCO2    = bootstrap_estimate(xBoot_va1_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_va1_4xCO2    = bootstrap_estimate(xBoot_va1_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_va2_4xCO2    = bootstrap_estimate(xBoot_va2_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_va2_4xCO2    = bootstrap_estimate(xBoot_va2_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_va2_4xCO2    = bootstrap_estimate(xBoot_va2_4xCO2, 0.975, False)   ; 97.5% upper confidence bound
;---------------------------------------------------------------------------------------------------------
stdarr1 = new((/3/),"float") ; SSH
stdarr1(0) =  xBoot_vaavg1_PD
stdarr1(1) =  xBoot_vaavg1_2xCO2 
stdarr1(2) =  xBoot_vaavg1_4xCO2

stdarr2 = new((/3/),"float") ; SSH
stdarr2(0) =  xBoot_vaavg2_PD
stdarr2(1) =  xBoot_vaavg2_2xCO2 
stdarr2(2) =  xBoot_vaavg2_4xCO2

outdir = "/Users/xueaoyun/Documents/work/TIW_future/revise_figures/"
wks = gsn_open_wks("pdf",outdir+"plot_CESM_PD_2xCO2_4xCO2_TIW_VVEL_mean_intensity_bar")
wks_type = "pdf"

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
res@tmXBLabelFontHeightF  = 0.010
res@tmYLLabelFontHeightF  = 0.010
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
; res@lbOrientation        = "Vertical"
res@lbBoxEndCapStyle     = "TriangleBothEnds"
res@pmLabelBarWidthF     = 0.35
res@pmLabelBarHeightF    = 0.05
res@pmLabelBarOrthogonalPosF = 0.19
 res@pmLabelBarParallelPosF   = 0.62
res@lbLabelFontHeightF       =0.008
res@gsnLeftString        = ""
res@gsnCenterStringFontHeightF  = 0.014

res@gsnSpreadColorEnd  =230

res@vpHeightF            = 0.15
res@vpWidthF             = 0.25

res@lbTitlePosition  = "Right"
res@lbTitleDirection = "Across" 
res@lbTitleFontHeightF= 0.012
res@lbTitleOffsetF    = 0.01

plot=new(5,graphic)

panels = (/"~F22~a", "~F22~b", "~F22~c", "~F22~d", "~F22~e"/)
text=(/"PD","2xCO2","4xCO2"/)
yspace=0.05
yf=0.9
res@tiDeltaF=3
res@gsnCenterStringOrthogonalPosF=0.13

; --------------> plot va---------------------

res@lbTitleString    = "m/s"
res@cnLevels         = fspan(0.10,0.20,21)
res@vpYF    = yf
do i=0,2
res@vpXF    = 0.08 + i*res@vpWidthF + i*yspace
if i.eq.0 .or. i.eq.2 then
res@lbLabelBarOn    = False
else
res@lbLabelBarOn    = True
end if
res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W",""/)  
res@gsnLeftString  = panels(i)
res@gsnLeftStringFontHeightF = 0.015
res@gsnLeftStringOrthogonalPosF = 0.05
res@gsnCenterString   = text(i)
plot(i)  = gsn_csm_contour_map(wks,va0(i,:,:),res)
end do

plres=True
plres@gsLineThicknessF  = 2.5
; plres@gsFillColor     = -1
plres@gsLineDashPattern= 2
plres@gsLineColor="blue"
plres@gsLineOpacityF=100
w1=360-160
e1=360-100

w2=360-160
e2=360-110

s1=0
n1=5

s2=2
n2=6
dum = new(3, graphic)

dum(0)=gsn_add_polyline(wks, plot(0), (/w2,e2,e2,w2,w1/), (/s2,s2,n2,n2,s2/), plres)
dum(1)=gsn_add_polyline(wks, plot(1), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
dum(2)=gsn_add_polyline(wks, plot(2), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)

plres=True
plres@gsLineThicknessF=2.5
; plres@gsFillColor     = -1
plres@gsLineDashPattern= 2
plres@gsLineColor="red"
plres@gsLineOpacityF=100
w1=360-140
e1=360-90

w2=360-140
e2=360-100

w3=360-150
e3=360-110

s1=-2
n1=1

dum1=new(3, graphic)

dum1(0)=gsn_add_polyline(wks, plot(0), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum1(1)=gsn_add_polyline(wks, plot(1), (/w2,e2,e2,w2,w2/), (/s1,s1,n1,n1,s1/), plres)
dum1(2) =gsn_add_polyline(wks,plot(2), (/w3,e3,e3,w3,w3/), (/s1,s1,n1,n1,s1/), plres)

 draw(plot(0))
 draw(plot(1))
 draw(plot(2))
; plot bar 
 xyres                       = True                 ; plot mods desired
 xyres@gsnFrame              = False                ; don't advance frame yet
 xyres@gsnDraw               = False
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.2                 ; bring bars down to zero
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
 xyres@tmYLValues        = (/0,0.1,0.2,0.3,0.4,0.5/)
 xyres@tmYLLabels        = (/"0","0.1","0.2","0.3","0.4","0.5"/)
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

 xyres@vpXF    = 0.35
 yspace=0.05
 yf1=0.65
 xyres@vpYF    = yf1
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.2                 ; bring bars down to zero
                       ; turn off right axis
 xyres@gsnLeftString  = panels(3)
 xyres@gsnLeftStringFontHeightF = 0.015
 xyres@gsnXYBarChartColors = (/"cyan3"/)
 xyres@tmYLValues        = (/0,0.05,0.10,0.15,0.20/)
 xyres@tmYLLabels        = (/"0","0.05","0.10","0.15","0.20"/)
 xyres@gsnCenterString       = "   VVEL_NEQ"
 xyres@gsnCenterStringParallelPosF = 0.45
 xyres@gsnCenterStringOrthogonalPosF  = -0.10
 plot(3)  = gsn_csm_xy(wks,ispan(1,3,1),stdarr1(:),xyres)   ; advancing the frame
 xyres@gsnLeftString  = panels(4)
 xyres@tmYLOn        = True  
 xyres@tmYRBorderOn  = False                        ; turn off right axis
 xyres@tmYLBorderOn  = True 
 xyres@vpXF = xyres@vpWidthF+xyres@vpXF+yspace 
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.2                 ; bring bars down to zero
 xyres@gsnXYBarChartColors = (/"orange3"/)
 xyres@tmYLValues        = (/0,0.05,0.10,0.15,0.20/)
 xyres@tmYLLabels        = (/"0","0.05","0.10","0.15","0.20"/)
 xyres@gsnCenterString       = "    VVEL_EQ"
 xyres@gsnCenterStringParallelPosF = 0.45
 plot(4)  = gsn_csm_xy(wks,ispan(1,3,1),stdarr2(:),xyres)   ; advancing the frame

  polyres=True
  polyres@gsLineDashPattern=0
  polyres@gsLineColor="black"
  polyres@gsLineThicknessF=1.0
  polyres@tfPolyDrawOrder = "postDraw"


  aumm1    = gsn_add_polyline(wks,plot(3),(/1,1/),(/xBootLow_va1_PD,xBootHi_va1_PD/), polyres)
  aumm01   = gsn_add_polyline(wks,plot(3),(/2,2/),(/xBootLow_va1_2xCO2,xBootHi_va1_2xCO2/), polyres)
  aumm001  = gsn_add_polyline(wks,plot(3),(/3,3/),(/xBootLow_va1_4xCO2,xBootHi_va1_4xCO2/), polyres)

  bumm1    = gsn_add_polyline(wks,plot(4),(/1,1/),(/xBootLow_va2_PD,xBootHi_va2_PD/), polyres)
  bumm01   = gsn_add_polyline(wks,plot(4),(/2,2/),(/xBootLow_va2_2xCO2,xBootHi_va2_2xCO2/), polyres)
  bumm001  = gsn_add_polyline(wks,plot(4),(/3,3/),(/xBootLow_va2_4xCO2,xBootHi_va2_4xCO2/), polyres)

  wid = 0.24/2
  ;----------------------------------------------------------------------------------------------------
  cumm1   = gsn_add_polyline(wks,plot(3),(/1-wid,1+wid/),(/xBootLow_va1_PD,xBootLow_va1_PD/), polyres)
  cumm01  = gsn_add_polyline(wks,plot(3),(/1-wid,1+wid/),(/xBootHi_va1_PD,xBootHi_va1_PD/), polyres)

  cumm2   = gsn_add_polyline(wks,plot(3),(/2-wid,2+wid/),(/xBootLow_va1_2xCO2,xBootLow_va1_2xCO2/), polyres)
  cumm02  = gsn_add_polyline(wks,plot(3),(/2-wid,2+wid/),(/xBootHi_va1_2xCO2,xBootHi_va1_2xCO2/), polyres)

  cumm3   = gsn_add_polyline(wks,plot(3),(/3-wid,3+wid/),(/xBootLow_va1_4xCO2,xBootLow_va1_4xCO2/), polyres)
  cumm03  = gsn_add_polyline(wks,plot(3),(/3-wid,3+wid/),(/xBootHi_va1_4xCO2,xBootHi_va1_4xCO2/), polyres)
  ;----------------------------------------------------------------------------------------------------
  dumm1   = gsn_add_polyline(wks,plot(4),(/1-wid,1+wid/),(/xBootLow_va2_PD,xBootLow_va2_PD/), polyres)
  dumm01  = gsn_add_polyline(wks,plot(4),(/1-wid,1+wid/),(/xBootHi_va2_PD,xBootHi_va2_PD/), polyres)

  dumm2   = gsn_add_polyline(wks,plot(4),(/2-wid,2+wid/),(/xBootLow_va2_2xCO2,xBootLow_va2_2xCO2/), polyres)
  dumm02  = gsn_add_polyline(wks,plot(4),(/2-wid,2+wid/),(/xBootHi_va2_2xCO2,xBootHi_va2_2xCO2/), polyres)

  dumm3   = gsn_add_polyline(wks,plot(4),(/3-wid,3+wid/),(/xBootLow_va2_4xCO2,xBootLow_va2_4xCO2/), polyres)
  dumm03  = gsn_add_polyline(wks,plot(4),(/3-wid,3+wid/),(/xBootHi_va2_4xCO2,xBootHi_va2_4xCO2/), polyres)
;-----------------------------------------------------------------------------------------------------

  draw(plot(0))
  draw(plot(1))
  draw(plot(2))
  draw(plot(3))
  draw(plot(4))
end 



