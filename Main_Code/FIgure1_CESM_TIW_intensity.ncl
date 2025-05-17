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

undef("readdata2")
function readdata2(ncdir, filelist, varname, scale)
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
ncdir  = "/Volumes/oly2/CESM/mean/"
scale=1. ;units degC
filelist=(/"CESM2_HR_SST_std_var_10-50_fft_PD_mean.nc","CESM2_HR_SST_std_var_10-50_fft_2xCO2_mean.nc","CESM2_HR_SST_std_var_10-50_fft_4xCO2_mean.nc"/)
sst=readdata1(ncdir, filelist, varname, scale)
sst@units="degC"

scale=100. ; units cm to m by scale
filelist=(/"CESM2_HR_SSH_std_var_10-50_fft_PD_mean.nc","CESM2_HR_SSH_std_var_10-50_fft_2xCO2_mean.nc","CESM2_HR_SSH_std_var_10-50_fft_4xCO2_mean.nc"/)
ssh=readdata1(ncdir, filelist, varname, scale)
ssh@units="meter"

scale=100. ; units cm to m by scale
filelist=(/"CESM2_HR_TIW_rms_var_va_10-50_fft_surface_PD_mean.nc","CESM2_HR_TIW_rms_var_va_10-50_fft_surface_2xCO2_mean.nc","CESM2_HR_TIW_rms_var_va_10-50_fft_surface_4xCO2_mean.nc"/)
va=readdata1(ncdir, filelist, varname, scale)
va@units="meter"
printMinMax(va, 0)
scale=10000.
varname="EKE" ;units cm2/s2 to m2/s2 by scale
filelist:=(/"CESM2_TIW_EKE_10-50_fft_surface_PD_mean.nc","CESM2_TIW_EKE_10-50_fft_surface_2xCO2_mean.nc","CESM2_TIW_EKE_10-50_fft_surface_4xCO2_mean.nc"/)
eke=readdata1(ncdir, filelist, varname, scale)
eke@units="J/m3"
printMinMax(eke(0,:,:), 0)
;=============================================================
; plot mean TIW EKE and RMS
;=============================================================
outdir = "/Users/xueaoyun/Documents/work/TIW_future/revision1/main_figure/"
wks = gsn_open_wks("pdf",outdir+"Figure1_CESM_PD_2xCO2_4xCO2_TIW_mean_intensity")
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
res@mpMinLatF            = -10
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
res@tmYLLabels = (/"","","","",""/)
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
; res@lbOrientation        = "Vertical"
res@lbBoxEndCapStyle     = "TriangleBothEnds"
res@pmLabelBarWidthF     = 0.20
; res@pmLabelBarHeightF    = 0.55
res@pmLabelBarOrthogonalPosF = 0.19
; res@pmLabelBarParallelPosF   = 0.6
res@lbLabelFontHeightF       =0.008
res@gsnLeftString        = ""
res@gsnCenterStringFontHeightF  = 0.014

res@gsnSpreadColorEnd  =230

res@vpHeightF            = 0.12
res@vpWidthF             = 0.20

res@lbTitlePosition  = "Bottom"
res@lbTitleDirection = "Across" 
res@lbTitleFontHeightF= 0.010
res@lbTitleOffsetF    = 0.05

plot=new(17,graphic)
dum=plot

; ----->1 plot sst
res@lbTitleString    = "degC"

;res@cnLevels             =(/0.15,0.2,0.25,0.3,0.35,0.4/)
res@cnLevels             = fspan(0.2,0.5,31)
; res@cnFillColors         = (/0,13,14,15,16,19,20/)
;res@cnFillColors         = (/0,138,168,178,188,204,230,246/)
; res@cnFillColors         = (/0,133,150,168,178,188,204,230,246/)
panels = (/"~F22~a", "~F22~d", "~F22~g", "~F22~j", "~F22~b", "~F22~e", "~F22~h", "~F22~k", "~F22~c","~F22~f","~F22~i","~F22~l", "~F22~m","~F22~n","","~F22~o","~F22~p"/)
text=(/"PD","2xCO2","4xCO2"/)
yspace=0.018
yf=0.9
res@tiDeltaF=3
res@gsnCenterStringOrthogonalPosF=0.13
;**********************************************************************************************************
; ----->1 plot sst
res@vpXF    = 0.08 + 2*res@vpWidthF + 2*yspace
do i=0,2
res@vpYF    = yf-i*yspace-i*res@vpHeightF
res@tiYAxisString=text(i)
res@tiYAxisOffsetXF=-0.02
res@tiYAxisFontHeightF=0.014
if i.eq.0 then
res@gsnCenterString = "SST"
res@tiMainString=""
else
res@gsnCenterString = ""
res@tiMainString=""
end if

if i.eq.2 then
res@lbLabelBarOn    = True
res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W",""/)  
else
res@lbLabelBarOn    = False
res@tmXBLabels = (/"","","","","",""/)
end if
res@gsnLeftString  = panels(2+i*4)
res@gsnLeftStringFontHeightF = 0.012
res@gsnLeftStringOrthogonalPosF = -0.13
res@gsnLeftStringParallelPosF  = 0.03
plot(i*4)  = gsn_csm_contour_map(wks,sst(i,:,:),res)
end do
delete([/res@tiYAxisString,res@tiYAxisOffsetXF,res@tiYAxisFontHeightF/])
;**********************************************************************************************************
; ----->2 plot ssh
res@vpXF    = 0.08 + 3*res@vpWidthF + 3*yspace
ssh0=ssh
ssh01=ssh*100
ssh=(/ssh01/)
copy_VarMeta(ssh0, ssh)
delete(res@cnLevels)
res@lbTitleString    = "10~S~-2~N~ m"
;res@cnLevels             :=(/0.15,0.2,0.25,0.3,0.35,0.4/)*10
res@cnLevels               := fspan(1.5,4.5,31)
res@cnFillPalette         = "WhiteBlueGreenYellowOrangeRed"
; res@cnFillColors         = (/0,12,13,14,15,17,19/)
res@tmYLLabels = (/"","","","",""/)
do i=0,2
res@vpYF    = yf-i*yspace-i*res@vpHeightF
if i.eq.0 then
res@gsnCenterString = "SSH"
else
res@gsnCenterString = ""
end if 

if i.eq.2 then
res@lbLabelBarOn    = True
res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W",""/)  
else
res@lbLabelBarOn    = False
res@tmXBLabels = (/"","","","","",""/)
end if
res@gsnLeftString  = panels(3+i*4)
res@gsnLeftStringFontHeightF = 0.012
res@gsnLeftStringOrthogonalPosF = -0.13
res@gsnLeftStringParallelPosF  = 0.03
plot(1+i*4)  = gsn_csm_contour_map(wks,ssh(i,:,:),res)
end do
;**********************************************************************************************************

; ----->3 plot va
res@vpXF    = 0.08 + res@vpWidthF + yspace
res@lbTitleString    = "m/s"
;res@cnLevels             :=(/15,20,25,30,35,40/)
delete(res@cnLevels)
res@cnLevels                = fspan(0.10,0.20, 21)
; res@cnFillColors         := (/0,12,13,14,15,17,19/)
res@tmYLLabels = (/"","","","",""/)
do i=0,2
res@vpYF    = yf-i*yspace-i*res@vpHeightF
if i.eq.0 then
res@gsnCenterString = "VVEL"
else
res@gsnCenterString = ""
end if
if i.eq.2 then
res@lbLabelBarOn    = True
res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W",""/)  
else
res@lbLabelBarOn    = False
res@tmXBLabels = (/"","","","","",""/)
end if
res@gsnLeftString  = panels(1+i*4)
res@gsnLeftStringFontHeightF = 0.012
res@gsnLeftStringOrthogonalPosF = -0.13
plot(2+i*4)  = gsn_csm_contour_map(wks,va(i,:,:),res)
end do
;**********************************************************************************************************

; ----->4 plot eke
res@vpXF    = 0.08
res@lbTitleString    = "J/m~S~3~N~"
;res@cnLevels             :=(/15,20,25,30,35,40/)
delete(res@cnLevels)
res@cnLevels                = fspan(15,45,31)
; res@cnFillColors         := (/0,12,13,14,15,17,19/)
res@tmYLLabels = (/-10,-5,0,5,10/)
res@tmYLLabels = (/"10~S~o~N~S","5~S~o~N~S","EQ","5~S~o~N~N","10~S~o~N~N"/)
do i=0,2
res@vpYF    = yf-i*yspace-i*res@vpHeightF
if i.eq.0 then
res@gsnCenterString = "EKE"
else
res@gsnCenterString = ""
end if
if i.eq.2 then
res@lbLabelBarOn    = True
res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W",""/)  
else
res@lbLabelBarOn    = False
res@tmXBLabels = (/"","","","","",""/)
end if
res@gsnLeftString  = panels(i*4)
res@gsnLeftStringFontHeightF = 0.012
res@gsnLeftStringOrthogonalPosF = -0.13
plot(3+i*4)  = gsn_csm_contour_map(wks,eke(i,:,:),res)
end do

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

s1=0
n1=5

s2 = 3
n2 = 8

s3 = 3
n3 = 6


dum(0)=gsn_add_polyline(wks, plot(0), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum(1)=gsn_add_polyline(wks, plot(1), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
dum(2)=gsn_add_polyline(wks, plot(2), (/w2,e2,e2,w2,w2/), (/s3,s3,n3,n3,s3/), plres)
dum(3)=gsn_add_polyline(wks, plot(3), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum(4)=gsn_add_polyline(wks, plot(4), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum(5)=gsn_add_polyline(wks, plot(5), (/w2,e2,e2,w2,w1/), (/s2,s2,n2,n2,s2/), plres)
dum(6)=gsn_add_polyline(wks, plot(6), (/w2,e2,e2,w2,w1/), (/s3,s3,n3,n3,s3/), plres)
dum(7)=gsn_add_polyline(wks, plot(7), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum(8)=gsn_add_polyline(wks, plot(8), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum(9) =gsn_add_polyline(wks,plot(9), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
dum(10)=gsn_add_polyline(wks, plot(10), (/w2,e2,e2,w2,w2/), (/s3,s3,n3,n3,s3/), plres)
dum(11)=gsn_add_polyline(wks, plot(11), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)

plres=True
plres@gsLineThicknessF=2.5
; plres@gsFillColor     = -1
plres@gsLineDashPattern= 2
plres@gsLineColor="black"
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

dum1(0)=gsn_add_polyline(wks, plot(2), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum1(1)=gsn_add_polyline(wks, plot(6), (/w2,e2,e2,w2,w2/), (/s1,s1,n1,n1,s1/), plres)
dum1(2) =gsn_add_polyline(wks,plot(10), (/w3,e3,e3,w3,w3/), (/s1,s1,n1,n1,s1/), plres)


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
pos=new((/3,2/), float)
do i=0,2
pos(i,:) = get_title_position(plot(i*4))
end do

;=============================================================
; read data
;=============================================================
varname="tiw_std"
ncdir  = "/Volumes/oly2/CESM/mean/"
scale=1. ;units degC
filelist=(/"CESM2_HR_SST_std_var_10-50_fft_PD.nc","CESM2_HR_SST_std_var_10-50_fft_2xCO2.nc","CESM2_HR_SST_std_var_10-50_fft_4xCO2.nc"/)
sst1=readdata2(ncdir, filelist, varname, scale)
sst1@units="degC"
printVarSummary(sst1)
scale=1. ; units cm to m by scale
filelist=(/"CESM2_HR_SSH_std_var_10-50_fft_PD.nc","CESM2_HR_SSH_std_var_10-50_fft_2xCO2.nc","CESM2_HR_SSH_std_var_10-50_fft_4xCO2.nc"/)
ssh1=readdata2(ncdir, filelist, varname, scale)
ssh1@units="meter"

scale=100. ; units cm to m by scale
filelist=(/"CESM2_HR_TIW_rms_var_va_10-50_fft_surface_PD.nc","CESM2_HR_TIW_rms_var_va_10-50_fft_surface_2xCO2.nc","CESM2_HR_TIW_rms_var_va_10-50_fft_surface_4xCO2.nc"/)
va1=readdata2(ncdir, filelist, varname, scale)
va1@units="meter"

scale=10000.
varname="EKE" ;units cm2/s2 to m2/s2 by scale
filelist:=(/"CESM2_TIW_EKE_10-50_fft_surface_PD.nc","CESM2_TIW_EKE_10-50_fft_surface_2xCO2.nc","CESM2_TIW_EKE_10-50_fft_surface_4xCO2.nc"/)
eke1=readdata2(ncdir, filelist, varname, scale)
eke1@units="J/m3"
sst_PD0 = dim_avg_n_Wrap(sst1(0,:,{0:5},{200:250}), 0)
sst_4xCO20 = dim_avg_n_Wrap(sst1(2,:,{0:5},{200:250}), 0)
print(max(sst_PD0))
print(max(sst_4xCO20))
sst_PD    = wgt_areaave_Wrap(sst1(0,:,{0:5},{200:250}), 1.0, 1.0, 0)
sst_2xCO2 = wgt_areaave_Wrap(sst1(1,:,{0:5},{200:250}), 1.0, 1.0, 0)
sst_4xCO2 = wgt_areaave_Wrap(sst1(2,:,{0:5},{200:250}), 1.0, 1.0, 0)

ssh_PD    = wgt_areaave_Wrap(ssh1(0,:,{3:8},{210:250}), 1.0, 1.0, 0)
ssh_2xCO2 = wgt_areaave_Wrap(ssh1(1,:,{3:8},{210:250}), 1.0, 1.0, 0)
ssh_4xCO2 = wgt_areaave_Wrap(ssh1(2,:,{3:8},{210:250}), 1.0, 1.0, 0)

va1_PD    = wgt_areaave_Wrap(va1(0,:,{3:6},{210:250}), 1.0, 1.0, 0)
va1_2xCO2 = wgt_areaave_Wrap(va1(1,:,{3:6},{210:250}), 1.0, 1.0, 0)
va1_4xCO2 = wgt_areaave_Wrap(va1(2,:,{3:6},{210:250}), 1.0, 1.0, 0)
va1_PD_avg    = dim_avg_n_Wrap(va1_PD, 0)
va1_2xCO2_avg = dim_avg_n_Wrap(va1_2xCO2, 0)
va1_4xCO2_avg = dim_avg_n_Wrap(va1_4xCO2, 0)

va2_PD    = wgt_areaave_Wrap(va1(0,:,{-2:1},{220:260}), 1.0, 1.0, 0)
va2_2xCO2 = wgt_areaave_Wrap(va1(1,:,{-2:1},{220:250}), 1.0, 1.0, 0)
va2_4xCO2 = wgt_areaave_Wrap(va1(2,:,{-2:1},{210:250}), 1.0, 1.0, 0)
va2_PD_avg    = dim_avg_n_Wrap(va2_PD, 0)
va2_2xCO2_avg = dim_avg_n_Wrap(va2_2xCO2, 0)
va2_4xCO2_avg = dim_avg_n_Wrap(va2_4xCO2, 0)
eke_PD    = wgt_areaave_Wrap(eke1(0,:,{0:5},{200:260}), 1.0, 1.0, 0)
eke_2xCO2 = wgt_areaave_Wrap(eke1(1,:,{0:5},{200:260}), 1.0, 1.0, 0)
eke_4xCO2 = wgt_areaave_Wrap(eke1(2,:,{0:5},{200:260}), 1.0, 1.0, 0)

eke_PD_avg = dim_avg_n_Wrap(eke_PD, 0)
eke_2xCO2_avg = dim_avg_n_Wrap(eke_2xCO2, 0)
eke_4xCO2_avg = dim_avg_n_Wrap(eke_4xCO2, 0)

stat  = 0                         ; 0=mean; 1=variance; 2=std. dev; ...
n     = 30                        ; sub-sample size;  default is n=N
nBoot = 10000                      ; bootstrap replications
N     = 30
opt   = False                     ; use defaults
if (isvar("n") .and. (n.lt.N)) then
    opt   = True
    opt@sample_size = n           ; default N
end if

 BootStrap_sst_PD     = bootstrap_stat(sst_PD, stat, nBoot, 0, opt)
 BootStrap_sst_2xCO2  = bootstrap_stat(sst_2xCO2, stat, nBoot, 0, opt)
 BootStrap_sst_4xCO2  = bootstrap_stat(sst_4xCO2, stat, nBoot, 0, opt)

 BootStrap_ssh_PD     = bootstrap_stat(ssh_PD, stat, nBoot, 0, opt)
 BootStrap_ssh_2xCO2  = bootstrap_stat(ssh_2xCO2, stat, nBoot, 0, opt)
 BootStrap_ssh_4xCO2  = bootstrap_stat(ssh_4xCO2, stat, nBoot, 0, opt)

 BootStrap_va1_PD     = bootstrap_stat(va1_PD, stat, nBoot, 0, opt)
 BootStrap_va1_2xCO2  = bootstrap_stat(va1_2xCO2, stat, nBoot, 0, opt)
 BootStrap_va1_4xCO2  = bootstrap_stat(va1_4xCO2, stat, nBoot, 0, opt)

 BootStrap_va2_PD     = bootstrap_stat(va2_PD, stat, nBoot, 0, opt)
 BootStrap_va2_2xCO2  = bootstrap_stat(va2_2xCO2, stat, nBoot, 0, opt)
 BootStrap_va2_4xCO2  = bootstrap_stat(va2_4xCO2, stat, nBoot, 0, opt)

 BootStrap_eke_PD     = bootstrap_stat(eke_PD, stat, nBoot, 0, opt)
 BootStrap_eke_2xCO2  = bootstrap_stat(eke_2xCO2, stat, nBoot, 0, opt)
 BootStrap_eke_4xCO2  = bootstrap_stat(eke_4xCO2, stat, nBoot, 0, opt)

 xBoot_sst_PD       = BootStrap_sst_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sstavg_PD    = BootStrap_sst_PD[1]       ; Mean of bootstrapped estimates
 xBoot_sststd_PD    = BootStrap_sst_PD[2]       ; Std dev of bootstrapped estimates

 xBoot_ssh_PD       = BootStrap_ssh_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sshavg_PD    = BootStrap_ssh_PD[1]       ; Mean of bootstrapped estimates
 xBoot_sshstd_PD    = BootStrap_ssh_PD[2]       ; Std dev of bootstrapped estimates

 xBoot_va1_PD       = BootStrap_va1_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_vaavg1_PD    = BootStrap_va1_PD[1]       ; Mean of bootstrapped estimates
 xBoot_vastd1_PD    = BootStrap_va1_PD[2]       ; Std dev of bootstrapped estimates

 xBoot_va2_PD       = BootStrap_va2_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_vaavg2_PD    = BootStrap_va2_PD[1]       ; Mean of bootstrapped estimates
 xBoot_vastd2_PD    = BootStrap_va2_PD[2]       ; Std dev of bootstrapped estimates

 xBoot_eke_PD       = BootStrap_eke_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_ekeavg_PD    = BootStrap_eke_PD[1]       ; Mean of bootstrapped estimates
 xBoot_ekestd_PD    = BootStrap_eke_PD[2]       ; Std dev of bootstrapped estimates
;=====================================================
    xBoot_sst_2xCO2    = BootStrap_sst_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sstavg_2xCO2    = BootStrap_sst_2xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_sststd_2xCO2    = BootStrap_sst_2xCO2[2]       ; Std dev of bootstrapped estimates

    xBoot_ssh_2xCO2    = BootStrap_ssh_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sshavg_2xCO2    = BootStrap_ssh_2xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_sshstd_2xCO2    = BootStrap_ssh_2xCO2[2]       ; Std dev of bootstrapped estimates

    xBoot_va1_2xCO2    = BootStrap_va1_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_vaavg1_2xCO2    = BootStrap_va1_2xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_vastd1_2xCO2    = BootStrap_va1_2xCO2[2]       ; Std dev of bootstrapped estimates

    xBoot_va2_2xCO2    = BootStrap_va2_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_vaavg2_2xCO2    = BootStrap_va2_2xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_vastd2_2xCO2    = BootStrap_va2_2xCO2[2]       ; Std dev of bootstrapped estimates

    xBoot_eke_2xCO2    = BootStrap_eke_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_ekeavg_2xCO2    = BootStrap_eke_2xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_ekestd_2xCO2    = BootStrap_eke_2xCO2[2]       ; Std dev of bootstrapped estimates
;-----------------------------------------------------
    xBoot_sst_4xCO2    = BootStrap_sst_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sstavg_4xCO2    = BootStrap_sst_4xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_sststd_4xCO2    = BootStrap_sst_4xCO2[2]       ; Std dev of bootstrapped estimates

    xBoot_ssh_4xCO2    = BootStrap_ssh_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sshavg_4xCO2    = BootStrap_ssh_4xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_sshstd_4xCO2    = BootStrap_ssh_4xCO2[2]       ; Std dev of bootstrapped estimates

    xBoot_va1_4xCO2    = BootStrap_va1_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_vaavg1_4xCO2    = BootStrap_va1_4xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_vastd1_4xCO2    = BootStrap_va1_4xCO2[2]       ; Std dev of bootstrapped estimates

    xBoot_va2_4xCO2    = BootStrap_va2_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_vaavg2_4xCO2    = BootStrap_va2_4xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_vastd2_4xCO2    = BootStrap_va2_4xCO2[2]       ; Std dev of bootstrapped estimates

    xBoot_eke_4xCO2    = BootStrap_eke_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_ekeavg_4xCO2    = BootStrap_eke_4xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_ekestd_4xCO2    = BootStrap_eke_4xCO2[2]       ; Std dev of bootstrapped estimates
;=====================================================
 xBootLow_sst_PD    = bootstrap_estimate(xBoot_sst_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_sst_PD    = bootstrap_estimate(xBoot_sst_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_sst_PD    = bootstrap_estimate(xBoot_sst_PD, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_ssh_PD    = bootstrap_estimate(xBoot_ssh_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_ssh_PD    = bootstrap_estimate(xBoot_ssh_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_ssh_PD    = bootstrap_estimate(xBoot_ssh_PD, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_va1_PD    = bootstrap_estimate(xBoot_va1_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_va1_PD    = bootstrap_estimate(xBoot_va1_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_va1_PD    = bootstrap_estimate(xBoot_va1_PD, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_va2_PD    = bootstrap_estimate(xBoot_va2_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_va2_PD    = bootstrap_estimate(xBoot_va2_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_va2_PD    = bootstrap_estimate(xBoot_va2_PD, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_eke_PD    = bootstrap_estimate(xBoot_eke_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_eke_PD    = bootstrap_estimate(xBoot_eke_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_eke_PD    = bootstrap_estimate(xBoot_eke_PD, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------
 xBootLow_sst_2xCO2    = bootstrap_estimate(xBoot_sst_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_sst_2xCO2    = bootstrap_estimate(xBoot_sst_2xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_sst_2xCO2    = bootstrap_estimate(xBoot_sst_2xCO2, 0.975, False)   ; 97.5% upper confidence bound
 
 xBootLow_ssh_2xCO2    = bootstrap_estimate(xBoot_ssh_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_ssh_2xCO2    = bootstrap_estimate(xBoot_ssh_2xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_ssh_2xCO2    = bootstrap_estimate(xBoot_ssh_2xCO2, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_va1_2xCO2    = bootstrap_estimate(xBoot_va1_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_va1_2xCO2    = bootstrap_estimate(xBoot_va1_2xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_va1_2xCO2    = bootstrap_estimate(xBoot_va1_2xCO2, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_va2_2xCO2    = bootstrap_estimate(xBoot_va2_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_va2_2xCO2    = bootstrap_estimate(xBoot_va2_2xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_va2_2xCO2    = bootstrap_estimate(xBoot_va2_2xCO2, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_eke_2xCO2    = bootstrap_estimate(xBoot_eke_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_eke_2xCO2    = bootstrap_estimate(xBoot_eke_2xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_eke_2xCO2    = bootstrap_estimate(xBoot_eke_2xCO2, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------
xBootLow_sst_4xCO2    = bootstrap_estimate(xBoot_sst_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_sst_4xCO2    = bootstrap_estimate(xBoot_sst_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_sst_4xCO2    = bootstrap_estimate(xBoot_sst_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_ssh_4xCO2    = bootstrap_estimate(xBoot_ssh_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_ssh_4xCO2    = bootstrap_estimate(xBoot_ssh_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_ssh_4xCO2    = bootstrap_estimate(xBoot_ssh_4xCO2, 0.975, False)   ; 97.5% upper confidence bound


 xBootLow_va1_4xCO2    = bootstrap_estimate(xBoot_va1_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_va1_4xCO2    = bootstrap_estimate(xBoot_va1_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_va1_4xCO2    = bootstrap_estimate(xBoot_va1_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_va2_4xCO2    = bootstrap_estimate(xBoot_va2_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_va2_4xCO2    = bootstrap_estimate(xBoot_va2_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_va2_4xCO2    = bootstrap_estimate(xBoot_va2_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_eke_4xCO2    = bootstrap_estimate(xBoot_eke_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_eke_4xCO2    = bootstrap_estimate(xBoot_eke_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_eke_4xCO2    = bootstrap_estimate(xBoot_eke_4xCO2, 0.975, False)   ; 97.5% upper confidence bound
;---------------------------------------------------------------------------------------------------------
stdarr1 = new((/3/),"float") ;  SST
stdarr1(0) =  xBoot_sstavg_PD
stdarr1(1) =  xBoot_sstavg_2xCO2 
stdarr1(2) =  xBoot_sstavg_4xCO2

stdarr2 = new((/3/),"float") ; SSH
stdarr2(0) =  xBoot_sshavg_PD
stdarr2(1) =  xBoot_sshavg_2xCO2 
stdarr2(2) =  xBoot_sshavg_4xCO2

stdarr3 = new((/3/),"float") ; SSH
stdarr3(0) = va1_PD_avg
stdarr3(1) = va1_2xCO2_avg
stdarr3(2) = va1_4xCO2_avg

stdarr4 = new((/3/),"float") ; SSH
stdarr4(0) =  va2_PD_avg
stdarr4(1) =  va2_2xCO2_avg
stdarr4(2) =  va2_4xCO2_avg

stdarr5 = new((/3/),"float") ; EKE 
stdarr5(0) =  eke_PD_avg
stdarr5(1) =  eke_2xCO2_avg
stdarr5(2) =  eke_4xCO2_avg 

 xyres                       = True                 ; plot mods desired
 xyres@gsnFrame              = False                ; don't advance frame yet
 xyres@gsnDraw               = False
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.4                 ; bring bars down to zero
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 4.                   ; of the 1st and last bars
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
;****************************************************
 xspace = 0.02
 yspace=0.05
 yf=0.9
 xyres@vpXF    = 0.08 + 2*res@vpWidthF + 2*0.035
 xyres@vpYF    = yf-3*yspace-3*res@vpHeightF
 xyres@gsnXYBarChartColors = (/"cyan3"/)
 xyres@gsnCenterString       = "    SST"
 xyres@gsnCenterStringOrthogonalPosF = -0.10
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString  = panels(15)
 xyres@gsnLeftStringFontHeightF = 0.012
 xyres@gsnLeftStringOrthogonalPosF = 0.
 xyres@gsnLeftStringParallelPosF = -0.05
 plot(12)  = gsn_csm_xy(wks,ispan(1,3,1),stdarr1(:),xyres)   ; advancing the frame
 ;****************************************************
 xyres@vpXF = 0.08 + 3*res@vpWidthF + 3*0.032
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 4                 ; bring bars down to zero
 xyres@gsnXYBarChartColors = (/"pink"/)
 xyres@tmYLValues        = (/0,1,2,3,4,5/)
 xyres@tmYLLabels        = (/"0","1","2","3","4","5"/)
 xyres@gsnCenterString       = "    SSH"
 xyres@gsnCenterStringParallelPosF = 0.45
 xyres@gsnLeftString  = panels(16)
 xyres@gsnLeftStringFontHeightF = 0.012
 xyres@gsnLeftStringOrthogonalPosF = 0.
 xyres@gsnLeftStringParallelPosF = -0.05
 plot(13)  = gsn_csm_xy(wks,ispan(1,3,1),stdarr2(:),xyres)   ; advancing the frame
  ;****************************************************
 xyres@vpXF = 0.08 + res@vpWidthF+0.02 
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@vpWidthF              = 0.10
 xyres@vpHeightF             = 0.15
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.2                 ; bring bars down to zero
                       ; turn off right axis
 xyres@gsnXYBarChartColors = (/"chartreuse2"/)
 xyres@tmYLValues        = (/0,0.05,0.10,0.15,0.20/)
 xyres@tmYLLabels        = (/"0","0.05","0.10","0.15","0.20"/)
 xyres@gsnCenterString       = "   VVEL_NEQ"
 xyres@gsnCenterStringOrthogonalPosF = -0.12
 xyres@gsnCenterStringParallelPosF = 0.40
 xyres@trXMinF               = 0.5                    ; adds space on either end
 xyres@trXMaxF               = 3.75                    ; of the 1st and last bars
 xyres@gsnLeftString  = panels(13)
 xyres@gsnLeftStringFontHeightF = 0.012
 xyres@gsnLeftStringOrthogonalPosF = 0.
 plot(14)  = gsn_csm_xy(wks,ispan(1,3,1),stdarr3(:),xyres)   ; advancing the frame
 xyres@tmYLOn        = False 
 xyres@tmYRBorderOn  = True                        ; turn off right axis
 xyres@tmYLBorderOn  = False
 xyres@vpXF = xyres@vpWidthF+xyres@vpXF
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.2                 ; bring bars down to zero
 xyres@trXMinF               = 0.25                    ; adds space on either end
 xyres@trXMaxF               = 3.5 
 xyres@gsnXYBarChartColors = (/"chartreuse2"/)
 xyres@tmYLValues        = (/0,0.05,0.10,0.15,0.20/)
 xyres@tmYLLabels        = (/"0","0.05","0.10","0.15","0.20"/)
 xyres@gsnCenterString       = "    VVEL_EQ"
 xyres@gsnCenterStringOrthogonalPosF = -0.10
 xyres@gsnCenterStringParallelPosF = 0.45
 xyres@gsnLeftString  = panels(14)
 xyres@gsnLeftStringFontHeightF = 0.012
 xyres@gsnLeftStringOrthogonalPosF = 0.
 plot(15)  = gsn_csm_xy(wks,ispan(1,3,1),stdarr4(:),xyres)   ; advancing the frame
  ;****************************************************
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@vpWidthF              = 0.15
 xyres@vpHeightF             = 0.15
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 4.                   ; of the 1st and last bars
 xyres@tmYLOn        = True
 xyres@tmYRBorderOn  = False                        ; turn off right axis
 xyres@tmYLBorderOn  = True 
 xyres@vpXF    = 0.08+0.03
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 40                ; bring bars down to zero
 xyres@gsnXYBarChartColors = (/"orange3"/)
 xyres@tmYLValues        = (/0,10,20,30,40,50/)
 xyres@tmYLLabels        = (/"0","10","20","30","40","50"/)
 xyres@gsnCenterString       = "    EKE"
 xyres@gsnCenterStringParallelPosF = 0.45
 xyres@gsnLeftString  = panels(12)
 xyres@gsnLeftStringFontHeightF = 0.012
 xyres@gsnLeftStringOrthogonalPosF = 0
 plot(16)  = gsn_csm_xy(wks,ispan(1,3,1),stdarr5(:),xyres)   ; advancing the frame
;=========================================================
  polyres=True
  polyres@gsLineDashPattern=0
  polyres@gsLineColor="black"
  polyres@gsLineThicknessF=1.0
  polyres@tfPolyDrawOrder = "postDraw"

  aumm1    = gsn_add_polyline(wks,plot(12),(/1,1/),(/xBootLow_sst_PD,xBootHi_sst_PD/), polyres)
  aumm01   = gsn_add_polyline(wks,plot(12),(/2,2/),(/xBootLow_sst_2xCO2,xBootHi_sst_2xCO2/), polyres)
  aumm001  = gsn_add_polyline(wks,plot(12),(/3,3/),(/xBootLow_sst_4xCO2,xBootHi_sst_4xCO2/), polyres)

  bumm1    = gsn_add_polyline(wks,plot(13),(/1,1/),(/xBootLow_ssh_PD,xBootHi_ssh_PD/), polyres)
  bumm01   = gsn_add_polyline(wks,plot(13),(/2,2/),(/xBootLow_ssh_2xCO2,xBootHi_ssh_2xCO2/), polyres)
  bumm001  = gsn_add_polyline(wks,plot(13),(/3,3/),(/xBootLow_ssh_4xCO2,xBootHi_ssh_4xCO2/), polyres)

  cumm1    = gsn_add_polyline(wks,plot(14),(/1,1/),(/xBootLow_va1_PD,xBootHi_va1_PD/), polyres)
  cumm01   = gsn_add_polyline(wks,plot(14),(/2,2/),(/xBootLow_va1_2xCO2,xBootHi_va1_2xCO2/), polyres)
  cumm001  = gsn_add_polyline(wks,plot(14),(/3,3/),(/xBootLow_va1_4xCO2,xBootHi_va1_4xCO2/), polyres)

  dumm1    = gsn_add_polyline(wks,plot(15),(/1,1/),(/xBootLow_va2_PD,xBootHi_va2_PD/), polyres)
  dumm01   = gsn_add_polyline(wks,plot(15),(/2,2/),(/xBootLow_va2_2xCO2,xBootHi_va2_2xCO2/), polyres)
  dumm001  = gsn_add_polyline(wks,plot(15),(/3,3/),(/xBootLow_va2_4xCO2,xBootHi_va2_4xCO2/), polyres)

  eumm1    = gsn_add_polyline(wks,plot(16),(/1,1/),(/xBootLow_eke_PD,xBootHi_eke_PD/), polyres)
  eumm01   = gsn_add_polyline(wks,plot(16),(/2,2/),(/xBootLow_eke_2xCO2,xBootHi_eke_2xCO2/), polyres)
  eumm001  = gsn_add_polyline(wks,plot(16),(/3,3/),(/xBootLow_eke_4xCO2,xBootHi_eke_4xCO2/), polyres)

  wid = 0.24/2
  fumm1   = gsn_add_polyline(wks,plot(12),(/1-wid,1+wid/),(/xBootLow_sst_PD,xBootLow_sst_PD/), polyres)
  fumm01  = gsn_add_polyline(wks,plot(12),(/1-wid,1+wid/),(/xBootHi_sst_PD,xBootHi_sst_PD/), polyres)

  gumm1   = gsn_add_polyline(wks,plot(12),(/2-wid,2+wid/),(/xBootLow_sst_2xCO2,xBootLow_sst_2xCO2/), polyres)
  gumm01  = gsn_add_polyline(wks,plot(12),(/2-wid,2+wid/),(/xBootHi_sst_2xCO2,xBootHi_sst_2xCO2/), polyres)

  humm1   = gsn_add_polyline(wks,plot(12),(/3-wid,3+wid/),(/xBootLow_sst_4xCO2,xBootLow_sst_4xCO2/), polyres)
  humm01  = gsn_add_polyline(wks,plot(12),(/3-wid,3+wid/),(/xBootHi_sst_4xCO2,xBootHi_sst_4xCO2/), polyres)
;----------------------------------------------------------------------------------------------------
  fumm2   = gsn_add_polyline(wks,plot(13),(/1-wid,1+wid/),(/xBootLow_ssh_PD,xBootLow_ssh_PD/), polyres)
  fumm02  = gsn_add_polyline(wks,plot(13),(/1-wid,1+wid/),(/xBootHi_ssh_PD,xBootHi_ssh_PD/), polyres)

  gumm2   = gsn_add_polyline(wks,plot(13),(/2-wid,2+wid/),(/xBootLow_ssh_2xCO2,xBootLow_ssh_2xCO2/), polyres)
  gumm02  = gsn_add_polyline(wks,plot(13),(/2-wid,2+wid/),(/xBootHi_ssh_2xCO2,xBootHi_ssh_2xCO2/), polyres)

  humm2   = gsn_add_polyline(wks,plot(13),(/3-wid,3+wid/),(/xBootLow_ssh_4xCO2,xBootLow_ssh_4xCO2/), polyres)
  humm02  = gsn_add_polyline(wks,plot(13),(/3-wid,3+wid/),(/xBootHi_ssh_4xCO2,xBootHi_ssh_4xCO2/), polyres)
  ;----------------------------------------------------------------------------------------------------
  fumm3   = gsn_add_polyline(wks,plot(14),(/1-wid,1+wid/),(/xBootLow_va1_PD,xBootLow_va1_PD/), polyres)
  fumm03  = gsn_add_polyline(wks,plot(14),(/1-wid,1+wid/),(/xBootHi_va1_PD,xBootHi_va1_PD/), polyres)

  gumm3   = gsn_add_polyline(wks,plot(14),(/2-wid,2+wid/),(/xBootLow_va1_2xCO2,xBootLow_va1_2xCO2/), polyres)
  gumm03  = gsn_add_polyline(wks,plot(14),(/2-wid,2+wid/),(/xBootHi_va1_2xCO2,xBootHi_va1_2xCO2/), polyres)

  humm3   = gsn_add_polyline(wks,plot(14),(/3-wid,3+wid/),(/xBootLow_va1_4xCO2,xBootLow_va1_4xCO2/), polyres)
  humm03  = gsn_add_polyline(wks,plot(14),(/3-wid,3+wid/),(/xBootHi_va1_4xCO2,xBootHi_va1_4xCO2/), polyres)
  ;----------------------------------------------------------------------------------------------------
  fumm4   = gsn_add_polyline(wks,plot(15),(/1-wid,1+wid/),(/xBootLow_va2_PD,xBootLow_va2_PD/), polyres)
  fumm04  = gsn_add_polyline(wks,plot(15),(/1-wid,1+wid/),(/xBootHi_va2_PD,xBootHi_va2_PD/), polyres)

  gumm4   = gsn_add_polyline(wks,plot(15),(/2-wid,2+wid/),(/xBootLow_va2_2xCO2,xBootLow_va2_2xCO2/), polyres)
  gumm04  = gsn_add_polyline(wks,plot(15),(/2-wid,2+wid/),(/xBootHi_va2_2xCO2,xBootHi_va2_2xCO2/), polyres)

  humm4   = gsn_add_polyline(wks,plot(15),(/3-wid,3+wid/),(/xBootLow_va2_4xCO2,xBootLow_va2_4xCO2/), polyres)
  humm04  = gsn_add_polyline(wks,plot(15),(/3-wid,3+wid/),(/xBootHi_va2_4xCO2,xBootHi_va2_4xCO2/), polyres)
;-----------------------------------------------------------------------------------------------------
  fumm5   = gsn_add_polyline(wks,plot(16),(/1-wid,1+wid/),(/xBootLow_eke_PD,xBootLow_eke_PD/), polyres)
  fumm05  = gsn_add_polyline(wks,plot(16),(/1-wid,1+wid/),(/xBootHi_eke_PD,xBootHi_eke_PD/), polyres)

  gumm5   = gsn_add_polyline(wks,plot(16),(/2-wid,2+wid/),(/xBootLow_eke_2xCO2,xBootLow_eke_2xCO2/), polyres)
  gumm05  = gsn_add_polyline(wks,plot(16),(/2-wid,2+wid/),(/xBootHi_eke_2xCO2,xBootHi_eke_2xCO2/), polyres)

  humm5   = gsn_add_polyline(wks,plot(16),(/3-wid,3+wid/),(/xBootLow_eke_4xCO2,xBootLow_eke_4xCO2/), polyres)
  humm05  = gsn_add_polyline(wks,plot(16),(/3-wid,3+wid/),(/xBootHi_eke_4xCO2,xBootHi_eke_4xCO2/), polyres)
;------------------------------------------------------------------------------------------------------

  draw(plot(12))
  draw(plot(13))
  draw(plot(14))
  draw(plot(15))
  draw(plot(16))
end 