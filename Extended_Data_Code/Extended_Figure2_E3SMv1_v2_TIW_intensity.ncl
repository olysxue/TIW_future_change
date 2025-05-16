;=================================================================
undef("readdata")
function readdata(ncdir, filelist, varname, scale)
local dims,temp,f
begin
   f= addfile(ncdir + filelist(0), "r")
   temp=f->$varname$(0,:,:)
   dims=getfilevardimsizes(f, varname)
   out=new((/3,dims(1),dims(2)/), typeof(temp))
   
   do i=0,1
   f:= addfile(ncdir + filelist(i), "r")
   out(i,:,:)=f->$varname$(0,:,:)
   out = out/scale
   end do 
   out(2,:,:)=out(1,:,:)-out(0,:,:) ; 4xCO2 - PD
   copy_VarMeta(temp,out(0,:,:))
   copy_VarMeta(temp,out(1,:,:))
   copy_VarMeta(temp,out(2,:,:))

  do i=0,2
  out(i,:,:)  = smth9_Wrap(out(i,:,:), 0.50, 0.25, False)
  end do
   return(out)
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
;-----------------------------------------------------------
; read E3SMv1 
ncdir  = "/Volumes/oly2/E3SM/mean/"
scale=1.
varname="tos" ;celsius
filelist1=(/"tos_day_E3SM-1-0_picontrol_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mean_new.nc",\
           "tos_day_E3SM-1-0_4xCO2_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mean_new.nc"/)
sst=readdata(ncdir, filelist1, varname, scale)
sst@units="celsius"

varname="zos" ;units m
filelist2=(/"zos_day_E3SM-1-0_picontrol_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mean_new.nc",\
           "zos_day_E3SM-1-0_4xCO2_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mean_new.nc"/)
ssh01=readdata(ncdir, filelist2, varname, scale)
ssh01@units="m"
ssh = ssh01 *100.
copy_VarMeta(ssh01, ssh)
printVarSummary(ssh)
delete(ssh01)

ncdir1  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/E3SM/E3SMv1//"
varname="KE" ;units m2/s2
filelist3=(/"EKE_picontrol_0001-0500_remap_equator_runave_mean.nc",\
           "EKE_abrupt_4xCO2_0001-0150_remap_equator_runave_mean.nc"/)
eke=readdata(ncdir1, filelist3, varname, scale)
eke@units="J/m2"
;=============================================================
; read E3SMv2 
ncdir1  = "/Volumes/oly2/E3SM/E3SMv2/"
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

sst1 = new((/2, dimsizes(tos_pattern1&lat), dimsizes(tos_pattern1&lon)/), typeof(tos_pattern1))
sst1(0,:,:) = tos_pattern1 
sst1(1,:,:) = tos_pattern2

ssh1 = new((/2, dimsizes(zos_pattern1&lat), dimsizes(zos_pattern1&lon)/), typeof(zos_pattern1))
ssh1(0,:,:) = zos_pattern1 
ssh1(1,:,:) = zos_pattern2

eke1 = new((/2, dimsizes(KE_pattern1&lat), dimsizes(KE_pattern1&lon)/), typeof(KE_pattern1))
eke1(0,:,:) = KE_pattern1 
eke1(1,:,:) = KE_pattern2
;=============================================================
; plot mean TIW EKE and RMS
;=============================================================
outdir = "/Users/xueaoyun/Documents/work/TIW_future/revision1/Extend_figure/"
wks = gsn_open_wks("pdf",outdir+"Extend_Figure2_E3SMv1_v2_PD_4xCO2_TIW_mean_intensity")
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

plot=new(18,graphic)
dum=plot

; ----->1 plot eke 
panels = (/"~F22~a", "~F22~c", "~F22~e", "~F22~b", "~F22~d", "~F22~f", "~F22~g", "~F22~i","~F22~k","~F22~h","~F22~j", "~F22~l","~F22~n","~F22~o","","~F22~p"/)
text=(/"PiControl","4xCO2"/)
yspace=0.02
yspace1 = 0.03
yf=0.9
res@tiDeltaF=3
res@gsnCenterStringOrthogonalPosF=0.13
res@cnLevels        =fspan(60,100,41)
res@vpXF    = 0.08
res@lbTitleString    = "J/m~S~3~N~"
do i=0,1
res@vpYF    = yf-i*yspace-i*res@vpHeightF
res@tiYAxisString=text(i)
res@tiYAxisOffsetXF=-0.02
res@tiYAxisFontHeightF=0.014
if i.eq.0 then
res@gsnCenterString = "E3SMv1 EKE"
res@tiMainString=""
else
res@gsnCenterString = ""
res@tiMainString=""
end if

if i.eq.1 then
res@lbLabelBarOn    = True
res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W",""/)  
else
res@lbLabelBarOn    = False
res@tmXBLabels = (/"","","","","",""/)
end if
res@gsnLeftString  = panels(i*3)
res@gsnLeftStringFontHeightF = 0.014
res@gsnLeftStringOrthogonalPosF = 0.
res@gsnLeftStringParallelPosF  = 0.03
plot(i*3)  = gsn_csm_contour_map(wks,eke(i,:,:),res)
end do
delete([/res@tiYAxisString,res@tiYAxisOffsetXF,res@tiYAxisFontHeightF/])

; ----->2 plot sst
res@lbTitleString    = "degC"
delete(res@cnLevels)
res@vpXF    = res@vpWidthF+res@vpXF+yspace1
res@cnLevels             = fspan(0.2,0.8,41)
res@cnFillPalette         = "WhiteBlueGreenYellowOrangeRed"
; res@cnFillColors         = (/0,12,13,14,15,17,19/)
res@tmYLLabels = (/"","","","",""/)
do i=0,1
res@vpYF    = yf-i*yspace-i*res@vpHeightF
if i.eq.0 then
res@gsnCenterString = "E3SMv1 SST"
else
res@gsnCenterString = ""
end if 
if i.eq.1 then
res@lbLabelBarOn    = True
res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W",""/)  
else
res@lbLabelBarOn    = False
res@tmXBLabels = (/"","","","","",""/)
end if
res@gsnLeftString  = panels(1+i*3)
res@gsnLeftStringFontHeightF = 0.014
res@gsnLeftStringOrthogonalPosF = 0.0
res@gsnLeftStringParallelPosF  = 0.03
plot(1+i*3)  = gsn_csm_contour_map(wks,sst(i,:,:),res)
end do
; ----->3 plot ssh 
delete(res@cnLevels)
res@vpXF    = res@vpWidthF+res@vpXF+yspace1
res@lbTitleString    = "10~S~-2~N~ m"
res@cnLevels              = fspan(2,4,21)
res@tmYLLabels = (/"","","","",""/)
do i=0,1
res@vpYF    = yf-i*yspace-i*res@vpHeightF
if i.eq.0 then
res@gsnCenterString = "E3SMv1 SSH"
else
res@gsnCenterString = ""
end if
if i.eq.1 then
res@lbLabelBarOn    = True
res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W",""/)  
else
res@lbLabelBarOn    = False
res@tmXBLabels = (/"","","","","",""/)
end if
res@gsnLeftString  = panels(2+i*3)
res@gsnLeftStringFontHeightF = 0.014
res@gsnLeftStringOrthogonalPosF = 0.0
plot(2+i*3)  = gsn_csm_contour_map(wks,ssh(i,:,:),res)
end do
;--------------------------------------------------------------------------------------------------------------
;==========================================E3SMv2================================================================
delete(res@cnLevels)
yf1 = res@vpYF-res@vpHeightF-0.12 
res@tmYLValues = (/-10,-5,0,5,10/)
res@tmYLLabels = (/"10~S~o~N~S","5~S~o~N~S","EQ","5~S~o~N~N","10~S~o~N~N"/)
; ----->1 plot eke
res@vpXF    = 0.08
res@gsnSpreadColorStart= 0
res@gsnSpreadColorEnd  =230
res@cnLevels        =fspan(40,80,41)
do i=0,1
res@vpYF    = yf1-i*yspace-i*res@vpHeightF
print(res@vpYF)
res@tiYAxisString=text(i)
res@gsnLeftString  = panels(3+i*3+3)
res@tiYAxisOffsetXF=-0.02
res@tiYAxisFontHeightF=0.014
if i.eq.0 then
res@gsnCenterString = "E3SMv2 EKE"
res@tiMainString=""
res@tmXBLabels = (/"","","","","",""/) 
else
res@gsnCenterString = ""
res@tiMainString=""
end if

if i.eq.1 then
res@lbLabelBarOn    = True
res@lbTitleString    = "J/m~S~3~N~"

res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W",""/)  
else
res@lbLabelBarOn    = False
end if

plot(3+i*3+3)  = gsn_csm_contour_map(wks,eke1(i,:,:),res)
end do
delete([/res@cnLevels,res@cnFillColors/])
delete([/res@tiYAxisString,res@tiYAxisOffsetXF,res@tiYAxisFontHeightF/])
;============================================================================
; ----->2 plot sst
res@vpXF    = res@vpWidthF+res@vpXF+yspace1
res@cnLevels                = fspan(0.2,0.6,41)

;res@cnFillColors         =  (/0,144,168,178,188,204,230,246/)
res@lbTitleString    = ""
; res@cnFillColors         = (/0,12,13,14,15,17,19/)
res@tmYLLabels = (/"","","","",""/)
do i=0,1
res@vpYF    = yf1-i*yspace-i*res@vpHeightF
res@gsnLeftString  = panels(4+i*3+3)
if i.eq.0 then
res@gsnCenterString = "E3SMv2 SST"
res@tmXBLabels = (/"","","","","",""/) 
else
res@gsnCenterString = ""
end if
if i.eq.1 then
res@lbLabelBarOn    = True
res@lbTitleString    = "degC"

res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W",""/)  
else
res@lbLabelBarOn    = False
end if
plot(4+i*3+3)  = gsn_csm_contour_map(wks,sst1(i,:,:),res)
end do
delete([/res@cnLevels,res@cnFillColors/])
; ----->3 plot ssh
res@vpXF    = res@vpWidthF+res@vpXF+yspace1
res@tmYLLabels = (/"","","","",""/)
res@cnLevels        =fspan(2,3,21)
do i=0,1
res@vpYF    = yf1-i*yspace-i*res@vpHeightF
res@gsnLeftString  = panels(5+i*3+3)
if i.eq.0 then
res@gsnCenterString = "E3SMv2 SSH"
res@tmXBLabels = (/"","","","","",""/) 
else
res@gsnCenterString = ""
end if
if i.eq.1 then
res@lbLabelBarOn    = True
res@lbTitleString    = "10~S~-2~N~ m"
res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W",""/)  
else
res@lbLabelBarOn    = False
res@tmXBLabels = (/"","","","","",""/)
end if
plot(5+i*3+3)  = gsn_csm_contour_map(wks,ssh1(i,:,:),res)
end do
delete([/res@cnLevels,res@cnFillColors/])
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
n2 = 8

s3 = 3
n3 = 6

dum(0)=gsn_add_polyline(wks, plot(0), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum(1)=gsn_add_polyline(wks, plot(2), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
dum(2)=gsn_add_polyline(wks, plot(1), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum(3)=gsn_add_polyline(wks, plot(3), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum(4)=gsn_add_polyline(wks, plot(5), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
dum(5)=gsn_add_polyline(wks, plot(4), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)

dum(6)=gsn_add_polyline(wks, plot(6), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum(7)=gsn_add_polyline(wks, plot(8), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
dum(8)=gsn_add_polyline(wks, plot(7), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum(9) =gsn_add_polyline(wks,plot(9), (/w3,e3,e3,w3,w3/), (/s1,s1,n1,n1,s1/), plres)
dum(10)=gsn_add_polyline(wks, plot(11), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
dum(11)=gsn_add_polyline(wks, plot(10), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)


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
delete([/sst, ssh, eke/])
delete([/sst1, ssh1, eke1/])
;=======================================================
; read data
;=======================================================
;read E3SMv1 
varname="tos"
ncdir  = "/Volumes/oly2/E3SM/mean/"
scale=1. ;units degC
filelist=(/"tos_day_E3SM-1-0_picontrol_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mon.nc","tos_day_E3SM-1-0_4xCO2_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mon.nc"/)
sst0=readdata1(ncdir, filelist, varname, scale)
sst0@units="degC"
sst1 = sqrt(sst0)
copy_VarMeta(sst0, sst1)
varname="zos"
scale=1. ; units cm to m by scale
filelist=(/"zos_day_E3SM-1-0_picontrol_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mon.nc","zos_day_E3SM-1-0_4xCO2_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mon.nc"/)
ssh0=readdata1(ncdir, filelist, varname, scale)
ssh0@units="meter"
printVarSummary(ssh0)
ssh1 = sqrt(ssh0)*100.
copy_VarMeta(ssh0, ssh1)
scale=1.

ncdir1  = "/Volumes/oly2/E3SM/E3SMv1/"
filename1   = "EKE_mon_picontrol_0001-0500_remap_equator_runave.nc"
filename2   = "EKE_mon_abrupt_4xCO2_0001-0155_remap_equator_runave.nc"
fil1        = systemfunc("ls " + ncdir1 + filename1)
fil2        = systemfunc("ls " + ncdir1 + filename2)
f1          =  addfile(fil1, "r")
f2          =  addfile(fil2, "r")
eke01       = f1->KE(0:140*12-1,:,:)
eke02       = f2->KE(0:140*12-1,:,:)

eke1 = new((/2, dimsizes(eke01&time),dimsizes(eke01&lat), dimsizes(eke01&lon)/), typeof(eke01))
eke1(0,:,:,:) = eke01
eke1(1,:,:,:) = eke02
copy_VarCoords(eke01, eke1(0,:,:,:))
copy_VarCoords(eke01, eke1(1,:,:,:))

sst1_PD    = wgt_areaave_Wrap(sst1(0,:,{0:5},{200:260}), 1.0, 1.0, 0)
sst1_4xCO2 = wgt_areaave_Wrap(sst1(1,:,{0:5},{200:260}), 1.0, 1.0, 0)

ssh1_PD    = wgt_areaave_Wrap(ssh1(0,:,{3:8},{210:250}), 1.0, 1.0, 0)
ssh1_4xCO2 = wgt_areaave_Wrap(ssh1(1,:,{3:8},{200:240}), 1.0, 1.0, 0)

eke1_PD    = wgt_areaave_Wrap(eke1(0,:,{0:5},{210:260}), 1.0, 1.0, 0)
eke1_4xCO2 = wgt_areaave_Wrap(eke1(1,:,{0:5},{200:250}), 1.0, 1.0, 0)
;=======================================================================
stat  = 0                         ; 0=mean; 1=variance; 2=std. dev; ...
n     = 1800                        ; sub-sample size;  default is n=N
nBoot = 10000                      ; bootstrap replications
N     = 1800
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
;---------------------------------------------------------------------------------------------------------
stdarr1 = new((/2/),typeof(sst1_PD)) ;  SST
stdarr1(0) =  xBoot_sstavg1_PD
stdarr1(1) =  xBoot_sstavg1_4xCO2

stdarr2 = new((/2/),typeof(ssh1_PD)) ; SSH
stdarr2(0) =  xBoot_sshavg1_PD
stdarr2(1) =  xBoot_sshavg1_4xCO2

stdarr3 = new((/2/),typeof(eke1_PD)) ; EKE 
stdarr3(0) =  xBoot_ekeavg1_PD
stdarr3(1) =  xBoot_ekeavg1_4xCO2 
;===================================================================================
; read E3SMv2 
;====================================================================================
ncdir1  = "/Volumes/oly2/E3SM/E3SMv2/"
fil1 = systemfunc("ls " + ncdir1+"tos_mon_E3SM-2-0_0001-0500_remap_equator_10-50fft_new_1mon_std_1degree.nc")
f1   = addfile(fil1,"r")
sst01 = f1->sst(:,:,:)
printVarSummary(sst01)
;---------------------------------------------------------
fil2 = systemfunc("ls " + ncdir1 + "sst_mon_abrupt_4xCO2_0001-0160_remap_equator_moving_20yr_detrend_10-50fft_1mon_runstd_1degree.nc")
f2   = addfile(fil2,"r")
sst02 = f2->sst(:,:,:)
printVarSummary(sst02)
;===========================================================================
fil01 = systemfunc("ls " + ncdir1+"zos_mon_E3SM-2-0_0001-0500_remap_equator_10-50fft_new_1mon_std_1degree.nc")
f01   = addfile(fil01,"r")
ssh001 = f01->ssh(:,:,:)
ssh01  = ssh001 *100.
copy_VarCoords(ssh001, ssh01)
printVarSummary(ssh01)
;---------------------------------------------------------
fil02 = systemfunc("ls " + ncdir1 +"ssh_mon_abrupt_4xCO2_0001-0160_remap_equator_moving_20yr_detrend_10-50fft_1mon_runstd_1degree.nc")
f02   = addfile(fil02,"r")
ssh002 = f02->ssh(:,:,:)
ssh02  = ssh002 *100.
copy_VarCoords(ssh002, ssh02)
printVarSummary(ssh02)
;===========================================================================
fil001 = systemfunc("ls " + ncdir1+"EKE_mon_picontrol_0001-0500_remap_equator_runave.nc")
f001   = addfile(fil001,"r")
KE01 = f001->KE(:,:,:)
printVarSummary(KE01)
;---------------------------------------------------------
fil002 = systemfunc("ls " + ncdir1 +"EKE_mon_abrupt_4xCO2_0001-0160_remap_equator_runave.nc")
f002   = addfile(fil002,"r")
KE02 = f002->KE(:,:,:)
printVarSummary(KE02)

sst01_PD    = wgt_areaave_Wrap(sst01(:,{0:5},{200:260}), 1.0, 1.0, 0)
sst01_4xCO2 = wgt_areaave_Wrap(sst02(:,{0:5},{200:260}), 1.0, 1.0, 0)

ssh01_PD    = wgt_areaave_Wrap(ssh01(:,{3:6},{210:250}), 1.0, 1.0, 0)
ssh01_4xCO2 = wgt_areaave_Wrap(ssh02(:,{3:6},{200:240}), 1.0, 1.0, 0)

eke01_PD    = wgt_areaave_Wrap(KE01(:,{0:5},{210:260}), 1.0, 1.0, 0)
eke01_4xCO2 = wgt_areaave_Wrap(KE02(:,{0:5},{200:250}), 1.0, 1.0, 0)

stat  = 0                         ; 0=mean; 1=variance; 2=std. dev; ...
n     = 400                       ; sub-sample size;  default is n=N
nBoot = 10000                      ; bootstrap replications
N     = 10000
opt   = False                     ; use defaults
if (isvar("n") .and. (n.lt.N)) then
    opt   = True
    opt@sample_size = n           ; default N
end if

 BootStrap_sst01_PD     = bootstrap_stat(sst01_PD, stat, nBoot, 0, opt)
 BootStrap_sst01_4xCO2  = bootstrap_stat(sst01_4xCO2, stat, nBoot, 0, opt)

 BootStrap_ssh01_PD     = bootstrap_stat(ssh01_PD, stat, nBoot, 0, opt)
 BootStrap_ssh01_4xCO2  = bootstrap_stat(ssh01_4xCO2, stat, nBoot, 0, opt)

 BootStrap_eke01_PD     = bootstrap_stat(eke01_PD, stat, nBoot, 0, opt)
 BootStrap_eke01_4xCO2  = bootstrap_stat(eke01_4xCO2, stat, nBoot, 0, opt)

    xBoot_sst01_PD    = BootStrap_sst01_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sstavg01_PD    = BootStrap_sst01_PD[1]       ; Mean of bootstrapped estimates
 xBoot_sststd01_PD    = BootStrap_sst01_PD[2]       ; Std dev of bootstrapped estimates

    xBoot_ssh01_PD    = BootStrap_ssh01_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sshavg01_PD    = BootStrap_ssh01_PD[1]       ; Mean of bootstrapped estimates
 xBoot_sshstd01_PD    = BootStrap_ssh01_PD[2]       ; Std dev of bootstrapped estimates

    xBoot_eke01_PD    = BootStrap_eke01_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_ekeavg01_PD    = BootStrap_eke01_PD[1]       ; Mean of bootstrapped estimates
 xBoot_ekestd01_PD    = BootStrap_eke01_PD[2]       ; Std dev of bootstrapped estimates
;=====================================================
    xBoot_sst01_4xCO2    = BootStrap_sst01_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sstavg01_4xCO2    = BootStrap_sst01_4xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_sststd01_4xCO2    = BootStrap_sst01_4xCO2[2]       ; Std dev of bootstrapped estimates

    xBoot_ssh01_4xCO2    = BootStrap_ssh01_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sshavg01_4xCO2    = BootStrap_ssh01_4xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_sshstd01_4xCO2    = BootStrap_ssh01_4xCO2[2]       ; Std dev of bootstrapped estimates

    xBoot_eke01_4xCO2    = BootStrap_eke01_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_ekeavg01_4xCO2    = BootStrap_eke01_4xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_ekestd01_4xCO2    = BootStrap_eke01_4xCO2[2]       ; Std dev of bootstrapped estimates
;=====================================================
 xBootLow_sst01_PD    = bootstrap_estimate(xBoot_sst01_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_sst01_PD    = bootstrap_estimate(xBoot_sst01_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_sst01_PD    = bootstrap_estimate(xBoot_sst01_PD, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_ssh01_PD    = bootstrap_estimate(xBoot_ssh01_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_ssh01_PD    = bootstrap_estimate(xBoot_ssh01_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_ssh01_PD    = bootstrap_estimate(xBoot_ssh01_PD, 0.975, False)   ; 97.5% upper confidence bound
 
 xBootLow_eke01_PD    = bootstrap_estimate(xBoot_eke01_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_eke01_PD    = bootstrap_estimate(xBoot_eke01_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_eke01_PD    = bootstrap_estimate(xBoot_eke01_PD, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------
xBootLow_sst01_4xCO2    = bootstrap_estimate(xBoot_sst01_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_sst01_4xCO2    = bootstrap_estimate(xBoot_sst01_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_sst01_4xCO2    = bootstrap_estimate(xBoot_sst01_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_ssh01_4xCO2    = bootstrap_estimate(xBoot_ssh01_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_ssh01_4xCO2    = bootstrap_estimate(xBoot_ssh01_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_ssh01_4xCO2    = bootstrap_estimate(xBoot_ssh01_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_eke01_4xCO2    = bootstrap_estimate(xBoot_eke01_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_eke01_4xCO2    = bootstrap_estimate(xBoot_eke01_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_eke01_4xCO2    = bootstrap_estimate(xBoot_eke01_4xCO2, 0.975, False)   ; 97.5% upper confidence bound
;=============================================================================
stdarr01 = new((/2/),typeof(sst01_PD)) ;  SST
stdarr01(0) =  xBoot_sstavg01_PD
stdarr01(1) =  xBoot_sstavg01_4xCO2

stdarr02 = new((/2/),typeof(ssh01_PD)) ; SSH
stdarr02(0) =  xBoot_sshavg01_PD
stdarr02(1) =  xBoot_sshavg01_4xCO2

stdarr03 = new((/2/),typeof(eke01_PD)) ; EKE 
stdarr03(0) =  xBoot_ekeavg01_PD
stdarr03(1) =  xBoot_ekeavg01_4xCO2 
;===================================================================================
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
 xyres@tmXBValues        = (/1,2/)
 xyres@tmXBLabels = (/"PI","4xCO2"/)
 xyres@tmXBLabelFontHeightF = 0.01
 xyres@tmYLMode          = "Explicit"         ; explicit labels
 xyres@tmYLValues        = (/0,20,40,60,80,100/)
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
 xyres@vpXF    = 0.08
 yspace=0.03
 yf2= res@vpYF- res@vpWidthF-0.02
 xyres@vpYF    = yf2
 xyres@gsnXYBarChartColors = (/"orange3"/)
 xyres@gsnCenterString       = "E3SMv1_EKE"
 xyres@gsnCenterStringOrthogonalPosF = -0.125
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString  = "~F22~m"
 xyres@gsnLeftStringFontHeightF = 0.014
 xyres@gsnLeftStringOrthogonalPosF = 0.
 xyres@gsnLeftStringParallelPosF = -0.05
 plot(12)  = gsn_csm_xy(wks,ispan(1,2,1),stdarr3(:),xyres)   ; advancing the frame
 delete(xyres@tmXBValues)
 delete(xyres@tmXBLabels)
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 3.                   ; of the 1st and last bars
 xyres@tmXBValues        = (/1,2/)
 xyres@tmXBLabels = (/"PI","4xCO2"/)
 xyres@tmYLOn        = False 
 xyres@tmYRBorderOn  = True                        ; turn off right axis
 xyres@tmYLBorderOn  = False
 xyres@vpXF = xyres@vpWidthF+xyres@vpXF
 xyres@gsnLeftString  = ""
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 100.                ; bring bars down to zero
 xyres@gsnCenterString       = "E3SMv2_EKE"
 xyres@gsnCenterStringOrthogonalPosF = -0.10
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftStringFontHeightF = 0.014
 xyres@gsnLeftStringOrthogonalPosF = 0.
 plot(13)  = gsn_csm_xy(wks,ispan(1,2,1),stdarr03(:),xyres)   ; advancing the frame
 xyres@vpXF = xyres@vpXF+xyres@vpWidthF+yspace
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 delete(xyres@tmXBValues)
 delete(xyres@tmXBLabels)
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 3.                   ; of the 1st and last bars
 xyres@tmXBValues        = (/1,2/)
 xyres@tmXBLabels = (/"PI","4xCO2"/)
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
 xyres@gsnCenterString       = "E3SMv1_SST"
 xyres@gsnCenterStringOrthogonalPosF = -0.125
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString  = "~F22~n"
 xyres@gsnLeftStringFontHeightF = 0.014
 xyres@gsnLeftStringOrthogonalPosF = 0.
 xyres@gsnLeftStringParallelPosF = -0.05
 plot(14)  = gsn_csm_xy(wks,ispan(1,2,1),stdarr1(:),xyres)   ; advancing the frame
 xyres@tmYLOn        = False 
 xyres@tmYRBorderOn  = True                        ; turn off right axis
 xyres@tmYLBorderOn  = False
 delete(xyres@tmXBValues)
 delete(xyres@tmXBLabels)
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 3.                   ; of the 1st and last bars
 xyres@tmXBValues        = (/1,2/)
 xyres@tmXBLabels = (/"PI","4CO2"/)
 xyres@vpXF = xyres@vpWidthF+xyres@vpXF
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.5                 ; bring bars down to zero
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 3.                   ; of the 1st and last bars
 xyres@gsnCenterString       = "E3SMv2_SST"
 xyres@gsnLeftString  = " "
 xyres@gsnCenterStringOrthogonalPosF = -0.10
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftStringFontHeightF = 0.014
 xyres@gsnLeftStringOrthogonalPosF = 0.
 plot(15)  = gsn_csm_xy(wks,ispan(1,2,1),stdarr01(:),xyres)   ; advancing the frame
 delete(xyres@tmXBValues)
 delete(xyres@tmXBLabels)
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 3.                   ; of the 1st and last bars
 xyres@tmXBValues        = (/1,2/)
 xyres@tmXBLabels = (/"PI","4xCO2"/)
 xyres@vpXF = xyres@vpXF+xyres@vpWidthF+yspace
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@tmYLOn        = True
 xyres@tmYRBorderOn  = False                        ; turn off right axis
 xyres@tmYLBorderOn  = True 
 xyres@vpWidthF              = 0.10
 xyres@vpHeightF             = 0.12
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 4.                 ; bring bars down to zero                   
 xyres@tmYLOn        = True 
 xyres@tmYLMode          = "Explicit"         ; explicit labels
 xyres@tmYLValues        = (/0,1,2,3,4,5/)
 xyres@tmYLLabels        = (/"0","1","2","3","4","5"/)
 xyres@gsnXYBarChartColors = (/"pink"/)
 xyres@gsnCenterString       = "E3SMv1_SSH"
 xyres@gsnCenterStringOrthogonalPosF = -0.125
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString  = "~F22~o"
 xyres@gsnLeftStringFontHeightF = 0.014
 xyres@gsnLeftStringOrthogonalPosF = 0.
 xyres@gsnLeftStringParallelPosF = -0.05
 plot(16)  = gsn_csm_xy(wks,ispan(1,2,1),stdarr2(:),xyres)   ; advancing the frame
 delete(xyres@tmXBValues)
 delete(xyres@tmXBLabels)
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 3.                   ; of the 1st and last bars
 xyres@tmXBValues        = (/1,2/)
 xyres@tmXBLabels = (/"PI","4xCO2"/)
 xyres@tmYLOn        = False 
 xyres@tmYRBorderOn  = True                        ; turn off right axis
 xyres@tmYLBorderOn  = False
 xyres@vpXF = xyres@vpWidthF+xyres@vpXF
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 4.                 ; bring bars down to zero
 xyres@gsnCenterString       = "E3SMv2_SSH"
 xyres@gsnLeftString  = ""
 xyres@gsnCenterStringOrthogonalPosF = -0.10
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftStringFontHeightF = 0.014
 xyres@gsnLeftStringOrthogonalPosF = 0.
 plot(17)  = gsn_csm_xy(wks,ispan(1,2,1),stdarr02(:),xyres)   ; advancing the frame


  polyres=True
  polyres@gsLineDashPattern=0
  polyres@gsLineColor="black"
  polyres@gsLineThicknessF=1.0
  polyres@tfPolyDrawOrder = "postDraw"

  
  eumm1    = gsn_add_polyline(wks,plot(12),(/1,1/),(/xBootLow_eke1_PD,xBootHi_eke1_PD/), polyres)
  eumm01   = gsn_add_polyline(wks,plot(12),(/2,2/),(/xBootLow_eke1_4xCO2,xBootHi_eke1_4xCO2/), polyres)

  fumm1    = gsn_add_polyline(wks,plot(13),(/1,1/),(/xBootLow_eke01_PD,xBootHi_eke01_PD/), polyres)
  fumm01   = gsn_add_polyline(wks,plot(13),(/2,2/),(/xBootLow_eke01_4xCO2,xBootHi_eke01_4xCO2/), polyres)


  aumm1    = gsn_add_polyline(wks,plot(14),(/1,1/),(/xBootLow_sst1_PD,xBootHi_sst1_PD/), polyres)
  aumm01   = gsn_add_polyline(wks,plot(14),(/2,2/),(/xBootLow_sst1_4xCO2,xBootHi_sst1_4xCO2/), polyres)

  bumm1    = gsn_add_polyline(wks,plot(15),(/1,1/),(/xBootLow_sst01_PD,xBootHi_sst01_PD/), polyres)
  bumm01   = gsn_add_polyline(wks,plot(15),(/2,2/),(/xBootLow_sst01_4xCO2,xBootHi_sst01_4xCO2/), polyres)

  cumm1    = gsn_add_polyline(wks,plot(16),(/1,1/),(/xBootLow_ssh1_PD,xBootHi_ssh1_PD/), polyres)
  cumm01   = gsn_add_polyline(wks,plot(16),(/2,2/),(/xBootLow_ssh1_4xCO2,xBootHi_ssh1_4xCO2/), polyres)

  dumm1    = gsn_add_polyline(wks,plot(17),(/1,1/),(/xBootLow_ssh01_PD,xBootHi_ssh01_PD/), polyres)
  dumm01   = gsn_add_polyline(wks,plot(17),(/2,2/),(/xBootLow_ssh01_4xCO2,xBootHi_ssh01_4xCO2/), polyres)


  wid = 0.24/2

  gumm3   = gsn_add_polyline(wks,plot(12),(/1-wid,1+wid/),(/xBootLow_eke1_PD,xBootLow_eke1_PD/), polyres)
  gumm03  = gsn_add_polyline(wks,plot(12),(/1-wid,1+wid/),(/xBootHi_eke1_PD,xBootHi_eke1_PD/), polyres)

  humm3   = gsn_add_polyline(wks,plot(12),(/2-wid,2+wid/),(/xBootLow_eke1_4xCO2,xBootLow_eke1_4xCO2/), polyres)
  humm03  = gsn_add_polyline(wks,plot(12),(/2-wid,2+wid/),(/xBootHi_eke1_4xCO2,xBootHi_eke1_4xCO2/), polyres)

  gumm1   = gsn_add_polyline(wks,plot(14),(/1-wid,1+wid/),(/xBootLow_sst1_PD,xBootLow_sst1_PD/), polyres)
  gumm01  = gsn_add_polyline(wks,plot(14),(/1-wid,1+wid/),(/xBootHi_sst1_PD,xBootHi_sst1_PD/), polyres)

  humm1   = gsn_add_polyline(wks,plot(14),(/2-wid,2+wid/),(/xBootLow_sst1_4xCO2,xBootLow_sst1_4xCO2/), polyres)
  humm01  = gsn_add_polyline(wks,plot(14),(/2-wid,2+wid/),(/xBootHi_sst1_4xCO2,xBootHi_sst1_4xCO2/), polyres)
;----------------------------------------------------------------------------------------------------
  gumm2   = gsn_add_polyline(wks,plot(16),(/1-wid,1+wid/),(/xBootLow_ssh1_PD,xBootLow_ssh1_PD/), polyres)
  gumm02  = gsn_add_polyline(wks,plot(16),(/1-wid,1+wid/),(/xBootHi_ssh1_PD,xBootHi_ssh1_PD/), polyres)
  
  humm2   = gsn_add_polyline(wks,plot(16),(/2-wid,2+wid/),(/xBootLow_ssh1_4xCO2,xBootLow_ssh1_4xCO2/), polyres)
  humm02  = gsn_add_polyline(wks,plot(16),(/2-wid,2+wid/),(/xBootHi_ssh1_4xCO2,xBootHi_ssh1_4xCO2/), polyres)
;----------------------------------------------------------------------------------------------------

  gumm6   = gsn_add_polyline(wks,plot(13),(/1-wid,1+wid/),(/xBootLow_eke01_PD,xBootLow_eke01_PD/), polyres)
  gumm06  = gsn_add_polyline(wks,plot(13),(/1-wid,1+wid/),(/xBootHi_eke01_PD,xBootHi_eke01_PD/), polyres)
  
  humm6   = gsn_add_polyline(wks,plot(13),(/2-wid,2+wid/),(/xBootLow_eke01_4xCO2,xBootLow_eke01_4xCO2/), polyres)
  humm06  = gsn_add_polyline(wks,plot(13),(/2-wid,2+wid/),(/xBootHi_eke01_4xCO2,xBootHi_eke01_4xCO2/), polyres)
;------------------------------------------------------------------------------------------------------
  gumm4   = gsn_add_polyline(wks,plot(15),(/1-wid,1+wid/),(/xBootLow_sst01_PD,xBootLow_sst01_PD/), polyres)
  gumm04  = gsn_add_polyline(wks,plot(15),(/1-wid,1+wid/),(/xBootHi_sst01_PD,xBootHi_sst01_PD/), polyres)
  
  humm4   = gsn_add_polyline(wks,plot(15),(/2-wid,2+wid/),(/xBootLow_sst01_4xCO2,xBootLow_sst01_4xCO2/), polyres)
  humm04  = gsn_add_polyline(wks,plot(15),(/2-wid,2+wid/),(/xBootHi_sst01_4xCO2,xBootHi_sst01_4xCO2/), polyres)
;----------------------------------------------------------------------------------------------------
  gumm5  = gsn_add_polyline(wks,plot(17),(/1-wid,1+wid/),(/xBootLow_ssh01_PD,xBootLow_ssh01_PD/), polyres)
  gumm05 = gsn_add_polyline(wks,plot(17),(/1-wid,1+wid/),(/xBootHi_ssh01_PD,xBootHi_ssh01_PD/), polyres)
  
  humm5  = gsn_add_polyline(wks,plot(17),(/2-wid,2+wid/),(/xBootLow_ssh01_4xCO2,xBootLow_ssh01_4xCO2/), polyres)
  humm05 = gsn_add_polyline(wks,plot(17),(/2-wid,2+wid/),(/xBootHi_ssh01_4xCO2,xBootHi_ssh01_4xCO2/), polyres)
;----------------------------------------------------------------------------------------------------
;======================================================================================================

  draw(plot(12))
  draw(plot(13))
  draw(plot(14))
  draw(plot(15))
  draw(plot(16))
  draw(plot(17))

end 