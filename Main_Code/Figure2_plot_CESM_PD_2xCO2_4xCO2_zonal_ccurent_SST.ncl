undef("readdata")
function readdata(ncdir, filelist, varname, scale, yrStrt, yrLast, ilev1, ilev2, lonw, lone)
local dims,temp,f
begin
f= addfile(ncdir + filelist(0), "r")
temp=f->$varname$(0,:,:,:) ; nlev,nlat,nlon
printVarSummary(temp)
dims=dimsizes(temp)
varave_xy=new((/5,dims(1),dims(2)/), float)
dady_xy  =new((/5,dims(1),dims(2)/), float) ;5,nlat,nlon

varave_yz=new((/5,dims(0),dims(1)/), float)
dady_yz  =new((/5,dims(0),dims(1)/), float) ;5,nlev,nlat

do i=0,2
f   = addfile(ncdir + filelist(i), "r")
YYYY= cd_calendar(f->time,-1)/100
iyr = ind(YYYY.ge.yrStrt(i) .and. YYYY.le.yrLast(i))
var0= f->$varname$(iyr,:,:,:)
var = var0/scale
copy_VarMeta(var0, var)

;---->cal dy
; parameter constant 
PI  = 3.14159265
R   = 6378388.
rad = PI/180.
dens= 1025 ; reference density of seawater
sec = 3600*24*30.42 ;one month

dlat   = (var&lat(2)-var&lat(1))*rad
dY     = R*dlat
clm = clmMonTLLL(var) ; 12month,lev,lat,lon
; varanom= calcMonAnomTLLL(var, varclm)
varave = dim_avg_n_Wrap(clm,0) ;lev,lat,lon
dady   = center_finite_diff_n(varave, dY, False, 0, 1)*1000000
copy_VarMeta(var(0,:,:,:), dady)

;---->vertical mean 
varave_xy(i,:,:)= dim_avg_n_Wrap(varave(ilev1:ilev2,:,:),0)
dady_xy(i,:,:)  = dim_avg_n_Wrap(dady(ilev1:ilev2,:,:),0)

;---->longitude mean 
varave_yz(i,:,:) = dim_avg_n_Wrap(varave(:,:,{lonw:lone}),2)
dady_yz(i,:,:)= dim_avg_n_Wrap(dady(:,:,{lonw:lone}),2)

delete([/YYYY,iyr,var0,var,varave,dady/])
print("i="+i)
end do

varave_xy(3,:,:)=varave_xy(2,:,:)-varave_xy(0,:,:)
varave_xy(4,:,:)=varave_xy(1,:,:)-varave_xy(0,:,:)

  dady_xy(3,:,:)=dady_xy(2,:,:)-dady_xy(0,:,:)
  dady_xy(4,:,:)=dady_xy(1,:,:)-dady_xy(0,:,:)

varave_yz(3,:,:)=varave_yz(2,:,:)-varave_yz(0,:,:)
varave_yz(4,:,:)=varave_yz(1,:,:)-varave_yz(0,:,:)

  dady_yz(3,:,:)=dady_yz(2,:,:)-dady_yz(0,:,:)
  dady_yz(4,:,:)=dady_yz(1,:,:)-dady_yz(0,:,:)

; varave_xy(1:2,:,:)=varave_xy(3:4,:,:)
;   dady_xy(1:2,:,:)=dady_xy(3:4,:,:)
; varave_yz(1:2,:,:)=varave_yz(3:4,:,:)
;   dady_yz(1:2,:,:)=dady_yz(3:4,:,:)

copy_VarMeta(temp(0,:,:), varave_xy(0,:,:))
copy_VarMeta(temp(0,:,:), dady_xy(0,:,:))
copy_VarMeta(temp(:,:,0), varave_yz(0,:,:))
copy_VarMeta(temp(:,:,0), dady_yz(0,:,:))
print("ok readdata")
return([/varave_xy,dady_xy,varave_yz,dady_yz/])
end
;--------------------------------------------------------------------------
undef("readdata1")
function readdata1(ncdir, filelist, varname, scale, yrStrt, yrLast, lonw, lone)
local dims1,temp1,f
begin
f= addfile(ncdir + filelist(0), "r")
temp1=f->$varname$(0,:,:) ; nlat,nlon
printVarSummary(temp1)
dims1=dimsizes(temp1)
varave_xy1=new((/5,dims1(0),dims1(1)/), float)

do i=0,2
f   = addfile(ncdir + filelist(i), "r")
YYYY= cd_calendar(f->time,-1)/100
iyr = ind(YYYY.ge.yrStrt(i) .and. YYYY.le.yrLast(i))
var00= f->$varname$(iyr,:,:)
printVarSummary(var00)
var0 = var00/scale
copy_VarMeta(var00, var0)

clm = clmMonTLL(var0) ; 12month,lev,lat,lon
; varanom= calcMonAnomTLLL(var, varclm)
varave_xy1(i,:,:) = dim_avg_n_Wrap(clm,0) ;lev,lat,lon

delete([/YYYY,iyr,var00,var0/])
end do

varave_xy1(3,:,:)=varave_xy1(2,:,:)-varave_xy1(0,:,:)
varave_xy1(4,:,:)=varave_xy1(1,:,:)-varave_xy1(0,:,:)

copy_VarMeta(temp1(:,:), varave_xy1(0,:,:))
print("ok readdata")
return([/varave_xy1/])
end

begin
;=========================================================================================
;=============================================================
; read SST data
;=============================================================
yrStrt  =  (/111,141,141/)
yrLast  =  (/140,170,170/)
ilev1   = 0
ilev2   = 5
lonw    = 200
lone    = 260
varname="TEMP"
ncdir  = "/Volumes/oly2/CESM/monthly/"
scale=1. ;degC
filelist=(/"TEMP.mon.0071-0140.20S_20N.120E_80W.levidx01_21_1degree.nc",\
           "TEMP.mon.2xCO2.0071-0170.20S_20N.120E_80W.levidx01_21_1degree.nc",\
           "TEMP.mon.4xCO2.0071-0170.20S_20N.120E_80W.levidx01_21_1degree.nc"/)
out=readdata(ncdir, filelist, varname, scale, yrStrt, yrLast, ilev1, ilev2, lonw, lone)
temp_xy=out[0]
dtdy_xy=out[1]
temp_yz=out[2]
dtdy_yz=out[3]
printVarSummary(temp_yz)
;=============================================================
; read UVEL data
;=============================================================
yrStrt  =  (/111,141,141/)
yrLast  =  (/140,170,170/)
ilev1   = 0
ilev2   =5
lonw    =210
lone    =250
varname="UVEL"
ncdir  = "/Volumes/oly2/CESM/monthly/"
scale=100. ;units cm/s to m/s
filelist=(/"UVEL.mon.0071-0140.20S_20N.120E_80W.levidx01_21_1degree.nc",\
           "UVEL.mon.2xCO2.0071-0170.20S_20N.120E_80W.levidx01_21_1degree.nc",\
           "UVEL.mon.4xCO2.0071-0170.20S_20N.120E_80W.levidx01_21_1degree.nc"/)
out:=readdata(ncdir, filelist, varname, scale, yrStrt, yrLast, ilev1, ilev2, lonw, lone)
uvel_xy=out[0]
dudy_xy=out[1]
uvel_yz=out[2]
dudy_yz=out[3]
uvel_xy@units="m/s"
dudy_xy@units="m/s"
uvel_yz@units="m/s"
dudy_yz@units="m/s"
printVarSummary(dudy_xy)
dudy_xy_new = dudy_xy
dudy_xy_new(:,{-10:-2},{220:280})= dudy_xy@missing_value
dudy_xy_new(:,{7:10},{180:280})= dudy_xy@missing_value
dudy_xy_new(1:2,:,{180:195})= dudy_xy@missing_value
dudy_xy_new(3,{-10:-2},{180:280})= dudy_xy@missing_value

copy_VarCoords(dudy_xy, dudy_xy_new)

;=============================================================
; read TAUX data
;=============================================================
yrStrt  =  (/41,141,141/)
yrLast  =  (/70,170,170/)
ilev1   = 0
ilev2   =5
lonw    =210
lone    =250
varname="TAUX"
ncdir  = "/Volumes/oly2/CESM/monthly/"
scale=1. ;units cm/s to m/s
filelist=(/"TAUX.mon.0071-0140.20S_20N.120E_80W_1degree.nc",\
           "TAUX.2xCO2.mon.0071-0170.20S_20N.120E_80W_1degree.nc",\
           "TAUX.4xCO2.mon.0071-0170.20S_20N.120E_80W_1degree.nc"/)
out1 =readdata1(ncdir, filelist, varname, scale, yrStrt, yrLast, lonw, lone)
TAUX_xy=out1[0]


yrStrt  =  (/41,141,141/)
yrLast  =  (/70,170,170/)
ilev1   = 0
ilev2   =10
lonw    =210
lone    =250
varname="TAUY"
ncdir  = "/Volumes/oly2/CESM/monthly/"
scale=1. ;units cm/s to m/s
filelist=(/"TAUY.mon.0071-0140.20S_20N.120E_80W_1degree.nc",\
           "TAUY.2xCO2.mon.0071-0170.20S_20N.120E_80W_1degree.nc",\
           "TAUY.4xCO2.mon.0071-0170.20S_20N.120E_80W_1degree.nc"/)
out2=readdata1(ncdir, filelist, varname, scale, yrStrt, yrLast, lonw, lone)
TAUY_xy=out2[0]
;==========================================================================================
leftstring  = (/"~F22~a  PD SST & wind stress","~F22~b  2xCO2 SST & wind stress","~F22~c  4xCO2 SST & wind stress","~F22~d  4xCO2-PD"/) 
leftstring1 = (/"~F22~e   PD UVEL","~F22~f   2xCO2 UVEL","~F22~g   4xCO2 UVEL ","~F22~h   4xCO2-PD"/) 
leftstring2 = (/"~F22~i","~F22~j","~F22~k","~F22~l"/) 
text=(/"~F22~PD","~F22~2xCO2","~F22~4xCO2","~F22~4xCO2 -PD"/)
outdir = "/Users/xueaoyun/Documents/work/TIW_future/revision1/main_figure/"
wks = gsn_open_wks("pdf",outdir+"CESM_PD_2xCO2_4xCO2_zonal_current_gradient")
; gsn_define_colormap(wks,"GMT_wysiwygcont")

latS   =  -10.
latN   =   10.
lonL   =  180
lonR   =  280
;------------------>
res                      = True         
res@gsnDraw              = False
res@gsnFrame             = False
res@gsnAddCyclic         = False      
res@mpFillOn             = False
res@mpLandFillColor      = "gray"
res@mpGeophysicalLineColor= "gray30"
res@mpGeophysicalLineThicknessF = 0.7
res@mpMinLatF            = latS
res@mpMaxLatF            = latN
res@mpMinLonF            = lonL  
res@mpMaxLonF            = lonR
res@mpCenterLonF         = 200
res@mpShapeMode          = "FreeAspect" 

res@tmXTOn   = False
res@tmYROn   = False
res@tmXBMode              = "Explicit"
res@tmXBValues            = (/120,160,200,240,280/)
res@tmXBLabels            = (/"","160~S~o~N~E","160~S~o~N~W","120~S~o~N~W","80~S~o~N~W"/)
res@tmXBMinorValues       = ispan(120,280,20)
res@tmYLMode              = "Explicit"
res@tmYLValues            = (/-10,0,10/)
res@tmYLLabels            = (/"10~S~o~N~S","0","10~S~o~N~N"/)
res@tmYLMinorValues       = (/-5,5/)
res@tmXBLabelFontHeightF     = 0.010
res@tmYLLabelFontHeightF     = 0.010
; res@tmXBLabelFontThicknessF  = 1.5
; res@tmYLLabelFontThicknessF  = 1.5 
res@tmBorderThicknessF       = 2
res@tmXBMajorLengthF         = -0.006
res@tmYLMajorLengthF         = -0.006
res@tmXBMinorLengthF         = -0.003
res@tmYLMinorLengthF         = -0.003 
res@tmXBMajorThicknessF      = 2.0
res@tmYLMajorThicknessF      = 2.0
res@tmXBMinorThicknessF      = 1.5
res@tmYLMinorThicknessF      = 1.5 
; res@tmXBLabelDeltaF          = 0.01
; res@tmYLLabelDeltaF          = 0.2
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
res@cnFillDrawOrder      = "PreDraw"
res@cnLinesOn            = False
res@cnLineLabelsOn       = False
res@cnInfoLabelOn        = False
res@cnSmoothingOn        = True
res@cnSmoothingTensionF  = 0.05
res@cnLevelSelectionMode = "ExplicitLevels"
res@lbLabelBarOn         = True
res@lbOrientation        = "Vertical"
res@lbBoxEndCapStyle     = "TriangleBothEnds"
res@lbLabelFontHeightF   = 0.010
res@lbBoxSeparatorLinesOn= False  
res@lbBoxLinesOn         = False
res@gsnRightString       = ""
res@gsnStringFontHeightF = 0.010

;------------------>
sres                      = True
sres@gsnDraw              = False
sres@gsnFrame             = False


sres@cnFillOn             = False
sres@cnInfoLabelOn        = False
sres@cnLinesOn            = True
sres@cnSmoothingOn        = True
sres@cnSmoothingTensionF  = 0.5
sres@cnLineColor          = "black"  
sres@cnLevelSelectionMode = "ExplicitLevels"
sres@cnLevels             = (/-1.5,1.5/)
sres@cnLineThicknessF     = 2
 sres@cnLineLabelDensityF=-0.2
 sres@cnLineLabelInterval  = 1
sres@cnLineLabelsOn       = True
sres@cnLineLabelFontColor = "black"
sres@cnLineLabelFontHeightF = 0.015
sres@cnLineLabelFontThicknessF = 2.0
sres@cnLineLabelPlacementMode = "computed"
sres@cnLineLabelFormat    = "@^sg"
; sres@cnLineLabelFontAspectF = 1.5
sres@cnLineLabelPerimColor = "transparent"
sres@cnLineLabelBackgroundColor = "white"

sres@cnSmoothingOn        = True
sres@cnSmoothingTensionF  = 3
sres@gsnContourNegLineDashPattern = 16
sres@cnLineDashSegLenF=0.1
sres@gsnLeftString       = ""
sres@gsnRightString      = ""
sres@gsnStringFontHeightF = 0.010

vcres                         = True
vcres@gsnDraw                 = False
vcres@gsnFrame                = False
vcres@vcGlyphStyle         = "LineArrow"  
vcres@vcFillArrowsOn            = True
vcres@vcFillArrowEdgeThicknessF = 1.5
vcres@vcLineArrowColor          = "black"
vcres@vcFillArrowFillColor      = "black"
vcres@vcRefAnnoBackgroundColor  = -1
vcres@vcFillArrowWidthF         = 0.05
vcres@vcFillArrowHeadInteriorXF = 0.2
vcres@vcFillArrowHeadXF         = 0.3
vcres@vcFillArrowHeadYF         = 0.15  
vcres@vcRefMagnitudeF           = 0.05                ; define vector ref mag
vcres@vcRefLengthF              = 0.015                ; define length of vec ref
vcres@vcMinDistanceF            = 0.03                ; thin out vectors
vcres@vcRefAnnoOrthogonalPosF = -1.0            ; move ref vector
vcres@vcRefAnnoParallelPosF   = 1.00
vcres@vcRefAnnoArrowSpaceF    = 1.5
vcres@vcRefAnnoFontColor      = "black"
;vcres@vcMonoLineArrowColor    = True                ; default 
vcres@vcRefAnnoFontHeightF    = 0.012
vcres@vcRefAnnoFontThicknessF = 1.0
vcres@vcRefAnnoOn          = True
vcres@vcRefAnnoPerimOn     = True
vcres@vcRefAnnoBackgroundColor = 0
vcres@vcRefAnnoPerimSpaceF = 0.05
vcres@vcRefAnnoArrowSpaceF = 0.1
vcres@vcLineArrowThicknessF= 1.0
vcres@vcRefAnnoString2On   = True
vcres@vcRefAnnoString2     = "dyne/cm~S~2~N"
vcres@vcVectorDrawOrder     = "PostDraw"
vcres@gsnLeftString = ""
vcres@gsnRightString=""

resp=True
resp@gsnFrame   =False
resp@gsnMaximize=False
resp@gsnPanelLabelBar  =True
resp@lbBoxLinesOn      =False
resp@lbOrientation     ="vertical"
resp@lbLabelFontHeightF=0.007
;=============================================================
; plot
;=============================================================
plot =new(4, graphic)
contour = new(4, graphic)
dumb =plot
plot1=plot
dumb1=plot
vector = plot
;res@cnFillPalette="GMT_wysiwygcont"
res@cnFillPalette="BlueWhiteOrangeRed"
res@gsnSpreadColorStart=15
res@gsnSpreadColorEnd  =245
res@vpHeightF            = 0.1
res@vpWidthF             = 0.23

; res@cnFillColors       = (/20,40,60,80,110,170,190,210,220,230/)
; res@cnLevels             = fspan(-0.6,0.6,31)
;--------------->
;--------------->
res@lbLabelBarOn   = False
res@tiYAxisOffsetXF=-0.03
res@tiYAxisFontHeightF=0.013
res@tiYAxisString=text(0)
res@vpXF           = 0.08
res@vpYF           = 0.90
res@gsnLeftString  = ""
res@cnLevels       = fspan(24,34,101)
vcres@vcRefMagnitudeF           = 0.5 
plot(0) = gsn_csm_contour_map_ce(wks,temp_xy(0,:,:),res)
vector(0) = gsn_csm_vector(wks, TAUX_xy(0,:,:), TAUY_xy(0,:,:),vcres)
overlay(plot(0), vector(0))

res@tiYAxisString=text(1)
res@vpXF           = 0.08
res@vpYF           = 0.75
res@gsnLeftString  = ""
res@cnLevels       = fspan(24,34,101)
vcres@vcRefMagnitudeF           = 0.5 
res@lbLabelBarOn         = True
res@pmLabelBarWidthF     = 0.03
res@pmLabelBarHeightF    = 0.38
res@pmLabelBarParallelPosF   = 0.5     ; 0.5 This is the right edge of the plot.
res@pmLabelBarOrthogonalPosF = 0.065    ; 0.5 This is the bottom edge of the plot.
res@lbLabelStride     = 10
vcres@vcRefAnnoOrthogonalPosF = -2.41            ; move ref vector
plot(1) = gsn_csm_contour_map_ce(wks,temp_xy(1,:,:),res)
vector(1) = gsn_csm_vector(wks, TAUX_xy(1,:,:), TAUY_xy(1,:,:),vcres)
overlay(plot(1), vector(1))

res@tiYAxisString=text(2)
res@vpXF           = 0.08
res@vpYF           = 0.60
res@gsnLeftString  = ""
res@cnLevels       = fspan(24,34,101)
vcres@vcRefMagnitudeF           = 0.5 
res@lbLabelBarOn         = False
res@pmLabelBarWidthF     = 0.03
res@pmLabelBarHeightF    = 0.38
res@pmLabelBarParallelPosF   = 1.5     ; 0.5 This is the right edge of the plot.
res@pmLabelBarOrthogonalPosF = 0.05    ; 0.5 This is the bottom edge of the plot.
vcres@vcRefAnnoOrthogonalPosF = -1.0            ; move ref vector
plot(2) = gsn_csm_contour_map_ce(wks,temp_xy(2,:,:),res)
vector(2) = gsn_csm_vector(wks, TAUX_xy(2,:,:), TAUY_xy(2,:,:),vcres)
overlay(plot(2), vector(2))
delete(res@cnLevels)
res@tiYAxisString=text(3)
res@vpXF           = 0.08
res@vpYF           = 0.43
res@gsnLeftString  = ""
res@cnLevels       = fspan(0,6,61)
vcres@vcRefMagnitudeF           = 0.3
res@lbLabelBarOn         = True
res@pmLabelBarWidthF     = 0.03
res@pmLabelBarHeightF    = 0.10
res@pmLabelBarParallelPosF   = 0.5    ; 0.5 This is the right edge of the plot.
res@pmLabelBarOrthogonalPosF = 0.06    ; 0.5 This is the bottom edge of the plot.
vcres@vcRefAnnoOrthogonalPosF = -1.03            ; move ref vector
plot(3) = gsn_csm_contour_map_ce(wks,temp_xy(3,:,:),res)
vector(3) = gsn_csm_vector(wks, TAUX_xy(3,:,:), TAUY_xy(3,:,:),vcres)
overlay(plot(3), vector(3))

txres                       = True
txres@txPerimOn             = False  
txres@txBackgroundFillColor = "White"
txres@txFontHeightF         = 0.012
txres@txFontThicknessF      = 1.5
txid1    = gsn_create_text(wks, leftstring(0), txres)
txid2    = gsn_create_text(wks, leftstring(1), txres)
txid3    = gsn_create_text(wks, leftstring(2), txres)
txid4    = gsn_create_text(wks, leftstring(3), txres)

amres                  = True
 amres@amParallelPosF   = -0.485     ; 0.5 This is the right edge of the plot.
 amres@amOrthogonalPosF = -0.74     ; 0.5 This is the bottom edge of the plot.
 amres@amJust           = "TopLeft"

 annoid1 = gsn_add_annotation(plot(0), txid1, amres)
 annoid2 = gsn_add_annotation(plot(1), txid2, amres)
 annoid3 = gsn_add_annotation(plot(2), txid3, amres)
 annoid4 = gsn_add_annotation(plot(3), txid4, amres)


;do i=0,3
;res@tiYAxisString=text(i)
;res@gsnLeftString  = leftstring(i)
;if i.le.2 then
;res@cnLevels       = fspan(24,34,101)
;vcres@vcRefMagnitudeF           = 0.5 
;else
;delete(res@cnLevels)
;res@cnLevels       = fspan(0,6,61)
;vcres@vcRefMagnitudeF           = 0.3   
;end if
;
;plot(i) = gsn_csm_contour_map_ce(wks,temp_xy(i,:,:),res)
;; dumb(i) = gsn_csm_contour(wks, dtdy_xy(i,:,:), sres)
;; overlay(plot(i), dumb(i))
;vector(i) = gsn_csm_vector(wks, TAUX_xy(i,:,:), TAUY_xy(i,:,:),vcres)
;overlay(plot(i), vector(i))
;print("i="+i)
;end do
draw(plot(0))
draw(plot(1))
draw(plot(2))
draw(plot(3))
delete(res@cnLevels)
;--------------->
res@lbLabelBarOn   = False
res@tiYAxisString=""
res@vpXF           = 0.40
res@vpYF           = 0.90
res@gsnLeftString  = ""
res@cnLevels       = fspan(-0.5,0.5,101)
plot1(0) = gsn_csm_contour_map_ce(wks,uvel_xy(0,:,:),res)
contour(0)  = gsn_csm_contour(wks, dudy_xy_new(0,:,:), sres)
overlay(plot1(0), contour(0))

res@vpXF           = 0.40
res@vpYF           = 0.75
res@gsnLeftString  =  ""
res@cnLevels       = fspan(-0.5,0.5,101)
res@lbLabelBarOn         = True
res@pmLabelBarWidthF     = 0.03
res@pmLabelBarHeightF    = 0.38
res@pmLabelBarParallelPosF   = 0.5     ; 0.5 This is the right edge of the plot.
res@pmLabelBarOrthogonalPosF = 0.065    ; 0.5 This is the bottom edge of the plot.
plot1(1) = gsn_csm_contour_map_ce(wks,uvel_xy(1,:,:),res)
contour(1)  = gsn_csm_contour(wks, dudy_xy_new(1,:,:), sres)
overlay(plot1(1), contour(1))

res@vpXF           = 0.40
res@vpYF           = 0.60
res@gsnLeftString  = ""
res@cnLevels       = fspan(-0.5,0.5,101)
res@lbLabelBarOn         = False 
res@pmLabelBarWidthF     = 0.03
res@pmLabelBarHeightF    = 0.38
res@pmLabelBarParallelPosF   = 1.5     ; 0.5 This is the right edge of the plot.
res@pmLabelBarOrthogonalPosF = 0.06    ; 0.5 This is the bottom edge of the plot.
plot1(2) = gsn_csm_contour_map_ce(wks,uvel_xy(2,:,:),res)
contour(2)  = gsn_csm_contour(wks, dudy_xy_new(2,:,:), sres)
overlay(plot1(2), contour(2))

delete(res@cnLevels)
delete(sres@cnLevels)
res@vpXF           = 0.40
res@vpYF           = 0.43
res@gsnLeftString  = ""
res@cnLevels       = fspan(-0.3,0.3,61)
sres@cnLevels             = (/-0.4/)
res@lbLabelBarOn         = True
res@pmLabelBarWidthF     = 0.03
res@pmLabelBarHeightF    = 0.10
res@pmLabelBarParallelPosF   = 0.5     ; 0.5 This is the right edge of the plot.
res@pmLabelBarOrthogonalPosF = 0.065    ; 0.5 This is the bottom edge of the plot.
plot1(3) = gsn_csm_contour_map_ce(wks,uvel_xy(3,:,:),res)
contour(3)  = gsn_csm_contour(wks, dudy_xy_new(3,:,:), sres)
overlay(plot1(3), contour(3))

txres                       = True
txres@txPerimOn             = False  
txres@txBackgroundFillColor = "White"
txres@txFontHeightF         = 0.012
txres@txFontThicknessF      = 1.5
txid01    = gsn_create_text(wks, leftstring1(0), txres)
txid02    = gsn_create_text(wks, leftstring1(1), txres)
txid03    = gsn_create_text(wks, leftstring1(2), txres)
txid04    = gsn_create_text(wks, leftstring1(3), txres)

amres                  = True
 amres@amParallelPosF   = -0.485     ; 0.5 This is the right edge of the plot.
 amres@amOrthogonalPosF = -0.74     ; 0.5 This is the bottom edge of the plot.
 amres@amJust           = "TopLeft"

 annoid01 = gsn_add_annotation(plot1(0), txid01, amres)
 annoid02 = gsn_add_annotation(plot1(1), txid02, amres)
 annoid03 = gsn_add_annotation(plot1(2), txid03, amres)
 annoid04 = gsn_add_annotation(plot1(3), txid04, amres)
draw(plot1(0))
draw(plot1(1))
draw(plot1(2))
draw(plot1(3))
;--------------->
;res@tiYAxisString=""
;res@vpXF           = 0.40
;res@lbLabelBarOn   = False
;do i=0,2
;printMinMax(uvel_xy(i,:,:), False)
;res@gsnLeftString  = leftstring1(i)
;res@cnLevels       = fspan(-0.5,0.5,101)
;plot1(i) = gsn_csm_contour_map_ce(wks,uvel_xy(i,:,:),res)
;contour(i)  = gsn_csm_contour(wks, dudy_xy_new(i,:,:), sres)
;overlay(plot1(i), contour(i))
;print("i="+i)
;end do
;delete(res@cnLevels)
;res@cnLevels       = fspan(-0.3,0.3,61)
;delete(sres@cnLevels)
;sres@cnLevels             = (/-0.4/)
;plot1(3) = gsn_csm_contour_map_ce(wks,uvel_xy(3,:,:),res)
;contour(3)  = gsn_csm_contour(wks, dudy_xy_new(3,:,:), sres)
;overlay(plot1(3), contour(3))
;;--------------->
;resp@pmLabelBarWidthF=0.04
;resp@pmLabelBarOrthogonalPosF=-0.002
;resp@gsnPanelYWhiteSpacePercent=4
;resp@gsnPanelMainString =""
;
;;;--------------->
;;; panel sst
;;;--------------->
;;resp@gsnPanelLeft =0.01
;;resp@gsnPanelRight=resp@gsnPanelLeft+0.34
;;resp@gsnPanelTop  =0.95
;;resp@gsnPanelBottom     =0.95-0.35*2
;;resp@pmLabelBarHeightF  =0.3
;;gsn_panel(wks, plot(0:2), (/3,1/), resp)
;;resp@gsnPanelTop        =0.95-0.24
;;resp@gsnPanelBottom     =0.01
;;resp@pmLabelBarHeightF  =0.08
;;gsn_panel(wks, plot(3), (/1,1/), resp)
; txres                       = True
; txres@txPerimOn             = False  
; txres@txBackgroundFillColor = "transparent"
; txres@txFontHeightF         = 0.015
; txres@txFontThicknessF      = 4
; txres@txFont                = 22
; txid1     = gsn_create_text(wks, "Mean States Change", txres)
; 
; amres                  = True
; amres@amParallelPosF   = -0.80    ; 0.5 This is the right edge of the plot.
; amres@amOrthogonalPosF = -1.1    ; 0.5 This is the bottom edge of the plot.
; amres@amJust           = "TopLeft"
; 
; ;gsn_text_ndc(wks, "Mean States Change", 0.50, 0.80, txres)
; ;--------------->
; ; panel uvel
; ;--------------->
; resp@gsnPanelLeft =0.34
; resp@gsnPanelRight=resp@gsnPanelLeft+0.34
; resp@gsnPanelTop  =0.95
; resp@gsnPanelBottom     =0.95-0.35*2
; resp@pmLabelBarHeightF  =0.3
; gsn_panel(wks, plot1(0:2), (/3,1/), resp)
; resp@gsnPanelTop        =0.95-0.24
; resp@gsnPanelBottom     =0.01
; resp@pmLabelBarHeightF  =0.08
; gsn_panel(wks, plot1(3), (/1,1/), resp)
;========================================================================
ares                      = True         
ares@gsnDraw              = False
ares@gsnFrame             = False
ares@gsnAddCyclic         = False     
;ares@gsnScalarContour     = True

ares@tmXTOn   = False
ares@tmYROn   = False
ares@tmYLOn   = True
ares@tmYLMode = "Explicit"
ares@tmYLValues = (/500,5000,10000,15000,20000/)
ares@tmYLLabels = (/"5","50","100","150","200"/)
ares@tmXBMode = "Explicit"
ares@tmXBValues = (/-10,-5,0,5,10/)
ares@tmXBLabels = (/"10~S~o~N~S","","Eq.","","10~S~o~N~N"/)  

ares@tmXBLabelFontHeightF  = 0.010
ares@tmYLLabelFontHeightF  = 0.010
ares@tmXBLabelFontThicknessF  = 1.5
ares@tmYLLabelFontThicknessF  = 1.5 
ares@tmBorderThicknessF       = 2.0
ares@tmXBMajorLengthF        = -0.006
ares@tmYLMajorLengthF        = -0.006
ares@tmXBMinorLengthF        = -0.003
ares@tmYLMinorLengthF        = -0.003 
ares@tmXBMajorThicknessF     = 1.5
ares@tmYLMajorThicknessF     = 1.5
ares@tmXBMinorThicknessF     = 1.0
ares@tmYLMinorThicknessF     = 1.0 
; ares@tmXBLabelDeltaF         = 0.01
; ares@tmYLLabelDeltaF         = 0.2
ares@tmXMajorGrid                = False
ares@tmXMajorGridLineDashPattern =  2
ares@tmXMajorGridLineColor       = "grey30"
ares@tmXMajorGridThicknessF      = 0.8
ares@tmYMajorGrid                = True
ares@tmYMajorGridLineDashPattern =  2
ares@tmYMajorGridLineColor       = "grey30"
ares@tmYMajorGridThicknessF      = 0.8 
ares@tmYMinorGrid                = False   
ares@tmYMinorGridLineDashPattern =  2
ares@tmYMinorGridLineColor       = "grey30"
ares@tmYMinorGridThicknessF      = 0.8 
ares@tmXMinorGrid                = False
ares@tmXMinorGridLineDashPattern =  2
ares@tmXMinorGridLineColor       = "grey30"
ares@tmXMinorGridThicknessF      = 0.8 
ares@trYReverse               = True             ; reverses y-axis
ares@gsnYAxisIrregular2Linear = True
ares@trXMinF = -10
ares@trXMaxF =  10

ares@cnFillOn             = True
ares@cnLinesOn            = False
ares@cnInfoLabelOn        = False
ares@cnLineLabelsOn       = False
;ares@cnMonoLineColor      = True
ares@cnLevelSelectionMode = "ExplicitLevels"  
;ares@cnLevels             =(/-0.10,-0.08,-0.06,-0.04,-0.02,0.02,0.04,0.06,0.08,0.10/)
;ares@cnFillColors         = (/4,5,6,7,8,0,11,12,13,14,15/)
ares@cnSmoothingOn        = True
ares@cnSmoothingTensionF  = -0.2
ares@gsnContourNegLineDashPattern = 4
ares@gsnContourZeroLineThicknessF = 2
ares@gsnRightString       = ""
ares@gsnStringFontHeightF  = 0.015
ares@cnFillPalette="BlueWhiteOrangeRed"

ares@lbLabelBarOn          = True
ares@lbOrientation      = "Vertical"
ares@lbBoxEndCapStyle     = "TriangleBothEnds"
ares@pmLabelBarWidthF     = 0.3
ares@pmLabelBarHeightF    = 0.05
ares@pmLabelBarOrthogonalPosF = 0.01
ares@pmLabelBarParallelPosF   = 0.5
ares@lbLabelFontHeightF       = 0.010
ares@lbBoxLinesOn             = False 
ares@lbLabelStride     = 20
;------------------>
bres                      = True
bres@gsnDraw              = False
bres@gsnFrame             = False
bres@gsnLeftString  = ""
bres@gsnRightString = ""

bres@cnInfoLabelOn        = False
bres@cnLinesOn            = True
bres@cnFillOn             = False
bres@cnLineColor          = "grey20"  
bres@cnLevelSelectionMode = "ExplicitLevels"
bres@cnLineLabelsOn       = True
bres@cnLineLabelFontColor = "grey20"
; bres@cnLineLabelCount=0.1
bres@cnLineLabelPlacementMode = "constant"
bres@cnLineLabelFormat    = "@^sg"
bres@cnLineLabelFontHeightF = 0.024
bres@cnLineLabelFontThicknessF = 2.0
bres@cnLineLabelFontAspectF = 1.5
bres@cnLineLabelInterval  = 2
bres@cnLineLabelPerimColor = "white"
bres@cnLineLabelPerimOn=True
bres@cnLineLabelBackgroundColor = "white"
; bres@cnLineDashSegLenF=0.1
bres@cnSmoothingOn        = True
bres@cnSmoothingTensionF  = -0.2

bres@cnLineThicknessF     = 1
bres@gsnContourNegLineDashPattern = 1
; bres@gsnContourZeroLineThicknessF = 2
bres@gsnContourPosLineDashPattern = 0

;=============================================================
; plot
;=============================================================
plot2 =new(4, graphic)
dumb2 =new(4, graphic)

; ares@cnFillColors         =(/20,40,60,80,110,170,190,210,220,230/)

ares@vpHeightF            = 0.1
ares@vpWidthF             = 0.15
ares@lbLabelBarOn   = False
ares@cnLevels             = fspan(-1.0, 1.0, 201)
bres@cnLevels             = fspan(14,28,8)
;--------------->
res@tiYAxisString=""
ares@vpXF           = 0.74
ares@vpYF           = 0.90
ares@gsnLeftString  = ""
plot2(0) = gsn_csm_contour(wks, uvel_yz(0,:,:),ares) 
dumb2(0) = gsn_csm_contour(wks, temp_yz(0,:,:),bres)
overlay(plot2(0), dumb2(0))

ares@vpXF           = 0.74
ares@vpYF           = 0.75
ares@gsnLeftString  = ""
ares@lbLabelBarOn         = True
ares@pmLabelBarWidthF     = 0.03
ares@pmLabelBarHeightF    = 0.38
ares@pmLabelBarParallelPosF   = 0.5     ; 0.5 This is the right edge of the plot.
ares@pmLabelBarOrthogonalPosF = 0.06    ; 0.5 This is the bottom edge of the plot.
plot2(1) = gsn_csm_contour(wks, uvel_yz(1,:,:),ares) 
dumb2(1) = gsn_csm_contour(wks, temp_yz(1,:,:),bres)
overlay(plot2(1), dumb2(1))
delete(ares@cnLevels)
delete(bres@cnLevels)
ares@cnLevels             = fspan(-1.0, 1.0, 201)
bres@cnLevels             = fspan(18,32,8)
ares@vpXF           = 0.74
ares@vpYF           = 0.60
ares@gsnLeftString  = ""
ares@lbLabelBarOn         = False 
ares@pmLabelBarWidthF     = 0.03
ares@pmLabelBarHeightF    = 0.38
ares@pmLabelBarParallelPosF   = 1.5     ; 0.5 This is the right edge of the plot.·
ares@pmLabelBarOrthogonalPosF = 0.06    ; 0.5 This is the bottom edge of the plot.
plot2(2) = gsn_csm_contour(wks, uvel_yz(2,:,:),ares) 
dumb2(2) = gsn_csm_contour(wks, temp_yz(2,:,:),bres)
overlay(plot2(2), dumb2(2))
delete(ares@cnLevels)
delete(bres@cnLevels)
ares@cnLevels             = fspan(-0.3, 0.3, 61)
bres@cnLevels             = (/-2,-1,0,1,2,3,4/)
ares@vpXF           = 0.74
ares@vpYF           = 0.43
ares@gsnLeftString  = ""
ares@lbLabelBarOn         = True
ares@pmLabelBarWidthF     = 0.03
ares@pmLabelBarHeightF    = 0.10
ares@pmLabelBarParallelPosF   = 0.5     ; 0.5 This is the right edge of the plot.
ares@pmLabelBarOrthogonalPosF = 0.06    ; 0.5 This is the bottom edge of the plot.
plot2(3) = gsn_csm_contour(wks, uvel_yz(3,:,:),ares) 
dumb2(3) = gsn_csm_contour(wks, temp_yz(3,:,:),bres)
overlay(plot2(3), dumb2(3))

txres                       = True
txres@txPerimOn             = False  
txres@txBackgroundFillColor = "White"
txres@txFontHeightF         = 0.0120
txres@txFontThicknessF      = 1.5
txid001    = gsn_create_text(wks, leftstring2(0), txres)
txid002    = gsn_create_text(wks, leftstring2(1), txres)
txid003    = gsn_create_text(wks, leftstring2(2), txres)
txid004    = gsn_create_text(wks, leftstring2(3), txres)

amres                  = True
 amres@amParallelPosF   = -0.485     ; 0.5 This is the right edge of the plot.
 amres@amOrthogonalPosF = -0.74     ; 0.5 This is the bottom edge of the plot.
 amres@amJust           = "TopLeft"

 annoid001 = gsn_add_annotation(plot2(0), txid001, amres)
 annoid002 = gsn_add_annotation(plot2(1), txid002, amres)
 annoid003 = gsn_add_annotation(plot2(2), txid003, amres)
 annoid004 = gsn_add_annotation(plot2(3), txid004, amres)


draw(plot2(0))
draw(plot2(1))
draw(plot2(2))
draw(plot2(3))
;====================================================
;do i=0,3
;ares@gsnLeftString =leftstring2(i)
;printMinMax(temp_yz(i,:,:), False)
;if i.le.1 then
;ares@cnLevels             = fspan(-1.0, 1.0, 201)
;bres@cnLevels             = fspan(14,28,8)
;else if i.eq. 2 then
;delete(ares@cnLevels)
;delete(bres@cnLevels)
;ares@cnLevels             = fspan(-1.0, 1.0, 201)
;bres@cnLevels             = fspan(18,32,8)
;else 
;delete(ares@cnLevels)
;delete(bres@cnLevels)
;ares@cnLevels             = fspan(-0.3, 0.3, 61)
;bres@cnLevels             = (/-2,-1,0,1,2,3,4/)
;end if 
;end if
;plot2(i) = gsn_csm_contour(wks, uvel_yz(i,:,:),ares) 
;dumb2(i) = gsn_csm_contour(wks, temp_yz(i,:,:),bres)
;overlay(plot2(i), dumb2(i))
;end do
;;--------------->
;; panel uvel
;;--------------->
;resp@lbLabelStride=10
;resp@gsnPanelLeft =0.66
;resp@gsnPanelRight=resp@gsnPanelLeft+0.24
;resp@gsnPanelTop  =0.95
;resp@gsnPanelBottom     =0.95-0.35*2
;resp@pmLabelBarHeightF  =0.3
;;resp@cnExplicitLabelBarLabelsOn = True
;;resp@lbLabelStrings = (/"-1.0","-0.9","-0.8","-0.7","-0.6","-0.5","-0.4","-0.3","-0.2","-0.1","0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0"/)
;gsn_panel(wks, plot2(0:2), (/3,1/), resp)
;resp@gsnPanelTop        =0.95-0.24
;resp@gsnPanelBottom     =0.01
;resp@pmLabelBarHeightF  =0.08
;gsn_panel(wks, plot2(3), (/1,1/), resp)


end