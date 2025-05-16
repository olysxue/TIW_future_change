undef("readdata")
function readdata(ncdir, filelist, varname, scale, yrStrt, yrLast, ilev1, ilev2, lonw, lone)
local dims,temp,f
begin
f= addfile(ncdir + filelist(0), "r")
temp=f->$varname$(0,:,:,:) ; nlev,nlat,nlon
printVarSummary(temp)
dims=dimsizes(temp)
varave_xy=new((/3,dims(1),dims(2)/), float)
dady_xy  =new((/3,dims(1),dims(2)/), double) ;5,nlat,nlon

varave_yz=new((/3,dims(0),dims(1)/), float)
dady_yz  =new((/3,dims(0),dims(1)/), double) ;5,nlev,nlat


do i=0,1
f   = addfile(ncdir + filelist(i), "r")
YYYY= cd_calendar(f->time,-1)/100
iyr = ind(YYYY.ge.yrStrt(i) .and. YYYY.le.yrLast(i))
var0= f->$varname$(iyr,:,:,:)
var = var0/scale
copy_VarMeta(var0, var)
printVarSummary(var)
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
dady   = center_finite_diff_n(varave, dY, False, 0, 1)*1000000.
copy_VarMeta(var(0,:,:,:), dady)
printVarSummary(dady)
;---->vertical mean 
varave_xy(i,:,:)= dim_avg_n_Wrap(varave(ilev1:ilev2,:,:),0)
dady_xy(i,:,:)  = dim_avg_n_Wrap(dady(ilev1:ilev2,:,:),0)

;---->longitude mean 
varave_yz(i,:,:) = dim_avg_n_Wrap(varave(:,:,{lonw:lone}),2)
dady_yz(i,:,:)= dim_avg_n_Wrap(dady(:,:,{lonw:lone}),2)

delete([/YYYY,iyr,var0,var,varave,dady/])
end do
; mean  
  varave_xy(2,:,:)=varave_xy(1,:,:)-varave_xy(0,:,:)
    dady_xy(2,:,:)=dady_xy(1,:,:)-dady_xy(0,:,:)
  varave_yz(2,:,:)=varave_yz(1,:,:)-varave_yz(0,:,:) 
    dady_yz(2,:,:)=dady_yz(1,:,:)-dady_yz(0,:,:)

;----------------------------------------------------
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
temp1=f->$varname$(0,:,0:179) ; nlat,nlon
printVarSummary(temp1)
dims1=dimsizes(temp1)
varave_xy01=new((/3,dims1(0),dims1(1)/), float)   ; mean 
do i=0,1
f   = addfile(ncdir + filelist(i), "r")
YYYY= cd_calendar(f->time,-1)/100
iyr = ind(YYYY.ge.yrStrt(i) .and. YYYY.le.yrLast(i))
var01 = f->$varname$(iyr,:,0:179)
printVarSummary(var01)
var1 = var01/scale
copy_VarMeta(var01, var1)
printVarSummary(var1)
varave_xy01(i,:,:) = dim_avg_n_Wrap(var1,0) ;lev,lat,lon
delete([/YYYY,iyr,var01,var1/])
end do
; mean 
varave_xy01(2,:,:)=varave_xy01(1,:,:)-varave_xy01(0,:,:)
copy_VarMeta(temp1(:,:), varave_xy01(0,:,:))
print("ok readdata")
return([/varave_xy01/])
end

begin
;=========================================================================================
;=============================================================
; read E3SMv1 SST data
;=============================================================
yrStrt  =  (/0001,0001,0001/)
yrLast  =  (/0080,0080,0080/)
ilev1   = 0
ilev2   = 4
lonw    = 200
lone    = 260
varname="thetao"
ncdir  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/E3SM/monthly/"
scale=1. ;degC
filelist=(/"thetao_Omon_E3SM-1-0_piControl_r1i1p1f1_gr_0001-0080_equator_0-200m.nc",\
         "thetao_Omon_E3SM-1-0_abrupt-4xCO2_r1i1p1f1_gr_0001-0080_equator_0-200m.nc"/)
out=readdata(ncdir, filelist, varname, scale, yrStrt, yrLast, ilev1, ilev2, lonw, lone)
; mean 
temp_xy1=out[0]
dtdy_xy1=out[1]
temp_yz1=out[2]
dtdy_yz1=out[3]
;=============================================================
; read UVEL data
;=============================================================
yrStrt  =  (/0001,0001,0001/)
yrLast  =  (/0080,0080,0080/)
ilev1   = 0
ilev2   =5
lonw    =210
lone    =250
varname="uo"
ncdir  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/E3SM/monthly/"
scale=1. ;units cm/s to m/s
filelist=(/"uo_Omon_E3SM-1-0_piControl_r1i1p1f1_gr_0001-0080_equator_0-200m.nc",\
           "uo_Omon_E3SM-1-0_abrupt-4xCO2_r1i1p1f1_gr_0001-0080_equator_0-200m.nc"/)
out:=readdata(ncdir, filelist, varname, scale, yrStrt, yrLast, ilev1, ilev2, lonw, lone)
uvel_xy1=out[0]
dudy_xy1=out[1]
uvel_yz1=out[2]
dudy_yz1=out[3]

uvel_xy1@units="m/s"
dudy_xy1@units="m/s"
uvel_yz1@units="m/s"
dudy_yz1@units="m/s"

dudy_xy1(0,{0:10},{120:155})= dudy_xy1@missing_value
;=============================================================
; read TAUX data
;=============================================================

yrStrt  =  (/0001,0001,0001/)
yrLast  =  (/0080,0080,0080/)
ilev1   = 0
ilev2   =5
lonw    =210
lone    =250
varname="tauu"
scale=1. ;units cm/s to m/s
ncdir  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/E3SM/monthly/"
filelist=(/"tauu_Amon_E3SM-1-0_piControl_r1i1p1f1_gr_0001-0080_equator.nc",\
           "tauu_Amon_E3SM-1-0_abrupt-4xCO2_r1i1p1f1_gr_0001-0080_equator.nc"/)
out1 =readdata1(ncdir, filelist, varname, scale, yrStrt, yrLast, lonw, lone)
TAUX_xy01=out1[0]
TAUX_xy1 = TAUX_xy01 *10.
copy_VarCoords(TAUX_xy01, TAUX_xy1)
yrStrt  =  (/0001,0001,0001/)
yrLast  =  (/0080,0080,0080/)
ilev1   = 0
ilev2   =5
lonw    =210
lone    =250
varname="tauv"
ncdir  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/E3SM/monthly/"
scale=1. ;units cm/s to m/s
filelist=(/"tauv_Amon_E3SM-1-0_piControl_r1i1p1f1_gr_0001-0080_equator.nc",\
           "tauv_Amon_E3SM-1-0_abrupt-4xCO2_r1i1p1f1_gr_0001-0080_equator.nc"/)
out2=readdata1(ncdir, filelist, varname, scale, yrStrt, yrLast, lonw, lone)
TAUY_xy01=out2[0]
TAUY_xy1 = TAUY_xy01 *10.
copy_VarCoords(TAUY_xy01, TAUY_xy1)
;==========================================================================================
;E3SMv2
yrStrt  =  (/00101,00101,00101/)
yrLast  =  (/15012,15012,15012/)
ilev1   = 0
ilev2   = 5
lonw    = 200
lone    = 260
varname="thetao"
ncdir  = "/Users/xueaoyun/Dropbox/TIW_future/data/E3SM/E3SMv2/monthly/"
scale=1. ;degC
filelist=(/"v2_piControl_temperature.0001-0150_equator_0-200m_new.nc",\
           "thetao_Omon_E3SM-2-0_abrupt-4xCO2_r2i1p1f1_gr_0001-0160_equator.nc"/)
out=readdata(ncdir, filelist, varname, scale, yrStrt, yrLast, ilev1, ilev2, lonw, lone)
; mean 
temp_xy2=out[0]
dtdy_xy2=out[1]
temp_yz2=out[2]
dtdy_yz2=out[3]
;=============================================================
; read UVEL data
;=============================================================
yrStrt  =  (/00101,00101,00101/)
yrLast  =  (/15012,15012,15012/)
ilev1   = 0
ilev2   =5
lonw    =210
lone    =250
varname="uo"
ncdir  = "/Users/xueaoyun/Dropbox/TIW_future/data/E3SM/E3SMv2/monthly/"
scale=1. ;units cm/s to m/s
filelist=(/"v2_piControl_uo_0001-0150_equator_0-200m_new.nc",\
           "uo_Omon_E3SM-2-0_abrupt-4xCO2_r2i1p1f1_gr_0001-0150_equator_0-200m.nc"/)
out:=readdata(ncdir, filelist, varname, scale, yrStrt, yrLast, ilev1, ilev2, lonw, lone)
uvel_xy2=out[0]
dudy_xy2=out[1]
uvel_yz2=out[2]
dudy_yz2=out[3]

uvel_xy2@units="m/s"
dudy_xy2@units="m/s"
uvel_yz2@units="m/s"
dudy_yz2@units="m/s"

dudy_xy2(0,{0:10},{120:155})= dudy_xy2@missing_value
;=============================================================
; read TAUX data
;=============================================================

yrStrt  =  (/00101,00101,00101/)
yrLast  =  (/15012,15012,15012/)
ilev1   = 0
ilev2   =5
lonw    =210
lone    =250
varname="tauu"
scale=1. ;units cm/s to m/s
ncdir  = "/Users/xueaoyun/Dropbox/TIW_future/data/E3SM/E3SMv2/monthly/"
filelist=(/"v2.FV1.piControl.eam.h0.TAUX.0001-0500_1degree_equator_new.nc",\
           "tauu_Amon_E3SM-2-0_abrupt-4xCO2_r2i1p1f1_gr_0001-0160_equator.nc"/)
out1 =readdata1(ncdir, filelist, varname, scale, yrStrt, yrLast, lonw, lone)
TAUX_xy02=out1[0]
TAUX_xy2 = TAUX_xy02 *10.
copy_VarCoords(TAUX_xy02, TAUX_xy2)

yrStrt  =  (/00101,00101,00101/)
yrLast  =  (/15012,15012,15012/)
ilev1   = 0
ilev2   =5
lonw    =210
lone    =250
varname="tauv"
ncdir  = "/Users/xueaoyun/Dropbox/TIW_future/data/E3SM/E3SMv2/monthly/"
scale=1. ;units cm/s to m/s
filelist=(/"v2.FV1.piControl.eam.h0.TAUY.0001-0500_remap_equator_new.nc",\
           "tauv_Amon_E3SM-2-0_abrupt-4xCO2_r2i1p1f1_gr_0001-0160_equator.nc"/)
out2=readdata1(ncdir, filelist, varname, scale, yrStrt, yrLast, lonw, lone)
TAUY_xy02=out2[0]
TAUY_xy2 = TAUY_xy02 *10.
copy_VarCoords(TAUY_xy02, TAUY_xy2)
TAUX_xy2(2,:,:) = TAUX_xy2(1,:,:) + TAUX_xy2(0,:,:)
TAUY_xy2(2,:,:) = TAUY_xy2(1,:,:) + TAUY_xy2(0,:,:)
copy_VarCoords(TAUY_xy2(0,:,:), TAUY_xy2(2,:,:))
copy_VarCoords(TAUX_xy2(0,:,:), TAUX_xy2(2,:,:))
TAUX_xy2(0,:,:) = -1*TAUX_xy2(0,:,:) 
TAUY_xy2(0,:,:) = -1*TAUY_xy2(0,:,:) 
copy_VarCoords(TAUY_xy2(2,:,:), TAUY_xy2(0,:,:))
copy_VarCoords(TAUX_xy2(2,:,:), TAUX_xy2(0,:,:))
;====================================================================================
leftstring0  = (/"~F22~a","~F22~b"/) 
leftstring1  = (/"~F22~c","~F22~d"/) 
leftstring2  = (/"~F22~e","~F22~f"/) 

leftstring00  = (/"~F22~g","~F22~h"/) 
leftstring01  = (/"~F22~i","~F22~j"/) 
leftstring02  = (/"~F22~k","~F22~l"/) 

Centerstring0  = (/"PI SST & Wind Stress","4xCO2-PI"/) 
Centerstring1 = (/"PI UVEL","4xCO2-PI"/) 

Centerstring00  = (/"PI SST & Wind Stress","4xCO2-PI"/) 
Centerstring01 = (/"PI UVEL","4xCO2-PI"/) 


text=(/"~F22~PI","~F22~4xCO2 -PI"/)
outdir = "/Users/xueaoyun/Documents/work/TIW_future/revision1/Extend_figure/"
wks = gsn_open_wks("pdf",outdir+"Extend_Figure5_E3SM1_E3SM2_PD_4xCO2_zonal_current_gradient")
; gsn_define_colormap(wks,"GMT_wysiwygcont")

latS   =  -10.
latN   =   10.
lonL   =  120
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
res@gsnRightStringFontHeightF  = 0.015

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
res@gsnStringFontHeightF = 0.012
res@lbBoxLinesOn      =False
res@lbLabelFontHeightF=0.007
res@lbLabelStride     = 10
res@pmLabelBarWidthF     = 0.03
res@pmLabelBarHeightF    = 0.10
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
sres@cnLevels             = (/-1.0,1.0/)
sres@cnLineThicknessF     = 2
sres@cnLineLabelDensityF=-0.3
sres@cnLineLabelInterval  = 1
;sres@cnLineDashSegLenF    = 0.5
sres@cnLineLabelsOn       = True
sres@cnLineLabelFontColor = "black"
sres@cnLineLabelFontHeightF = 0.018
sres@cnLineLabelFontThicknessF = 2.0
sres@cnLineLabelPlacementMode = "computed"
sres@cnLineLabelFormat    = "@^sg"
; sres@cnLineLabelFontAspectF = 1.5
sres@cnLineLabelPerimColor = "white"
sres@cnLineLabelBackgroundColor = "white"

sres@cnSmoothingOn        = True
sres@cnSmoothingTensionF  = 3
sres@gsnContourNegLineDashPattern = 16
sres@cnLineDashSegLenF=0.1
sres@gsnLeftString       = ""
sres@gsnRightString      = ""
sres@gsnStringFontHeightF = 0.010
;=====================================================
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
vcres@vcRefAnnoOrthogonalPosF = -1.30            ; move ref vector
vcres@vcRefAnnoParallelPosF   = 1.00
vcres@vcRefAnnoArrowSpaceF    = 1.5
vcres@vcRefAnnoFontColor      = "black"
;vcres@vcMonoLineArrowColor    = True                ; default 
vcres@vcRefAnnoFontHeightF    = 0.015
vcres@vcRefAnnoFontThicknessF = 1.0
vcres@vcRefAnnoOn          = True
vcres@vcRefAnnoPerimOn     = False
vcres@vcRefAnnoPerimSpaceF = 0.1
vcres@vcRefAnnoArrowSpaceF = 0.1
vcres@vcLineArrowThicknessF= 1.0
vcres@vcRefAnnoString2On   = True
vcres@vcRefAnnoString2     = "dyne/cm~S~2~N"
vcres@vcVectorDrawOrder     = "PostDraw"
vcres@gsnLeftString = ""
vcres@gsnRightString=""
;=============================================================
; plot
;=============================================================
;res@cnFillPalette="GMT_wysiwygcont"
res@cnFillPalette="BlueWhiteOrangeRed"
;res@gsnSpreadColorStart=15
;res@gsnSpreadColorEnd  =245
res@vpHeightF            = 0.1
res@vpWidthF             = 0.25
;--------------->
; plot SST 
res@vpXF           = 0.10
res@vpYF           = 0.90
res@lbLabelBarOn   = True
res@tiYAxisOffsetXF=-0.03
res@tiYAxisFontHeightF=0.012
res@tiYAxisString=text(0)
res@gsnLeftString  = leftstring0(0)
res@gsnCenterString = Centerstring0(0)
res@cnLevels       = fspan(18,32,71)
vcres@vcRefMagnitudeF           = 0.5 
plot0 = gsn_csm_contour_map_ce(wks,temp_xy1(0,:,:),res)
vector0 = gsn_csm_vector(wks, TAUX_xy1(0,:,:), TAUY_xy1(0,:,:),vcres)
overlay(plot0, vector0)
delete(res@cnLevels)
res@vpXF           = 0.10
res@vpYF           = 0.75
res@lbLabelBarOn   = True
res@tiYAxisOffsetXF=-0.03
res@tiYAxisFontHeightF=0.012
res@tiYAxisString= text(1)
res@gsnLeftString  = leftstring0(1)
res@gsnCenterString = Centerstring0(1)
res@cnLevels       = fspan(1,7,61)
vcres@vcRefMagnitudeF     = 0.3   
res@lbLabelStride     = 5
plot1 = gsn_csm_contour_map_ce(wks,temp_xy1(2,:,:),res)
vector1 = gsn_csm_vector(wks, TAUX_xy1(2,:,:), TAUY_xy1(2,:,:),vcres)
overlay(plot1, vector1)
;=========================================
; plot uvel
delete(res@cnLevels)
res@tmYLLabels            = (/"","",""/)
res@tiYAxisString=""
res@vpXF           = 0.40
res@vpYF           = 0.90
res@lbLabelBarOn   = True
res@tiYAxisString  = ""
res@gsnLeftString  = leftstring1(0)
res@gsnCenterString = Centerstring1(0)
res@cnLevels       = fspan(-0.5,0.5,101)
plot2  = gsn_csm_contour_map_ce(wks,uvel_xy1(0,:,:),res)
contour2  = gsn_csm_contour(wks, dudy_xy1(0,:,:), sres)
overlay(plot2, contour2)
delete(res@cnLevels)
res@tiYAxisString=""
res@vpXF           = 0.40
res@vpYF           = 0.75
res@lbLabelBarOn   = True
res@tiYAxisString  = ""
res@gsnLeftString  = leftstring1(1)
res@gsnCenterString = Centerstring1(1)
res@cnLevels       = fspan(-0.3,0.3,61)
plot3  = gsn_csm_contour_map_ce(wks,uvel_xy1(2,:,:),res)
delete(sres@cnLevels)           
sres@cnLevels             = (/-0.4/)   ; set the contour levels
sres@cnLineLabelInterval      = 1
sres@cnLineLabelDensityF      = 0.3
sres@cnLineLabelPlacementMode = "constant"
contour3  = gsn_csm_contour(wks, dudy_xy1(2,:,:), sres)
overlay(plot3, contour3)
;===================================================================
; plot E3SMV2
delete(res@cnLevels)
res@vpXF           = 0.10
res@vpYF           = 0.54
res@lbLabelBarOn   = True
res@tiYAxisOffsetXF=-0.03
res@tiYAxisFontHeightF=0.012
res@tiYAxisString=text(0)
res@gsnLeftString  = leftstring00(0)
res@gsnCenterString  = Centerstring00(0)
res@cnLevels       = fspan(18,32,71)
res@lbLabelStride  = 10
vcres@vcRefMagnitudeF           = 0.5 
res@tmYLValues            = (/-10,0,10/)
res@tmYLLabels            = (/"10~S~o~N~S","0","10~S~o~N~N"/)
plot00 = gsn_csm_contour_map_ce(wks,temp_xy2(0,:,:),res)
vector00 = gsn_csm_vector(wks, TAUX_xy2(0,:,:), TAUY_xy2(0,:,:),vcres)
overlay(plot00, vector00)
delete(res@cnLevels) 
res@vpXF           = 0.10
res@vpYF           = 0.39
res@lbLabelBarOn   = True
res@tiYAxisOffsetXF=-0.03
res@tiYAxisFontHeightF=0.012
res@tiYAxisString= text(1)
res@gsnLeftString  = leftstring00(1)
res@gsnCenterString  = Centerstring00(1)
res@cnLevels       = fspan(1,7,61)
vcres@vcRefMagnitudeF     = 0.3   
res@lbLabelStride     = 10
plot01 = gsn_csm_contour_map_ce(wks,temp_xy2(2,:,:),res)
vector01 = gsn_csm_vector(wks, TAUX_xy2(2,:,:), TAUY_xy2(2,:,:),vcres)
overlay(plot01, vector01)
;=========================================
; plot uvel
delete(res@cnLevels)
res@tmYLLabels            = (/"","",""/)
res@tiYAxisString=""
res@vpXF           = 0.40
res@vpYF           = 0.54
res@lbLabelBarOn   = True
res@tiYAxisString  = ""
res@gsnLeftString  = leftstring01(0)
res@gsnCenterString  = Centerstring01(0)
res@cnLevels       = fspan(-0.5,0.5,101)
res@lbLabelStride  = 5
plot02  = gsn_csm_contour_map_ce(wks,uvel_xy2(0,:,:),res)
delete(sres@cnLevels)
sres@cnLevels             = (/-1.0,1.0/)   ; set the contour levels
sres@cnLineLabelDensityF= -0.3
contour02 = gsn_csm_contour(wks, dudy_xy2(0,:,:), sres)
overlay(plot02, contour02)
delete(res@cnLevels)
res@tiYAxisString=""
res@vpXF           = 0.40
res@vpYF           = 0.39
res@lbLabelBarOn   = True
res@tiYAxisString  = ""
res@gsnLeftString  = leftstring01(1)
res@gsnCenterString  = Centerstring01(1)
res@cnLevels       = fspan(-0.3,0.3,61)
plot03  = gsn_csm_contour_map_ce(wks,uvel_xy2(2,:,:),res)
delete(sres@cnLevels)
sres@cnLevels             = (/-0.4/)   ; set the contour levels
sres@cnLineLabelDensityF  = 0.1
contour03 = gsn_csm_contour(wks, dudy_xy2(2,:,:), sres)
overlay(plot03, contour03)
;========================================================================
; plot  vertical 
ares                      = True         
ares@gsnDraw              = False
ares@gsnFrame             = False
ares@gsnAddCyclic         = False     
;ares@gsnScalarContour     = True

ares@tmXTOn   = False
ares@tmYROn   = False
ares@tmYLOn   = True
ares@tmYLMode = "Explicit"
ares@tmYLValues = (/5,50,100,150,200/)
ares@tmYLLabels = (/"5","50","100","150","200"/)
ares@tmXBMode = "Explicit"
ares@tmXBValues = (/-10,-5,0,5,10/)
ares@tmXBLabels = (/"10~S~o~N~S","","Eq","","10~S~o~N~N"/)  

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
ares@trYMinF =   0
ares@trYMaxF =  200
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

ares@pmLabelBarOrthogonalPosF = 0.01
ares@pmLabelBarParallelPosF   = 0.5
ares@lbLabelFontHeightF       = 0.007
ares@lbBoxLinesOn      =False
ares@lbLabelFontHeightF=0.007
ares@lbLabelStride     = 10
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

; ares@cnFillColors         =(/20,40,60,80,110,170,190,210,220,230/)

ares@vpHeightF            = 0.1
ares@vpWidthF             = 0.2
ares@vpXF                  = 0.73
ares@vpYF                  = 0.90
ares@lbLabelBarOn         = True
ares@gsnLeftString =leftstring2(0)
ares@gsnLeftStringFontHeightF = 0.012
ares@cnLevels             = fspan(-1.0, 1.0, 201)
bres@cnLevels             = fspan(14,28,8)

plot4 = gsn_csm_contour(wks, uvel_yz1(0,:,:),ares) 
dumb4 = gsn_csm_contour(wks, temp_yz1(0,:,:),bres)
overlay(plot4, dumb4)
delete(ares@cnLevels)
delete(bres@cnLevels)
ares@cnLevels             = fspan(-0.3, 0.3, 61)
bres@cnLevels             = (/-2,-1,0,1,2,3,4/)
ares@vpHeightF            = 0.1
ares@vpWidthF             = 0.2
ares@vpXF                  = 0.73
ares@vpYF                  = 0.75
ares@lbLabelBarOn         = True
ares@gsnLeftString =leftstring2(1)
ares@lbLabelStride     = 5
plot5 = gsn_csm_contour(wks, uvel_yz1(2,:,:),ares) 
dumb5 = gsn_csm_contour(wks, temp_yz1(2,:,:),bres)
overlay(plot5, dumb5)

; plot seasonal cycle
delete(ares@cnLevels)
delete(bres@cnLevels)
ares@vpHeightF            = 0.1
ares@vpWidthF             = 0.2
ares@vpXF                  = 0.73
ares@vpYF                  = 0.54
ares@lbLabelBarOn         = True
ares@gsnLeftString =leftstring02(0)
ares@gsnLeftStringFontHeightF = 0.012
ares@cnLevels             = fspan(-1.0, 1.0, 201)
bres@cnLevels             = fspan(14,28,8)
bres@cnLineLabelInterval  = 1
plot04 = gsn_csm_contour(wks, uvel_yz2(0,:,:),ares) 
dumb04 = gsn_csm_contour(wks, temp_yz2(0,:,:),bres)
overlay(plot04, dumb04)
delete(ares@cnLevels)
delete(bres@cnLevels)
ares@cnLevels             = fspan(-0.3, 0.3, 61)
bres@cnLevels             = (/-2,-1,0,1,2,3,4/)
ares@vpHeightF            = 0.1
ares@vpWidthF             = 0.2
ares@vpXF                  = 0.73
ares@vpYF                  = 0.39
ares@lbLabelBarOn         = True
ares@gsnLeftString =leftstring02(1)
plot05 = gsn_csm_contour(wks, uvel_yz2(2,:,:),ares) 
dumb05 = gsn_csm_contour(wks, temp_yz2(2,:,:),bres)
overlay(plot05, dumb05)

txres                       = True
txres@txPerimOn             = False  
txres@txBackgroundFillColor = "transparent"
txres@txFontHeightF         = 0.02
txres@txFontThicknessF      = 4
txres@txFont                = 22
txid1     = gsn_create_text(wks, "~F22~E3SMv1", txres)
txid2     = gsn_create_text(wks, "~F22~E3SMv2", txres)

amres                  = True
amres@amParallelPosF   = 0.0    ; 0.5 This is the right edge of the plot.
amres@amOrthogonalPosF = -1.2    ; 0.5 This is the bottom edge of the plot.
amres@amJust           = "TopCenter"

annoid1 = gsn_add_annotation(plot2, txid1, amres)
annoid2 = gsn_add_annotation(plot02, txid2, amres)

draw(plot0)
draw(plot1)
draw(plot2)
draw(plot3)
draw(plot4)
draw(plot5)

draw(plot00)
draw(plot01)
draw(plot02)
draw(plot03)
draw(plot04)
draw(plot05)

end 