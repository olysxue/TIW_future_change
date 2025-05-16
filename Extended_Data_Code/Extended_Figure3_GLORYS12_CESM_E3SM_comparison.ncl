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

delete([/YYYY,iyr,var0,var,varave,dady/])
print("i="+i)
end do

varave_xy(3,:,:)=varave_xy(2,:,:)-varave_xy(0,:,:)
varave_xy(4,:,:)=varave_xy(1,:,:)-varave_xy(0,:,:)

  dady_xy(3,:,:)=dady_xy(2,:,:)-dady_xy(0,:,:)
  dady_xy(4,:,:)=dady_xy(1,:,:)-dady_xy(0,:,:)


copy_VarMeta(temp(0,:,:), varave_xy(0,:,:))
copy_VarMeta(temp(0,:,:), dady_xy(0,:,:))
print("ok readdata")
return([/varave_xy,dady_xy/])
end
;-----------------------------------------------------------------------------------
undef("readdata1")
function readdata1(ncdir, filelist, varname, scale, yrStrt, yrLast, ilev1, ilev2, lonw, lone)
local dims,temp,f
begin
f= addfile(ncdir + filelist(0), "r")
temp=f->$varname$(0,:,:,:) ; nlev,nlat,nlon
printVarSummary(temp)
dims=dimsizes(temp)
varave_xy=new((/3,dims(1),dims(2)/), float)
dady_xy  =new((/3,dims(1),dims(2)/), double) ;5,nlat,nlon

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

delete([/YYYY,iyr,var0,var,varave,dady/])
end do
; mean  
  varave_xy(2,:,:)=varave_xy(1,:,:)-varave_xy(0,:,:)
    dady_xy(2,:,:)=dady_xy(1,:,:)-dady_xy(0,:,:)
;----------------------------------------------------
copy_VarMeta(temp(0,:,:), varave_xy(0,:,:))
copy_VarMeta(temp(0,:,:), dady_xy(0,:,:))
print("ok readdata")
return([/varave_xy,dady_xy/])
end
;--------------------------------------------------------------------------
begin 
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
;=============================================================
CESM_temp = temp_xy(0,:,:)
CESM_uvel = uvel_xy(0,:,:)
delete(filelist)
delete(out)
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
out=readdata1(ncdir, filelist, varname, scale, yrStrt, yrLast, ilev1, ilev2, lonw, lone)
; mean 
temp_xy1=out[0]
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
out:=readdata1(ncdir, filelist, varname, scale, yrStrt, yrLast, ilev1, ilev2, lonw, lone)
uvel_xy1=out[0]
;==========================================================================================
E3SMv1_temp = temp_xy1(0,:,:)
E3SMv1_uvel = uvel_xy1(0,:,:)
;---------------------------------------------------
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
out=readdata1(ncdir, filelist, varname, scale, yrStrt, yrLast, ilev1, ilev2, lonw, lone)
; mean 
temp_xy2=out[0]
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
out:=readdata1(ncdir, filelist, varname, scale, yrStrt, yrLast, ilev1, ilev2, lonw, lone)
uvel_xy2=out[0]

E3SMv2_temp = temp_xy2(0,:,:)
E3SMv2_uvel = uvel_xy2(0,:,:)
;=============================================================

datadir1 = "/Volumes/oly2/CMEMS/mon/"
fil1    = systemfunc("ls " + datadir1 + "GLOBAL-REANALYSIS-PHY-001-030-monthly_uo_0-50m_mean_1degree.nc") 
f1      = addfile(fil1, "r")
CMEMS_uvel = f1->uo(0,:,:)

fil2    = systemfunc("ls " + datadir1 + "GLOBAL-REANALYSIS-PHY-001-030-monthly_thetao_0-50m_mean_1degree.nc") 
f2      = addfile(fil2, "r")
CMEMS_temp = f2->thetao(0,:,:)
;=================================================================
pattern_temp_CMEMS = smth9_Wrap(CMEMS_temp, 0.50, 0.25, False)
pattern_temp_CESM  = smth9_Wrap(CESM_temp, 0.50, 0.25, False)
pattern_temp_E3SMv1 = smth9_Wrap(E3SMv1_temp, 0.50, 0.25, False)
pattern_temp_E3SMv2 = smth9_Wrap(E3SMv2_temp, 0.50, 0.25, False)

pattern_uvel_CMEMS = smth9_Wrap(CMEMS_uvel, 0.50, 0.25, False)
pattern_uvel_CESM  = smth9_Wrap(CESM_uvel, 0.50, 0.25, False)
pattern_uvel_E3SMv1 = smth9_Wrap(E3SMv1_uvel, 0.50, 0.25,False)
pattern_uvel_E3SMv2 = smth9_Wrap(E3SMv2_uvel, 0.50, 0.25, False)

leftstring  = (/"~F22~a","~F22~b","~F22~c","~F22~d"/) 
Centerstring = (/"~F22~GLORYS12 SST","~F22~CESM-UHR SST","~F22~E3SMv1 SST","~F22~E3SMv2 SST"/)

leftstring1  = (/"~F22~e","~F22~f","~F22~g","~F22~h"/) 
Centerstring1 = (/"~F22~GLORYS12 UVEL","~F22~CESM-UHR UVEL","~F22~E3SMv1 UVEL","~F22~E3SMv2 UVEL"/)

outdir = "/Users/xueaoyun/Documents/work/TIW_future/revision1/Review_figure/"
wks = gsn_open_wks("pdf",outdir+"Review_GLORYS12_CESM_E3SM_mean_state")
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
res@gsnRightStringFontHeightF  = 0.012

res@cnFillOn             = True
res@cnFillDrawOrder      = "PreDraw"
res@cnLinesOn            = False
res@cnLineLabelsOn       = False
res@cnInfoLabelOn        = False
res@cnSmoothingOn        = True
res@cnSmoothingTensionF  = 0.05
res@cnLevelSelectionMode = "ExplicitLevels"
res@lbLabelBarOn         = True
res@lbOrientation        = "Horizontal"
res@lbBoxEndCapStyle     = "TriangleBothEnds"
res@lbLabelFontHeightF   = 0.010
res@lbBoxSeparatorLinesOn= False  
res@lbBoxLinesOn         = False
res@gsnRightString       = ""
res@gsnStringFontHeightF = 0.014
;=============================================================
; plot
;=============================================================
xspace = 0.05
plot =new(8, graphic)
;res@cnFillPalette="GMT_wysiwygcont"
res@cnFillPalette="BlueWhiteOrangeRed"
res@gsnSpreadColorStart=15
res@gsnSpreadColorEnd  =245
res@vpHeightF            = 0.1
res@vpWidthF             = 0.18

; res@cnFillColors       = (/20,40,60,80,110,170,190,210,220,230/)
; res@cnLevels             = fspan(-0.6,0.6,31)
;--------------->
;--------------->
res@lbLabelBarOn   = False
res@tiYAxisOffsetXF=-0.03
res@tiYAxisFontHeightF=0.013
res@tiYAxisString=""
res@vpXF           = 0.08
res@vpYF           = 0.90
res@gsnLeftString  = ""
res@cnLevels       = fspan(20,30,101)
plot(0) = gsn_csm_contour_map_ce(wks,pattern_temp_CMEMS,res)

res@vpXF           = 0.08 +res@vpWidthF + xspace
res@vpYF           = 0.90
res@gsnLeftString  = ""
res@cnLevels       = fspan(20,30,101)
res@lbLabelBarOn         = True
res@pmLabelBarWidthF     = 0.50
res@pmLabelBarHeightF    = 0.05
res@pmLabelBarParallelPosF   = 1.17     ; 0.5 This is the right edge of the plot.
res@pmLabelBarOrthogonalPosF = 0.20    ; 0.5 This is the bottom edge of the plot.
res@lbLabelStride     = 10
plot(1) = gsn_csm_contour_map_ce(wks,pattern_temp_CESM,res)

res@vpXF           = 0.08 + 2*res@vpWidthF + 2*xspace
res@vpYF           = 0.90
res@gsnLeftString  = ""
res@cnLevels       = fspan(20,30,101)
res@lbLabelBarOn         = False
plot(2) = gsn_csm_contour_map_ce(wks,pattern_temp_E3SMv1,res)

res@vpXF           = 0.08 + 3*res@vpWidthF + 3*xspace
res@vpYF           = 0.90
res@gsnLeftString  = ""
res@cnLevels       = fspan(20,30,101)
res@lbLabelBarOn         = False
plot(3) = gsn_csm_contour_map_ce(wks,pattern_temp_E3SMv2,res)


txres                       = True
txres@txPerimOn             = False  
txres@txBackgroundFillColor = "White"
txres@txFontHeightF         = 0.012
txres@txFontThicknessF      = 1.5
txid1    = gsn_create_text(wks, leftstring(0), txres)
txid2    = gsn_create_text(wks, leftstring(1), txres)
txid3    = gsn_create_text(wks, leftstring(2), txres)
txid4    = gsn_create_text(wks, leftstring(3), txres)

txid01    = gsn_create_text(wks, Centerstring(0), txres)
txid02    = gsn_create_text(wks, Centerstring(1), txres)
txid03    = gsn_create_text(wks, Centerstring(2), txres)
txid04    = gsn_create_text(wks, Centerstring(3), txres)

 amres                  = True
 amres@amParallelPosF   = -0.485     ; 0.5 This is the right edge of the plot.
 amres@amOrthogonalPosF = -0.74     ; 0.5 This is the bottom edge of the plot.
 amres@amJust           = "TopLeft"

 annoid1 = gsn_add_annotation(plot(0), txid1, amres)
 annoid2 = gsn_add_annotation(plot(1), txid2, amres)
 annoid3 = gsn_add_annotation(plot(2), txid3, amres)
 annoid4 = gsn_add_annotation(plot(3), txid4, amres)
 amres                  = True
 amres@amParallelPosF   = -0.3     ; 0.5 This is the right edge of the plot.
 amres@amOrthogonalPosF = -0.74     ; 0.5 This is the bottom edge of the plot.
 amres@amJust           = "TopLeft"
 annoid01 = gsn_add_annotation(plot(0), txid01, amres)
 annoid02 = gsn_add_annotation(plot(1), txid02, amres)
 annoid03 = gsn_add_annotation(plot(2), txid03, amres)
 annoid04 = gsn_add_annotation(plot(3), txid04, amres)

 draw(plot(0))
 draw(plot(1))
 draw(plot(2))
 draw(plot(3))
;-----------------------------------------------------------
 ; plot uvel
delete(res@cnLevels)
res@tiYAxisString=""
res@vpXF           = 0.08
res@vpYF           = 0.70
res@lbLabelBarOn   = False
res@tiYAxisString  = ""
res@cnLevels       = fspan(-0.5,0.5,101)
plot(4)  = gsn_csm_contour_map_ce(wks,pattern_uvel_CMEMS,res)

delete(res@cnLevels)
res@tiYAxisString=""
res@vpXF           =0.08 + res@vpWidthF + xspace
res@vpYF           = 0.70
res@lbLabelBarOn   = True
res@tiYAxisString  = ""
res@cnLevels       =  fspan(-0.5,0.5,101)
plot(5)  = gsn_csm_contour_map_ce(wks,pattern_uvel_CESM,res)

delete(res@cnLevels)
res@tiYAxisString=""
res@vpXF           = 0.08 + 2*res@vpWidthF + 2*xspace
res@vpYF           = 0.70
res@lbLabelBarOn   = False
res@tiYAxisString  = ""
res@cnLevels       = fspan(-0.5,0.5,101)
plot(6)  = gsn_csm_contour_map_ce(wks,pattern_uvel_E3SMv1,res)

delete(res@cnLevels)
res@tiYAxisString=""
res@vpXF           =0.08 + 3*res@vpWidthF + 3*xspace
res@vpYF           = 0.70
res@lbLabelBarOn   = False
res@tiYAxisString  = ""
res@cnLevels       =  fspan(-0.5,0.5,101)
plot(7)  = gsn_csm_contour_map_ce(wks,pattern_uvel_E3SMv2,res)

txres                       = True
txres@txPerimOn             = False  
txres@txBackgroundFillColor = "White"
txres@txFontHeightF         = 0.012
txres@txFontThicknessF      = 1.5
txid001    = gsn_create_text(wks, leftstring1(0), txres)
txid002    = gsn_create_text(wks, leftstring1(1), txres)
txid003    = gsn_create_text(wks, leftstring1(2), txres)
txid004    = gsn_create_text(wks, leftstring1(3), txres)

txid0001    = gsn_create_text(wks, Centerstring1(0), txres)
txid0002    = gsn_create_text(wks, Centerstring1(1), txres)
txid0003    = gsn_create_text(wks, Centerstring1(2), txres)
txid0004    = gsn_create_text(wks, Centerstring1(3), txres)

amres                  = True
amres@amParallelPosF   = -0.485     ; 0.5 This is the right edge of the plot.
amres@amOrthogonalPosF = -0.74     ; 0.5 This is the bottom edge of the plot.
amres@amJust           = "TopLeft"

bnnoid1 = gsn_add_annotation(plot(4), txid001, amres)
bnnoid2 = gsn_add_annotation(plot(5), txid002, amres)
bnnoid3 = gsn_add_annotation(plot(6), txid003, amres)
bnnoid4 = gsn_add_annotation(plot(7), txid004, amres)
amres                  = True
amres@amParallelPosF   = -0.3     ; 0.5 This is the right edge of the plot.
amres@amOrthogonalPosF = -0.74     ; 0.5 This is the bottom edge of the plot.
amres@amJust           = "TopLeft"
bnnoid01 = gsn_add_annotation(plot(4), txid0001, amres)
bnnoid02 = gsn_add_annotation(plot(5), txid0002, amres)
bnnoid03 = gsn_add_annotation(plot(6), txid0003, amres)
bnnoid04 = gsn_add_annotation(plot(7), txid0004, amres)

draw(plot(4))
draw(plot(5))
draw(plot(6))
draw(plot(7))
end 