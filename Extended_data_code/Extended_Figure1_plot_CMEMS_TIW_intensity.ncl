load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
begin 
;======================================================================================== 
 ; read GLORYS12
 datadir1 = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CMEMS/"
 ; read TIW variance with bandpass
 fil1    = systemfunc("ls " + datadir1 + "CMEMS_tiw_rms_var_ta_10-50_fft_surface_1degree.nc")
 f1 = addfile (fil1, "r")
 tiw_ssta_std_CMEMS = f1->tiw_std
 tiw_ssta_var_CMEMS = f1->tiw_var

 tiw_meanstd_CMEMS1 = dim_avg_n_Wrap(tiw_ssta_std_CMEMS, 0)
 tiw_meanvar_CMEMS1 = dim_avg_n_Wrap(tiw_ssta_var_CMEMS, 0)
;==========================================================================================
 ; read GLORYS12
 datadir2 = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CMEMS/"
 ; read TIW variance with bandpass
 fil2    = systemfunc("ls " + datadir2 + "CMEMS_tiw_rms_var_zos_10-50_fft_surface_1degree.nc")
 f2 = addfile (fil2, "r")
 tiw_ssha_std_CMEMS = f2->tiw_std
 tiw_ssha_var_CMEMS = f2->tiw_var

 tiw_meanstd_CMEMS2 = dim_avg_n_Wrap(tiw_ssha_std_CMEMS, 0)
 tiw_meanvar_CMEMS2 = dim_avg_n_Wrap(tiw_ssha_var_CMEMS, 0)
;=========================================================================================
 ; read GLORYS12
 datadir3 = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CMEMS/"
 ; read TIW variance with bandpass
 fil3    = systemfunc("ls " + datadir3 + "CMEMS_tiw_rms_var_va_10-50_fft_surface_1degree.nc")
 f3 = addfile (fil3, "r")
 tiw_va_std_CMEMS = f3->tiw_std
 tiw_va_var_CMEMS = f3->tiw_var

 tiw_meanstd_CMEMS3 = dim_avg_n_Wrap(tiw_va_std_CMEMS, 0)
 tiw_meanvar_CMEMS3 = dim_avg_n_Wrap(tiw_va_var_CMEMS, 0)
 printMinMax(tiw_meanstd_CMEMS3, 0)
;=========================================================================================
; read GLORYS12
datadir4 = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CMEMS/"
; read TIW variance with bandpass
fil4    = systemfunc("ls " + datadir4 + "CMEMS_EKE_10-50_fft_surface_mon_1degree.nc")
f4 = addfile (fil4, "r")
tiw_std_CMEMS = f4->EKE
tiw_meanstd_CMEMS4 = dim_avg_n_Wrap(tiw_std_CMEMS, 0)
;=============================================================
; plot mean TIW EKE and RMS
;=============================================================
outdir = "/Users/xueaoyun/Documents/work/TIW_future/final_figures/"
wks = gsn_open_wks("pdf",outdir+"plot_CMEMS_TIW_mean_intensity")
wks_type = "eps"
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
res@tmYLValues = (/-10,-5,0,5,10/)
res@tmYLLabels = (/"10~S~o~N~S","5~S~o~N~S","EQ","5~S~o~N~N","10~S~o~N~N"/)
res@tmYLMinorValues         = ispan(-10,10,5) 
res@tmXBMode = "Explicit"
res@tmXBValues = (/180,200,220,240,260,280/)
res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W","80~S~o~N~W"/)  
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
res@lbOrientation        = "Vertical"
res@lbBoxEndCapStyle     = "TriangleBothEnds"
res@pmLabelBarWidthF     = 0.02
res@pmLabelBarHeightF    = 0.18
res@pmLabelBarOrthogonalPosF = 0.05
res@pmLabelBarParallelPosF   = 0.40
res@lbLabelFontHeightF       =0.008
res@gsnLeftString        = ""
res@gsnLeftStringOrthogonalPosF =0.02
res@gsnCenterStringFontHeightF  = 0.015

res@gsnSpreadColorEnd  =230

res@vpHeightF            = 0.15
res@vpWidthF             = 0.30

res@lbTitlePosition  = "Bottom"
res@lbTitleDirection = "Across" 
res@lbTitleFontHeightF= 0.008

plot=new(4,graphic)
dum=plot
dum1 = plot

; ----->1 plot sst
res@lbTitleString    = ""

res@cnLevels             = fspan(0.2,0.6,21)

text=(/"","",""/)
yspace1=0.05
yspace2 = 0.1
yf=0.9
res@tiDeltaF=3
res@gsnCenterStringOrthogonalPosF=0.03
; ----->1 plot sst
res@vpXF    = 0.08
res@vpYF    = yf
res@tiYAxisString=text(0)
res@tiYAxisOffsetXF=-0.02
res@tiYAxisFontHeightF=0.014

res@gsnCenterString = "SST"
res@tiMainString=""
res@lbLabelBarOn    = True
res@gsnLeftString  = "~F22~a "
plot(0)  = gsn_csm_contour_map(wks, tiw_meanstd_CMEMS1(:,:),res)
delete([/res@tiYAxisString,res@tiYAxisOffsetXF,res@tiYAxisFontHeightF/])

; ----->2 plot ssh
res@vpXF    = res@vpWidthF+res@vpXF+yspace2

tiw_meanstd_CMEMS02 = tiw_meanstd_CMEMS2 *100
copy_VarMeta(tiw_meanstd_CMEMS2, tiw_meanstd_CMEMS02)
delete(res@cnLevels)
res@gsnCenterString = "SSH"
res@lbTitleString = ""
;res@lbTitleString    = "10~S~-2~N~m"
res@cnLevels               := fspan(1.5,4,26)
res@cnFillPalette         = "WhiteBlueGreenYellowOrangeRed"
res@lbLabelBarOn    = True
res@gsnLeftString  = "~F22~b"
plot(1)  = gsn_csm_contour_map(wks,tiw_meanstd_CMEMS02(:,:),res)

; ----->2 plot ssh
res@gsnCenterString = "VVEL"
res@vpXF    = 0.08
res@vpYF    = yf-res@vpHeightF -yspace1
delete(res@cnLevels)
res@lbTitleString    = ""
res@cnLevels               := fspan(0.1,0.25,16)
res@cnFillPalette         = "WhiteBlueGreenYellowOrangeRed"
res@lbLabelBarOn    = True
res@gsnLeftString  = "~F22~c"
plot(2)  = gsn_csm_contour_map(wks,tiw_meanstd_CMEMS3(:,:),res)

; ----->3 plot eke
res@gsnCenterString = "EKE"
res@vpXF    = res@vpWidthF+res@vpXF+yspace2
res@vpYF    = yf-res@vpHeightF -yspace1
res@lbTitleString    = ""

delete(res@cnLevels)
res@cnLevels                = fspan(30,60,31)

res@gsnCenterString = "EKE"
res@lbLabelBarOn    = True
res@gsnLeftString  = "~F22~d"  
plot(3)  = gsn_csm_contour_map(wks,tiw_meanstd_CMEMS4(:,:),res)


;=========================================================================
plres=True
plres@gsLineThicknessF=2.5
; plres@gsFillColor     = -1
plres@gsLineDashPattern=0
plres@gsLineColor="blue"
plres@gsLineOpacityF=100
w1=360-140
e1=360-90
s1=0
n1=2

w2=360-160
e2=360-110
s2 = 3
n2 = 8

w3=360-160
e3=360-110
s3 = 2
n3 = 6

w03=360-150
e03=360-100
s03= -2
n03= 1

w4=360-150
e4=360-110
s4 = 0
n4 = 5

dum(0) =gsn_add_polyline(wks, plot(0), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum(1) =gsn_add_polyline(wks, plot(1), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
dum(2) =gsn_add_polyline(wks, plot(2), (/w3,e3,e3,w3,w3/), (/s3,s3,n3,n3,s3/), plres)
dum1(2)=gsn_add_polyline(wks, plot(2), (/w03,e03,e03,w03,w03/), (/s03,s03,n03,n03,s03/), plres)
dum(3) =gsn_add_polyline(wks, plot(3), (/w4,e4,e4,w4,w4/), (/s4,s4,n4,n4,s4/), plres)


draw(plot(0))
draw(plot(1))
draw(plot(2))
draw(plot(3))
end 
