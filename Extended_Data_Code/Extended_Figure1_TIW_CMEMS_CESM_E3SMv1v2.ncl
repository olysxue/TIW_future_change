load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

;=================================================================
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
;-----------------------------------------------------------
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
;===============================================================
begin 
;======================================================================================== 
 ; read GLORYS12
 datadir1 = "/Volumes/oly2/CMEMS/"
 ; read TIW variance with bandpass
 fil1    = systemfunc("ls " + datadir1 + "CMEMS_tiw_rms_var_ta_10-50_fft_surface_1degree.nc")
 f1 = addfile (fil1, "r")
 tiw_ssta_std_CMEMS = f1->tiw_std
 tiw_ssta_var_CMEMS = f1->tiw_var

 tiw_meanstd_CMEMS1 = dim_avg_n_Wrap(tiw_ssta_std_CMEMS, 0)
 tiw_meanvar_CMEMS1 = dim_avg_n_Wrap(tiw_ssta_var_CMEMS, 0)
 printVarSummary(tiw_meanstd_CMEMS1)
;==========================================================================================
 ; read GLORYS12
 datadir2 = "/Volumes/oly2/CMEMS/"
 ; read TIW variance with bandpass
 fil2    = systemfunc("ls " + datadir2 + "CMEMS_tiw_rms_var_zos_10-50_fft_surface_1degree.nc")
 f2 = addfile (fil2, "r")
 tiw_ssha_std_CMEMS = f2->tiw_std
 tiw_ssha_var_CMEMS = f2->tiw_var

 tiw_meanstd_CMEMS2 = dim_avg_n_Wrap(tiw_ssha_std_CMEMS, 0)
 tiw_meanvar_CMEMS2 = dim_avg_n_Wrap(tiw_ssha_var_CMEMS, 0)
;=========================================================================================
 ; read GLORYS12
 datadir3 = "/Volumes/oly2/CMEMS/"
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
datadir4 ="/Volumes/oly2/CMEMS/"
; read TIW variance with bandpass
fil4    = systemfunc("ls " + datadir4 + "CMEMS_EKE_10-50_fft_surface_mon_1degree.nc")
f4 = addfile (fil4, "r")
tiw_std_CMEMS = f4->EKE
tiw_meanstd_CMEMS4 = dim_avg_n_Wrap(tiw_std_CMEMS, 0)
;========================================================================================
; calculate GLORYS12  mean
sst_PD    = wgt_areaave_Wrap(tiw_meanstd_CMEMS1({0:5},{-160:-100}), 1.0, 1.0, 0)
ssh_PD    = wgt_areaave_Wrap(tiw_meanstd_CMEMS2({3:8},{-150:-110}), 1.0, 1.0, 0)
va1_PD    = wgt_areaave_Wrap(tiw_meanstd_CMEMS3({3:6},{-160:-110}), 1.0, 1.0, 0)
va2_PD    = wgt_areaave_Wrap(tiw_meanstd_CMEMS3({-2:1},{-150:-100}), 1.0, 1.0, 0)
eke_PD    = wgt_areaave_Wrap(tiw_meanstd_CMEMS4({0:5},{-160:-100}), 1.0, 1.0, 0)
;================================================
sst_PD1 = wgt_areaave_Wrap(tiw_ssta_std_CMEMS(:,{0:5},{-160:-100}), 1.0, 1.0, 0)
ssh_PD1 = wgt_areaave_Wrap(tiw_ssha_std_CMEMS(:,{3:8},{-150:-110}), 1.0, 1.0, 0)
va1_PD1 = wgt_areaave_Wrap(tiw_va_std_CMEMS(:,{3:6},{-160:-110}), 1.0, 1.0, 0)
va2_PD1 = wgt_areaave_Wrap(tiw_va_std_CMEMS(:,{-2:1},{-140:-90}), 1.0, 1.0, 0)
eke_PD1 = wgt_areaave_Wrap(tiw_std_CMEMS(:,{0:5},{-160:-100}), 1.0, 1.0, 0)

ntim1  = dimsizes(sst_PD1&time)
ntim2  = dimsizes(sst_PD1&time)
ntim3  = dimsizes(sst_PD1&time)

nyears1 = ntim1/12
nyears2 = ntim2/12
nyears3 = ntim3/12

; calculate annul cycle
sst_PD_ann  = reshape(sst_PD1,(/nyears1,12/))
ssh_PD_ann  = reshape(ssh_PD1,(/nyears2,12/))
va1_PD_ann  = reshape(va1_PD1,(/nyears3,12/))
va2_PD_ann  = reshape(va2_PD1,(/nyears3,12/))
eke_PD_ann  = reshape(eke_PD1,(/nyears3,12/))

; annual cycle
sst_PD_ann_ave = dim_avg_n_Wrap(sst_PD_ann, 0)
ssh_PD_ann_ave = dim_avg_n_Wrap(ssh_PD_ann, 0)
va1_PD_ann_ave = dim_avg_n_Wrap(va1_PD_ann, 0)
va2_PD_ann_ave = dim_avg_n_Wrap(va2_PD_ann, 0)
eke_PD_ann_ave = dim_avg_n_Wrap(eke_PD_ann, 0)

;=============================================================
; read CESM data
;=============================================================
; read data
;=============================================================
varname="tiw_std"
ncdir  = "/Volumes/oly2/CESM/mean/"
scale=1. ;units degC
filelist=(/"CESM2_HR_SST_std_var_10-50_fft_PD.nc","CESM2_HR_SST_std_var_10-50_fft_2xCO2.nc","CESM2_HR_SST_std_var_10-50_fft_4xCO2.nc"/)
sst=readdata(ncdir, filelist, varname, scale)
sst@units="degC"
printVarSummary(sst)
scale=1. ; units cm to m by scale
filelist=(/"CESM2_HR_SSH_std_var_10-50_fft_PD.nc","CESM2_HR_SSH_std_var_10-50_fft_2xCO2.nc","CESM2_HR_SSH_std_var_10-50_fft_4xCO2.nc"/)
ssh=readdata(ncdir, filelist, varname, scale)
ssh@units="meter"

scale=100. ; units cm to m by scale
filelist=(/"CESM2_HR_TIW_rms_var_va_10-50_fft_surface_PD.nc","CESM2_HR_TIW_rms_var_va_10-50_fft_surface_2xCO2.nc","CESM2_HR_TIW_rms_var_va_10-50_fft_surface_4xCO2.nc"/)
va=readdata(ncdir, filelist, varname, scale)
va@units="meter"

scale=10000.
varname="EKE" ;units cm2/s2 to m2/s2 by scale
filelist:=(/"CESM2_TIW_EKE_10-50_fft_surface_PD.nc","CESM2_TIW_EKE_10-50_fft_surface_2xCO2.nc","CESM2_TIW_EKE_10-50_fft_surface_4xCO2.nc"/)
eke=readdata(ncdir, filelist, varname, scale)
eke@units="J/m3"
; calculate index PD simulation
cesm_sst    = wgt_areaave_Wrap(sst(0,:,{2:5},{200:260}), 1.0, 1.0, 0)
cesm_ssh    = wgt_areaave_Wrap(ssh(0,:,{3:8},{210:250}), 1.0, 1.0, 0)
cesm_va1    = wgt_areaave_Wrap(va(0,:,{3:6},{200:250}), 1.0, 1.0, 0)
cesm_va2    = wgt_areaave_Wrap(va(0,:,{-2:1},{220:270}), 1.0, 1.0, 0)
cesm_eke    = wgt_areaave_Wrap(eke(0,:,{0:5},{200:260}), 1.0, 1.0, 0)

; calculate mean
cesm_sst_PD    = dim_avg_n_Wrap(cesm_sst,0)
cesm_ssh_PD    = dim_avg_n_Wrap(cesm_ssh,0)
cesm_va1_PD    = dim_avg_n_Wrap(cesm_va1,0)
cesm_va2_PD    = dim_avg_n_Wrap(cesm_va2,0)
cesm_eke_PD    = dim_avg_n_Wrap(cesm_eke,0)

ntim1  = dimsizes(cesm_sst&time)
ntim2  = dimsizes(cesm_sst&time)
ntim3  = dimsizes(cesm_sst&time)
nyears1 = ntim1/12
nyears2 = ntim2/12
nyears3 = ntim3/12

; calculate annual cycle
cesm_sst_PD_ann  = reshape(cesm_sst,(/nyears1,12/))
cesm_ssh_PD_ann  = reshape(cesm_ssh,(/nyears2,12/))
cesm_va1_PD_ann  = reshape(cesm_va1,(/nyears3,12/))
cesm_va2_PD_ann  = reshape(cesm_va2,(/nyears3,12/))
cesm_eke_PD_ann  = reshape(cesm_eke,(/nyears3,12/))

cesm_sst_PD_ann_ave  = dim_avg_n_Wrap(cesm_sst_PD_ann, 0)
cesm_ssh_PD_ann_ave  = dim_avg_n_Wrap(cesm_ssh_PD_ann, 0)
cesm_va1_PD_ann_ave  = dim_avg_n_Wrap(cesm_va1_PD_ann, 0)
cesm_va2_PD_ann_ave  = dim_avg_n_Wrap(cesm_va2_PD_ann, 0)
cesm_eke_PD_ann_ave  = dim_avg_n_Wrap(cesm_eke_PD_ann, 0)
print(cesm_ssh_PD_ann_ave)
;-----------------------------------------------------------
delete(filelist)
;read E3SM 
varname="tos"
ncdir  = "/Volumes/oly2/E3SM/mean/"
scale=1. ;units degC
filelist=(/"tos_day_E3SM-1-0_picontrol_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mon.nc","tos_day_E3SM-1-0_4xCO2_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mon.nc"/)
sst0=readdata1(ncdir, filelist, varname, scale)
sst0@units="degC"
sst1 = tofloat(sqrt(sst0))
copy_VarMeta(sst0, sst1)
printVarSummary(sst1)
varname="zos"
scale=1. ; units cm to m by scale
filelist=(/"zos_day_E3SM-1-0_picontrol_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mon.nc","zos_day_E3SM-1-0_4xCO2_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mon.nc"/)
ssh0=readdata1(ncdir, filelist, varname, scale)
ssh0@units="meter"
printVarSummary(ssh0)
ssh1 = tofloat(sqrt(ssh0)*100.)
copy_VarMeta(ssh0, ssh1)
scale=1.

ncdir1  = "/Volumes/oly2/E3SM/mean/"
filename1   = "TIW_KE_mon_E3SM-1-0_picontrol_gr_00010101-00831231_regrid_1deg_equator_new.nc"
filename2   = "TIW_KE_mon_E3SM-1-0_4xCO2_gr_00010101-00831231_regrid_1deg_equator_new.nc"
fil1        = systemfunc("ls " + ncdir1 + filename1)
fil2        = systemfunc("ls " + ncdir1 + filename2)
f1          =  addfile(fil1, "r")
f2          =  addfile(fil2, "r")
eke01       = f1->KE(:,:,:)
eke02       = f2->KE(:,:,:)

eke1 = new((/2, dimsizes(eke01&time),dimsizes(eke01&lat), dimsizes(eke01&lon)/), float)
eke1(0,:,:,:) = tofloat(eke01)
eke1(1,:,:,:) = tofloat(eke02)
copy_VarCoords(eke01, eke1(0,:,:,:))
copy_VarCoords(eke01, eke1(1,:,:,:))

; calculate E3SMv1 index
e3sm1_sst    = wgt_areaave_Wrap(sst1(0,0:959,{0:5},{200:260}), 1.0, 1.0, 0)
e3sm1_ssh    = wgt_areaave_Wrap(ssh1(0,0:959,{3:8},{210:250}), 1.0, 1.0, 0)
e3sm1_eke    = wgt_areaave_Wrap(eke1(0,0:959,{0:5},{200:260}), 1.0, 1.0, 0)
ntim  = 960
;======================================================================================
delete([/sst1, ssh1, eke01, eke02/]);
;calculate mean
e3sm1_sst_PD    = dim_avg_n_Wrap(e3sm1_sst,0)
e3sm1_ssh_PD    = dim_avg_n_Wrap(e3sm1_ssh,0)
e3sm1_eke_PD    = dim_avg_n_Wrap(e3sm1_eke,0)
; calculate annual cycle
nyears = ntim/12
e3sm1_sst_PD_ann  = reshape(e3sm1_sst,(/nyears,12/))
e3sm1_ssh_PD_ann  = reshape(e3sm1_ssh,(/nyears,12/))
e3sm1_eke_PD_ann  = reshape(e3sm1_eke,(/nyears,12/))

e3sm1_sst_PD_ann_ave = dim_avg_n_Wrap(e3sm1_sst_PD_ann,0)
e3sm1_ssh_PD_ann_ave = dim_avg_n_Wrap(e3sm1_ssh_PD_ann,0)
e3sm1_eke_PD_ann_ave = dim_avg_n_Wrap(e3sm1_eke_PD_ann,0)
;================================================================================
;==========================================read E3SMv1======================================================
ncdir1  = "/Volumes/oly2/E3SM/E3SMv2/"

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
;============================================================================
; calculate index 
e3sm2_sst    = wgt_areaave_Wrap(sst1(0:1799,{0:5},{200:260}), 1.0, 1.0, 0)
e3sm2_ssh    = wgt_areaave_Wrap(ssh1(0:1799,{3:8},{210:250}), 1.0, 1.0, 0)
e3sm2_eke    = wgt_areaave_Wrap(KE1(0:1799,{0:5},{200:260}), 1.0, 1.0, 0)
delete([/sst1, ssh1, KE1/]);
delete([/sst2, ssh2, KE2/]);
;calculate mean
e3sm2_sst_PD    = dim_avg_n_Wrap(e3sm2_sst,0)
e3sm2_ssh_PD    = dim_avg_n_Wrap(e3sm2_ssh,0)
e3sm2_eke_PD    = dim_avg_n_Wrap(e3sm2_eke,0)

ntim1  = dimsizes(e3sm2_sst&time)
nyears1 = ntim1/12

; calculate annual cycle
nyears = ntim/12
e3sm2_sst_PD_ann  = reshape(e3sm2_sst,(/nyears1,12/))
e3sm2_ssh_PD_ann  = reshape(e3sm2_ssh,(/nyears1,12/))
e3sm2_eke_PD_ann  = reshape(e3sm2_eke,(/nyears1,12/))

e3sm2_sst_PD_ann_ave = dim_avg_n_Wrap(e3sm2_sst_PD_ann,0)
e3sm2_ssh_PD_ann_ave = dim_avg_n_Wrap(e3sm2_ssh_PD_ann,0)
e3sm2_eke_PD_ann_ave = dim_avg_n_Wrap(e3sm2_eke_PD_ann,0)
;---------------------------------------------------------------------------------------------------------
stdarr1 = new((/4/),"float") ;  SST
stdarr1(0) =  tofloat(sst_PD)
stdarr1(1) =  tofloat(cesm_sst_PD)
stdarr1(2) =  tofloat(e3sm1_sst_PD)
stdarr1(3) =  tofloat(e3sm2_sst_PD)

stdarr2 = new((/4/),"float") ; SSH
stdarr2(0) =  tofloat(ssh_PD)*100
stdarr2(1) =  tofloat(cesm_ssh_PD)
stdarr2(2) =  tofloat(e3sm1_ssh_PD)
stdarr2(3) =  tofloat(e3sm2_ssh_PD)
print(stdarr2)
stdarr3 = new((/2/),"float") ; VA1
stdarr3(0) =  tofloat(va1_PD)
stdarr3(1) =  tofloat(cesm_va1_PD)

stdarr4 = new((/2/),"float") ; VA2
stdarr4(0) =  tofloat(va2_PD)
stdarr4(1) =  tofloat(cesm_va2_PD)

stdarr5 = new((/4/),"float") ; EKE 
stdarr5(0) =  tofloat(eke_PD)
stdarr5(1) =  tofloat(cesm_eke_PD)
stdarr5(2) =  tofloat(e3sm1_eke_PD)
stdarr5(3) =  tofloat(e3sm2_eke_PD)
;------------------------------------------------------
; eke annual cycle
data1 = new((/4,12/), float)
data1(3,:) = tofloat(eke_PD_ann_ave)
data1(2,:) = tofloat(cesm_eke_PD_ann_ave)
data1(1,:) = tofloat(e3sm1_eke_PD_ann_ave)
data1(0,:) = tofloat(e3sm2_eke_PD_ann_ave)

data2 = new((/4,12/), float)
data2(3,:) = tofloat(sst_PD_ann_ave)
data2(2,:) = tofloat(cesm_sst_PD_ann_ave)
data2(1,:) = tofloat(e3sm1_sst_PD_ann_ave)
data2(0,:) = tofloat(e3sm2_sst_PD_ann_ave)

data3 = new((/4,12/), float)
data3(3,:) = tofloat(ssh_PD_ann_ave)*100
data3(2,:) = tofloat(cesm_ssh_PD_ann_ave)
data3(1,:) = tofloat(e3sm1_ssh_PD_ann_ave)
data3(0,:) = tofloat(e3sm2_ssh_PD_ann_ave)

data4 = new((/2,12/), float)
data4(1,:) = tofloat(va1_PD_ann_ave)
data4(0,:) = tofloat(cesm_va1_PD_ann_ave)

data5 = new((/2,12/), float)
data5(1,:) = tofloat(va2_PD_ann_ave)
data5(0,:) = tofloat(cesm_va2_PD_ann_ave)
;=========================================================
; plot mean TIW EKE and RMS
;=============================================================
outdir = "/Users/xueaoyun/Documents/work/TIW_future/revision1/Extend_figure/"
wks = gsn_open_wks("pdf",outdir+"Extended_Figure1_plot_CMEMS_CESM_E3SM_TIW_mean_intensity")
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
res@gsnRightStringFontHeightF  = 0.010

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
res@lbLabelStride        = 5
res@pmLabelBarWidthF     = 0.02
res@pmLabelBarHeightF    = 0.10
res@pmLabelBarOrthogonalPosF = 0.05
res@pmLabelBarParallelPosF   = 0.50
res@lbLabelFontHeightF       =0.008
res@gsnLeftString        = ""
res@gsnLeftStringOrthogonalPosF =0.012
res@gsnCenterStringFontHeightF  = 0.010

res@gsnSpreadColorEnd  =230

res@vpHeightF            = 0.10
res@vpWidthF             = 0.18


; ----->1 plot sst


res@cnLevels             = fspan(0.2,0.6,21)

text=(/"","",""/)
yspace1=0.05
yspace2 = 0.05
yf=0.9
res@tiDeltaF=3
res@gsnCenterStringOrthogonalPosF=0.03
; ----->1 plot sst
res@vpXF    = 0.08 + 2*res@vpWidthF+2*yspace2
res@vpYF    = yf
res@tiYAxisString=text(0)
res@tiYAxisOffsetXF=-0.02
res@tiYAxisFontHeightF=0.014

res@gsnCenterString = "~F22~TIW-SST"
res@tiMainString=""
res@lbLabelBarOn    = True
res@gsnLeftString  = "~F22~c "
res@gsnLeftStringFontHeightF = 0.012
res@gsnCenterStringFontHeightF = 0.010
res@tmYLLabels = (/"","","","",""/)
plot0  = gsn_csm_contour_map(wks, tiw_meanstd_CMEMS1(:,:),res)
delete([/res@tiYAxisString,res@tiYAxisOffsetXF,res@tiYAxisFontHeightF/])

; ----->2 plot ssh
tiw_meanstd_CMEMS02 = tiw_meanstd_CMEMS2 *100
copy_VarMeta(tiw_meanstd_CMEMS2, tiw_meanstd_CMEMS02)
delete(res@cnLevels)
res@vpXF    = 0.08 + 3*res@vpWidthF+3*yspace2
res@vpYF    = yf
res@gsnCenterString = "~F22~TIW-SSH"
;res@lbTitleString    = "10~S~-2~N~m"
res@cnLevels               := fspan(1.5,4,26)
res@cnFillPalette         = "WhiteBlueGreenYellowOrangeRed"
res@lbLabelBarOn    = True
res@gsnLeftString  = "~F22~d"
plot1  = gsn_csm_contour_map(wks,tiw_meanstd_CMEMS02(:,:),res)

; ----->2 plot ssh
res@gsnCenterString = "~F22~TIW-VVEL"

res@vpXF    = 0.08 + res@vpWidthF +yspace2
res@vpYF    = yf
delete(res@cnLevels)
res@cnLevels               := fspan(0.1,0.25,16)
res@cnFillPalette         = "WhiteBlueGreenYellowOrangeRed"
res@lbLabelBarOn    = True
res@gsnLeftString  = "~F22~b"
plot2  = gsn_csm_contour_map(wks,tiw_meanstd_CMEMS3(:,:),res)

; ----->3 plot eke
res@gsnCenterString = "~F22~TIW-EKE"

res@vpXF    = 0.08
res@vpYF    = yf
delete(res@cnLevels)
res@cnLevels                = fspan(30,60,31)

res@lbLabelBarOn    = True
res@gsnLeftString  = "~F22~a"  
res@tmYLLabels = (/"10~S~o~N~S","5~S~o~N~S","EQ","5~S~o~N~N","10~S~o~N~N"/)
plot3  = gsn_csm_contour_map(wks,tiw_meanstd_CMEMS4(:,:),res)

;=========================================================================
plres=True
plres@gsLineThicknessF=2.5
; plres@gsFillColor     = -1
plres@gsLineDashPattern=0
plres@gsLineColor="blue"
plres@gsLineOpacityF=100
w1=360-160
e1=360-100
s1=0
n1=5

w2=360-160
e2=360-110
s2 = 3
n2 = 8

w3=360-160
e3=360-110
s3 = 2
n3 = 6

w03=360-140
e03=360-90
s03= -2
n03= 1

w4=360-160
e4=360-100
s4 = 0
n4 = 5

dum0 =gsn_add_polyline(wks, plot0, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum1 =gsn_add_polyline(wks, plot1, (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
dum2 =gsn_add_polyline(wks, plot2, (/w3,e3,e3,w3,w3/), (/s3,s3,n3,n3,s3/), plres)
dum3 =gsn_add_polyline(wks, plot2, (/w03,e03,e03,w03,w03/), (/s03,s03,n03,n03,s03/), plres)
dum4 =gsn_add_polyline(wks, plot3, (/w4,e4,e4,w4,w4/), (/s4,s4,n4,n4,s4/), plres)


draw(plot0)
draw(plot1)
draw(plot2)
draw(plot3)

bres               = True 
bres@gsnFrame      = False
bres@gsnDraw       = False
bres@trXMinF       = 0.
bres@trXMaxF       = 13.
bres@trYMinF       = 0.
bres@trYMaxF       = 180.
bres@vpXF          = 0.08                  ; Move slightly to left
bres@vpYF          = 0.75
bres@vpHeightF     = 0.1    
bres@vpWidthF      = 0.18 
bres@gsnYRefLine              = 0. 
bres@gsnYRefLineColors        = (/"black"/)
bres@gsnYRefLineThicknessF    = 2.0
bres@gsnYRefLineDashPattern   = 0 
;bres@gsnXRefLine              = 11.0 
;bres@gsnXRefLineColors        = (/"black"/)
;bres@gsnXRefLineThicknessF    = 2.0
;bres@gsnXRefLineDashPattern   = 1
bres@tfPolyDrawOrder          = "PreDraw" 
bres@tmXBMode                 = "Explicit"                 
bres@tmXBValues               = ispan(0,13,1)
bres@tmXBLabels               = (/"","J","F","M","A","M","J","J","A","S","O","N","D",""/)  
bres@tmYLMode                 = "Explicit"   
bres@tmYLValues               = (/0,30.,60.,90.,120.,150.,180./)
bres@tmYLLabels               = (/"0","30","60","90","120","150","180"/)  
bres@tmXBLabelFontHeightF     = 0.001
bres@xyMarkLineModes          = (/"MarkLines","MarkLines","MarkLines","MarkLines"/)  
bres@xyMarkers                = (/16,16,16,16/)
bres@xyMarkerColors           = (/"green","red","blue","black"/)
bres@xyLineColors             = (/"green","red","blue","black"/)
bres@xyDashPatterns           = (/0/)
bres@xyLineThicknesses        = (/2.0,2.0,2.0,2.0/)
bres@xyMarkerSizes            = (/0.003,0.003,0.003,0.003/) 
bres@tmXBLabelFontHeightF     = 0.008
bres@tmYLLabelFontHeightF     = 0.008
bres@tmXBLabelFont            = 21
bres@tmYLLabelFont            = 21
bres@tmXBLabelFontThicknessF  = 2.0
bres@tmYLLabelFontThicknessF  = 2.0 
bres@tmBorderThicknessF       = 2.5
bres@tmXBMajorLengthF         = 0.005
bres@tmYLMajorLengthF         = 0.005
bres@tmXBMinorLengthF         = 0.006
bres@tmYLMinorLengthF         = 0.006 
bres@tmXBMajorThicknessF      = 2.0
bres@tmYLMajorThicknessF      = 2.0
bres@tmXBMinorThicknessF      = 1.5
bres@tmYLMinorThicknessF      = 1.5 
bres@tmXBLabelDeltaF          = -0.2
bres@tmYLLabelDeltaF          = -0.2 
bres@tiYAxisString            = " "
bres@tiXAxisString            = " "
bres@tiMainString             = " "
bres@tmXBLabelAngleF          = 0
bres@pmLegendDisplayMode    = "Always"            ; turn on legend
bres@pmLegendSide           = "Top"               ; Change location of 
bres@pmLegendZone           = 0        
bres@pmLegendParallelPosF   = -0.30                ; move units right
bres@pmLegendOrthogonalPosF = 0.30                ; move units down
bres@pmLegendWidthF         = 0.04                ; Change width and
bres@pmLegendHeightF        = 0.04                ; height of legend.
bres@lgPerimOn              = False               ; turn off box around
bres@lgLabelFontHeightF     = 0.008              ; label font height
bres@lgLabelFont            = 22                  ; label font height
bres@lgLabelFontThicknessF  = 2                   ; label font height
bres@xyExplicitLegendLabels =(/"E3SMv2", "E3SMv1","CESM","GLORYS"/)
bres@tmXTOn                      = False 
bres@tmYROn                      = False 
bres@tmXMajorGrid                = False 
bres@tmXMajorGridLineDashPattern =  2
bres@tmXMajorGridLineColor       = "grey40"
bres@tmXMajorGridThicknessF      = 0.8
bres@tmYMajorGrid                = True
bres@tmYMajorGridLineDashPattern =  2
bres@tmYMajorGridLineColor       = "grey40"
bres@tmYMajorGridThicknessF      = 0.8 
bres@tmYMinorGrid                = False  
bres@tmYMinorGridLineDashPattern =  2
bres@tmYMinorGridLineColor       = "grey40"
bres@tmYMinorGridThicknessF      = 0.3 
bres@tmXMinorGrid                = False 
bres@tmXMinorGridLineDashPattern =  2
bres@tmXMinorGridLineColor       = "grey40"
bres@tmXMinorGridThicknessF      = 0.3 
; bres@tiMainString                ="northward jump of WNPSH"
bres@tiMainFontHeightF =0.02 
bres@gsnLeftString      ="~F22~e"
bres@gsnCenterString    ="~F22~TIW-EKE"
bres@gsnLeftStringFontHeightF =0.012
bres@gsnCenterStringFontHeightF = 0.010
x = ispan(1,12,1)
plot4 = gsn_csm_xy(wks,x,data1,bres)
delete(bres@tmYLValues)
delete(bres@tmYLLabels)
delete(bres@xyExplicitLegendLabels)
delete(bres@xyMarkerColors)
delete(bres@xyLineColors)
bres@pmLegendOrthogonalPosF = 0.35                ; move units down
bres@xyExplicitLegendLabels = (/"CESM","GLORYS"/)  
bres@xyMarkerColors           = (/"blue","black"/)
bres@xyLineColors             = (/"blue","black"/)
bres@pmLegendWidthF         = 0.04                ; Change width and
bres@pmLegendHeightF        = 0.03                ; height of legend.
bres@trYMinF       = 0.
bres@trYMaxF       = 0.3
bres@vpXF          = bres@vpXF+bres@vpWidthF+yspace1                      ; Move slightly to left
bres@vpYF          = bres@vpYF
bres@vpHeightF     = 0.1    
bres@vpWidthF      = 0.18 
bres@tmYLValues               = (/0,0.1,0.2,0.3,0.4,0.5,0.6/)
bres@tmYLLabels               = (/"0","0.1","0.2","0.3","0.4","0.5","0.6"/)  
bres@gsnLeftString      ="~F22~f"
bres@gsnCenterString    ="~F22~TIW-VVEL (Yanai mode)" 
plot5 = gsn_csm_xy(wks,x,data4,bres)
delete(bres@tmYLValues)
delete(bres@tmYLLabels)
delete(bres@xyExplicitLegendLabels)
delete(bres@xyMarkerColors)
delete(bres@xyLineColors)
bres@pmLegendWidthF         = 0.04                ; Change width and
bres@pmLegendHeightF        = 0.04                ; height of legend.
bres@pmLegendOrthogonalPosF = 0.30                ; move units down
bres@xyExplicitLegendLabels = (/"E3SMv2", "E3SMv1","CESM","GLORYS"/)
bres@xyMarkerColors           = (/"green","red","blue","black"/)
bres@xyLineColors             = (/"green","red","blue","black"/)
bres@trYMinF       = 0.
bres@trYMaxF       = 0.8
bres@vpXF          = bres@vpXF+bres@vpWidthF+yspace1                     ; Move slightly to left
bres@vpYF          = bres@vpYF
bres@vpHeightF     = 0.1    
bres@vpWidthF      = 0.18 
bres@tmYLValues               = (/0,0.2,0.4,0.6,0.8/)
bres@tmYLLabels               = (/"0","0.2","0.4","0.6","0.8"/)  
bres@gsnLeftString      ="~F22~g"
bres@gsnCenterString    ="~F22~TIW-SST"
plot6 = gsn_csm_xy(wks,x,data2,bres)

delete(bres@tmYLValues)
delete(bres@tmYLLabels)
bres@trYMinF       = 0.
bres@trYMaxF       = 6.
bres@vpXF          = bres@vpXF+bres@vpWidthF+yspace1                 ; Move slightly to left
bres@vpYF          = bres@vpYF
bres@vpHeightF     = 0.1    
bres@vpWidthF      = 0.18 
bres@tmYLValues               = (/0,1,2,3,4,5,6/)
bres@tmYLLabels               = (/"0","1","2","3","4","5","6"/) 
bres@gsnLeftString      ="~F22~h"
bres@gsnCenterString    ="~F22~TIW-SSH"
plot7 = gsn_csm_xy(wks,x,data3,bres)
draw(plot4)
draw(plot5)
draw(plot6)
draw(plot7)
; plot labelbar
;EKE
 xyres                       = True                 ; plot mods desired
 xyres@gsnFrame              = False                ; don't advance frame yet
 xyres@gsnDraw               = False
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 5.                   ; of the 1st and last bars
 xyres@vpWidthF              = 0.18
 xyres@vpHeightF             = 0.10

 xyres@gsnXYBarChart         = True                 ; turn on bar chart
 xyres@gsnXYBarChartBarWidth = 0.24                  ; change bar widths
 xyres@tmXBOn                = True                 ; turn off tickmarks at bot
 xyres@tmYLOn                = True
 xyres@tmYROn                = False                 ; turn off tickmarks at bot
 xyres@tmXTOn                = False                 ; turn off tickmarks at bot
 xyres@tmXBMode          = "Explicit"         ; explicit labels
 xyres@tmXBValues        = (/1,2,3,4/)
 xyres@tmXBLabels = (/"GLORYS","CESM","E3SMv1","E3SMv2"/)
 xyres@tmXBLabelFontHeightF = 0.008
 xyres@tmXBLabelAngleF   = 45
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 100.                ; bring bars down to zero
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
;****************************************************
 xspace = 0.05
 yspace=0.05
 yf=0.9
 xyres@vpXF    = 0.08
 xyres@vpYF    = yf-2*yspace-2*res@vpHeightF+0.01
 xyres@gsnXYBarChartColors = (/"orange3"/)
 xyres@gsnCenterString       = "EKE"
 xyres@gsnCenterStringOrthogonalPosF = -0.10
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString  = "~F22~i"
 xyres@gsnLeftStringFontHeightF = 0.012
 xyres@gsnLeftStringOrthogonalPosF = 0.
 xyres@gsnLeftStringParallelPosF = -0.05
 plot8  = gsn_csm_xy(wks,ispan(1,4,1),stdarr5(:),xyres)   ; advancing the frame
 ;****************************************************
 xyres@vpXF = xyres@vpXF + xyres@vpWidthF + xspace
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@vpWidthF              = 0.09
 xyres@vpHeightF             = 0.10
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.25                 ; bring bars down to zero
 xyres@gsnXYBarChartColors = (/"chartreuse2"/)
 xyres@tmYLValues        = (/0,0.05,0.10,0.15,0.20,0.25/)
 xyres@tmYLLabels        = (/"0","0.05","0.10","0.15","0.20","0.25"/)
 xyres@gsnCenterString       = "   VVEL_NEQ"
 xyres@gsnCenterStringOrthogonalPosF = -0.12
 xyres@gsnCenterStringParallelPosF = 0.40
 xyres@trXMinF               = 0.5                    ; adds space on either end
 xyres@trXMaxF               = 2.75                    ; of the 1st and last bars
 xyres@gsnLeftString  = "~F22~j"
 xyres@gsnLeftStringFontHeightF = 0.012
 xyres@gsnLeftStringOrthogonalPosF = 0.
 plot9  = gsn_csm_xy(wks,ispan(1,2,1),stdarr3(:),xyres)   ; advancing the frame

 xyres@tmYLOn        = False 
 xyres@tmYRBorderOn  = False                        ; turn off right axis
 xyres@tmYLBorderOn  = False
 xyres@vpXF = xyres@vpWidthF+xyres@vpXF
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.2                 ; bring bars down to zero
 xyres@trXMinF               = 0.25                    ; adds space on either end
 xyres@trXMaxF               = 2.5 
 xyres@gsnXYBarChartColors = (/"chartreuse2"/)
 xyres@tmYLValues        = (/0,0.05,0.10,0.15,0.20,0.25/)
 xyres@tmYLLabels        = (/"0","0.05","0.10","0.15","0.20","0.25"/)
 xyres@gsnCenterString       = "    VVEL_EQ"
 xyres@gsnCenterStringOrthogonalPosF = -0.10
 xyres@gsnCenterStringParallelPosF = 0.45
 xyres@gsnLeftString  = ""
 xyres@gsnLeftStringFontHeightF = 0.012
 xyres@gsnLeftStringOrthogonalPosF = 0.
 plot10  = gsn_csm_xy(wks,ispan(1,2,1),stdarr4(:),xyres)   ; advancing the frame
  ;****************************************************
 xyres@vpXF = xyres@vpXF + xyres@vpWidthF + xspace+0.01
  ;****************************************************
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@vpWidthF              = 0.18
 xyres@vpHeightF             = 0.10
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 5.                   ; of the 1st and last bars
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.6                 ; bring bars down to zero
 xyres@tmYLOn        = True
 xyres@tmYRBorderOn  = False                        ; turn off right axis
 xyres@tmYLBorderOn  = True 
 xyres@gsnXYBarChartColors = (/"cyan3"/)
 xyres@gsnCenterString       = "SST"
 xyres@gsnCenterStringOrthogonalPosF = -0.10
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString  = "~F22~k"
 xyres@gsnLeftStringFontHeightF = 0.012
 xyres@gsnLeftStringOrthogonalPosF = 0.
 xyres@gsnLeftStringParallelPosF = -0.05
 xyres@tmYLMode          = "Explicit"         ; explicit labels
 xyres@tmYLValues        = (/0,0.1,0.2,0.3,0.4,0.5,0.6/)
 xyres@tmYLLabels        = (/"0","0.1","0.2","0.3","0.4","0.5","0.6"/)
 plot11 = gsn_csm_xy(wks,ispan(1,4,1),stdarr1(:),xyres)   ; advancing the frame
;=========================================================
delete(xyres@tmYLValues)
delete(xyres@tmYLLabels)
xyres@vpXF = xyres@vpXF + xyres@vpWidthF + xspace
xyres@trYMinF               = 0                    ; bring bars down to zero
xyres@trYMaxF               = 5.                 ; bring bars down to zero
xyres@gsnXYBarChartColors = (/"pink"/)
xyres@tmYLValues        = (/0,1,2,3,4,5/)
xyres@tmYLLabels        = (/"0","1","2","3","4","5"/)
xyres@gsnCenterString       = "SSH"
xyres@gsnCenterStringParallelPosF = 0.45
xyres@gsnLeftString  = "~F22~l"
xyres@gsnLeftStringFontHeightF = 0.012
xyres@gsnLeftStringOrthogonalPosF = 0.
xyres@gsnLeftStringParallelPosF = -0.05
plot12  = gsn_csm_xy(wks,ispan(1,4,1),stdarr2(:),xyres)   ; advancing the frame
 ;****************************************************
 draw(plot8)
 draw(plot9)
 draw(plot10)
 draw(plot11)
 draw(plot12)

end 
