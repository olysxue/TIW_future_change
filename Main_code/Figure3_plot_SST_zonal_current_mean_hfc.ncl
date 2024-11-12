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

; seasonal cycle

varave_xy1 = new((/5,dims(1),dims(2)/), float)
dady_xy1   = new((/5,dims(1),dims(2)/), float) ;5,nlat,nlon

varave_yz1 = new((/5,dims(0),dims(1)/), float)
dady_yz1   = new((/5,dims(0),dims(1)/), float) ;5,nlev,nlat

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
varave1= (clm(8,:,:,:)+clm(9,:,:,:)+clm(10,:,:,:) - clm(2,:,:,:)-clm(3,:,:,:)-clm(4,:,:,:))/6.
copy_VarCoords(varave, varave1)
dady   = center_finite_diff_n(varave, dY, False, 0, 1)*1000000
dady1  = center_finite_diff_n(varave1, dY, False, 0, 1)*1000000
copy_VarMeta(var(0,:,:,:), dady)
copy_VarMeta(var(0,:,:,:), dady1)

;---->vertical mean 
varave_xy(i,:,:)= dim_avg_n_Wrap(varave(ilev1:ilev2,:,:),0)
dady_xy(i,:,:)  = dim_avg_n_Wrap(dady(ilev1:ilev2,:,:),0)

varave_xy1(i,:,:)= dim_avg_n_Wrap(varave1(ilev1:ilev2,:,:),0)
dady_xy1(i,:,:)  = dim_avg_n_Wrap(dady1(ilev1:ilev2,:,:),0)

;---->longitude mean 
varave_yz(i,:,:) = dim_avg_n_Wrap(varave(:,:,{lonw:lone}),2)
dady_yz(i,:,:)= dim_avg_n_Wrap(dady(:,:,{lonw:lone}),2)

varave_yz1(i,:,:) = dim_avg_n_Wrap(varave1(:,:,{lonw:lone}),2)
dady_yz1(i,:,:)= dim_avg_n_Wrap(dady1(:,:,{lonw:lone}),2)

delete([/YYYY,iyr,var0,var,varave,dady,varave1,dady1/])
end do
; mean  
  varave_xy(3,:,:)=varave_xy(2,:,:)-varave_xy(0,:,:)
  varave_xy(4,:,:)=varave_xy(1,:,:)-varave_xy(0,:,:)
  
    dady_xy(3,:,:)=dady_xy(2,:,:)-dady_xy(0,:,:)
    dady_xy(4,:,:)=dady_xy(1,:,:)-dady_xy(0,:,:)
  
  varave_yz(3,:,:)=varave_yz(2,:,:)-varave_yz(0,:,:)
  varave_yz(4,:,:)=varave_yz(1,:,:)-varave_yz(0,:,:)
  
    dady_yz(3,:,:)=dady_yz(2,:,:)-dady_yz(0,:,:)
    dady_yz(4,:,:)=dady_yz(1,:,:)-dady_yz(0,:,:)

; seasonal cycle
  varave_xy1(3,:,:)=varave_xy1(2,:,:)-varave_xy1(0,:,:)
  varave_xy1(4,:,:)=varave_xy1(1,:,:)-varave_xy1(0,:,:)
  
  dady_xy1(3,:,:)=dady_xy1(2,:,:)-dady_xy1(0,:,:)
  dady_xy1(4,:,:)=dady_xy1(1,:,:)-dady_xy1(0,:,:)
  
  varave_yz1(3,:,:)=varave_yz1(2,:,:)-varave_yz1(0,:,:)
  varave_yz1(4,:,:)=varave_yz1(1,:,:)-varave_yz1(0,:,:)
  
  dady_yz1(3,:,:)=dady_yz1(2,:,:)-dady_yz1(0,:,:)
  dady_yz1(4,:,:)=dady_yz1(1,:,:)-dady_yz1(0,:,:)
;----------------------------------------------------
copy_VarMeta(temp(0,:,:), varave_xy(0,:,:))
copy_VarMeta(temp(0,:,:), dady_xy(0,:,:))
copy_VarMeta(temp(:,:,0), varave_yz(0,:,:))
copy_VarMeta(temp(:,:,0), dady_yz(0,:,:))

copy_VarMeta(temp(0,:,:), varave_xy1(0,:,:))
copy_VarMeta(temp(0,:,:), dady_xy1(0,:,:))
copy_VarMeta(temp(:,:,0), varave_yz1(0,:,:))
copy_VarMeta(temp(:,:,0), dady_yz1(0,:,:))
print("ok readdata")
return([/varave_xy,dady_xy,varave_yz,dady_yz,varave_xy1,dady_xy1,varave_yz1,dady_yz1/])
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
varave_xy01=new((/5,dims1(0),dims1(1)/), float)   ; mean 
varave_xy02=new((/5,dims1(0),dims1(1)/), float)   ; seasonal cycle

do i=0,2
f   = addfile(ncdir + filelist(i), "r")
YYYY= cd_calendar(f->time,-1)/100

iyr = ind(YYYY.ge.yrStrt(i) .and. YYYY.le.yrLast(i))

var01 = f->$varname$(iyr,:,:)
printVarSummary(var01)
var1 = var01/scale
copy_VarMeta(var01, var1)

clm = clmMonTLL(var1) ; 12month,lev,lat,lon
varave_xy01(i,:,:) = dim_avg_n_Wrap(clm,0) ;lev,lat,lon
varave_xy02(i,:,:) = (clm(8,:,:)+clm(9,:,:)+clm(10,:,:) - clm(2,:,:)-clm(3,:,:)-clm(4,:,:))/6.
copy_VarCoords(varave_xy01, varave_xy02)
delete([/YYYY,iyr,var01,var1/])
end do
; mean 
varave_xy01(3,:,:)=varave_xy01(2,:,:)-varave_xy01(0,:,:)
varave_xy01(4,:,:)=varave_xy01(1,:,:)-varave_xy01(0,:,:)
; seasonal cycle
varave_xy02(3,:,:)=varave_xy02(2,:,:)-varave_xy02(0,:,:)
varave_xy02(4,:,:)=varave_xy02(1,:,:)-varave_xy02(0,:,:)
copy_VarMeta(temp1(:,:), varave_xy01(0,:,:))
copy_VarMeta(temp1(:,:), varave_xy02(0,:,:))

print("ok readdata")
return([/varave_xy01, varave_xy02/])
end

begin
; read SST data
;=============================================================
yrStrt  =  (/111,141,141/)
yrLast  =  (/140,170,170/)
ilev1   = 0
ilev2   = 5
lonw    = 200
lone    = 260
varname="TEMP"
ncdir  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/monthly/"
scale=1. ;degC
filelist=(/"TEMP.mon.0071-0140.20S_20N.120E_80W.levidx01_21_1degree.nc",\
           "TEMP.mon.2xCO2.0071-0170.20S_20N.120E_80W.levidx01_21_1degree.nc",\
           "TEMP.mon.4xCO2.0071-0170.20S_20N.120E_80W.levidx01_21_1degree.nc"/)
out=readdata(ncdir, filelist, varname, scale, yrStrt, yrLast, ilev1, ilev2, lonw, lone)
; mean 
temp_xy1=out[0]
dtdy_xy1=out[1]
temp_yz1=out[2]
dtdy_yz1=out[3]
;==============================
; mean 
temp1_PD    = wgt_areaave_Wrap(temp_xy1(0,{4:6},{210:250}), 1.0, 1.0, 0)
temp1_2xCO2 = wgt_areaave_Wrap(temp_xy1(1,{4:6},{210:250}), 1.0, 1.0, 0)
temp1_4xCO2 = wgt_areaave_Wrap(temp_xy1(2,{4:6},{210:250}), 1.0, 1.0, 0)

temp2_PD    = wgt_areaave_Wrap(temp_xy1(0,{-1:1},{200:250}), 1.0, 1.0, 0)
temp2_2xCO2 = wgt_areaave_Wrap(temp_xy1(1,{-1:1},{200:250}), 1.0, 1.0, 0)
temp2_4xCO2 = wgt_areaave_Wrap(temp_xy1(2,{-1:1},{200:250}), 1.0, 1.0, 0)
;==============================
; calculate the difference
temp1_diff = temp1_4xCO2 - temp1_PD
temp2_diff = temp2_4xCO2 - temp2_PD
print(temp2_diff)

; seasonal cycle
temp_xy2=out[4]
dtdy_xy2=out[5]
temp_yz2=out[6]
dtdy_yz2=out[7]

; seasonal cycle 
temp3_PD    = wgt_areaave_Wrap(temp_xy2(0,{4:6},{210:250}), 1.0, 1.0, 0)
temp3_2xCO2 = wgt_areaave_Wrap(temp_xy2(1,{4:6},{210:250}), 1.0, 1.0, 0)
temp3_4xCO2 = wgt_areaave_Wrap(temp_xy2(2,{4:6},{210:250}), 1.0, 1.0, 0)

temp4_PD    = wgt_areaave_Wrap(temp_xy2(0,{-1:1},{200:250}), 1.0, 1.0, 0)
temp4_2xCO2 = wgt_areaave_Wrap(temp_xy2(1,{-1:1},{200:250}), 1.0, 1.0, 0)
temp4_4xCO2 = wgt_areaave_Wrap(temp_xy2(2,{-1:1},{200:250}), 1.0, 1.0, 0)

; calculate the difference
temp3_diff = temp3_4xCO2 - temp3_PD
temp4_diff = temp4_4xCO2 - temp4_PD

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
ncdir  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/monthly/"
scale=100. ;units cm/s to m/s
filelist=(/"UVEL.mon.0071-0140.20S_20N.120E_80W.levidx01_21_1degree.nc",\
           "UVEL.mon.2xCO2.0071-0170.20S_20N.120E_80W.levidx01_21_1degree.nc",\
           "UVEL.mon.4xCO2.0071-0170.20S_20N.120E_80W.levidx01_21_1degree.nc"/)
out:=readdata(ncdir, filelist, varname, scale, yrStrt, yrLast, ilev1, ilev2, lonw, lone)
;mean
uvel_xy1=out[0]
dudy_xy1=out[1]
uvel_yz1=out[2]
dudy_yz1=out[3]
; seasonal cycle 
uvel_xy2=out[4]
dudy_xy2=out[5]
uvel_yz2=out[6]
dudy_yz2=out[7]

uvel_xy1@units="m/s"
dudy_xy1@units="m/s"
uvel_yz1@units="m/s"
dudy_yz1@units="m/s"

uvel_xy2@units="m/s"
dudy_xy2@units="m/s"
uvel_yz2@units="m/s"
dudy_yz2@units="m/s"

uvel1_PD    = wgt_areaave_Wrap(uvel_xy1(0,{1:3},{200:250}), 1.0, 1.0, 0)
uvel1_2xCO2 = wgt_areaave_Wrap(uvel_xy1(1,{1:3},{200:250}), 1.0, 1.0, 0)
uvel1_4xCO2 = wgt_areaave_Wrap(uvel_xy1(2,{1:3},{200:250}), 1.0, 1.0, 0)

uvel2_PD    = wgt_areaave_Wrap(uvel_xy1(0,{-1:0.5},{200:260}), 1.0, 1.0, 0)
uvel2_2xCO2 = wgt_areaave_Wrap(uvel_xy1(1,{-1:0.5},{200:260}), 1.0, 1.0, 0)
uvel2_4xCO2 = wgt_areaave_Wrap(uvel_xy1(2,{-1:0.5},{200:260}), 1.0, 1.0, 0)

; calculate the difference
uvel1_diff = uvel1_4xCO2 - uvel1_PD
uvel2_diff = uvel2_4xCO2 - uvel2_PD

print(uvel2_diff)
;calculate seasonal cycle 
uvel3_PD    = wgt_areaave_Wrap(uvel_xy2(0,{1:3},{200:250}), 1.0, 1.0, 0)
uvel3_2xCO2 = wgt_areaave_Wrap(uvel_xy2(1,{1:3},{200:250}), 1.0, 1.0, 0)
uvel3_4xCO2 = wgt_areaave_Wrap(uvel_xy2(2,{1:3},{200:250}), 1.0, 1.0, 0)

uvel4_PD    = wgt_areaave_Wrap(uvel_xy2(0,{-1:0.5},{200:260}), 1.0, 1.0, 0)
uvel4_2xCO2 = wgt_areaave_Wrap(uvel_xy2(1,{-1:0.5},{200:260}), 1.0, 1.0, 0)
uvel4_4xCO2 = wgt_areaave_Wrap(uvel_xy2(2,{-1:0.5},{200:260}), 1.0, 1.0, 0)

; calculate the difference
uvel3_diff = uvel3_4xCO2 - uvel3_PD
uvel4_diff = uvel4_4xCO2 - uvel4_PD
;=============================================

dudy_xy_new = dudy_xy1 
dudy_xy_new(:,{-10:-2},{240:280})= dudy_xy1@missing_value
dudy_xy_new(3,{6:10},{220:280})= dudy_xy1@missing_value
dudy_xy_new(0,{6:10},{240:280})= dudy_xy1@missing_value
copy_VarCoords(dudy_xy1, dudy_xy_new)

dudy_xy_new2 = dudy_xy2 
dudy_xy_new2(0,{-10:-1},{240:280})= dudy_xy2@missing_value
dudy_xy_new2(0,{-10:-1},{120:180})= dudy_xy2@missing_value
dudy_xy_new2(0,{6:10},{120:180})= dudy_xy2@missing_value
copy_VarCoords(dudy_xy2, dudy_xy_new2)
printMinMax(dudy_xy_new2(3,:,:),0)
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
ncdir  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/monthly/"
scale=1. ;units cm/s to m/s
filelist=(/"TAUX.mon.0071-0140.20S_20N.120E_80W_1degree.nc",\
           "TAUX.2xCO2.mon.0071-0170.20S_20N.120E_80W_1degree.nc",\
           "TAUX.4xCO2.mon.0071-0170.20S_20N.120E_80W_1degree.nc"/)
out1 =readdata1(ncdir, filelist, varname, scale, yrStrt, yrLast, lonw, lone)
TAUX_xy1=out1[0]
TAUX_xy2=out1[1]

yrStrt  =  (/41,141,141/)
yrLast  =  (/70,170,170/)
ilev1   = 0
ilev2   =10
lonw    =210
lone    =250
varname="TAUY"
ncdir  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/monthly/"
scale=1. ;units cm/s to m/s
filelist=(/"TAUY.mon.0071-0140.20S_20N.120E_80W_1degree.nc",\
           "TAUY.2xCO2.mon.0071-0170.20S_20N.120E_80W_1degree.nc",\
           "TAUY.4xCO2.mon.0071-0170.20S_20N.120E_80W_1degree.nc"/)
out2=readdata1(ncdir, filelist, varname, scale, yrStrt, yrLast, lonw, lone)
TAUY_xy1=out2[0]
TAUY_xy2=out2[1]
;==========================================================================================
 ncdir1  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/mean/flux/"
 ; read mean patterns
 fil01    = systemfunc("ls " + ncdir1 + "CESM2_HR_monthly_heat_flux_convergence_10-50_mean.nc")
 fil02    = systemfunc("ls " + ncdir1 + "CESM2_HR_monthly_momentum_flux_convergence_10-50_mean.nc")
   f01    = addfile(fil01, "r")
   f02    = addfile(fil02, "r")
  
      hfc0_PD     = f01->hfc_PD(0,0:1,:,:)
      hfc0_2xCO2  = f01->hfc_2xCO2(0,0:1,:,:)
      hfc0_4xCO2  = f01->hfc_4xCO2(0,0:1,:,:)
      
      hfc_PD = dim_avg_n_Wrap(hfc0_PD,0)
      hfc_2xCO2 = dim_avg_n_Wrap(hfc0_2xCO2, 0)
      hfc_4xCO2 = dim_avg_n_Wrap(hfc0_4xCO2, 0)

      printVarSummary(hfc_PD)
      printVarSummary(hfc_2xCO2)
      printVarSummary(hfc_4xCO2)
      hfc_diff   = hfc_4xCO2 - hfc_PD
      copy_VarMeta(hfc_PD, hfc_diff)  
      mfc0_PD     = f02->mfc_PD(0,0:1,:,:)
      mfc0_2xCO2  = f02->mfc_2xCO2(0,0:1,:,:)
      mfc0_4xCO2  = f02->mfc_4xCO2(0,0:1,:,:)
      mfc_PD = dim_avg_n_Wrap(mfc0_PD,0)
      mfc_2xCO2 = dim_avg_n_Wrap(mfc0_2xCO2, 0)
      mfc_4xCO2 = dim_avg_n_Wrap(mfc0_4xCO2, 0)
      printVarSummary(mfc_PD)
      printVarSummary(mfc_2xCO2)
      printVarSummary(mfc_4xCO2)
      mfc_diff   = mfc_4xCO2 - mfc_PD
      copy_VarMeta(mfc_PD, mfc_diff)  

      pattern1 = smth9_Wrap(hfc_PD, 0.50, 0.25, False)
      pattern2 = smth9_Wrap(hfc_2xCO2, 0.50, 0.25, False)
      pattern3 = smth9_Wrap(hfc_4xCO2, 0.50, 0.25, False)
      pattern4 = smth9_Wrap(hfc_diff, 0.50, 0.25, False)
      pattern_hfc = new((/2, dimsizes(pattern1&lat), dimsizes(pattern1&lon)/), typeof(pattern1))
      pattern_hfc(0,:,:) = pattern1 
      pattern_hfc(1,:,:) = pattern3

      pattern01 = smth9_Wrap(mfc_PD, 0.50, 0.25, False)
      pattern02 = smth9_Wrap(mfc_2xCO2, 0.50, 0.25, False)
      pattern03 = smth9_Wrap(mfc_4xCO2, 0.50, 0.25, False)
      pattern04 = smth9_Wrap(mfc_diff, 0.50, 0.25, False)
      pattern_mfc = new((/2, dimsizes(pattern01&lat), dimsizes(pattern01&lon)/), typeof(pattern01))
      pattern_mfc(0,:,:) = pattern01 
      pattern_mfc(1,:,:) = pattern03

      delete([/hfc0_PD, hfc0_2xCO2, hfc0_4xCO2/])
      delete([/hfc_PD,  hfc_2xCO2,  hfc_4xCO2/])
      
      delete([/mfc0_PD,mfc0_2xCO2,mfc0_4xCO2/])
      delete([/mfc_PD, mfc_2xCO2, mfc_4xCO2/])
;================================================================================
ncdir01  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/mean/flux/mon/"
fil001    = systemfunc("ls " + ncdir01 + "CESM2_HR_monthly_heat_flux_convergence_10-50_1degree.nc")
fil002    = systemfunc("ls " + ncdir01 + "CESM2_HR_monthly_momentum_flux_convergence_10-50_1degree.nc")
  f001    = addfile(fil001, "r")
  f002    = addfile(fil002, "r")
 
   hfc0_PD     = f001->hfc_PD(:,0:1,:,:)
   hfc0_2xCO2  = f001->hfc_2xCO2(:,0:1,:,:)
   hfc0_4xCO2  = f001->hfc_4xCO2(:,0:1,:,:)
   
   mfc0_PD     = f002->mfc_PD(:,0:1,:,:)
   mfc0_2xCO2  = f002->mfc_2xCO2(:,0:1,:,:)
   mfc0_4xCO2  = f002->mfc_4xCO2(:,0:1,:,:)

   hfc_PD    = dim_avg_n_Wrap(hfc0_PD,1)
   hfc_2xCO2 = dim_avg_n_Wrap(hfc0_2xCO2, 1)
   hfc_4xCO2 = dim_avg_n_Wrap(hfc0_4xCO2, 1)

  
   mfc_PD    = dim_avg_n_Wrap(mfc0_PD,1)
   mfc_2xCO2 = dim_avg_n_Wrap(mfc0_2xCO2, 1)
   mfc_4xCO2 = dim_avg_n_Wrap(mfc0_4xCO2, 1)

  
;===================================================================================================
leftstring0  = (/"~F22~a","~F22~b"/) 
leftstring1 = (/"~F22~c","~F22~d"/) 
leftstring2 = (/"~F22~e","~F22~f"/) 

leftstring00  = (/"~F22~g","~F22~h"/) 
leftstring01  = (/"~F22~i","~F22~j"/) 
leftstring02  = (/"~F22~k","~F22~l"/) 

text=(/"PD","4xCO2 -PD"/)
outdir = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/revise_figures/"
wks = gsn_open_wks("pdf",outdir+"new_Figure3_PD_4xCO2_zonal_current_Mean_HFC_MFC")
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
sres@cnLevels             = (/-1.5,1.5/)
sres@cnLineThicknessF     = 1.5
sres@cnLineLabelInterval  = 1
;sres@cnLineDashSegLenF    = 0.2
sres@cnLineLabelsOn       = True
sres@cnLineLabelFontColor = "black"
sres@cnLineLabelFontHeightF = 0.012
sres@cnLineLabelFontThicknessF = 1.5
sres@cnLineLabelPlacementMode = "computed"
sres@cnLineLabelFormat    = "@^sg"
;sres@cnLineLabelFontAspectF = 0.5
sres@cnLineLabelPerimColor = "transparent"
sres@cnLineLabelBackgroundColor = "white"

sres@cnSmoothingOn        = True
sres@cnSmoothingTensionF  = 3
sres@gsnContourNegLineDashPattern = 16
sres@cnLineDashSegLenF=0.1
sres@gsnLeftString       = ""
sres@gsnRightString      = ""
sres@gsnStringFontHeightF = 0.015
;------------------------------------------
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
res@cnLevels       = fspan(24,34,101)
res@gsnCenterString = " PD SST & Wind Stress"
res@gsnCenterStringFontHeightF = 0.012
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
res@cnLevels       = fspan(0,6,61)
res@gsnCenterString = " 4xCO2-PD"
res@gsnCenterStringFontHeightF = 0.012
vcres@vcRefMagnitudeF     = 0.3   
res@lbLabelStride     = 5
plot1 = gsn_csm_contour_map_ce(wks,temp_xy1(3,:,:),res)
vector1 = gsn_csm_vector(wks, TAUX_xy1(3,:,:), TAUY_xy1(3,:,:),vcres)
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
res@cnLevels       = fspan(-0.5,0.5,101)
res@gsnCenterString = " PD UVEL"
res@gsnCenterStringFontHeightF = 0.012
plot2  = gsn_csm_contour_map_ce(wks,uvel_xy1(0,:,:),res)
contour2  = gsn_csm_contour(wks, dudy_xy_new(0,:,:), sres)
overlay(plot2, contour2)
delete(res@cnLevels)
res@tiYAxisString=""
res@vpXF           = 0.40
res@vpYF           = 0.75
res@lbLabelBarOn   = True
res@tiYAxisString  = ""
res@gsnLeftString  = leftstring1(1)
res@cnLevels       = fspan(-0.3,0.3,61)
res@gsnCenterString = " 4xCO2-PD"
res@gsnCenterStringFontHeightF = 0.012
plot3  = gsn_csm_contour_map_ce(wks,uvel_xy1(3,:,:),res)
delete(sres@cnLevels)           
sres@cnLevels             = (/-0.5/)   ; set the contour levels
sres@cnLineLabelDensityF= 1.5
contour3  = gsn_csm_contour(wks, dudy_xy_new(3,:,:), sres)
overlay(plot3, contour3)
;===================================================================
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
ares@tmYLValues = (/500,5000,10000,15000,20000/)
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
ares@gsnCenterString = " PD"
ares@gsnCenterStringFontHeightF = 0.012
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
ares@gsnCenterString = " 4xCO2-PD"
ares@gsnCenterStringFontHeightF = 0.012
plot5 = gsn_csm_contour(wks, uvel_yz1(3,:,:),ares) 
dumb5 = gsn_csm_contour(wks, temp_yz1(3,:,:),bres)
overlay(plot5, dumb5)
;==========================================================
; plot hfc and mfc

latS   =  -5.
latN   =   10. 
lonL   =  180
lonR   =  280
  mres                      = True         
  mres@gsnDraw              = False        ; don't draw yet
  mres@gsnFrame             = False        ; don't advance frame yet
  mres@gsnAddCyclic         = False      
  mres@cnFillOn             = True               ; turn on color for contours
  mres@cnLinesOn            = False              ; turn off contour lines
  mres@mpFillOn             = False
  mres@cnLineLabelsOn       = False              ; turn off contour line labels
  mres@cnInfoLabelOn        = False                ; turn off info 
  mres@mpLandFillColor      = "gray"             ; color of land
  mres@mpGeophysicalLineColor      = "gray30"
  mres@cnFillPalette="BlueWhiteOrangeRed"
  mres@mpGeophysicalLineThicknessF = 0.7
  mres@mpMaxLatF            = latN
  mres@mpMinLatF            = latS
  mres@mpMaxLonF            = lonR
  mres@mpMinLonF            = lonL  
  mres@mpCenterLonF         = 200
  mres@tmXBOn               = False 
  mres@tmXBMode              = "Explicit"
  mres@tmXBValues            = (/180,200,220,240,260,280/)
  mres@tmXBLabels            = (/" ","","","","",""/)
  mres@tmXBMinorValues       = ispan(120,280,20)
  mres@tmXTOn   = False
  mres@tmYROn   = False
  mres@tmYLMode = "Explicit"
  mres@tmYLValues = (/-20,-10,-5,0,5,10,20/)
  mres@tmYLLabels = (/"20~S~o~N~S","10~S~o~N~S","5~S~o~N~S","0","5~S~o~N~N","10~S~o~N~N","20~S~o~N~N"/)
  mres@tmYLMinorValues =ispan(-20,20,5) 
  mres@cnSmoothingOn         = True
  mres@cnSmoothingTensionF   = 0.05
  mres@tmXBLabelFontHeightF     = 0.010
  mres@tmYLLabelFontHeightF     = 0.010
  mres@tmXBLabelFontThicknessF  = 1.5
  mres@tmYLLabelFontThicknessF  = 1.5 
  mres@tmBorderThicknessF       = 2.5
  mres@tmXBMajorLengthF        = -0.006
  mres@tmYLMajorLengthF        = -0.006
  mres@tmXBMinorLengthF        = -0.003
  mres@tmYLMinorLengthF        = -0.003 
  mres@tmXBMajorThicknessF     = 1.5
  mres@tmYLMajorThicknessF     = 1.5
  mres@tmXBMinorThicknessF     = 1
  mres@tmYLMinorThicknessF     = 1
  mres@tmXBLabelDeltaF          = -0.4
  mres@tmYLLabelDeltaF          = -0.4
  mres@cnFillOn              = True                 ; turn on color fill

  mres@tmXBOn   = True
  mres@tmXTOn   = False
  mres@tmYROn   = False
  mres@tmXMajorGrid                = True
  mres@tmXMajorGridLineDashPattern =  2
  mres@tmXMajorGridLineColor       = "grey30"
  mres@tmXMajorGridThicknessF      = 0.8
  mres@tmYMajorGrid                = True
  mres@tmYMajorGridLineDashPattern =  2
  mres@tmYMajorGridLineColor       = "grey30"
  mres@tmYMajorGridThicknessF      = 0.8 
  mres@tmYMinorGrid                = True  
  mres@tmYMinorGridLineDashPattern =  2
  mres@tmYMinorGridLineColor       = "grey30"
  mres@tmYMinorGridThicknessF      = 0.8 
  mres@tmXMinorGrid                = True     
  mres@tmXMinorGridLineDashPattern =  2
  mres@tmXMinorGridLineColor       = "grey30"
  mres@tmXMinorGridThicknessF      = 0.8 

  mres@cnLevelSelectionMode = "ExplicitLevels"
  mres@lbLabelBarOn         = True
  mres@lbOrientation        = "Vertical"
  mres@lbBoxEndCapStyle     = "TriangleBothEnds"
  mres@gsnStringFontHeightF = 0.012
  mres@lbBoxLinesOn      =False
  mres@lbLabelFontHeightF=0.007
  mres@lbLabelStride     = 5
  mres@pmLabelBarWidthF     = 0.03
  mres@pmLabelBarHeightF    = 0.10

 ;mres@cnFillDrawOrder       = "PreDraw"
  mres@gsnRightString          = " "
  mres@tiMainString            = " "
  mres@mpShapeMode             = "FreeAspect" 

; ----->1 plot mfc
mres@vpHeightF            = 0.1
mres@vpWidthF             = 0.25
mres@vpXF                  = 0.10
mres@vpYF                  = 0.55
mres@cnLevelSelectionMode  = "ExplicitLevels"
mres@cnLevels               := fspan(-0.5,0.5,21)

mres@tiYAxisString= "PD"
mres@tiYAxisOffsetXF=-0.03
mres@tiYAxisFontHeightF=0.012

mres@gsnLeftString  = leftstring00(0)
mres@gsnLeftStringFontHeightF = 0.015
mres@gsnCenterString = "PD MFC_mean"
mres@gsnCenterStringFontHeightF = 0.012
mres@tiMainString =""
mres@lbLabelBarOn    = True
mres@tmXBLabels = (/"","","","","",""/)
mplot0  = gsn_csm_contour_map(wks, pattern_mfc(0,:,:),mres)

mres@vpHeightF            = 0.1
mres@vpWidthF             = 0.25
mres@vpXF                  = 0.10
mres@vpYF                  = 0.40

mres@tiYAxisString= "4xCO2"
mres@tiYAxisOffsetXF=-0.03
mres@tiYAxisFontHeightF=0.012

mres@gsnLeftString  = leftstring00(1)
mres@gsnLeftStringFontHeightF = 0.015
mres@gsnCenterString = "4xCO2 MFC_mean"
mres@gsnCenterStringFontHeightF = 0.012
mres@tiMainString =""
mres@lbLabelBarOn    = True 
mres@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W","80~S~o~N~W"/)
mplot1  = gsn_csm_contour_map(wks, pattern_mfc(1,:,:),mres)
delete([/res@tiYAxisString,res@tiYAxisOffsetXF,res@tiYAxisFontHeightF/])
;---------------------------------------------------------------------------
; ----->2 plot hfc
delete(mres@cnLevels)
mres@cnLevels                = fspan(-0.8,0.8,17)
mres@tmYLLabels = (/"","","","","","",""/)

mres@vpHeightF            = 0.1
mres@vpWidthF             = 0.25
mres@vpXF                  = 0.40
mres@vpYF                  = 0.55
mres@gsnLeftString  = leftstring01(0)
mres@gsnLeftStringFontHeightF = 0.015
mres@gsnCenterString = "PD HFC_mean"
mres@gsnCenterStringFontHeightF = 0.012

mres@lbLabelBarOn    = True 
mres@tiYAxisString   = ""
mres@tmXBLabels = (/"","","","","",""/)
mplot2  = gsn_csm_contour_map(wks,pattern_hfc(0,:,:),mres)

mres@vpHeightF            = 0.1
mres@vpWidthF             = 0.25
mres@vpXF                  = 0.40
mres@vpYF                  = 0.40
mres@gsnLeftString  = leftstring01(1)
mres@gsnLeftStringFontHeightF = 0.015
mres@gsnCenterString = "4xCO2 HFC_mean"
mres@gsnCenterStringFontHeightF = 0.012

mres@lbLabelBarOn    = True
mres@tiYAxisString   = ""
mres@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W","80~S~o~N~W"/)
mplot3  = gsn_csm_contour_map(wks,pattern_hfc(1,:,:),mres)

txres                       = True
txres@txPerimOn             = False  
txres@txBackgroundFillColor = "transparent"
txres@txFontHeightF         = 0.02
txres@txFontThicknessF      = 4
txres@txFont                = 22
txid1     = gsn_create_text(wks, "Mean State: Background States", txres)
txid2     = gsn_create_text(wks, "Mean State: Nonlinear feedbacks", txres)

amres                  = True
amres@amParallelPosF   = -0.80    ; 0.5 This is the right edge of the plot.
amres@amOrthogonalPosF = -1.1    ; 0.5 This is the bottom edge of the plot.
amres@amJust           = "TopLeft"

annoid1 = gsn_add_annotation(plot2, txid1, amres)
annoid2 = gsn_add_annotation(mplot2, txid2, amres)


draw(plot0)
draw(plot1)
draw(plot2)
draw(plot3)
draw(plot4)
draw(plot5)
;==================================================================
w1=360-160
e1=360-110

s1= 1
n1= 3 

plres=True
plres@gsLineThicknessF=2
plres@gsLineDashPattern= 2
plres@gsLineColor="blue"
plres@gsLineOpacityF=100

ddum1 = gsn_add_polyline(wks, mplot0, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum2 = gsn_add_polyline(wks, mplot1, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)

w2=360-160
e2=360-100

s2= -1
n2= 0.5

plres@gsLineThicknessF=2
plres@gsLineDashPattern= 2
plres@gsLineColor="red"
plres@gsLineOpacityF=100


eum1=gsn_add_polyline(wks, mplot0, (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
eum2=gsn_add_polyline(wks, mplot1, (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)

w3=360-160
e3=360-100


w03=360-160
e03=360-110

s03 = -1
n03 = 1 

s3= 0
n3= 2

w4=360-150
e4=360-110

s4= 4
n4= 6

plres@gsLineThicknessF=2
plres@gsLineDashPattern= 2
plres@gsLineColor="blue"
plres@gsLineOpacityF=100

ddum01=gsn_add_polyline(wks, mplot2, (/w3,e3,e3,w3,w3/), (/s03,s03,n03,n03,s03/), plres)
ddum02=gsn_add_polyline(wks, mplot3, (/w3,e3,e3,w3,w3/), (/s03,s03,n03,n03,s03/), plres)

plres=True
plres@gsLineThicknessF=2
plres@gsLineDashPattern= 2
plres@gsLineColor="red"
plres@gsLineOpacityF=100

eum01=gsn_add_polyline(wks, mplot2, (/w4,e4,e4,w4,w4/), (/s4,s4,n4,n4,s4/), plres)
eum02=gsn_add_polyline(wks, mplot3, (/w4,e4,e4,w4,w4/), (/s4,s4,n4,n4,s4/), plres)
;=======================================================================
hfc1_PD    = wgt_areaave_Wrap(hfc_PD(:,{4:6},{210:250}), 1.0, 1.0, 0)
hfc1_2xCO2 = wgt_areaave_Wrap(hfc_2xCO2(:,{4:6},{210:250}), 1.0, 1.0, 0)
hfc1_4xCO2 = wgt_areaave_Wrap(hfc_4xCO2(:,{4:6},{210:250}), 1.0, 1.0, 0)

hfc2_PD    = wgt_areaave_Wrap(hfc_PD(:,{-1:1},{200:250}), 1.0, 1.0, 0)
hfc2_2xCO2 = wgt_areaave_Wrap(hfc_2xCO2(:,{-1:1},{200:250}), 1.0, 1.0, 0)
hfc2_4xCO2 = wgt_areaave_Wrap(hfc_4xCO2(:,{-1:1},{200:250}), 1.0, 1.0, 0)

mfc1_PD    = wgt_areaave_Wrap(mfc_PD(:,{1:3},{200:250}), 1.0, 1.0, 0)
mfc1_2xCO2 = wgt_areaave_Wrap(mfc_2xCO2(:,{1:3},{200:250}), 1.0, 1.0, 0)
mfc1_4xCO2 = wgt_areaave_Wrap(mfc_4xCO2(:,{1:3},{200:250}), 1.0, 1.0, 0)

mfc2_PD    = wgt_areaave_Wrap(mfc_PD(:,{-1:0.5},{200:260}), 1.0, 1.0, 0)
mfc2_2xCO2 = wgt_areaave_Wrap(mfc_2xCO2(:,{-1:0.5},{200:260}), 1.0, 1.0, 0)
mfc2_4xCO2 = wgt_areaave_Wrap(mfc_4xCO2(:,{-1:0.5},{200:260}), 1.0, 1.0, 0)

hfc01_PD    = dim_avg_n_Wrap(hfc1_PD, 0) 
hfc01_2xCO2 = dim_avg_n_Wrap(hfc1_2xCO2, 0) 
hfc01_4xCO2 = dim_avg_n_Wrap(hfc1_4xCO2, 0)  

hfc02_PD    = dim_avg_n_Wrap(hfc2_PD, 0) 
hfc02_2xCO2 = dim_avg_n_Wrap(hfc2_2xCO2, 0) 
hfc02_4xCO2 = dim_avg_n_Wrap(hfc2_4xCO2, 0) 

mfc01_PD    = dim_avg_n_Wrap(mfc1_PD, 0) 
mfc01_2xCO2 = dim_avg_n_Wrap(mfc1_2xCO2, 0) 
mfc01_4xCO2 = dim_avg_n_Wrap(mfc1_4xCO2, 0)  

mfc02_PD    = dim_avg_n_Wrap(mfc2_PD, 0) 
mfc02_2xCO2 = dim_avg_n_Wrap(mfc2_2xCO2, 0) 
mfc02_4xCO2 = dim_avg_n_Wrap(mfc2_4xCO2, 0) 

print("hfc02_4xCO2="+hfc02_4xCO2)
print("mfc02_4xCO2="+mfc02_4xCO2)

; calculate ratio
; mean 
ratio_temp_n = hfc01_4xCO2/temp1_diff
ratio_temp_e = hfc02_4xCO2/temp2_diff

ratio_uvel_n = mfc01_4xCO2/uvel1_diff
ratio_uvel_e = mfc02_4xCO2/uvel2_diff

print("ratio_temp_north="+ratio_temp_n)
print("ratio_temp_equator="+ratio_temp_e)

print("ratio_uvel_north="  +ratio_uvel_n)
print("ratio_uvel_equator="+ratio_uvel_e)

stat  = 0                         ; 0=mean; 1=variance; 2=std. dev; ...
n     = 360                        ; sub-sample size;  default is n=N
nBoot = 10000                      ; bootstrap replications
N     = 360
opt   = False                     ; use defaults
if (isvar("n") .and. (n.lt.N)) then
    opt   = True
    opt@sample_size = n           ; default N
end if

 BootStrap_hfc1_PD     = bootstrap_stat(hfc1_PD, stat, nBoot, 0, opt)
 BootStrap_hfc1_2xCO2  = bootstrap_stat(hfc1_2xCO2, stat, nBoot, 0, opt)
 BootStrap_hfc1_4xCO2  = bootstrap_stat(hfc1_4xCO2, stat, nBoot, 0, opt)

 BootStrap_hfc2_PD     = bootstrap_stat(hfc2_PD, stat, nBoot, 0, opt)
 BootStrap_hfc2_2xCO2  = bootstrap_stat(hfc2_2xCO2, stat, nBoot, 0, opt)
 BootStrap_hfc2_4xCO2  = bootstrap_stat(hfc2_4xCO2, stat, nBoot, 0, opt)

 BootStrap_mfc1_PD     = bootstrap_stat(mfc1_PD, stat, nBoot, 0, opt)
 BootStrap_mfc1_2xCO2  = bootstrap_stat(mfc1_2xCO2, stat, nBoot, 0, opt)
 BootStrap_mfc1_4xCO2  = bootstrap_stat(mfc1_4xCO2, stat, nBoot, 0, opt)

 BootStrap_mfc2_PD     = bootstrap_stat(mfc2_PD, stat, nBoot, 0, opt)
 BootStrap_mfc2_2xCO2  = bootstrap_stat(mfc2_2xCO2, stat, nBoot, 0, opt)
 BootStrap_mfc2_4xCO2  = bootstrap_stat(mfc2_4xCO2, stat, nBoot, 0, opt)
;=============================================================
 xBoot_hfc1_PD      = BootStrap_hfc1_PD[1]   ; mean value
 xBoot_hfc2_PD      = BootStrap_hfc2_PD[1]   ; mean value
 xBoot_mfc1_PD      = BootStrap_mfc1_PD[1]   ; mean value
 xBoot_mfc2_PD      = BootStrap_mfc2_PD[1]   ; mean value  
;---------------------------------------------------
xBoot_hfc1_2xCO2      = BootStrap_hfc1_2xCO2[1]   ; mean value
xBoot_hfc2_2xCO2      = BootStrap_hfc2_2xCO2[1]   ; mean value
xBoot_mfc1_2xCO2      = BootStrap_mfc1_2xCO2[1]   ; mean value
xBoot_mfc2_2xCO2      = BootStrap_mfc2_2xCO2[1]   ; mean value 
;----------------------------------------------------   
xBoot_hfc1_4xCO2      = BootStrap_hfc1_4xCO2[1]   ; mean value 
xBoot_hfc2_4xCO2      = BootStrap_hfc2_4xCO2[1]   ; mean value 
xBoot_mfc1_4xCO2      = BootStrap_mfc1_4xCO2[1]   ; mean value 
xBoot_mfc2_4xCO2      = BootStrap_mfc2_4xCO2[1]   ; mean value 
;------------------------------------------------
 xBoot_hfc01_PD      = BootStrap_hfc1_PD[0]      
 xBoot_hfc02_PD      = BootStrap_hfc2_PD[0]      
 xBoot_mfc01_PD      = BootStrap_mfc1_PD[0]      
 xBoot_mfc02_PD      = BootStrap_mfc2_PD[0]  
;------------------------------------------------
 xBoot_hfc01_2xCO2   = BootStrap_hfc1_2xCO2[0]      
 xBoot_hfc02_2xCO2   = BootStrap_hfc2_2xCO2[0]      
 xBoot_mfc01_2xCO2   = BootStrap_mfc1_2xCO2[0]      
 xBoot_mfc02_2xCO2   = BootStrap_mfc2_2xCO2[0]
;------------------------------------------------
xBoot_hfc01_4xCO2      = BootStrap_hfc1_4xCO2[0]   
xBoot_hfc02_4xCO2      = BootStrap_hfc2_4xCO2[0]   
xBoot_mfc01_4xCO2      = BootStrap_mfc1_4xCO2[0]   
xBoot_mfc02_4xCO2      = BootStrap_mfc2_4xCO2[0] 
;=======================================================================
xBootLow_hfc1_PD    = bootstrap_estimate(xBoot_hfc01_PD, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_hfc1_PD    = bootstrap_estimate(xBoot_hfc01_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_hfc1_PD    = bootstrap_estimate(xBoot_hfc01_PD, 0.95, False)   ; 97.5% upper confidence bound

xBootLow_hfc2_PD    = bootstrap_estimate(xBoot_hfc02_PD, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_hfc2_PD    = bootstrap_estimate(xBoot_hfc02_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_hfc2_PD    = bootstrap_estimate(xBoot_hfc02_PD, 0.95, False)   ; 97.5% upper confidence bound

xBootLow_mfc1_PD    = bootstrap_estimate(xBoot_mfc01_PD, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_mfc1_PD    = bootstrap_estimate(xBoot_mfc01_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_mfc1_PD    = bootstrap_estimate(xBoot_mfc01_PD, 0.95, False)   ; 97.5% upper confidence bound

xBootLow_mfc2_PD    = bootstrap_estimate(xBoot_mfc02_PD, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_mfc2_PD    = bootstrap_estimate(xBoot_mfc02_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_mfc2_PD    = bootstrap_estimate(xBoot_mfc02_PD, 0.95, False)   ; 97.5% upper confidence bound
;--------------------------------------------------------------------------------
xBootLow_hfc1_2xCO2    = bootstrap_estimate(xBoot_hfc01_2xCO2, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_hfc1_2xCO2    = bootstrap_estimate(xBoot_hfc01_2xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_hfc1_2xCO2    = bootstrap_estimate(xBoot_hfc01_2xCO2, 0.95, False)   ; 97.5% upper confidence bound

xBootLow_hfc2_2xCO2    = bootstrap_estimate(xBoot_hfc02_2xCO2, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_hfc2_2xCO2    = bootstrap_estimate(xBoot_hfc02_2xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_hfc2_2xCO2    = bootstrap_estimate(xBoot_hfc02_2xCO2, 0.95, False)   ; 97.5% upper confidence bound

xBootLow_mfc1_2xCO2    = bootstrap_estimate(xBoot_mfc01_2xCO2, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_mfc1_2xCO2    = bootstrap_estimate(xBoot_mfc01_2xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_mfc1_2xCO2    = bootstrap_estimate(xBoot_mfc01_2xCO2, 0.95, False)   ; 97.5% upper confidence bound

xBootLow_mfc2_2xCO2    = bootstrap_estimate(xBoot_mfc02_2xCO2, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_mfc2_2xCO2    = bootstrap_estimate(xBoot_mfc02_2xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_mfc2_2xCO2    = bootstrap_estimate(xBoot_mfc02_2xCO2, 0.95, False)   ; 97.5% upper confidence bound
;---------------------------------------------------------------------------------------------------------
xBootLow_hfc1_4xCO2    = bootstrap_estimate(xBoot_hfc01_4xCO2, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_hfc1_4xCO2    = bootstrap_estimate(xBoot_hfc01_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_hfc1_4xCO2    = bootstrap_estimate(xBoot_hfc01_4xCO2, 0.95, False)   ; 97.5% upper confidence bound

xBootLow_hfc2_4xCO2    = bootstrap_estimate(xBoot_hfc02_4xCO2, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_hfc2_4xCO2    = bootstrap_estimate(xBoot_hfc02_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_hfc2_4xCO2    = bootstrap_estimate(xBoot_hfc02_4xCO2, 0.95, False)   ; 97.5% upper confidence bound

xBootLow_mfc1_4xCO2    = bootstrap_estimate(xBoot_mfc01_4xCO2, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_mfc1_4xCO2    = bootstrap_estimate(xBoot_mfc01_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_mfc1_4xCO2    = bootstrap_estimate(xBoot_mfc01_4xCO2, 0.95, False)   ; 97.5% upper confidence bound

xBootLow_mfc2_4xCO2    = bootstrap_estimate(xBoot_mfc02_4xCO2, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_mfc2_4xCO2    = bootstrap_estimate(xBoot_mfc02_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_mfc2_4xCO2    = bootstrap_estimate(xBoot_mfc02_4xCO2, 0.95, False)   ; 97.5% upper confidence bound
;-------------------------------------------------------------------------------

stdarr1 = new((/3/),"double") ; mfc1
stdarr1(0) =  xBoot_mfc1_PD
stdarr1(1) =  xBoot_mfc1_2xCO2
stdarr1(2) =  xBoot_mfc1_4xCO2

stdarr2 = new((/3/),"double") ; mfc2
stdarr2(0) =  -1 * xBoot_mfc2_PD
stdarr2(1) =  -1 * xBoot_mfc2_2xCO2
stdarr2(2) =  -1 * xBoot_mfc2_4xCO2

stdarr3 = new((/3/),"double") ; hfc1
stdarr3(0) =  -1*xBoot_hfc1_PD
stdarr3(1) =  -1*xBoot_hfc1_2xCO2
stdarr3(2) =  -1*xBoot_hfc1_4xCO2

stdarr4 = new((/3/),"double") ; hfc2
stdarr4(0) =  xBoot_hfc2_PD
stdarr4(1) =  xBoot_hfc2_2xCO2
stdarr4(2) =  xBoot_hfc2_4xCO2
;===========================================

xyres                       = True                 ; plot mods desired
 xyres@gsnFrame              = False                ; don't advance frame yet
 xyres@gsnDraw               = False
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.30                 ; bring bars down to zero
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
 xyres@tmYLMode          = "Explicit"         ; explicit labels
 xyres@tmYLValues        = (/0,0.05,0.10,0.15,0.20,0.25,0.30/)
 xyres@tmYLLabels        = (/"0","0.05","0.10","0.15","0.20","0.25","0.30"/)
 xyres@tmYLLabelFontHeightF = 0.006
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

 xyres@vpXF    = 0.73
 yspace=0.06
 xyres@tmYLOn         = True
 xyres@tmYROn        = False ; turn off right tickmarks
 xyres@tmYRBorderOn  = False ; turn off right axis
 xyres@tmXTBorderOn  = False 
 xyres@vpYF    = 0.55
 xyres@gsnXYBarChartColors = (/"orange3"/)
 xyres@gsnCenterString       = "MFC (NEQ)"
 xyres@gsnCenterStringFontColor  = "red"
 xyres@gsnCenterStringOrthogonalPosF = -0.15
 xyres@gsnCenterStringParallelPosF = 0.45
 xyres@gsnLeftString  = leftstring02(0)
 xyres@gsnLeftStringFontHeightF = 0.015
 xyres@gsnLeftStringOrthogonalPosF = 0.
 xyres@gsnLeftStringParallelPosF = -0.05
 xyres@tmYLValues        = (/0,0.05,0.10,0.15,0.20,0.25,0.30/)
 xyres@tmYLLabels        = (/"0","0.05","0.10","0.15","0.20","0.25","0.30"/)
 mplot4  = gsn_csm_xy(wks,ispan(1,3,1),stdarr1(:),xyres)   ; advancing the frame
 xyres@tmYROn        = False ; turn off right tickmarks
 xyres@tmYLOn        = False
 xyres@tmYRBorderOn  = False                        ; turn off right axis
 xyres@tmYLBorderOn  = False  
 xyres@vpXF = 0.83
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.3               ; bring bars down to zero
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 3.5   
 xyres@gsnXYBarChartColors = (/"cyan3"/)
 xyres@gsnCenterString       =  " -1xMFC (EQ)"
 xyres@gsnCenterStringFontColor   = "blue"
 xyres@gsnCenterStringParallelPosF = 0.55
 xyres@gsnCenterStringOrthogonalPosF = -0.12
 xyres@gsnLeftString  = " "
 xyres@gsnLeftStringFontHeightF = 0.015
 xyres@gsnLeftStringOrthogonalPosF = 0.
 xyres@gsnLeftStringParallelPosF = -0.05
 mplot5  = gsn_csm_xy(wks,ispan(1,3,1),stdarr2(:),xyres)   ; advancing the frame
 xyres@vpXF = 0.73
 xyres@vpYF = 0.40
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@tmYLOn         = True 
 xyres@tmYROn        = False ; turn off right tickmarks
 xyres@tmYRBorderOn  = False    ; turn off right axis
 xyres@tmYLBorderOn  = True   ; turn off right axis
 xyres@vpWidthF              = 0.10
 xyres@vpHeightF             = 0.10
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.5                  ; bring bars down to zero
                       ; turn off right axis
 xyres@gsnXYBarChartColors = (/"orange3"/)
 xyres@tmYLValues        = (/0,0.1,0.2,0.3,0.4,0.5/)
 xyres@tmYLLabels        = (/"0","0.1","0.2","0.3","0.4","0.5"/)
 xyres@gsnCenterString       = " -1*HFC_NEQ"
 xyres@gsnCenterStringFontColor   = "blue"
 xyres@gsnCenterStringOrthogonalPosF = -0.15
 xyres@gsnCenterStringParallelPosF = 0.40
 xyres@trXMinF               = 0.5                    ; adds space on either end
 xyres@trXMaxF               = 4.0                   ; of the 1st and last bars
 xyres@gsnLeftString  = leftstring02(1)
 xyres@gsnLeftStringFontHeightF = 0.015
 xyres@gsnLeftStringOrthogonalPosF = 0.
 mplot6  = gsn_csm_xy(wks,ispan(1,3,1),stdarr3(:),xyres)   ; advancing the frame
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@tmYROn        = False  ; turn off right tickmarks
 xyres@tmYLOn        = False
 xyres@tmYRBorderOn  = False                       ; turn off right axis
 xyres@tmYLBorderOn  = False
 xyres@vpXF = 0.83
 xyres@vpYF = 0.40
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.5                ; bring bars down to zero
 xyres@trXMinF               = 0                    ; adds space on either end
 xyres@trXMaxF               = 3.5 
 xyres@gsnXYBarChartColors = (/"cyan3"/)
 xyres@tmYLValues        = (/0,0.05,0.10,0.15,0.20,0.25,0.30/)
 xyres@tmYLLabels        = (/"0","0.05","0.10","0.15","0.20","0.25","0.30"/)
 xyres@gsnCenterString       = "HFC_EQ"
 xyres@gsnCenterStringFontColor = "red"
 xyres@gsnCenterStringOrthogonalPosF = -0.12
 xyres@gsnCenterStringParallelPosF = 0.55
 xyres@gsnLeftString  = ""
 xyres@gsnLeftStringFontHeightF = 0.015
 xyres@gsnLeftStringOrthogonalPosF = 0.
 mplot7  = gsn_csm_xy(wks,ispan(1,3,1),stdarr4(:),xyres)   ; advancing the frame


 polyres=True
 polyres@gsLineDashPattern=0
 polyres@gsLineColor="black"
 polyres@gsLineThicknessF=1.0
 polyres@tfPolyDrawOrder = "postDraw"
;
 aumm1    = gsn_add_polyline(wks,mplot4,(/1,1/),(/xBootLow_mfc1_PD,   xBootHi_mfc1_PD/), polyres)
 aumm01   = gsn_add_polyline(wks,mplot4,(/2,2/),(/xBootLow_mfc1_2xCO2,xBootHi_mfc1_2xCO2/), polyres)
 aumm001  = gsn_add_polyline(wks,mplot4,(/3,3/),(/xBootLow_mfc1_4xCO2,xBootHi_mfc1_4xCO2/), polyres)
;
 bumm1    = gsn_add_polyline(wks,mplot5,(/1,1/),(/-1*xBootHi_mfc2_PD,   -1*xBootLow_mfc2_PD/), polyres)
 bumm01   = gsn_add_polyline(wks,mplot5,(/2,2/),(/-1*xBootHi_mfc2_2xCO2,-1*xBootLow_mfc2_2xCO2/), polyres)
 bumm001  = gsn_add_polyline(wks,mplot5,(/3,3/),(/-1*xBootHi_mfc2_4xCO2,-1*xBootLow_mfc2_4xCO2/), polyres)

 cumm1    = gsn_add_polyline(wks,mplot6,(/1,1/),(/-1*xBootHi_hfc1_PD,   -1*xBootLow_hfc1_PD/), polyres)
 cumm01   = gsn_add_polyline(wks,mplot6,(/2,2/),(/-1*xBootHi_hfc1_2xCO2,-1*xBootLow_hfc1_2xCO2/), polyres)
 cumm001  = gsn_add_polyline(wks,mplot6,(/3,3/),(/-1*xBootHi_hfc1_4xCO2,-1*xBootLow_hfc1_4xCO2/), polyres)

 dumm1    = gsn_add_polyline(wks,mplot7,(/1,1/),(/xBootLow_hfc2_PD,   xBootHi_hfc2_PD/), polyres)
 dumm01   = gsn_add_polyline(wks,mplot7,(/2,2/),(/xBootLow_hfc2_2xCO2,xBootHi_hfc2_2xCO2/), polyres)
 dumm001  = gsn_add_polyline(wks,mplot7,(/3,3/),(/xBootLow_hfc2_4xCO2,xBootHi_hfc2_4xCO2/), polyres)

 wid = 0.24/2
 ;-------------------------------------------------
 cum1   = gsn_add_polyline(wks,mplot4,(/1-wid,1+wid/),(/xBootLow_mfc1_PD,xBootLow_mfc1_PD/), polyres)
 cum01  = gsn_add_polyline(wks,mplot4,(/1-wid,1+wid/),(/xBootHi_mfc1_PD,xBootHi_mfc1_PD/), polyres)

 cum2   = gsn_add_polyline(wks,mplot4,(/2-wid,2+wid/),(/xBootLow_mfc1_2xCO2,xBootLow_mfc1_2xCO2/), polyres)
 cum02  = gsn_add_polyline(wks,mplot4,(/2-wid,2+wid/),(/xBootHi_mfc1_2xCO2,xBootHi_mfc1_2xCO2/), polyres)

 cum3   = gsn_add_polyline(wks,mplot4,(/3-wid,3+wid/),(/xBootLow_mfc1_4xCO2,xBootLow_mfc1_4xCO2/), polyres)
 cum03  = gsn_add_polyline(wks,mplot4,(/3-wid,3+wid/),(/xBootHi_mfc1_4xCO2,xBootHi_mfc1_4xCO2/), polyres)

 dum1   = gsn_add_polyline(wks,mplot5,(/1-wid,1+wid/),(/-1*xBootLow_mfc2_PD,-1*xBootLow_mfc2_PD/), polyres)
 dum01  = gsn_add_polyline(wks,mplot5,(/1-wid,1+wid/),(/-1*xBootHi_mfc2_PD, -1*xBootHi_mfc2_PD/), polyres)

 dum2   = gsn_add_polyline(wks,mplot5,(/2-wid,2+wid/),(/-1*xBootLow_mfc2_2xCO2,-1*xBootLow_mfc2_2xCO2/), polyres)
 dum02  = gsn_add_polyline(wks,mplot5,(/2-wid,2+wid/),(/-1*xBootHi_mfc2_2xCO2, -1*xBootHi_mfc2_2xCO2/), polyres)

 dum3   = gsn_add_polyline(wks,mplot5,(/3-wid,3+wid/),(/-1*xBootLow_mfc2_4xCO2,-1*xBootLow_mfc2_4xCO2/), polyres)
 dum03  = gsn_add_polyline(wks,mplot5,(/3-wid,3+wid/),(/-1*xBootHi_mfc2_4xCO2, -1*xBootHi_mfc2_4xCO2/), polyres)

 fum1   = gsn_add_polyline(wks,mplot6,(/1-wid,1+wid/),(/-1*xBootLow_hfc1_PD,-1*xBootLow_hfc1_PD/), polyres)
 fum01  = gsn_add_polyline(wks,mplot6,(/1-wid,1+wid/),(/-1*xBootHi_hfc1_PD, -1*xBootHi_hfc1_PD/), polyres)

 fum2   = gsn_add_polyline(wks,mplot6,(/2-wid,2+wid/),(/-1*xBootLow_hfc1_2xCO2,-1*xBootLow_hfc1_2xCO2/), polyres)
 fum02  = gsn_add_polyline(wks,mplot6,(/2-wid,2+wid/),(/-1*xBootHi_hfc1_2xCO2, -1*xBootHi_hfc1_2xCO2/), polyres)

 fum3   = gsn_add_polyline(wks,mplot6,(/3-wid,3+wid/),(/-1*xBootLow_hfc1_4xCO2,-1*xBootLow_hfc1_4xCO2/), polyres)
 fum03  = gsn_add_polyline(wks,mplot6,(/3-wid,3+wid/),(/-1*xBootHi_hfc1_4xCO2, -1*xBootHi_hfc1_4xCO2/), polyres)

 gum1   = gsn_add_polyline(wks,mplot7,(/1-wid,1+wid/),(/xBootLow_hfc2_PD,xBootLow_hfc2_PD/), polyres)
 gum01  = gsn_add_polyline(wks,mplot7,(/1-wid,1+wid/),(/xBootHi_hfc2_PD, xBootHi_hfc2_PD/), polyres)

 gum2   = gsn_add_polyline(wks,mplot7,(/2-wid,2+wid/),(/xBootLow_hfc2_2xCO2,xBootLow_hfc2_2xCO2/), polyres)
 gum02  = gsn_add_polyline(wks,mplot7,(/2-wid,2+wid/),(/xBootHi_hfc2_2xCO2, xBootHi_hfc2_2xCO2/), polyres)

 gum3   = gsn_add_polyline(wks,mplot7,(/3-wid,3+wid/),(/xBootLow_hfc2_4xCO2,xBootLow_hfc2_4xCO2/), polyres)
 gum03  = gsn_add_polyline(wks,mplot7,(/3-wid,3+wid/),(/xBootHi_hfc2_4xCO2, xBootHi_hfc2_4xCO2/), polyres)

draw(mplot0)
draw(mplot1)
draw(mplot2)
draw(mplot3)
draw(mplot4)
draw(mplot5)
draw(mplot6)
draw(mplot7)
end 