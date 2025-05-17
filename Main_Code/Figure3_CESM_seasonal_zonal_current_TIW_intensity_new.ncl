; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
; line 65: check units of EKE whether cm2/s2 in file
; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
; !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
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
;=============================================
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
;=================================================
begin 
; parameter constant 
PI  = 3.14159265
R   = 6378388.
rad = PI/180.
dens= 1025 ; reference density of seawater
sec = 3600*24*30.42 ;one month
yrStrt  =  (/111,141,141/)
yrLast  =  (/140,170,170/)
ilev1   = 0
ilev2   = 5
lonw    = 200
lone    = 260
ncdir1  = "/Volumes/oly2/CESM/monthly/"
fil1    = systemfunc("ls "+ ncdir1 + "TEMP.mon.0071-0140.20S_20N.120E_80W.levidx01_21_1degree.nc")
fil2    = systemfunc("ls "+ ncdir1 + "TEMP.mon.2xCO2.0071-0170.20S_20N.120E_80W.levidx01_21_1degree.nc")
fil3    = systemfunc("ls "+ ncdir1 + "TEMP.mon.4xCO2.0071-0170.20S_20N.120E_80W.levidx01_21_1degree.nc")
   
f1      = addfile(fil1, "r")
f2      = addfile(fil2, "r")
f3      = addfile(fil3, "r")

YYYY1= cd_calendar(f1->time,-1)/100
iyr1 = ind(YYYY1.ge.yrStrt(0) .and. YYYY1.le.yrLast(0))

YYYY2= cd_calendar(f2->time,-1)/100
iyr2 = ind(YYYY2.ge.yrStrt(1) .and. YYYY2.le.yrLast(1))

YYYY3= cd_calendar(f3->time,-1)/100
iyr3 = ind(YYYY3.ge.yrStrt(2) .and. YYYY3.le.yrLast(2))

temp1  = f1->TEMP(iyr1,:,:,:)
temp2  = f2->TEMP(iyr2,:,:,:)
temp3  = f3->TEMP(iyr3,:,:,:)

temp1_xy  = dim_avg_n_Wrap(temp1(:,0:1,:,:),1)
temp2_xy  = dim_avg_n_Wrap(temp2(:,0:1,:,:),1)
temp3_xy  = dim_avg_n_Wrap(temp3(:,0:1,:,:),1)

; calculate domain-avg
temp1_nino3 = wgt_areaave_Wrap(temp1_xy(:,{-2:2},{200:270}), 1.0, 1.0, 0)
temp2_nino3 = wgt_areaave_Wrap(temp2_xy(:,{-2:2},{200:270}), 1.0, 1.0, 0)
temp3_nino3 = wgt_areaave_Wrap(temp3_xy(:,{-2:2},{200:270}), 1.0, 1.0, 0)
printVarSummary(temp1_nino3)
dlat   = (temp1&lat(2)-temp1&lat(1))*rad
dY     = R*dlat

dTdy1_xy   = center_finite_diff_n(temp1_xy, dY, False, 0, 1)*1000000.
dTdy2_xy   = center_finite_diff_n(temp2_xy, dY, False, 0, 1)*1000000.
dTdy3_xy   = center_finite_diff_n(temp3_xy, dY, False, 0, 1)*1000000.

copy_VarMeta(temp1(:,0,:,:), dTdy1_xy)
copy_VarMeta(temp2(:,0,:,:), dTdy2_xy)
copy_VarMeta(temp3(:,0,:,:), dTdy3_xy)

dTdy1_xy_nino3  = runave(wgt_areaave_Wrap(dTdy1_xy(:,{2:5},{200:260}), 1.0, 1.0, 0),3,0)
dTdy2_xy_nino3  = runave(wgt_areaave_Wrap(dTdy2_xy(:,{2:5},{200:260}), 1.0, 1.0, 0),3,0)
dTdy3_xy_nino3  = runave(wgt_areaave_Wrap(dTdy3_xy(:,{2:5},{200:260}), 1.0, 1.0, 0),3,0)


ntim  = dimsizes(temp1&time)
nyears = ntim/12
print(nyears)
temp1_nino3_ann  = reshape(temp1_nino3,(/nyears,12/))
temp2_nino3_ann  = reshape(temp2_nino3,(/nyears,12/))
temp3_nino3_ann  = reshape(temp3_nino3,(/nyears,12/))

dTdy1_index_ann  = reshape(dTdy1_xy_nino3,(/nyears,12/))
dTdy2_index_ann  = reshape(dTdy2_xy_nino3,(/nyears,12/))
dTdy3_index_ann  = reshape(dTdy3_xy_nino3,(/nyears,12/))

; annual cycle
temp1_nino3_ave = dim_avg_n_Wrap(temp1_nino3_ann, 0) - avg(dim_avg_n_Wrap(temp1_nino3_ann, 0))
temp2_nino3_ave = dim_avg_n_Wrap(temp2_nino3_ann, 0) - avg(dim_avg_n_Wrap(temp2_nino3_ann, 0))
temp3_nino3_ave = dim_avg_n_Wrap(temp3_nino3_ann, 0) - avg(dim_avg_n_Wrap(temp3_nino3_ann, 0))

dTdy1_index_ave = dim_avg_n_Wrap(dTdy1_index_ann, 0)
dTdy2_index_ave = dim_avg_n_Wrap(dTdy2_index_ann, 0)
dTdy3_index_ave = dim_avg_n_Wrap(dTdy3_index_ann, 0)

; std 
temp1_nino3_std = dim_stddev_n_Wrap(temp1_nino3_ann, 0) 
temp2_nino3_std = dim_stddev_n_Wrap(temp2_nino3_ann, 0) 
temp3_nino3_std = dim_stddev_n_Wrap(temp3_nino3_ann, 0) 

dTdy1_index_std = dim_stddev_n_Wrap(dTdy1_index_ann, 0) 
dTdy2_index_std = dim_stddev_n_Wrap(dTdy2_index_ann, 0) 
dTdy3_index_std = dim_stddev_n_Wrap(dTdy3_index_ann, 0) 
;================================================================================
; read UVEL data
;=============================================================
yrStrt  =  (/111,141,141/)
yrLast  =  (/140,170,170/)
ilev1   = 0
ilev2   = 5
lonw    =210
lone    =250
varname="UVEL"
ncdir2  ="/Volumes/oly2/CESM/monthly/"
fil01    = systemfunc("ls "+ ncdir2 + "UVEL.mon.0071-0140.20S_20N.120E_80W.levidx01_21_1degree.nc")
fil02    = systemfunc("ls "+ ncdir2 + "UVEL.mon.2xCO2.0071-0170.20S_20N.120E_80W.levidx01_21_1degree.nc")
fil03    = systemfunc("ls "+ ncdir2 + "UVEL.mon.4xCO2.0071-0170.20S_20N.120E_80W.levidx01_21_1degree.nc")
   
f01      = addfile(fil01, "r")
f02      = addfile(fil02, "r")
f03      = addfile(fil03, "r")

YYYY01= cd_calendar(f01->time,-1)/100
iyr01 = ind(YYYY01.ge.yrStrt(0) .and. YYYY01.le.yrLast(0))

YYYY02= cd_calendar(f02->time,-1)/100
iyr02 = ind(YYYY02.ge.yrStrt(1) .and. YYYY02.le.yrLast(1))

YYYY03= cd_calendar(f03->time,-1)/100
iyr03 = ind(YYYY03.ge.yrStrt(2) .and. YYYY03.le.yrLast(2))

uvel01  = f01->UVEL(iyr01,:,:,:)
uvel02  = f02->UVEL(iyr02,:,:,:)
uvel03  = f03->UVEL(iyr03,:,:,:)
uvel1   = uvel01/100.
uvel2   = uvel02/100.
uvel3   = uvel03/100.
copy_VarMeta(uvel01, uvel1)
copy_VarMeta(uvel02, uvel2)
copy_VarMeta(uvel03, uvel3)

uvel1_xy  = dim_avg_n_Wrap(uvel1(:,ilev1:ilev2,:,:),1)
uvel2_xy  = dim_avg_n_Wrap(uvel2(:,ilev1:ilev2,:,:),1)
uvel3_xy  = dim_avg_n_Wrap(uvel3(:,ilev1:ilev2,:,:),1)


dlat   = (temp1&lat(2)-temp1&lat(1))*rad
dY     = R*dlat

dUdy1_xy   = center_finite_diff_n(uvel1_xy, dY, False, 0, 1)*1000000.
dUdy2_xy   = center_finite_diff_n(uvel2_xy, dY, False, 0, 1)*1000000.
dUdy3_xy   = center_finite_diff_n(uvel3_xy, dY, False, 0, 1)*1000000.

copy_VarMeta(uvel1(:,0,:,:), dUdy1_xy)
copy_VarMeta(uvel2(:,0,:,:), dUdy2_xy)
copy_VarMeta(uvel3(:,0,:,:), dUdy3_xy)

dUdy1_xy_index  = runave(wgt_areaave_Wrap(dUdy1_xy(:,{3:6},{200:250}), 1.0, 1.0, 0),3,0)
dUdy2_xy_index  = runave(wgt_areaave_Wrap(dUdy2_xy(:,{3:6},{200:250}), 1.0, 1.0, 0),3,0)
dUdy3_xy_index  = runave(wgt_areaave_Wrap(dUdy3_xy(:,{3:6},{200:250}), 1.0, 1.0, 0),3,0)

dUdy01_xy_index  = runave(wgt_areaave_Wrap(dUdy1_xy(:,{-1:1},{200:250}), 1.0, 1.0, 0),3,0)
dUdy02_xy_index  = runave(wgt_areaave_Wrap(dUdy2_xy(:,{-1:1},{200:250}), 1.0, 1.0, 0),3,0)
dUdy03_xy_index  = runave(wgt_areaave_Wrap(dUdy3_xy(:,{-1:1},{200:250}), 1.0, 1.0, 0),3,0)


ntim  = dimsizes(temp1&time)
nyears = ntim/12

dUdy1_index_ann  = reshape(dUdy1_xy_index,(/nyears,12/))
dUdy2_index_ann  = reshape(dUdy2_xy_index,(/nyears,12/))
dUdy3_index_ann  = reshape(dUdy3_xy_index,(/nyears,12/))

dUdy01_index_ann  = reshape(dUdy01_xy_index,(/nyears,12/))
dUdy02_index_ann  = reshape(dUdy02_xy_index,(/nyears,12/))
dUdy03_index_ann  = reshape(dUdy03_xy_index,(/nyears,12/))

; annual cycle
dUdy1_index_ave = dim_avg_n_Wrap(dUdy1_index_ann, 0) 
dUdy2_index_ave = dim_avg_n_Wrap(dUdy2_index_ann, 0) 
dUdy3_index_ave = dim_avg_n_Wrap(dUdy3_index_ann, 0) 

dUdy01_index_ave = dim_avg_n_Wrap(dUdy01_index_ann, 0)
dUdy02_index_ave = dim_avg_n_Wrap(dUdy02_index_ann, 0)
dUdy03_index_ave = dim_avg_n_Wrap(dUdy03_index_ann, 0)

; std 
dUdy1_index_std = dim_stddev_n_Wrap(dUdy1_index_ann, 0) 
dUdy2_index_std = dim_stddev_n_Wrap(dUdy2_index_ann, 0) 
dUdy3_index_std = dim_stddev_n_Wrap(dUdy3_index_ann, 0) 
print(dUdy1_index_std)
print(dUdy3_index_std)
dUdy01_index_std = dim_stddev_n_Wrap(dUdy01_index_ann, 0) 
dUdy02_index_std = dim_stddev_n_Wrap(dUdy02_index_ann, 0) 
dUdy03_index_std = dim_stddev_n_Wrap(dUdy03_index_ann, 0) 

;==========================================================================================

data1 = new((/2,12/), float)
data1(0,:) = dTdy1_index_ave
data1(1,:) = dTdy3_index_ave

data2 = new((/2,12/), float)
data2(0,:) = dUdy01_index_ave
data2(1,:) = dUdy03_index_ave

data3 = new((/2,12/), float)
data3(0,:) = dUdy1_index_ave
data3(1,:) = dUdy3_index_ave

xp1_temp_PD = new((/24/), float)
yp1_temp_PD = new((/24/), float)

do i=0,11
   xp1_temp_PD(i)          = i+1
   xp1_temp_PD(2*12-1-i)   = i+1
   yp1_temp_PD(i)          = temp1_nino3_ave(i) + temp1_nino3_std(i)/2.
   yp1_temp_PD(2*12-1-i)   = temp1_nino3_ave(i) - temp1_nino3_std(i)/2. 
 end do

 xp1_temp_4xCO2 = new((/24/), float)
 yp1_temp_4xCO2 = new((/24/), float)
 
 do i=0,11
    xp1_temp_4xCO2(i)          = i+1
    xp1_temp_4xCO2(2*12-1-i)   = i+1
    yp1_temp_4xCO2(i)          = temp3_nino3_ave(i) + temp3_nino3_std(i)/2.
    yp1_temp_4xCO2(2*12-1-i)   = temp3_nino3_ave(i) - temp3_nino3_std(i)/2. 
  end do

  ;==================================================
xp1_dTdy_PD = new((/24/), float)
yp1_dTdy_PD = new((/24/), float)

do i=0,11
   xp1_dTdy_PD(i)          = i+1
   xp1_dTdy_PD(2*12-1-i)   = i+1
   yp1_dTdy_PD(i)          = dTdy1_index_ave(i) + dTdy1_index_std(i)/2.
   yp1_dTdy_PD(2*12-1-i)   = dTdy1_index_ave(i) - dTdy1_index_std(i)/2. 
 end do

 xp1_dTdy_4xCO2 = new((/24/), float)
 yp1_dTdy_4xCO2 = new((/24/), float)
 
 do i=0,11
    xp1_dTdy_4xCO2(i)          = i+1
    xp1_dTdy_4xCO2(2*12-1-i)   = i+1
    yp1_dTdy_4xCO2(i)          = dTdy3_index_ave(i) + dTdy3_index_std(i)/2.
    yp1_dTdy_4xCO2(2*12-1-i)   = dTdy3_index_ave(i) - dTdy3_index_std(i)/2. 
  end do
;==================================================
xp1_dUdy_PD = new((/24/), float)
yp1_dUdy_PD = new((/24/), float)

do i=0,11
   xp1_dUdy_PD(i)          = i+1
   xp1_dUdy_PD(2*12-1-i)   = i+1
   yp1_dUdy_PD(i)          = dUdy1_index_ave(i) + dUdy1_index_std(i)/2.
   yp1_dUdy_PD(2*12-1-i)   = dUdy1_index_ave(i) - dUdy1_index_std(i)/2. 
 end do

 xp1_dUdy_4xCO2 = new((/24/), float)
 yp1_dUdy_4xCO2 = new((/24/), float)
 
 do i=0,11
    xp1_dUdy_4xCO2(i)          = i+1
    xp1_dUdy_4xCO2(2*12-1-i)   = i+1
    yp1_dUdy_4xCO2(i)          = dUdy3_index_ave(i) + dUdy3_index_std(i)/2.
    yp1_dUdy_4xCO2(2*12-1-i)   = dUdy3_index_ave(i) - dUdy3_index_std(i)/2. 
  end do
;--------------------------------------------------------
xp2_dUdy_PD = new((/24/), float)
yp2_dUdy_PD = new((/24/), float)

do i=0,11
   xp2_dUdy_PD(i)          = i+1
   xp2_dUdy_PD(2*12-1-i)   = i+1
   yp2_dUdy_PD(i)          = dUdy01_index_ave(i) + dUdy01_index_std(i)/2.
   yp2_dUdy_PD(2*12-1-i)   = dUdy01_index_ave(i) - dUdy01_index_std(i)/2. 
 end do

 xp2_dUdy_4xCO2 = new((/24/), float)
 yp2_dUdy_4xCO2 = new((/24/), float)
 
 do i=0,11
    xp2_dUdy_4xCO2(i)          = i+1
    xp2_dUdy_4xCO2(2*12-1-i)   = i+1
    yp2_dUdy_4xCO2(i)          = dUdy03_index_ave(i) + dUdy03_index_std(i)/2.
    yp2_dUdy_4xCO2(2*12-1-i)   = dUdy03_index_ave(i) - dUdy03_index_std(i)/2. 
  end do
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

delete(filelist)
sst_PD    = wgt_areaave_Wrap(sst(0,:,{0:5},{200:260}), 1.0, 1.0, 0)
sst_2xCO2 = wgt_areaave_Wrap(sst(1,:,{0:5},{200:260}), 1.0, 1.0, 0)
sst_4xCO2 = wgt_areaave_Wrap(sst(2,:,{0:5},{200:260}), 1.0, 1.0, 0)

ssh_PD    = wgt_areaave_Wrap(ssh(0,:,{3:8},{210:250}), 1.0, 1.0, 0)
ssh_2xCO2 = wgt_areaave_Wrap(ssh(1,:,{3:8},{210:250}), 1.0, 1.0, 0)
ssh_4xCO2 = wgt_areaave_Wrap(ssh(2,:,{3:8},{210:250}), 1.0, 1.0, 0)

va1_PD    = wgt_areaave_Wrap(va(0,:,{3:6},{210:250}), 1.0, 1.0, 0)
va1_2xCO2 = wgt_areaave_Wrap(va(1,:,{3:6},{210:250}), 1.0, 1.0, 0)
va1_4xCO2 = wgt_areaave_Wrap(va(2,:,{3:6},{210:250}), 1.0, 1.0, 0)

va2_PD    = wgt_areaave_Wrap(va(0,:,{-1:1},{220:260}), 1.0, 1.0, 0)
va2_2xCO2 = wgt_areaave_Wrap(va(1,:,{-1:1},{220:260}), 1.0, 1.0, 0)
va2_4xCO2 = wgt_areaave_Wrap(va(2,:,{-1:1},{210:250}), 1.0, 1.0, 0)

eke_PD    = wgt_areaave_Wrap(eke(0,:,{0:5},{210:250}), 1.0, 1.0, 0)
eke_2xCO2 = wgt_areaave_Wrap(eke(1,:,{0:5},{210:250}), 1.0, 1.0, 0)
eke_4xCO2 = wgt_areaave_Wrap(eke(2,:,{0:5},{210:250}), 1.0, 1.0, 0)

ntim  = dimsizes(sst&time)
delete([/sst, ssh, va, eke/]);
nyears = ntim/12
sst_PD_ann  = reshape(sst_PD,(/nyears,12/))
sst_4xCO2_ann  = reshape(sst_4xCO2,(/nyears,12/))

ssh_PD_ann  = reshape(ssh_PD,(/nyears,12/))
ssh_4xCO2_ann  = reshape(ssh_4xCO2,(/nyears,12/))

va1_PD_ann  = reshape(va1_PD,(/nyears,12/))
va1_4xCO2_ann  = reshape(va1_4xCO2,(/nyears,12/))

va2_PD_ann  = reshape(va2_PD,(/nyears,12/))
va2_4xCO2_ann  = reshape(va2_4xCO2,(/nyears,12/))

eke_PD_ann  = reshape(eke_PD,(/nyears,12/))
eke_4xCO2_ann  = reshape(eke_4xCO2,(/nyears,12/))

; annual cycle
sst_PD_ann_ave = dim_avg_n_Wrap(sst_PD_ann, 0)
sst_4xCO2_ann_ave = dim_avg_n_Wrap(sst_4xCO2_ann, 0)

ssh_PD_ann_ave = dim_avg_n_Wrap(ssh_PD_ann, 0)
ssh_4xCO2_ann_ave = dim_avg_n_Wrap(ssh_4xCO2_ann, 0)

va1_PD_ann_ave = dim_avg_n_Wrap(va1_PD_ann, 0)
va1_4xCO2_ann_ave = dim_avg_n_Wrap(va1_4xCO2_ann, 0)

va2_PD_ann_ave = dim_avg_n_Wrap(va2_PD_ann, 0)
va2_4xCO2_ann_ave = dim_avg_n_Wrap(va2_4xCO2_ann, 0)

eke_PD_ann_ave = dim_avg_n_Wrap(eke_PD_ann, 0)
eke_4xCO2_ann_ave = dim_avg_n_Wrap(eke_4xCO2_ann, 0)

; std 
sst_PD_ann_std = dim_stddev_n_Wrap(sst_PD_ann, 0)
sst_4xCO2_ann_std = dim_stddev_n_Wrap(sst_4xCO2_ann, 0)

ssh_PD_ann_std = dim_stddev_n_Wrap(ssh_PD_ann, 0)
ssh_4xCO2_ann_std = dim_stddev_n_Wrap(ssh_4xCO2_ann, 0)

va1_PD_ann_std = dim_stddev_n_Wrap(va1_PD_ann, 0)
va1_4xCO2_ann_std = dim_stddev_n_Wrap(va1_4xCO2_ann, 0)

va2_PD_ann_std = dim_stddev_n_Wrap(va2_PD_ann, 0)
va2_4xCO2_ann_std = dim_stddev_n_Wrap(va2_4xCO2_ann, 0)

eke_PD_ann_std = dim_stddev_n_Wrap(eke_PD_ann, 0)
eke_4xCO2_ann_std = dim_stddev_n_Wrap(eke_4xCO2_ann, 0)

data01 = new((/2,12/), float)
data01(0,:) = sst_PD_ann_ave
data01(1,:) = sst_4xCO2_ann_ave

data02 = new((/2,12/), float)
data02(0,:) = ssh_PD_ann_ave
data02(1,:) = ssh_4xCO2_ann_ave

data03 = new((/2,12/), float) ; Rossby-mode 
data03(0,:) = va1_PD_ann_ave
data03(1,:) = va1_4xCO2_ann_ave

data04 = new((/2,12/), float) ; Yanai-mode
data04(0,:) = va2_PD_ann_ave
data04(1,:) = va2_4xCO2_ann_ave

data05 = new((/2,12/), float)
data05(0,:) = eke_PD_ann_ave
data05(1,:) = eke_4xCO2_ann_ave

;==========================================================
xp1_sst_PD = new((/24/), float)
yp1_sst_PD = new((/24/), float)

do i=0,11
   xp1_sst_PD(i)          = i+1
   xp1_sst_PD(2*12-1-i)   = i+1
   yp1_sst_PD(i)          = sst_PD_ann_ave(i) + sst_PD_ann_std(i)/2.
   yp1_sst_PD(2*12-1-i)   = sst_PD_ann_ave(i) - sst_PD_ann_std(i)/2. 
 end do

 xp1_sst_4xCO2 = new((/24/), float)
 yp1_sst_4xCO2 = new((/24/), float)
 
 do i=0,11
    xp1_sst_4xCO2(i)          = i+1
    xp1_sst_4xCO2(2*12-1-i)   = i+1
    yp1_sst_4xCO2(i)          = sst_4xCO2_ann_ave(i) + sst_4xCO2_ann_std(i)/2.
    yp1_sst_4xCO2(2*12-1-i)   = sst_4xCO2_ann_ave(i) - sst_4xCO2_ann_std(i)/2. 
  end do
;=======================================================
xp1_ssh_PD = new((/24/), float)
yp1_ssh_PD = new((/24/), float)

do i=0,11
   xp1_ssh_PD(i)          = i+1
   xp1_ssh_PD(2*12-1-i)   = i+1
   yp1_ssh_PD(i)          = ssh_PD_ann_ave(i) + ssh_PD_ann_std(i)/2.
   yp1_ssh_PD(2*12-1-i)   = ssh_PD_ann_ave(i) - ssh_PD_ann_std(i)/2. 
 end do

 xp1_ssh_4xCO2 = new((/24/), float)
 yp1_ssh_4xCO2 = new((/24/), float)
 
 do i=0,11
    xp1_ssh_4xCO2(i)          = i+1
    xp1_ssh_4xCO2(2*12-1-i)   = i+1
    yp1_ssh_4xCO2(i)          = ssh_4xCO2_ann_ave(i) + ssh_4xCO2_ann_std(i)/2.
    yp1_ssh_4xCO2(2*12-1-i)   = ssh_4xCO2_ann_ave(i) - ssh_4xCO2_ann_std(i)/2. 
  end do
;-------------------------------------------------------

xp1_va1_PD = new((/24/), float)
yp1_va1_PD = new((/24/), float)

do i=0,11
   xp1_va1_PD(i)          = i+1
   xp1_va1_PD(2*12-1-i)   = i+1
   yp1_va1_PD(i)          = va1_PD_ann_ave(i) + va1_PD_ann_std(i)/2.
   yp1_va1_PD(2*12-1-i)   = va1_PD_ann_ave(i) - va1_PD_ann_std(i)/2. 
 end do

 xp1_va1_4xCO2 = new((/24/), float)
 yp1_va1_4xCO2 = new((/24/), float)
 
 do i=0,11
    xp1_va1_4xCO2(i)          = i+1
    xp1_va1_4xCO2(2*12-1-i)   = i+1
    yp1_va1_4xCO2(i)          = va1_4xCO2_ann_ave(i) + va1_4xCO2_ann_std(i)/2.
    yp1_va1_4xCO2(2*12-1-i)   = va1_4xCO2_ann_ave(i) - va1_4xCO2_ann_std(i)/2. 
  end do
;=========================================================
xp1_va2_PD = new((/24/), float)
yp1_va2_PD = new((/24/), float)

do i=0,11
   xp1_va2_PD(i)          = i+1
   xp1_va2_PD(2*12-1-i)   = i+1
   yp1_va2_PD(i)          = va2_PD_ann_ave(i) + va2_PD_ann_std(i)/2.
   yp1_va2_PD(2*12-1-i)   = va2_PD_ann_ave(i) - va2_PD_ann_std(i)/2. 
 end do

 xp1_va2_4xCO2 = new((/24/), float)
 yp1_va2_4xCO2 = new((/24/), float)
 
 do i=0,11
    xp1_va2_4xCO2(i)          = i+1
    xp1_va2_4xCO2(2*12-1-i)   = i+1
    yp1_va2_4xCO2(i)          = va2_4xCO2_ann_ave(i) + va2_4xCO2_ann_std(i)/2.
    yp1_va2_4xCO2(2*12-1-i)   = va2_4xCO2_ann_ave(i) - va2_4xCO2_ann_std(i)/2. 
  end do
;========================================================== 
xp1_eke_PD = new((/24/), float)
yp1_eke_PD = new((/24/), float)

do i=0,11
   xp1_eke_PD(i)          = i+1
   xp1_eke_PD(2*12-1-i)   = i+1
   yp1_eke_PD(i)          = eke_PD_ann_ave(i) + eke_PD_ann_std(i)/2.
   yp1_eke_PD(2*12-1-i)   = eke_PD_ann_ave(i) - eke_PD_ann_std(i)/2. 
 end do

 xp1_eke_4xCO2 = new((/24/), float)
 yp1_eke_4xCO2 = new((/24/), float)
 
 do i=0,11
    xp1_eke_4xCO2(i)          = i+1
    xp1_eke_4xCO2(2*12-1-i)   = i+1
    yp1_eke_4xCO2(i)          = eke_4xCO2_ann_ave(i) + eke_4xCO2_ann_std(i)/2.
    yp1_eke_4xCO2(2*12-1-i)   = eke_4xCO2_ann_ave(i) - eke_4xCO2_ann_std(i)/2. 
  end do
;=========================================================================================

;SON TIWs intensity 
sst_PD_SON    = (sst_PD(8::12)    + sst_PD(9::12)   + sst_PD(10::12))/3.
sst_2xCO2_SON = (sst_2xCO2(8::12) + sst_2xCO2(9::12)+ sst_2xCO2(10::12))/3.
sst_4xCO2_SON = (sst_4xCO2(8::12) + sst_4xCO2(9::12)+ sst_4xCO2(10::12))/3.

ssh_PD_SON    = (ssh_PD(8::12)    + ssh_PD(9::12)   + ssh_PD(10::12))/3.
ssh_2xCO2_SON = (ssh_2xCO2(8::12) + ssh_2xCO2(9::12)+ ssh_2xCO2(10::12))/3.
ssh_4xCO2_SON = (ssh_4xCO2(8::12) + ssh_4xCO2(9::12)+ ssh_4xCO2(10::12))/3.

va1_PD_SON    = (va1_PD(8::12)    + va1_PD(9::12)   + va1_PD(10::12))/3.
va1_2xCO2_SON = (va1_2xCO2(8::12) + va1_2xCO2(9::12)+ va1_2xCO2(10::12))/3.
va1_4xCO2_SON = (va1_4xCO2(8::12) + va1_4xCO2(9::12)+ va1_4xCO2(10::12))/3.

va2_PD_SON    = (va2_PD(8::12)    + va2_PD(9::12)   + va2_PD(10::12))/3.
va2_2xCO2_SON = (va2_2xCO2(8::12) + va2_2xCO2(9::12)+ va2_2xCO2(10::12))/3.
va2_4xCO2_SON = (va2_4xCO2(8::12) + va2_4xCO2(9::12)+ va2_4xCO2(10::12))/3.

eke_PD_SON    = (eke_PD(8::12)    + eke_PD(9::12)   + eke_PD(10::12))/3.
eke_2xCO2_SON = (eke_2xCO2(8::12) + eke_2xCO2(9::12)+ eke_2xCO2(10::12))/3.
eke_4xCO2_SON = (eke_4xCO2(8::12) + eke_4xCO2(9::12)+ eke_4xCO2(10::12))/3.

; MAM TIWs activity
sst_PD_MAM    = (sst_PD(2::12)    + sst_PD(3::12)   + sst_PD(4::12))/3.
sst_2xCO2_MAM = (sst_2xCO2(2::12) + sst_2xCO2(3::12)+ sst_2xCO2(4::12))/3.
sst_4xCO2_MAM = (sst_4xCO2(2::12) + sst_4xCO2(3::12)+ sst_4xCO2(4::12))/3.

ssh_PD_MAM    = (ssh_PD(2::12)    + ssh_PD(3::12)   + ssh_PD(4::12))/3.
ssh_2xCO2_MAM = (ssh_2xCO2(2::12) + ssh_2xCO2(3::12)+ ssh_2xCO2(4::12))/3.
ssh_4xCO2_MAM = (ssh_4xCO2(2::12) + ssh_4xCO2(3::12)+ ssh_4xCO2(4::12))/3.

va1_PD_MAM    = (va1_PD(2::12)    + va1_PD(3::12)   + va1_PD(4::12))/3.
va1_2xCO2_MAM = (va1_2xCO2(2::12) + va1_2xCO2(3::12)+ va1_2xCO2(4::12))/3.
va1_4xCO2_MAM = (va1_4xCO2(2::12) + va1_4xCO2(3::12)+ va1_4xCO2(4::12))/3.

va2_PD_MAM    = (va2_PD(2::12)    + va2_PD(3::12)   + va2_PD(4::12))/3.
va2_2xCO2_MAM = (va2_2xCO2(2::12) + va2_2xCO2(3::12)+ va2_2xCO2(4::12))/3.
va2_4xCO2_MAM = (va2_4xCO2(2::12) + va2_4xCO2(3::12)+ va2_4xCO2(4::12))/3.

eke_PD_MAM    = (eke_PD(2::12)    + eke_PD(3::12)   + eke_PD(4::12))/3.
eke_2xCO2_MAM = (eke_2xCO2(2::12) + eke_2xCO2(3::12)+ eke_2xCO2(4::12))/3.
eke_4xCO2_MAM = (eke_4xCO2(2::12) + eke_4xCO2(3::12)+ eke_4xCO2(4::12))/3.

;=========================================================================================
stat  = 0                         ; 0=mean; 1=variance; 2=std. dev; ...
n     = 30                        ; sub-sample size;  default is n=N
nBoot = 10000                      ; bootstrap replications
N     = 30
opt   = False                     ; use defaults
if (isvar("n") .and. (n.lt.N)) then
    opt   = True
    opt@sample_size = n           ; default N
end if

 BootStrap_sst1_PD     = bootstrap_stat(sst_PD_SON, stat, nBoot, 0, opt)
 BootStrap_sst1_2xCO2  = bootstrap_stat(sst_2xCO2_SON, stat, nBoot, 0, opt)
 BootStrap_sst1_4xCO2  = bootstrap_stat(sst_4xCO2_SON, stat, nBoot, 0, opt)

 BootStrap_ssh1_PD     = bootstrap_stat(ssh_PD_SON, stat, nBoot, 0, opt)
 BootStrap_ssh1_2xCO2  = bootstrap_stat(ssh_2xCO2_SON, stat, nBoot, 0, opt)
 BootStrap_ssh1_4xCO2  = bootstrap_stat(ssh_4xCO2_SON, stat, nBoot, 0, opt)

 BootStrap_va1_PD     = bootstrap_stat(va1_PD_SON, stat, nBoot, 0, opt)
 BootStrap_va1_2xCO2  = bootstrap_stat(va1_2xCO2_SON, stat, nBoot, 0, opt)
 BootStrap_va1_4xCO2  = bootstrap_stat(va1_4xCO2_SON, stat, nBoot, 0, opt)

 BootStrap_va2_PD     = bootstrap_stat(va2_PD_SON, stat, nBoot, 0, opt)
 BootStrap_va2_2xCO2  = bootstrap_stat(va2_2xCO2_SON, stat, nBoot, 0, opt)
 BootStrap_va2_4xCO2  = bootstrap_stat(va2_4xCO2_SON, stat, nBoot, 0, opt)

 BootStrap_eke1_PD     = bootstrap_stat(eke_PD_SON, stat, nBoot, 0, opt)
 BootStrap_eke1_2xCO2  = bootstrap_stat(eke_2xCO2_SON, stat, nBoot, 0, opt)
 BootStrap_eke1_4xCO2  = bootstrap_stat(eke_4xCO2_SON, stat, nBoot, 0, opt)
 ;=====================================================
 xBoot_sst1_PD       = BootStrap_sst1_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sstavg1_PD    = BootStrap_sst1_PD[1]       ; Mean of bootstrapped estimates

 xBoot_ssh1_PD       = BootStrap_ssh1_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sshavg1_PD    = BootStrap_ssh1_PD[1]       ; Mean of bootstrapped estimates

 xBoot_va1_PD       = BootStrap_va1_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_vaavg1_PD    = BootStrap_va1_PD[1]       ; Mean of bootstrapped estimates

 xBoot_va2_PD       = BootStrap_va2_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_vaavg2_PD    = BootStrap_va2_PD[1]       ; Mean of bootstrapped estimates

 xBoot_eke1_PD       = BootStrap_eke1_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_ekeavg1_PD    = BootStrap_eke1_PD[1]       ; Mean of bootstrapped estimates
;=====================================================
    xBoot_sst1_2xCO2    = BootStrap_sst1_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sstavg1_2xCO2    = BootStrap_sst1_2xCO2[1]       ; Mean of bootstrapped estimates

    xBoot_ssh1_2xCO2    = BootStrap_ssh1_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sshavg1_2xCO2    = BootStrap_ssh1_2xCO2[1]       ; Mean of bootstrapped estimates

    xBoot_va1_2xCO2    = BootStrap_va1_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_vaavg1_2xCO2    = BootStrap_va1_2xCO2[1]       ; Mean of bootstrapped estimates

    xBoot_va2_2xCO2    = BootStrap_va2_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_vaavg2_2xCO2    = BootStrap_va2_2xCO2[1]       ; Mean of bootstrapped estimates

    xBoot_eke1_2xCO2    = BootStrap_eke1_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_ekeavg1_2xCO2    = BootStrap_eke1_2xCO2[1]       ; Mean of bootstrapped estimates
;-----------------------------------------------------
    xBoot_sst1_4xCO2    = BootStrap_sst1_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sstavg1_4xCO2    = BootStrap_sst1_4xCO2[1]       ; Mean of bootstrapped estimates

    xBoot_ssh1_4xCO2    = BootStrap_ssh1_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sshavg1_4xCO2    = BootStrap_ssh1_4xCO2[1]       ; Mean of bootstrapped estimates
 
    xBoot_va1_4xCO2    = BootStrap_va1_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_vaavg1_4xCO2    = BootStrap_va1_4xCO2[1]       ; Mean of bootstrapped estimates

    xBoot_va2_4xCO2    = BootStrap_va2_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_vaavg2_4xCO2    = BootStrap_va2_4xCO2[1]       ; Mean of bootstrapped estimates

    xBoot_eke1_4xCO2    = BootStrap_eke1_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_ekeavg1_4xCO2    = BootStrap_eke1_4xCO2[1]       ; Mean of bootstrapped estimates
;=====================================================
 xBootLow_sst1_PD    = bootstrap_estimate(xBoot_sst1_PD, 0.025, False)   ;  2.5% lower confidence bound
  xBootHi_sst1_PD    = bootstrap_estimate(xBoot_sst1_PD, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_ssh1_PD    = bootstrap_estimate(xBoot_ssh1_PD, 0.025, False)   ;  2.5% lower confidence bound
  xBootHi_ssh1_PD    = bootstrap_estimate(xBoot_ssh1_PD, 0.975, False)   ; 97.5% upper confidence bound

  xBootLow_va1_PD    = bootstrap_estimate(xBoot_va1_PD, 0.025, False)   ;  2.5% lower confidence bound
   xBootHi_va1_PD    = bootstrap_estimate(xBoot_va1_PD, 0.975, False)   ; 97.5% upper confidence bound
 
  xBootLow_va2_PD    = bootstrap_estimate(xBoot_va2_PD, 0.025, False)   ;  2.5% lower confidence bound
   xBootHi_va2_PD    = bootstrap_estimate(xBoot_va2_PD, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_eke1_PD    = bootstrap_estimate(xBoot_eke1_PD, 0.025, False)   ;  2.5% lower confidence bound
  xBootHi_eke1_PD    = bootstrap_estimate(xBoot_eke1_PD, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------
 xBootLow_sst1_2xCO2    = bootstrap_estimate(xBoot_sst1_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
  xBootHi_sst1_2xCO2    = bootstrap_estimate(xBoot_sst1_2xCO2, 0.975, False)   ; 97.5% upper confidence bound
 
 xBootLow_ssh1_2xCO2    = bootstrap_estimate(xBoot_ssh1_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
  xBootHi_ssh1_2xCO2    = bootstrap_estimate(xBoot_ssh1_2xCO2, 0.975, False)   ; 97.5% upper confidence bound

  xBootLow_va1_2xCO2    = bootstrap_estimate(xBoot_va1_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
   xBootHi_va1_2xCO2    = bootstrap_estimate(xBoot_va1_2xCO2, 0.975, False)   ; 97.5% upper confidence bound
 
  xBootLow_va2_2xCO2    = bootstrap_estimate(xBoot_va2_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
   xBootHi_va2_2xCO2    = bootstrap_estimate(xBoot_va2_2xCO2, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_eke1_2xCO2    = bootstrap_estimate(xBoot_eke1_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
  xBootHi_eke1_2xCO2    = bootstrap_estimate(xBoot_eke1_2xCO2, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------
xBootLow_sst1_4xCO2    = bootstrap_estimate(xBoot_sst1_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_sst1_4xCO2    = bootstrap_estimate(xBoot_sst1_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_ssh1_4xCO2    = bootstrap_estimate(xBoot_ssh1_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_ssh1_4xCO2    = bootstrap_estimate(xBoot_ssh1_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_va1_4xCO2    = bootstrap_estimate(xBoot_va1_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
  xBootHi_va1_4xCO2    = bootstrap_estimate(xBoot_va1_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_va2_4xCO2    = bootstrap_estimate(xBoot_va2_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
  xBootHi_va2_4xCO2    = bootstrap_estimate(xBoot_va2_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_eke1_4xCO2    = bootstrap_estimate(xBoot_eke1_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_eke1_4xCO2    = bootstrap_estimate(xBoot_eke1_4xCO2, 0.975, False)   ; 97.5% upper confidence bound
;---------------------------------------------------------------------------------------------------------
;MAM
BootStrap_sst2_PD     = bootstrap_stat(sst_PD_MAM, stat, nBoot, 0, opt)
BootStrap_sst2_2xCO2  = bootstrap_stat(sst_2xCO2_MAM, stat, nBoot, 0, opt)
BootStrap_sst2_4xCO2  = bootstrap_stat(sst_4xCO2_MAM, stat, nBoot, 0, opt)

BootStrap_ssh2_PD     = bootstrap_stat(ssh_PD_MAM, stat, nBoot, 0, opt)
BootStrap_ssh2_2xCO2  = bootstrap_stat(ssh_2xCO2_MAM, stat, nBoot, 0, opt)
BootStrap_ssh2_4xCO2  = bootstrap_stat(ssh_4xCO2_MAM, stat, nBoot, 0, opt)

BootStrap_va01_PD     = bootstrap_stat(va1_PD_MAM, stat, nBoot, 0, opt)
BootStrap_va01_2xCO2  = bootstrap_stat(va1_2xCO2_MAM, stat, nBoot, 0, opt)
BootStrap_va01_4xCO2  = bootstrap_stat(va1_4xCO2_MAM, stat, nBoot, 0, opt)

BootStrap_va02_PD     = bootstrap_stat(va2_PD_MAM, stat, nBoot, 0, opt)
BootStrap_va02_2xCO2  = bootstrap_stat(va2_2xCO2_MAM, stat, nBoot, 0, opt)
BootStrap_va02_4xCO2  = bootstrap_stat(va2_4xCO2_MAM, stat, nBoot, 0, opt)

BootStrap_eke2_PD     = bootstrap_stat(eke_PD_MAM, stat, nBoot, 0, opt)
BootStrap_eke2_2xCO2  = bootstrap_stat(eke_2xCO2_MAM, stat, nBoot, 0, opt)
BootStrap_eke2_4xCO2  = bootstrap_stat(eke_4xCO2_MAM, stat, nBoot, 0, opt)
;=====================================================
xBoot_sst2_PD       = BootStrap_sst2_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_sstavg2_PD    = BootStrap_sst2_PD[1]       ; Mean of bootstrapped estimates

xBoot_ssh2_PD       = BootStrap_ssh2_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_sshavg2_PD    = BootStrap_ssh2_PD[1]       ; Mean of bootstrapped estimates

xBoot_va01_PD       = BootStrap_va01_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_vaavg01_PD    = BootStrap_va01_PD[1]       ; Mean of bootstrapped estimates

xBoot_va02_PD       = BootStrap_va02_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_vaavg02_PD    = BootStrap_va02_PD[1]       ; Mean of bootstrapped estimates

xBoot_eke2_PD       = BootStrap_eke2_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_ekeavg2_PD    = BootStrap_eke2_PD[1]       ; Mean of bootstrapped estimates
;=====================================================
   xBoot_sst2_2xCO2    = BootStrap_sst2_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_sstavg2_2xCO2    = BootStrap_sst2_2xCO2[1]       ; Mean of bootstrapped estimates

   xBoot_ssh2_2xCO2    = BootStrap_ssh2_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_sshavg2_2xCO2    = BootStrap_ssh2_2xCO2[1]       ; Mean of bootstrapped estimates

   xBoot_va01_2xCO2    = BootStrap_va01_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_vaavg01_2xCO2    = BootStrap_va01_2xCO2[1]       ; Mean of bootstrapped estimates

   xBoot_va02_2xCO2    = BootStrap_va02_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_vaavg02_2xCO2    = BootStrap_va02_2xCO2[1]       ; Mean of bootstrapped estimates

   xBoot_eke2_2xCO2    = BootStrap_eke2_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_ekeavg2_2xCO2    = BootStrap_eke2_2xCO2[1]       ; Mean of bootstrapped estimates
;-----------------------------------------------------
   xBoot_sst2_4xCO2    = BootStrap_sst2_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_sstavg2_4xCO2    = BootStrap_sst2_4xCO2[1]       ; Mean of bootstrapped estimates

   xBoot_ssh2_4xCO2    = BootStrap_ssh2_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_sshavg2_4xCO2    = BootStrap_ssh2_4xCO2[1]       ; Mean of bootstrapped estimates

   xBoot_va01_4xCO2    = BootStrap_va01_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_vaavg01_4xCO2    = BootStrap_va01_4xCO2[1]       ; Mean of bootstrapped estimates

   xBoot_va02_4xCO2    = BootStrap_va02_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_vaavg02_4xCO2    = BootStrap_va02_4xCO2[1]       ; Mean of bootstrapped estimates

   xBoot_eke2_4xCO2    = BootStrap_eke2_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_ekeavg2_4xCO2    = BootStrap_eke2_4xCO2[1]       ; Mean of bootstrapped estimates
;=====================================================
xBootLow_sst2_PD    = bootstrap_estimate(xBoot_sst2_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_sst2_PD    = bootstrap_estimate(xBoot_sst2_PD, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_ssh2_PD    = bootstrap_estimate(xBoot_ssh2_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_ssh2_PD    = bootstrap_estimate(xBoot_ssh2_PD, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_va01_PD    = bootstrap_estimate(xBoot_va01_PD, 0.025, False)   ;  2.5% lower confidence bound
  xBootHi_va01_PD    = bootstrap_estimate(xBoot_va01_PD, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_va02_PD    = bootstrap_estimate(xBoot_va02_PD, 0.025, False)   ;  2.5% lower confidence bound
  xBootHi_va02_PD    = bootstrap_estimate(xBoot_va02_PD, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_eke2_PD    = bootstrap_estimate(xBoot_eke2_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_eke2_PD    = bootstrap_estimate(xBoot_eke2_PD, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------
xBootLow_sst2_2xCO2    = bootstrap_estimate(xBoot_sst2_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_sst2_2xCO2    = bootstrap_estimate(xBoot_sst2_2xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_ssh2_2xCO2    = bootstrap_estimate(xBoot_ssh2_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_ssh2_2xCO2    = bootstrap_estimate(xBoot_ssh2_2xCO2, 0.975, False)   ; 97.5% upper confidence bound

    xBootLow_va01_2xCO2    = bootstrap_estimate(xBoot_va01_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
     xBootHi_va01_2xCO2    = bootstrap_estimate(xBoot_va01_2xCO2, 0.975, False)   ; 97.5% upper confidence bound

   xBootLow_va02_2xCO2    = bootstrap_estimate(xBoot_va02_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
    xBootHi_va02_2xCO2    = bootstrap_estimate(xBoot_va02_2xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_eke2_2xCO2    = bootstrap_estimate(xBoot_eke2_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_eke2_2xCO2    = bootstrap_estimate(xBoot_eke2_2xCO2, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------
xBootLow_sst2_4xCO2    = bootstrap_estimate(xBoot_sst2_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_sst2_4xCO2    = bootstrap_estimate(xBoot_sst2_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_ssh2_4xCO2    = bootstrap_estimate(xBoot_ssh2_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_ssh2_4xCO2    = bootstrap_estimate(xBoot_ssh2_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_va01_4xCO2    = bootstrap_estimate(xBoot_va01_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
  xBootHi_va01_4xCO2    = bootstrap_estimate(xBoot_va01_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

  xBootLow_va02_4xCO2    = bootstrap_estimate(xBoot_va02_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
   xBootHi_va02_4xCO2    = bootstrap_estimate(xBoot_va02_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_eke2_4xCO2    = bootstrap_estimate(xBoot_eke2_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_eke2_4xCO2    = bootstrap_estimate(xBoot_eke2_4xCO2, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------------------------------------------------------------
; boreal <MAM>
stdarr1 = new((/3/),"float") ;  SST
stdarr1(0) =  xBoot_sstavg2_PD
stdarr1(1) =  xBoot_sstavg2_2xCO2 
stdarr1(2) =  xBoot_sstavg2_4xCO2
printMinMax(stdarr1, 0)
stdarr2 = new((/3/),"float") ; SSH
stdarr2(0) =  xBoot_sshavg2_PD
stdarr2(1) =  xBoot_sshavg2_2xCO2 
stdarr2(2) =  xBoot_sshavg2_4xCO2

stdarr3 = new((/3/),"float") ; EKE 
stdarr3(0) =  xBoot_ekeavg2_PD
stdarr3(1) =  xBoot_ekeavg2_2xCO2 
stdarr3(2) =  xBoot_ekeavg2_4xCO2 

stdarr4 = new((/3/),"float") ; va1
stdarr4(0) =  xBoot_vaavg01_PD
stdarr4(1) =  xBoot_vaavg01_2xCO2 
stdarr4(2) =  xBoot_vaavg01_4xCO2 

stdarr5 = new((/3/),"float") ; va2 
stdarr5(0) =  xBoot_vaavg02_PD
stdarr5(1) =  xBoot_vaavg02_2xCO2 
stdarr5(2) =  xBoot_vaavg02_4xCO2 
;===========================================================
; boreal SON
stdarr01 = new((/3/),"float") ;  SST
stdarr01(0) =  xBoot_sstavg1_PD
stdarr01(1) =  xBoot_sstavg1_2xCO2 
stdarr01(2) =  xBoot_sstavg1_4xCO2

stdarr02 = new((/3/),"float") ; SSH
stdarr02(0) =  xBoot_sshavg1_PD
stdarr02(1) =  xBoot_sshavg1_2xCO2 
stdarr02(2) =  xBoot_sshavg1_4xCO2

stdarr03 = new((/3/),"float") ; EKE 
stdarr03(0) =  xBoot_ekeavg1_PD
stdarr03(1) =  xBoot_ekeavg1_2xCO2 
stdarr03(2) =  xBoot_ekeavg1_4xCO2 

stdarr04 = new((/3/),"float") ; va1
stdarr04(0) =  xBoot_vaavg1_PD
stdarr04(1) =  xBoot_vaavg1_2xCO2 
stdarr04(2) =  xBoot_vaavg1_4xCO2 

stdarr05 = new((/3/),"float") ; va2 
stdarr05(0) =  xBoot_vaavg2_PD
stdarr05(1) =  xBoot_vaavg2_2xCO2 
stdarr05(2) =  xBoot_vaavg2_4xCO2 

;===========================================================================================================
outdir = "/Users/xueaoyun/Documents/work/TIW_future/revision1/main_figure/"
wks = gsn_open_wks("pdf",outdir+"Figure3_plot_CESM_PD_2xCO2_4xCO2_zonal_current_TIW_annual_cycle_new")
wks_type = "pdf"
yspace   = 0.05
ares               = True 
ares@gsnFrame      = False
ares@gsnDraw       = False
ares@trXMinF       = 0.
ares@trXMaxF       = 13.
ares@trYMinF       = -1.
ares@trYMaxF       = 6.
ares@vpXF          = 0.10                  ; Move slightly to left
ares@vpYF          = 0.9
ares@vpHeightF     = 0.1    
ares@vpWidthF      = 0.25 

ares@gsnYRefLine              = 0. 
ares@gsnYRefLineColors        = (/"black"/)
ares@gsnYRefLineThicknessF    = 1.0
ares@gsnYRefLineDashPattern   = 2 
;ares@gsnXRefLine              = 11.0 
;ares@gsnXRefLineColors        = (/"black"/)
;ares@gsnXRefLineThicknessF    = 2.0
;ares@gsnXRefLineDashPattern   = 1
ares@tfPolyDrawOrder          = "PreDraw" 
ares@tmXBMode                 = "Explicit"                 
ares@tmXBValues               = ispan(0,13,1)
ares@tmXBLabels               = (/"","J","F","M","A","M","J","J","A","S","O","N","D",""/)  
ares@tmYLMode                 = "Explicit"   
ares@tmYLValues               = (/-1.0,0,1.0,2.0,3.0,4.0,5.0,6.0,7.0,8.0/)
ares@tmYLLabels               = (/"-1","0","1","2","3","4","5","6","7","8"/)  
ares@tmXBLabelFontHeightF     = 0.001
ares@xyMarkLineModes          = (/"MarkLines","MarkLines"/)  
ares@xyMarkers                = (/16,16/)
ares@xyLineColors             = (/"blue","red"/)
ares@xyDashPatterns           = (/0/)
ares@xyLineThicknesses        = (/2.0,2.0/)
ares@xyMarkerSizes            = (/0.005,0.005/) 
ares@tmXBLabelFontHeightF     = 0.008
ares@tmYLLabelFontHeightF     = 0.008
ares@tmXBLabelFont            = 21
ares@tmYLLabelFont            = 21
ares@tmXBLabelFontThicknessF  = 2.0
ares@tmYLLabelFontThicknessF  = 2.0 
ares@tmBorderThicknessF       = 2.5
ares@tmXBMajorLengthF         = 0.005
ares@tmYLMajorLengthF         = 0.005
ares@tmXBMinorLengthF         = 0.006
ares@tmYLMinorLengthF         = 0.006 
ares@tmXBMajorThicknessF      = 2.0
ares@tmYLMajorThicknessF      = 2.0
ares@tmXBMinorThicknessF      = 1.5
ares@tmYLMinorThicknessF      = 1.5 
ares@tmXBLabelDeltaF          = -0.2
ares@tmYLLabelDeltaF          = -0.2 
ares@tiYAxisString            = " "
ares@tiXAxisString            = " "
ares@tiMainString             = " "
ares@tmXBLabelAngleF          = 0
ares@pmLegendDisplayMode    = "Always"            ; turn on legend
ares@pmLegendSide           = "Top"               ; Change location of 
ares@pmLegendZone           = 0        
ares@pmLegendParallelPosF   = 0.30                ; move units right
ares@pmLegendOrthogonalPosF = -0.35                ; move units down
ares@pmLegendWidthF         = 0.04                ; Change width and
ares@pmLegendHeightF        = 0.03                ; height of legend.
ares@lgPerimOn              = False               ; turn off box around
ares@lgLabelFontHeightF     = 0.008              ; label font height
ares@lgLabelFont            = 21                  ; label font height
ares@lgLabelFontThicknessF  = 3                   ; label font height
ares@xyExplicitLegendLabels = (/"PD","4xCO2"/)  
ares@tmXTOn                      = False 
ares@tmYROn                      = False 
ares@tmXMajorGrid                = False 
ares@tmXMajorGridLineDashPattern =  2
ares@tmXMajorGridLineColor       = "grey40"
ares@tmXMajorGridThicknessF      = 0.8
ares@tmYMajorGrid                = True
ares@tmYMajorGridLineDashPattern =  2
ares@tmYMajorGridLineColor       = "grey40"
ares@tmYMajorGridThicknessF      = 0.8 
ares@tmYMinorGrid                = False  
ares@tmYMinorGridLineDashPattern =  2
ares@tmYMinorGridLineColor       = "grey40"
ares@tmYMinorGridThicknessF      = 0.3 
ares@tmXMinorGrid                = False 
ares@tmXMinorGridLineDashPattern =  2
ares@tmXMinorGridLineColor       = "grey40"
ares@tmXMinorGridThicknessF      = 0.3 
;ares@tiMainString                ="northward jump of WNPSH"
ares@tiMainFontHeightF =0.02 
ares@gsnLeftString      ="~F22~a "
ares@gsnLeftStringFontHeightF =0.012
ares@gsnCenterString   = "~F22~Meridional SST gradient"
ares@gsnCenterStringFontHeightF = 0.010
x = ispan(1,12,1)
plot0 = gsn_csm_xy(wks,x,data1,ares)
delete(ares@tmYLValues)
delete(ares@tmYLLabels)
ares@pmLegendParallelPosF   = -0.30                ; move units right
ares@trYMinF       = -2.
ares@trYMaxF       = 1.
ares@vpXF          = 0.40              ; Move slightly to left
ares@vpYF          = ares@vpYF
ares@vpHeightF     = 0.1    
ares@vpWidthF      = 0.25 
ares@tmYLValues               = (/-2.0,-1.0,0,1,2,3,4/)
ares@tmYLLabels               = (/"-2","-1","0","1","2","3","4"/) 
ares@gsnLeftString      ="~F22~b "
ares@gsnCenterString    = "~F22~EQ Shear (EUC & nSEC)"
plot1 = gsn_csm_xy(wks,x,data2,ares)
delete(ares@tmYLValues)
delete(ares@tmYLLabels)
ares@pmLegendParallelPosF   = 0.30                ; move units right
ares@trYMinF       = -2.
ares@trYMaxF       = 4.
ares@vpXF          = 0.70              ; Move slightly to left
ares@vpYF          = ares@vpYF
ares@vpHeightF     = 0.1    
ares@vpWidthF      = 0.25 
ares@tmYLValues               = (/-2.0,-1.0,0,1,2,3,4/)
ares@tmYLLabels               = (/"-2","-1","0","1","2","3","4"/)  
ares@gsnLeftString      ="~F22~c "
ares@gsnCenterString    = "~F22~NEQ Shear (NECC & nSEC)"
plot2 = gsn_csm_xy(wks,x,data3,ares)

gsres                   = True                        ; poly res
gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
gsres@gsFillColor       = "royalblue1"                 ; color chosen
ummy1 = gsn_add_polygon (wks,plot0,xp1_dTdy_PD,yp1_dTdy_PD,gsres)

gsres                   = True                        ; poly res
;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
gsres@gsFillColor       = "pink1"                 ; color chosen
gsres@gsFillOpacityF     = 0.5
ummy2 = gsn_add_polygon (wks,plot0,xp1_dTdy_4xCO2,yp1_dTdy_4xCO2,gsres)

gsres                   = True                        ; poly res
gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
gsres@gsFillColor       = "royalblue1"                 ; color chosen
ummy01 = gsn_add_polygon (wks,plot1,xp2_dUdy_PD,yp2_dUdy_PD,gsres)

gsres                   = True                        ; poly res
;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
gsres@gsFillColor       = "pink1"                 ; color chosen
gsres@gsFillOpacityF     = 0.5
ummy02 = gsn_add_polygon (wks,plot1,xp2_dUdy_4xCO2,yp2_dUdy_4xCO2,gsres)

gsres                   = True                        ; poly res
gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
gsres@gsFillColor       = "royalblue1"                 ; color chosen
ummy001 = gsn_add_polygon (wks,plot2,xp1_dUdy_PD,yp1_dUdy_PD,gsres)

gsres                   = True                        ; poly res
;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
gsres@gsFillColor       = "pink1"                 ; color chosen
gsres@gsFillOpacityF     = 0.5
ummy002 = gsn_add_polygon (wks,plot2,xp1_dUdy_4xCO2,yp1_dUdy_4xCO2,gsres)
;==================================================================================
yspace1 = 0.03
bres               = True 
bres@gsnFrame      = False
bres@gsnDraw       = False
bres@trXMinF       = 0.
bres@trXMaxF       = 13.
bres@trYMinF       = 0.
bres@trYMaxF       = 0.4
bres@vpXF          = 0.55                  ; Move slightly to left
bres@vpYF          = 0.73
bres@vpHeightF     = 0.1    
bres@vpWidthF      = 0.19 
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
bres@tmYLValues               = (/0,0.1,0.2,0.3,0.4,0.5,0.6/)
bres@tmYLLabels               = (/"0","0.1","0.2","0.3","0.4","0.5","0.6"/)  
bres@tmXBLabelFontHeightF     = 0.001
bres@xyMarkLineModes          = (/"MarkLines","MarkLines"/)  
bres@xyMarkers                = (/16,16/)
bres@xyLineColors             = (/"blue","red"/)
bres@xyDashPatterns           = (/0/)
bres@xyLineThicknesses        = (/2.0,2.0/)
bres@xyMarkerSizes            = (/0.005,0.005/) 
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
bres@pmLegendOrthogonalPosF = 0.35                ; move units down
bres@pmLegendWidthF         = 0.04                ; Change width and
bres@pmLegendHeightF        = 0.03                ; height of legend.
bres@lgPerimOn              = False               ; turn off box around
bres@lgLabelFontHeightF     = 0.008              ; label font height
bres@lgLabelFont            = 21                  ; label font height
bres@lgLabelFontThicknessF  = 2                   ; label font height
bres@xyExplicitLegendLabels = (/"PD","4xCO2"/)  
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
bres@gsnLeftString      ="~F22~f "
bres@gsnCenterString    ="~F22~TIW-VVEL (Rossby mode)"
bres@gsnLeftStringFontHeightF =0.012
bres@gsnCenterStringFontHeightF = 0.010
x = ispan(1,12,1)
plot3 = gsn_csm_xy(wks,x,data03,bres)
delete(bres@tmYLValues)
delete(bres@tmYLLabels)
bres@trYMinF       = 0.
bres@trYMaxF       = 0.3
bres@vpXF          = 0.32                   ; Move slightly to left
bres@vpYF          = bres@vpYF
bres@vpHeightF     = 0.1    
bres@vpWidthF      = 0.19 
bres@tmYLValues               = (/0,0.1,0.2,0.3,0.4,0.5,0.6/)
bres@tmYLLabels               = (/"0","0.1","0.2","0.3","0.4","0.5","0.6"/)  
bres@gsnLeftString      ="~F22~e"
bres@gsnCenterString    ="~F22~TIW-VVEL (Yanai mode)" 
plot4 = gsn_csm_xy(wks,x,data04,bres)
delete(bres@tmYLValues)
delete(bres@tmYLLabels)
bres@trYMinF       = 0.
bres@trYMaxF       = 0.6
bres@vpXF          = 0.10                ; Move slightly to left
bres@vpYF          = bres@vpYF
bres@vpHeightF     = 0.1    
bres@vpWidthF      = 0.19 
bres@tmYLValues               = (/0,0.1,0.2,0.3,0.4,0.5,0.6/)
bres@tmYLLabels               = (/"0","0.1","0.2","0.3","0.4","0.5","0.6"/)  
bres@gsnLeftString      ="~F22~d"
bres@gsnCenterString    ="~F22~TIW-SST"
plot5 = gsn_csm_xy(wks,x,data01,bres)

delete(bres@tmYLValues)
delete(bres@tmYLLabels)
bres@trYMinF       = 0.
bres@trYMaxF       = 5
bres@vpXF          = 0.77                ; Move slightly to left
bres@vpYF          = bres@vpYF
bres@vpHeightF     = 0.1    
bres@vpWidthF      = 0.19 
bres@tmYLValues               = (/0,1,2,3,4,5/)
bres@tmYLLabels               = (/"0","1","2","3","4","5"/) 
bres@gsnLeftString      ="~F22~g"
bres@gsnCenterString    ="~F22~TIW-SSH"
plot6 = gsn_csm_xy(wks,x,data02,bres)

;error_bar1 = new(12,graphic)
;error_bar2 = new(12,graphic)
;
;error_bar01 = new(12,graphic)
;error_bar02 = new(12,graphic)
;
;error_bar001 = new(12,graphic)
;error_bar002 = new(12,graphic)
;polyres                   = True                       ; marker resources
; do t = 0,11
;  error_bar1(t) = gsn_add_polyline(wks,plot0,(/t+1,t+1/),(/yp1_sst_PD(t),yp2_sst_PD(t)/),polyres)
;  error_bar2(t) = gsn_add_polyline(wks,plot0,(/t+1,t+1/),(/yp1_sst_4xCO2(t),yp2_sst_4xCO2(t)/),polyres)
;
;  error_bar01(t) = gsn_add_polyline(wks,plot1,(/t+1,t+1/),(/yp1_ssh_PD(t),yp2_ssh_PD(t)/),polyres)
;  error_bar02(t) = gsn_add_polyline(wks,plot1,(/t+1,t+1/),(/yp1_ssh_4xCO2(t),yp2_ssh_4xCO2(t)/),polyres)
;
;  error_bar001(t) = gsn_add_polyline(wks,plot2,(/t+1,t+1/),(/yp1_eke_PD(t),yp2_eke_PD(t)/),polyres)
;  error_bar002(t) = gsn_add_polyline(wks,plot2,(/t+1,t+1/),(/yp1_eke_4xCO2(t),yp2_eke_4xCO2(t)/),polyres)
; end do


 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "royalblue1"                 ; color chosen
 dummy1 = gsn_add_polygon (wks,plot5,xp1_sst_PD,yp1_sst_PD,gsres)

 gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "pink1"                 ; color chosen
 gsres@gsFillOpacityF     = 0.5
 dummy2 = gsn_add_polygon (wks,plot5,xp1_sst_4xCO2,yp1_sst_4xCO2,gsres)

 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "royalblue1"                 ; color chosen
 dummy01 = gsn_add_polygon (wks,plot4,xp1_va2_PD,yp1_va2_PD,gsres)

 gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "pink1"                 ; color chosen
 gsres@gsFillOpacityF     = 0.5
 dummy02 = gsn_add_polygon (wks,plot4,xp1_va2_4xCO2,yp1_va2_4xCO2,gsres)

 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "royalblue1"                 ; color chosen
 dummy001 = gsn_add_polygon (wks,plot3,xp1_va1_PD,yp1_va1_PD,gsres)

 gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "pink1"                 ; color chosen
 gsres@gsFillOpacityF     = 0.5
 dummy002 = gsn_add_polygon (wks,plot3,xp1_va1_4xCO2,yp1_va1_4xCO2,gsres)


 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "royalblue1"                 ; color chosen
 dummy0001 = gsn_add_polygon (wks,plot6,xp1_ssh_PD,yp1_ssh_PD,gsres)

 gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "pink1"                 ; color chosen
 gsres@gsFillOpacityF     = 0.5
 dummy0002 = gsn_add_polygon (wks,plot6,xp1_ssh_4xCO2,yp1_ssh_4xCO2,gsres)

;txres                       = True
;txres@txPerimOn             = False  
;txres@txBackgroundFillColor = "transparent"
;txres@txFontHeightF         = 0.02
;txres@txFontThicknessF      = 4
;txres@txFont                = 22
;txid1     = gsn_create_text(wks, "CESM-UHR", txres)
;
;amres                  = True
;amres@amParallelPosF   =  0.2    ; 0.5 This is the right edge of the plot.
;amres@amOrthogonalPosF = -1.2    ; 0.5 This is the bottom edge of the plot.
;amres@amJust           = "TopLeft"
;
;annoid1 = gsn_add_annotation(plot0, txid1, amres)

draw(plot0)
draw(plot1)
draw(plot2)
draw(plot3)
draw(plot4)
draw(plot5)
draw(plot6)

;==================================================================
 xyres                       = True                 ; plot mods desired
 xyres@gsnFrame              = False                ; don't advance frame yet
 xyres@gsnDraw               = False
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.3                  ; bring bars down to zero
 xyres@trXMinF               = 0                    ; adds space on either end
 xyres@trXMaxF               = 4                   ; of the 1st and last bars
 xyres@vpWidthF              = 0.095
 xyres@vpHeightF             = 0.1

 xyres@gsnXYBarChart         = True                 ; turn on bar chart
 xyres@gsnXYBarChartBarWidth = 0.24                  ; change bar widths
 xyres@tmXBOn                = True                 ; turn off tickmarks at bot
 xyres@tmYLOn                = True
 xyres@tmYROn                = False                 ; turn off tickmarks at bot
 xyres@tmXTOn                = False                 ; turn off tickmarks at bot
 xyres@tmXBMode          = "Explicit"         ; explicit labels
 xyres@tmXBValues        = (/1,2,3,4/)
 xyres@tmXBLabels = (/"PD","2CO2","4CO2","" /)
 xyres@tmXBLabelFontHeightF = 0.005
 xyres@tmYLMode          = "Explicit"         ; explicit labels
 xyres@tmYLValues        = (/0,0.1,0.2,0.3,0.4/)
 xyres@tmYLLabels        = (/"0","0.1","0.2","0.3","0.4"/)
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

 xyres@vpXF    = 0.10 + 4*xyres@vpWidthF + 2*yspace1
 yf=0.58
 xyres@vpYF    = yf
 xyres@gsnXYBarChartColors = (/"royalblue1"/)
 xyres@gsnCenterString       = "MAM"
 xyres@gsnCenterStringOrthogonalPosF = -0.13
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString      ="~F22~j"
 xyres@gsnLeftStringFontHeightF =0.012
 plot7  = gsn_csm_xy(wks,ispan(1,3,1),stdarr4(:),xyres)   ; advancing the frame
 xyres@tmYLOn                = False 
 xyres@tmYLBorderOn          = False 
 xyres@vpXF = 0.10 + 5*xyres@vpWidthF+ 2*yspace1
 xyres@gsnXYBarChartColors = (/"firebrick1"/)
 xyres@gsnCenterString       = "SON"
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.4               ; bring bars down to zero
 xyres@tmYLValues        = (/0,0.1,0.2,0.3,0.4/)
 xyres@tmYLLabels        = (/"0","0.1","0.2","0.3","0.4"/)
 xyres@gsnCenterStringOrthogonalPosF = -0.10
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString      = ""
 plot8  = gsn_csm_xy(wks,ispan(1,3,1),stdarr04(:),xyres)   ; advancing the frame
 ;===================================================================
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@tmYLOn                = True 
 xyres@tmYLBorderOn          = True 
 xyres@vpXF    = 0.10 + 2*xyres@vpWidthF +yspace1
 xyres@vpYF    = yf
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.3                ; bring bars down to zero
 xyres@gsnXYBarChartColors = (/"royalblue1"/)
 xyres@tmYLValues        = (/0,0.1,0.2,0.3,0.4/)
 xyres@tmYLLabels        = (/"0","0.1","0.2","0.3","0.4"/)
 xyres@gsnCenterString       = "MAM"
 xyres@gsnCenterStringOrthogonalPosF = -0.13
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString      ="~F22~i "
 xyres@gsnLeftStringFontHeightF =0.012
 plot9  = gsn_csm_xy(wks,ispan(1,3,1),stdarr5(:),xyres)   ; advancing the frame
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@tmYLOn                = False 
 xyres@tmYLBorderOn          = False 
 xyres@vpXF    = 0.10 + 3*xyres@vpWidthF+yspace1
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.4               ; bring bars down to zero
 xyres@gsnXYBarChartColors = (/"firebrick1"/)
 xyres@tmYLValues        = (/0,0.1,0.2,0.3,0.4/)
 xyres@tmYLLabels        = (/"0","0.1","0.2","0.3","0.4"/)
 xyres@gsnCenterString       = "SON"
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnCenterStringOrthogonalPosF = -0.10
  xyres@gsnLeftString      =""
 plot10  = gsn_csm_xy(wks,ispan(1,3,1),stdarr05(:),xyres)   ; advancing the frame
;==========================================================
;SST metric
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@tmYLOn                = True
 xyres@vpXF    =  0.10
 xyres@vpYF    = yf
 xyres@tmYLBorderOn  = True 
 xyres@gsnXYBarChartColors = (/"royalblue1"/)
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.8               ; bring bars down to zero
 xyres@tmYLValues        = (/0,0.2,0.4,0.6,0.8/)
 xyres@tmYLLabels        = (/"0","0.2","0.4","0.6","0.8"/)
 xyres@gsnCenterString       = "MAM"
 xyres@gsnCenterStringOrthogonalPosF = -0.13
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString      ="~F22~h"
 plot11  = gsn_csm_xy(wks,ispan(1,3,1),stdarr1(:),xyres)   ; advancing the frame
 xyres@vpXF = xyres@vpWidthF+xyres@vpXF
 xyres@tmYLOn        = False
 xyres@tmYLBorderOn  = False                        ; turn off right axis
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.6               ; bring bars down to zero
 xyres@gsnXYBarChartColors = (/"firebrick1"/)
 xyres@tmYLValues        = (/0,0.2,0.4,0.6,0.8/)
 xyres@tmYLLabels        = (/"0","0.2","0.4","0.6","0.8"/)
 xyres@gsnCenterString       = "SON"
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnCenterStringOrthogonalPosF = -0.10
 xyres@gsnLeftString      =""
 plot12  = gsn_csm_xy(wks,ispan(1,3,1),stdarr01(:),xyres)   ; advancing the frame
;-----------------------------------------------------------------------------------
; SSH metric
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@tmYLOn                = True
 xyres@vpXF    =  0.10 + 6*xyres@vpWidthF +3*yspace1
 xyres@vpYF    = yf
 xyres@tmYLBorderOn  = True 
 xyres@gsnXYBarChartColors = (/"royalblue1"/)
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 5              ; bring bars down to zero
 xyres@tmYLValues        = (/0,1,2,3,4,5/)
 xyres@tmYLLabels        = (/"0","1","2","3","4","5"/)
 xyres@gsnCenterString       = "MAM"
 xyres@gsnCenterStringOrthogonalPosF = -0.13
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString      ="~F22~k"
 plot13  = gsn_csm_xy(wks,ispan(1,3,1),stdarr2(:),xyres)   ; advancing the frame

 xyres@vpXF = 0.10 + 7*xyres@vpWidthF +3*yspace1
 xyres@tmYLOn        = False
 xyres@tmYLBorderOn  = False                        ; turn off right axis
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 5              ; bring bars down to zero
 xyres@gsnXYBarChartColors = (/"firebrick1"/)
 xyres@tmYLValues        = (/0,1,2,3,4,5/)
 xyres@tmYLLabels        = (/"0","1","2","3","4","5"/)
 xyres@gsnCenterString       = "SON"
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnCenterStringOrthogonalPosF = -0.10
 xyres@gsnLeftString      =""
 plot14  = gsn_csm_xy(wks,ispan(1,3,1),stdarr02(:),xyres)   ; advancing the frame

  polyres=True
  polyres@gsLineDashPattern=0
  polyres@gsLineColor="black"
  polyres@gsLineThicknessF=1.0
  polyres@tfPolyDrawOrder = "postDraw"
; SON
  amm1    = gsn_add_polyline(wks,plot12,(/1,1/),(/xBootLow_sst1_PD,xBootHi_sst1_PD/), polyres)
  amm01   = gsn_add_polyline(wks,plot12,(/2,2/),(/xBootLow_sst1_2xCO2,xBootHi_sst1_2xCO2/), polyres)
  amm001  = gsn_add_polyline(wks,plot12,(/3,3/),(/xBootLow_sst1_4xCO2,xBootHi_sst1_4xCO2/), polyres)

  bmm1    = gsn_add_polyline(wks,plot14,(/1,1/),(/xBootLow_ssh1_PD,xBootHi_ssh1_PD/), polyres)
  bmm01   = gsn_add_polyline(wks,plot14,(/2,2/),(/xBootLow_ssh1_2xCO2,xBootHi_ssh1_2xCO2/), polyres)
  bmm001  = gsn_add_polyline(wks,plot14,(/3,3/),(/xBootLow_ssh1_4xCO2,xBootHi_ssh1_4xCO2/), polyres)

 ; cumm1    = gsn_add_polyline(wks,plot11,(/1,1/),(/xBootLow_eke1_PD,xBootHi_eke1_PD/), polyres)
 ; cumm01   = gsn_add_polyline(wks,plot11,(/2,2/),(/xBootLow_eke1_2xCO2,xBootHi_eke1_2xCO2/), polyres)
 ; cumm001  = gsn_add_polyline(wks,plot11,(/3,3/),(/xBootLow_eke1_4xCO2,xBootHi_eke1_4xCO2/), polyres)

  dmm1    = gsn_add_polyline(wks,plot8,(/1,1/),(/xBootLow_va1_PD,xBootHi_va1_PD/), polyres)  ; Rossby mode
  dmm01   = gsn_add_polyline(wks,plot8,(/2,2/),(/xBootLow_va1_2xCO2,xBootHi_va1_2xCO2/), polyres)
  dmm001  = gsn_add_polyline(wks,plot8,(/3,3/),(/xBootLow_va1_4xCO2,xBootHi_va1_4xCO2/), polyres)

  emm1    = gsn_add_polyline(wks,plot10,(/1,1/),(/xBootLow_va2_PD,xBootHi_va2_PD/), polyres)  ; Yanai mode
  emm01   = gsn_add_polyline(wks,plot10,(/2,2/),(/xBootLow_va2_2xCO2,xBootHi_va2_2xCO2/), polyres)
  emm001  = gsn_add_polyline(wks,plot10,(/3,3/),(/xBootLow_va2_4xCO2,xBootHi_va2_4xCO2/), polyres)
;====================================================================================================
; MAM
  amm2    = gsn_add_polyline(wks,plot11,(/1,1/),(/xBootLow_sst2_PD,xBootHi_sst2_PD/), polyres)
  amm02   = gsn_add_polyline(wks,plot11,(/2,2/),(/xBootLow_sst2_2xCO2,xBootHi_sst2_2xCO2/), polyres)
  amm002  = gsn_add_polyline(wks,plot11,(/3,3/),(/xBootLow_sst2_4xCO2,xBootHi_sst2_4xCO2/), polyres)

  bmm2    = gsn_add_polyline(wks,plot13,(/1,1/),(/xBootLow_ssh2_PD,xBootHi_ssh2_PD/), polyres)
  bmm02   = gsn_add_polyline(wks,plot13,(/2,2/),(/xBootLow_ssh2_2xCO2,xBootHi_ssh2_2xCO2/), polyres)
  bmm002  = gsn_add_polyline(wks,plot13,(/3,3/),(/xBootLow_ssh2_4xCO2,xBootHi_ssh2_4xCO2/), polyres)

 ;cumm2    = gsn_add_polyline(wks,plot10,(/1,1/),(/xBootLow_eke2_PD,xBootHi_eke2_PD/), polyres)
 ;cumm02   = gsn_add_polyline(wks,plot10,(/2,2/),(/xBootLow_eke2_2xCO2,xBootHi_eke2_2xCO2/), polyres)
 ;cumm002  = gsn_add_polyline(wks,plot10,(/3,3/),(/xBootLow_eke2_4xCO2,xBootHi_eke2_4xCO2/), polyres)
  
  dmm2    = gsn_add_polyline(wks,plot7,(/1,1/),(/xBootLow_va01_PD,xBootHi_va01_PD/), polyres)
  dmm02   = gsn_add_polyline(wks,plot7,(/2,2/),(/xBootLow_va01_2xCO2,xBootHi_va01_2xCO2/), polyres)
  dmm002  = gsn_add_polyline(wks,plot7,(/3,3/),(/xBootLow_va01_4xCO2,xBootHi_va01_4xCO2/), polyres)

  emm2    = gsn_add_polyline(wks,plot9,(/1,1/),(/xBootLow_va02_PD,xBootHi_va02_PD/), polyres)
  emm02   = gsn_add_polyline(wks,plot9,(/2,2/),(/xBootLow_va02_2xCO2,xBootHi_va02_2xCO2/), polyres)
  emm002  = gsn_add_polyline(wks,plot9,(/3,3/),(/xBootLow_va02_4xCO2,xBootHi_va02_4xCO2/), polyres)
;====================================================================================================

  wid = 0.24/2
  dumm1   = gsn_add_polyline(wks,plot12,(/1-wid,1+wid/),(/xBootLow_sst1_PD,xBootLow_sst1_PD/), polyres)
  dumm01  = gsn_add_polyline(wks,plot12,(/1-wid,1+wid/),(/xBootHi_sst1_PD,xBootHi_sst1_PD/), polyres)

  eumm1   = gsn_add_polyline(wks,plot12,(/2-wid,2+wid/),(/xBootLow_sst1_2xCO2,xBootLow_sst1_2xCO2/), polyres)
  eumm01  = gsn_add_polyline(wks,plot12,(/2-wid,2+wid/),(/xBootHi_sst1_2xCO2,xBootHi_sst1_2xCO2/), polyres)

  fumm1   = gsn_add_polyline(wks,plot12,(/3-wid,3+wid/),(/xBootLow_sst1_4xCO2,xBootLow_sst1_4xCO2/), polyres)
  fumm01  = gsn_add_polyline(wks,plot12,(/3-wid,3+wid/),(/xBootHi_sst1_4xCO2,xBootHi_sst1_4xCO2/), polyres)
;----------------------------------------------------------------------------------------------------
  dumm2   = gsn_add_polyline(wks,plot14,(/1-wid,1+wid/),(/xBootLow_ssh1_PD,xBootLow_ssh1_PD/), polyres)
  dumm02  = gsn_add_polyline(wks,plot14,(/1-wid,1+wid/),(/xBootHi_ssh1_PD,xBootHi_ssh1_PD/), polyres)

  eumm2   = gsn_add_polyline(wks,plot14,(/2-wid,2+wid/),(/xBootLow_ssh1_2xCO2,xBootLow_ssh1_2xCO2/), polyres)
  eumm02  = gsn_add_polyline(wks,plot14,(/2-wid,2+wid/),(/xBootHi_ssh1_2xCO2,xBootHi_ssh1_2xCO2/), polyres)

  fumm2   = gsn_add_polyline(wks,plot14,(/3-wid,3+wid/),(/xBootLow_ssh1_4xCO2,xBootLow_ssh1_4xCO2/), polyres)
  fumm02  = gsn_add_polyline(wks,plot14,(/3-wid,3+wid/),(/xBootHi_ssh1_4xCO2,xBootHi_ssh1_4xCO2/), polyres)
;-----------------------------------------------------------------------------------------------------
; dumm3   = gsn_add_polyline(wks,plot11,(/1-wid,1+wid/),(/xBootLow_eke1_PD,xBootLow_eke1_PD/), polyres)
; dumm03  = gsn_add_polyline(wks,plot11,(/1-wid,1+wid/),(/xBootHi_eke1_PD,xBootHi_eke1_PD/), polyres)
;
; eumm3   = gsn_add_polyline(wks,plot11,(/2-wid,2+wid/),(/xBootLow_eke1_2xCO2,xBootLow_eke1_2xCO2/), polyres)
; eumm03  = gsn_add_polyline(wks,plot11,(/2-wid,2+wid/),(/xBootHi_eke1_2xCO2,xBootHi_eke1_2xCO2/), polyres)
;
; fumm3   = gsn_add_polyline(wks,plot11,(/3-wid,3+wid/),(/xBootLow_eke1_4xCO2,xBootLow_eke1_4xCO2/), polyres)
; fumm03  = gsn_add_polyline(wks,plot11,(/3-wid,3+wid/),(/xBootHi_eke1_4xCO2,xBootHi_eke1_4xCO2/), polyres)
;-----------------------------------------------------------------------------------------------------
  dumm4   = gsn_add_polyline(wks,plot8,(/1-wid,1+wid/),(/xBootLow_va1_PD,xBootLow_va1_PD/), polyres)
  dumm04  = gsn_add_polyline(wks,plot8,(/1-wid,1+wid/),(/xBootHi_va1_PD,xBootHi_va1_PD/), polyres)

  eumm4   = gsn_add_polyline(wks,plot8,(/2-wid,2+wid/),(/xBootLow_va1_2xCO2,xBootLow_va1_2xCO2/), polyres)
  eumm04  = gsn_add_polyline(wks,plot8,(/2-wid,2+wid/),(/xBootHi_va1_2xCO2,xBootHi_va1_2xCO2/), polyres)

  fumm4   = gsn_add_polyline(wks,plot8,(/3-wid,3+wid/),(/xBootLow_va1_4xCO2,xBootLow_va1_4xCO2/), polyres)
  fumm04  = gsn_add_polyline(wks,plot8,(/3-wid,3+wid/),(/xBootHi_va1_4xCO2,xBootHi_va1_4xCO2/), polyres)
  ;-----------------------------------------------------------------------------------------------------
  dumm5   = gsn_add_polyline(wks,plot10,(/1-wid,1+wid/),(/xBootLow_va2_PD,xBootLow_va2_PD/), polyres)
  dumm05  = gsn_add_polyline(wks,plot10,(/1-wid,1+wid/),(/xBootHi_va2_PD,xBootHi_va2_PD/), polyres)

  eumm5   = gsn_add_polyline(wks,plot10,(/2-wid,2+wid/),(/xBootLow_va2_2xCO2,xBootLow_va2_2xCO2/), polyres)
  eumm05  = gsn_add_polyline(wks,plot10,(/2-wid,2+wid/),(/xBootHi_va2_2xCO2,xBootHi_va2_2xCO2/), polyres)

  fumm5   = gsn_add_polyline(wks,plot10,(/3-wid,3+wid/),(/xBootLow_va2_4xCO2,xBootLow_va2_4xCO2/), polyres)
  fumm05  = gsn_add_polyline(wks,plot10,(/3-wid,3+wid/),(/xBootHi_va2_4xCO2,xBootHi_va2_4xCO2/), polyres)
;------------------------------------------------------------------------------------------------------
  ddumm1   = gsn_add_polyline(wks,plot11,(/1-wid,1+wid/),(/xBootLow_sst2_PD,xBootLow_sst2_PD/), polyres)
  ddumm01  = gsn_add_polyline(wks,plot11,(/1-wid,1+wid/),(/xBootHi_sst2_PD,xBootHi_sst2_PD/), polyres)
  
  eeumm1   = gsn_add_polyline(wks,plot11,(/2-wid,2+wid/),(/xBootLow_sst2_2xCO2,xBootLow_sst2_2xCO2/), polyres)
  eeumm01  = gsn_add_polyline(wks,plot11,(/2-wid,2+wid/),(/xBootHi_sst2_2xCO2,xBootHi_sst2_2xCO2/), polyres)
  
  ffumm1   = gsn_add_polyline(wks,plot11,(/3-wid,3+wid/),(/xBootLow_sst2_4xCO2,xBootLow_sst2_4xCO2/), polyres)
  ffumm01  = gsn_add_polyline(wks,plot11,(/3-wid,3+wid/),(/xBootHi_sst2_4xCO2,xBootHi_sst2_4xCO2/), polyres)
;----------------------------------------------------------------------------------------------------
  ddumm2   = gsn_add_polyline(wks,plot13,(/1-wid,1+wid/),(/xBootLow_ssh2_PD,xBootLow_ssh2_PD/), polyres)
  ddumm02  = gsn_add_polyline(wks,plot13,(/1-wid,1+wid/),(/xBootHi_ssh2_PD,xBootHi_ssh2_PD/), polyres)
  
  eeumm2   = gsn_add_polyline(wks,plot13,(/2-wid,2+wid/),(/xBootLow_ssh2_2xCO2,xBootLow_ssh2_2xCO2/), polyres)
  eeumm02  = gsn_add_polyline(wks,plot13,(/2-wid,2+wid/),(/xBootHi_ssh2_2xCO2,xBootHi_ssh2_2xCO2/), polyres)
  
  ffumm2   = gsn_add_polyline(wks,plot13,(/3-wid,3+wid/),(/xBootLow_ssh2_4xCO2,xBootLow_ssh2_4xCO2/), polyres)
  ffumm02  = gsn_add_polyline(wks,plot13,(/3-wid,3+wid/),(/xBootHi_ssh2_4xCO2,xBootHi_ssh2_4xCO2/), polyres)
;-----------------------------------------------------------------------------------------------------
; ddumm3   = gsn_add_polyline(wks,plot10,(/1-wid,1+wid/),(/xBootLow_eke2_PD,xBootLow_eke2_PD/), polyres)
; ddumm03  = gsn_add_polyline(wks,plot10,(/1-wid,1+wid/),(/xBootHi_eke2_PD,xBootHi_eke2_PD/), polyres)
; 
; eeumm3   = gsn_add_polyline(wks,plot10,(/2-wid,2+wid/),(/xBootLow_eke2_2xCO2,xBootLow_eke2_2xCO2/), polyres)
; eeumm03  = gsn_add_polyline(wks,plot10,(/2-wid,2+wid/),(/xBootHi_eke2_2xCO2,xBootHi_eke2_2xCO2/), polyres)
; 
; ffumm3   = gsn_add_polyline(wks,plot10,(/3-wid,3+wid/),(/xBootLow_eke2_4xCO2,xBootLow_eke2_4xCO2/), polyres)
; ffumm03  = gsn_add_polyline(wks,plot10,(/3-wid,3+wid/),(/xBootHi_eke2_4xCO2,xBootHi_eke2_4xCO2/), polyres)
  ;-----------------------------------------------------------------------------------------------------
  ddumm4   = gsn_add_polyline(wks,plot7,(/1-wid,1+wid/),(/xBootLow_va01_PD,xBootLow_va01_PD/), polyres)
  ddumm04  = gsn_add_polyline(wks,plot7,(/1-wid,1+wid/),(/xBootHi_va01_PD,xBootHi_va01_PD/), polyres)
  
  eeumm4   = gsn_add_polyline(wks,plot7,(/2-wid,2+wid/),(/xBootLow_va01_2xCO2,xBootLow_va01_2xCO2/), polyres)
  eeumm04  = gsn_add_polyline(wks,plot7,(/2-wid,2+wid/),(/xBootHi_va01_2xCO2,xBootHi_va01_2xCO2/), polyres)
  
  ffumm4   = gsn_add_polyline(wks,plot7,(/3-wid,3+wid/),(/xBootLow_va01_4xCO2,xBootLow_va01_4xCO2/), polyres)
  ffumm04  = gsn_add_polyline(wks,plot7,(/3-wid,3+wid/),(/xBootHi_va01_4xCO2,xBootHi_va01_4xCO2/), polyres)
;-----------------------------------------------------------------------------------------------------
ddumm5   = gsn_add_polyline(wks,plot9,(/1-wid,1+wid/),(/xBootLow_va02_PD,xBootLow_va02_PD/), polyres)
ddumm05  = gsn_add_polyline(wks,plot9,(/1-wid,1+wid/),(/xBootHi_va02_PD,xBootHi_va02_PD/), polyres)

eeumm5   = gsn_add_polyline(wks,plot9,(/2-wid,2+wid/),(/xBootLow_va02_2xCO2,xBootLow_va02_2xCO2/), polyres)
eeumm05  = gsn_add_polyline(wks,plot9,(/2-wid,2+wid/),(/xBootHi_va02_2xCO2,xBootHi_va02_2xCO2/), polyres)

ffumm5   = gsn_add_polyline(wks,plot9,(/3-wid,3+wid/),(/xBootLow_va02_4xCO2,xBootLow_va02_4xCO2/), polyres)
ffumm05  = gsn_add_polyline(wks,plot9,(/3-wid,3+wid/),(/xBootHi_va02_4xCO2,xBootHi_va02_4xCO2/), polyres)
;==================================================================================================================
draw(plot7)
draw(plot8)
draw(plot9)
draw(plot10)
draw(plot11)
draw(plot12)
draw(plot13)
draw(plot14)

end 
