begin
; parameter constant 
PI  = 3.14159265
R   = 6378388.
rad = PI/180.
dens= 1025 ; reference density of seawater
sec = 3600*24*30.42 ;one month
;=============================================================
; read SST data
;=============================================================
yrStrt  =  (/111,141,141/)
yrLast  =  (/140,170,170/)
ilev1   = 0
ilev2   = 5
lonw    = 200
lone    = 260
ncdir1  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/monthly/"
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

dTdy1_xy_nino3  = runave(wgt_areaave_Wrap(dTdy1_xy(:,{0:5},{200:270}), 1.0, 1.0, 0),3,0)
dTdy2_xy_nino3  = runave(wgt_areaave_Wrap(dTdy2_xy(:,{0:5},{200:270}), 1.0, 1.0, 0),3,0)
dTdy3_xy_nino3  = runave(wgt_areaave_Wrap(dTdy3_xy(:,{0:5},{200:270}), 1.0, 1.0, 0),3,0)


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

dTdy1_index_ave = dim_avg_n_Wrap(dTdy1_index_ann, 0) - avg(dim_avg_n_Wrap(dTdy1_index_ann, 0))
dTdy2_index_ave = dim_avg_n_Wrap(dTdy2_index_ann, 0) - avg(dim_avg_n_Wrap(dTdy2_index_ann, 0))
dTdy3_index_ave = dim_avg_n_Wrap(dTdy3_index_ann, 0) - avg(dim_avg_n_Wrap(dTdy3_index_ann, 0))

; std 
temp1_nino3_std = dim_stddev_n_Wrap(temp1_nino3_ann, 0) 
temp2_nino3_std = dim_stddev_n_Wrap(temp2_nino3_ann, 0) 
temp3_nino3_std = dim_stddev_n_Wrap(temp3_nino3_ann, 0) 

dTdy1_index_std = dim_stddev_n_Wrap(dTdy1_index_ann, 0) 
dTdy2_index_std = dim_stddev_n_Wrap(dTdy2_index_ann, 0) 
dTdy3_index_std = dim_stddev_n_Wrap(dTdy3_index_ann, 0) 
;=============================================================
; read UVEL data
;=============================================================
yrStrt  =  (/111,141,141/)
yrLast  =  (/140,170,170/)
ilev1   = 0
ilev2   = 5
lonw    =210
lone    =250
varname="UVEL"
ncdir2  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/monthly/"
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

dUdy1_xy_index  = runave(wgt_areaave_Wrap(dUdy1_xy(:,{3:5},{200:250}), 1.0, 1.0, 0),3,0)
dUdy2_xy_index  = runave(wgt_areaave_Wrap(dUdy2_xy(:,{3:5},{200:250}), 1.0, 1.0, 0),3,0)
dUdy3_xy_index  = runave(wgt_areaave_Wrap(dUdy3_xy(:,{3:5},{200:250}), 1.0, 1.0, 0),3,0)

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
dUdy1_index_ave = dim_avg_n_Wrap(dUdy1_index_ann, 0) - avg(dim_avg_n_Wrap(dUdy1_index_ann, 0))
dUdy2_index_ave = dim_avg_n_Wrap(dUdy2_index_ann, 0) - avg(dim_avg_n_Wrap(dUdy2_index_ann, 0))
dUdy3_index_ave = dim_avg_n_Wrap(dUdy3_index_ann, 0) - avg(dim_avg_n_Wrap(dUdy3_index_ann, 0))

dUdy01_index_ave = dim_avg_n_Wrap(dUdy01_index_ann, 0) - avg(dim_avg_n_Wrap(dUdy01_index_ann, 0))
dUdy02_index_ave = dim_avg_n_Wrap(dUdy02_index_ann, 0) - avg(dim_avg_n_Wrap(dUdy02_index_ann, 0))
dUdy03_index_ave = dim_avg_n_Wrap(dUdy03_index_ann, 0) - avg(dim_avg_n_Wrap(dUdy03_index_ann, 0))

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

printMinMax(data1, 0)
printMinMax(data2, 0)
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
;=============================================================================

ncdir01  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/mean/flux/mon/"
fil01    = systemfunc("ls " + ncdir01 + "CESM2_HR_monthly_heat_flux.10-50.1degree.nc")
fil02    = systemfunc("ls " + ncdir01 + "CESM2_HR_monthly_momentum_flux.10-50.1degree.nc")
  f01    = addfile(fil01, "r")
  f02    = addfile(fil02, "r")
 
   hf0_PD     = f01->hf_mon_PD(:,0:1,:,:)
   hf0_2xCO2  = f01->hf_mon_2xCO2(:,0:1,:,:)
   hf0_4xCO2  = f01->hf_mon_4xCO2(:,0:1,:,:)
   
   mf0_PD     = f02->mf_mon_PD(:,0:1,:,:)
   mf0_2xCO2  = f02->mf_mon_2xCO2(:,0:1,:,:)
   mf0_4xCO2  = f02->mf_mon_4xCO2(:,0:1,:,:)
            
   hf_PD    = dim_avg_n_Wrap(hf0_PD,1)
   hf_2xCO2 = dim_avg_n_Wrap(hf0_2xCO2, 1)
   hf_4xCO2 = dim_avg_n_Wrap(hf0_4xCO2, 1)
  
   mf_PD    = dim_avg_n_Wrap(mf0_PD,1)
   mf_2xCO2 = dim_avg_n_Wrap(mf0_2xCO2, 1)
   mf_4xCO2 = dim_avg_n_Wrap(mf0_4xCO2, 1)


hf1_PD    = runave(wgt_areaave_Wrap(hf_PD(:,{2:5},{210:250}), 1.0, 1.0, 0),3,0)
hf1_2xCO2 = runave(wgt_areaave_Wrap(hf_2xCO2(:,{2:5},{210:250}), 1.0, 1.0, 0),3,0)
hf1_4xCO2 = runave(wgt_areaave_Wrap(hf_4xCO2(:,{2:5},{210:250}), 1.0, 1.0, 0),3,0)

mf2_PD    = runave(wgt_areaave_Wrap(mf_PD(:,{3:5},{200:250}), 1.0, 1.0, 0),3,0)
mf2_2xCO2 = runave(wgt_areaave_Wrap(mf_2xCO2(:,{3:5},{200:250}), 1.0, 1.0, 0),3,0)
mf2_4xCO2 = runave(wgt_areaave_Wrap(mf_4xCO2(:,{3:5},{200:250}), 1.0, 1.0, 0),3,0)

mf1_PD    = runave(wgt_areaave_Wrap(mf_PD(:,{-1:1},{200:250}), 1.0, 1.0, 0),3,0)
mf1_2xCO2 = runave(wgt_areaave_Wrap(mf_2xCO2(:,{-1:1},{200:250}), 1.0, 1.0, 0),3,0)
mf1_4xCO2 = runave(wgt_areaave_Wrap(mf_4xCO2(:,{-1:1},{200:240}), 1.0, 1.0, 0),3,0)


hf1_PD_ann     = reshape(hf1_PD,(/nyears,12/))
hf1_2xCO2_ann  = reshape(hf1_2xCO2,(/nyears,12/))
hf1_4xCO2_ann  = reshape(hf1_4xCO2,(/nyears,12/))

mf1_PD_ann     = reshape(mf1_PD,(/nyears,12/))
mf1_2xCO2_ann  = reshape(mf1_2xCO2,(/nyears,12/))
mf1_4xCO2_ann  = reshape(mf1_4xCO2,(/nyears,12/))

mf2_PD_ann     = reshape(mf2_PD,(/nyears,12/))
mf2_2xCO2_ann  = reshape(mf2_2xCO2,(/nyears,12/))
mf2_4xCO2_ann  = reshape(mf2_4xCO2,(/nyears,12/))

; annual cycle
hf1_PD_ave    = dim_avg_n_Wrap(hf1_PD_ann, 0)    - avg(dim_avg_n_Wrap(hf1_PD_ann, 0))
hf1_2xCO2_ave = dim_avg_n_Wrap(hf1_2xCO2_ann, 0) - avg(dim_avg_n_Wrap(hf1_2xCO2_ann, 0))
hf1_4xCO2_ave = dim_avg_n_Wrap(hf1_4xCO2_ann, 0) - avg(dim_avg_n_Wrap(hf1_4xCO2_ann, 0))

mf1_PD_ave    = dim_avg_n_Wrap(mf1_PD_ann, 0)    - avg(dim_avg_n_Wrap(mf1_PD_ann, 0))
mf1_2xCO2_ave = dim_avg_n_Wrap(mf1_2xCO2_ann, 0) - avg(dim_avg_n_Wrap(mf1_2xCO2_ann, 0))
mf1_4xCO2_ave = dim_avg_n_Wrap(mf1_4xCO2_ann, 0) - avg(dim_avg_n_Wrap(mf1_4xCO2_ann, 0))

mf2_PD_ave    = dim_avg_n_Wrap(mf2_PD_ann, 0)    - avg(dim_avg_n_Wrap(mf2_PD_ann, 0))
mf2_2xCO2_ave = dim_avg_n_Wrap(mf2_2xCO2_ann, 0) - avg(dim_avg_n_Wrap(mf2_2xCO2_ann, 0))
mf2_4xCO2_ave = dim_avg_n_Wrap(mf2_4xCO2_ann, 0) - avg(dim_avg_n_Wrap(mf2_4xCO2_ann, 0))

; std 
hf1_PD_std    = dim_stddev_n_Wrap(hf1_PD_ann, 0)
hf1_2xCO2_std = dim_stddev_n_Wrap(hf1_2xCO2_ann, 0)
hf1_4xCO2_std = dim_stddev_n_Wrap(hf1_4xCO2_ann, 0) 

mf1_PD_std    = dim_stddev_n_Wrap(mf1_PD_ann, 0)/100.
mf1_2xCO2_std = dim_stddev_n_Wrap(mf1_2xCO2_ann, 0)/100. 
mf1_4xCO2_std = dim_stddev_n_Wrap(mf1_4xCO2_ann, 0)/100. 

mf2_PD_std    = dim_stddev_n_Wrap(mf2_PD_ann, 0)/100.
mf2_2xCO2_std = dim_stddev_n_Wrap(mf2_2xCO2_ann, 0)/100. 
mf2_4xCO2_std = dim_stddev_n_Wrap(mf2_4xCO2_ann, 0)/100. 


data01 = new((/2,12/), float)
data01(0,:) = (/hf1_PD_ave/)
data01(1,:) = (/hf1_4xCO2_ave/)

data02 = new((/2,12/), float)
data02(0,:) = (/mf1_PD_ave/)/100.
data02(1,:) = (/mf1_4xCO2_ave/)/100.

data03 = new((/2,12/), float)
data03(0,:) = (/mf2_PD_ave/)/100.
data03(1,:) = (/mf2_4xCO2_ave/)/100.
printMinMax(data01, 0)
printMinMax(data02, 0)
printMinMax(data03, 0)
;===============================================
xp1_hf_PD = new((/24/), float)
yp1_hf_PD = new((/24/), float)

do i=0,11
   xp1_hf_PD(i)          = i+1
   xp1_hf_PD(2*12-1-i)   = i+1
   yp1_hf_PD(i)          = hf1_PD_ave(i) + hf1_PD_std(i)/2.
   yp1_hf_PD(2*12-1-i)   = hf1_PD_ave(i) - hf1_PD_std(i)/2. 
 end do

 xp1_hf_4xCO2 = new((/24/), float)
 yp1_hf_4xCO2 = new((/24/), float)
 
 do i=0,11
    xp1_hf_4xCO2(i)          = i+1
    xp1_hf_4xCO2(2*12-1-i)   = i+1
    yp1_hf_4xCO2(i)          = hf1_4xCO2_ave(i) + hf1_4xCO2_std(i)/2.
    yp1_hf_4xCO2(2*12-1-i)   = hf1_4xCO2_ave(i) - hf1_4xCO2_std(i)/2. 
  end do
;--------------------------------------------------------
xp1_mf_PD = new((/24/), float)
yp1_mf_PD = new((/24/), float)

do i=0,11
   xp1_mf_PD(i)          = i+1
   xp1_mf_PD(2*12-1-i)   = i+1
   yp1_mf_PD(i)          = mf1_PD_ave(i)/100. + mf1_PD_std(i)/2.
   yp1_mf_PD(2*12-1-i)   = mf1_PD_ave(i)/100. - mf1_PD_std(i)/2. 
 end do

 xp1_mf_4xCO2 = new((/24/), float)
 yp1_mf_4xCO2 = new((/24/), float)
 
 do i=0,11
    xp1_mf_4xCO2(i)          = i+1
    xp1_mf_4xCO2(2*12-1-i)   = i+1
    yp1_mf_4xCO2(i)          = mf1_4xCO2_ave(i)/100. + mf1_4xCO2_std(i)/2.
    yp1_mf_4xCO2(2*12-1-i)   = mf1_4xCO2_ave(i)/100. - mf1_4xCO2_std(i)/2. 
  end do
;=========================================================
xp2_mf_PD = new((/24/), float)
yp2_mf_PD = new((/24/), float)

do i=0,11
   xp2_mf_PD(i)          = i+1
   xp2_mf_PD(2*12-1-i)   = i+1
   yp2_mf_PD(i)          = mf2_PD_ave(i)/100. + mf2_PD_std(i)/2.
   yp2_mf_PD(2*12-1-i)   = mf2_PD_ave(i)/100. - mf2_PD_std(i)/2. 
 end do

 xp2_mf_4xCO2 = new((/24/), float)
 yp2_mf_4xCO2 = new((/24/), float)
 
 do i=0,11
    xp2_mf_4xCO2(i)          = i+1
    xp2_mf_4xCO2(2*12-1-i)   = i+1
    yp2_mf_4xCO2(i)          = mf2_4xCO2_ave(i)/100. + mf2_4xCO2_std(i)/2.
    yp2_mf_4xCO2(2*12-1-i)   = mf2_4xCO2_ave(i)/100. - mf2_4xCO2_std(i)/2. 
  end do
;======================================================================
outdir = "/Users/xueaoyun/Dropbox/TIW_future/figures_annual/"
wks = gsn_open_wks("pdf",outdir+"plot_CESM_PD_2xCO2_4xCO2_Seasonal_cycle_background_time_series_MFC_HFC_ac")
wks_type = "eps"
yspace   = 0.05
xspace   = 0.05
bres               = True 
bres@gsnFrame      = False
bres@gsnDraw       = False
bres@trXMinF       = 0.
bres@trXMaxF       = 13.
bres@trYMinF       = -4.
bres@trYMaxF       = 3.
bres@vpXF          = 0.1                  ; Move slightly to left
bres@vpYF          = 0.9
bres@vpHeightF     = 0.1    
bres@vpWidthF      = 0.25 
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
bres@tmYLValues               = (/-4.0,-3.0,-2.0,-1.0,0,1.0,2.0,3.0,4.0/)
bres@tmYLLabels               = (/"-4","-3","-2","-1","0","1","2","3","4"/)  
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
bres@pmLegendParallelPosF   = 0.30                ; move units right
bres@pmLegendOrthogonalPosF = -0.35                ; move units down
bres@pmLegendWidthF         = 0.04                ; Change width and
bres@pmLegendHeightF        = 0.03                ; height of legend.
bres@lgPerimOn              = False               ; turn off box around
bres@lgLabelFontHeightF     = 0.008              ; label font height
bres@lgLabelFont            = 21                  ; label font height
bres@lgLabelFontThicknessF  = 3                   ; label font height
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
bres@gsnLeftString      ="~F22~a "
bres@gsnLeftStringFontHeightF =0.015
bres@gsnCenterString   = "Meridional SST gradient"
bres@gsnCenterStringFontHeightF = 0.012
x = ispan(1,12,1)
plot0 = gsn_csm_xy(wks,x,data1,bres)
delete(bres@tmYLValues)
delete(bres@tmYLLabels)
bres@pmLegendParallelPosF   = -0.30                ; move units right
bres@trYMinF       = -1.
bres@trYMaxF       = 1.
bres@vpXF          = bres@vpXF+bres@vpWidthF+xspace                  ; Move slightly to left
bres@vpYF          = bres@vpYF
bres@vpHeightF     = 0.1    
bres@vpWidthF      = 0.25 
bres@tmYLValues               = (/-1.0,-0.5,0,0.5,1/)
bres@tmYLLabels               = (/"-1","-0.5","0","0.5","1"/) 
bres@gsnLeftString      ="~F22~b "
bres@gsnCenterString    = "EQ Shear(EUC & nSEC)"
plot1 = gsn_csm_xy(wks,x,data2,bres)
delete(bres@tmYLValues)
delete(bres@tmYLLabels)
bres@pmLegendParallelPosF   = 0.30                ; move units right
bres@trYMinF       = -3.
bres@trYMaxF       = 3.
bres@vpXF          = bres@vpXF+bres@vpWidthF+xspace                   ; Move slightly to left
bres@vpYF          = bres@vpYF
bres@vpHeightF     = 0.1    
bres@vpWidthF      = 0.25 
bres@tmYLValues               = (/-3.0,-2.0,-1.0,0,1.0,2.0,3.0/)
bres@tmYLLabels               = (/"-3","-2","-1","0","1","2","3"/)  
bres@gsnLeftString      ="~F22~c "
bres@gsnCenterString    = "NEQ Shear(NECC & nSEC)"
plot2 = gsn_csm_xy(wks,x,data3,bres)
delete(bres@tmYLValues)
delete(bres@tmYLLabels)
bres@trYMinF       = -5.
bres@trYMaxF       = 5.
bres@vpXF          = 0.1                 ; Move slightly to left
bres@vpYF          = 0.75
bres@vpHeightF     = 0.1    
bres@vpWidthF      = 0.25 
bres@tmYLValues               = (/-1.0,-0.5,0,0.5,1/)*4
bres@tmYLLabels               = (/"-0.04","-0.02","0","0.02","0.04"/) 
bres@gsnLeftString      ="~F22~d "
bres@gsnCenterString    = "Heat flux"
plot3 = gsn_csm_xy(wks,x,data01,bres)
delete(bres@tmYLValues)
delete(bres@tmYLLabels)
bres@pmLegendParallelPosF   = -0.30                ; move units right
bres@trYMinF       = -1.
bres@trYMaxF       = 1.
bres@vpXF          = bres@vpXF+bres@vpWidthF+xspace                   ; Move slightly to left
bres@vpYF          = 0.75
bres@vpHeightF     = 0.1    
bres@vpWidthF      = 0.25 
bres@tmYLValues               = (/-1.0,-0.5,0,0.5,1.0/)
bres@tmYLLabels               = (/"-0.01","-0.005","0","0.005","0.01"/)  
bres@gsnLeftString      ="~F22~e "
bres@gsnCenterString    = "EQ momentum flux"
plot4 = gsn_csm_xy(wks,x,data02,bres)
delete(bres@tmYLValues)
delete(bres@tmYLLabels)
bres@pmLegendParallelPosF   = 0.30                ; move units right
bres@trYMinF       = -1.
bres@trYMaxF       = 2.
bres@vpXF          = bres@vpXF+bres@vpWidthF+xspace                   ; Move slightly to left
bres@vpYF          = 0.75
bres@vpHeightF     = 0.1    
bres@vpWidthF      = 0.25
bres@tmYLValues               = (/-3.0,-2.0,-1.0,0,1.0,2.0,3.0/)
bres@tmYLLabels               = (/"-0.03","-0.02","-0.01","0","0.01","0.02","0.03"/)  
bres@gsnLeftString      ="~F22~f "
bres@gsnCenterString    = "NEQ Momentum flux"
plot5 = gsn_csm_xy(wks,x,data03,bres)


 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "royalblue1"                 ; color chosen
 dummy1 = gsn_add_polygon (wks,plot0,xp1_dTdy_PD,yp1_dTdy_PD,gsres)

 gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "pink1"                 ; color chosen
 gsres@gsFillOpacityF     = 0.5
 dummy2 = gsn_add_polygon (wks,plot0,xp1_dTdy_4xCO2,yp1_dTdy_4xCO2,gsres)

 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "royalblue1"                 ; color chosen
 dummy01 = gsn_add_polygon (wks,plot1,xp2_dUdy_PD,yp2_dUdy_PD,gsres)

 gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "pink1"                 ; color chosen
 gsres@gsFillOpacityF     = 0.5
 dummy02 = gsn_add_polygon (wks,plot1,xp2_dUdy_4xCO2,yp2_dUdy_4xCO2,gsres)

 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "royalblue1"                 ; color chosen
 dummy001 = gsn_add_polygon (wks,plot2,xp1_dUdy_PD,yp1_dUdy_PD,gsres)

 gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "pink1"                 ; color chosen
 gsres@gsFillOpacityF     = 0.5
 dummy002 = gsn_add_polygon (wks,plot2,xp1_dUdy_4xCO2,yp1_dUdy_4xCO2,gsres)
;==================================================================================
 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "royalblue1"                 ; color chosen
 ddummy1 = gsn_add_polygon (wks,plot3,xp1_hf_PD,yp1_hf_PD,gsres)

 gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "pink1"                 ; color chosen
 gsres@gsFillOpacityF     = 0.5
 ddummy2 = gsn_add_polygon (wks,plot3,xp1_hf_4xCO2,yp1_hf_4xCO2,gsres)

 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "royalblue1"                 ; color chosen
 ddummy01 = gsn_add_polygon (wks,plot4,xp1_mf_PD,yp1_mf_PD,gsres)

 gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "pink1"                 ; color chosen
 gsres@gsFillOpacityF     = 0.5
 ddummy02 = gsn_add_polygon (wks,plot4,xp1_mf_4xCO2,yp1_mf_4xCO2,gsres)

 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "royalblue1"                 ; color chosen
 ddummy001 = gsn_add_polygon (wks,plot5,xp2_mf_PD,yp2_mf_PD,gsres)

 gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "pink1"                 ; color chosen
 gsres@gsFillOpacityF     = 0.5
 ddummy002 = gsn_add_polygon (wks,plot5,xp2_mf_4xCO2,yp2_mf_4xCO2,gsres)

txres                       = True
txres@txPerimOn             = False  
txres@txBackgroundFillColor = "transparent"
txres@txFontHeightF         = 0.02
txres@txFontThicknessF      = 4
txres@txFont                = 22
txid1     = gsn_create_text(wks, "Seasonality", txres)

amres                  = True
amres@amParallelPosF   =  -0.30    ; 0.5 This is the right edge of the plot.
amres@amOrthogonalPosF = -1.2    ; 0.5 This is the bottom edge of the plot.
amres@amJust           = "TopLeft"

annoid1 = gsn_add_annotation(plot1, txid1, amres)

draw(plot0)
draw(plot1)
draw(plot2)
draw(plot3)
draw(plot4)
draw(plot5)
;-----------------------------------------------------------------------------------------------
ncdir  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/mean/flux/mon/"
; read SST var
fil001    = systemfunc("ls " + ncdir + "CESM2_HR_monthly_momentum_flux_convergence_10-50_1degree.nc")
  f001    = addfile(fil001, "r")
  tiw_mfc_PD0     = f001->mfc_PD(:,0:1,:,:)
  tiw_mfc_2xCO20  = f001->mfc_2xCO2(:,0:1,:,:)
  tiw_mfc_4xCO20  = f001->mfc_4xCO2(:,0:1,:,:)
  tiw_mfc_PD      = dim_avg_n_Wrap(tiw_mfc_PD0, 1)
  tiw_mfc_2xCO2   = dim_avg_n_Wrap(tiw_mfc_2xCO20, 1)
  tiw_mfc_4xCO2   = dim_avg_n_Wrap(tiw_mfc_4xCO20, 1)
  printVarSummary(tiw_mfc_PD)

  tiw_mfc_PD_ac1    = clmMonTLL(tiw_mfc_PD)
  tiw_mfc_2xCO2_ac1 = clmMonTLL(tiw_mfc_2xCO2)
  tiw_mfc_4xCO2_ac1 = clmMonTLL(tiw_mfc_4xCO2)
  printVarSummary(tiw_mfc_PD_ac1)

  tiw_mfc_PD_ac      = (tiw_mfc_PD_ac1(8,:,:)   +tiw_mfc_PD_ac1(9,:,:)   +tiw_mfc_PD_ac1(10,:,:)   -tiw_mfc_PD_ac1(2,:,:)   -tiw_mfc_PD_ac1(3,:,:)   -tiw_mfc_PD_ac1(4,:,:))/6.
  tiw_mfc_2xCO2_ac   = (tiw_mfc_2xCO2_ac1(8,:,:)+tiw_mfc_2xCO2_ac1(9,:,:)+tiw_mfc_2xCO2_ac1(10,:,:)-tiw_mfc_2xCO2_ac1(2,:,:)-tiw_mfc_2xCO2_ac1(3,:,:)-tiw_mfc_2xCO2_ac1(4,:,:))/6.
  tiw_mfc_4xCO2_ac   = (tiw_mfc_4xCO2_ac1(8,:,:)+tiw_mfc_4xCO2_ac1(9,:,:)+tiw_mfc_4xCO2_ac1(10,:,:)-tiw_mfc_4xCO2_ac1(2,:,:)-tiw_mfc_4xCO2_ac1(3,:,:)-tiw_mfc_4xCO2_ac1(4,:,:))/6.
  
  copy_VarMeta(tiw_mfc_PD(0,:,:), tiw_mfc_PD_ac)
  copy_VarMeta(tiw_mfc_2xCO2(0,:,:), tiw_mfc_2xCO2_ac)
  copy_VarMeta(tiw_mfc_4xCO2(0,:,:), tiw_mfc_4xCO2_ac)

  pattern001 = smth9_Wrap(tiw_mfc_PD_ac, 0.50, 0.25, False)
  pattern002 = smth9_Wrap(tiw_mfc_2xCO2_ac, 0.50, 0.25, False)
  pattern003 = smth9_Wrap(tiw_mfc_4xCO2_ac, 0.50, 0.25, False)

  pattern_mfc_ac = new((/2, dimsizes(pattern001&lat), dimsizes(pattern001&lon)/), typeof(pattern001))
  pattern_mfc_ac(0,:,:) = pattern001
  pattern_mfc_ac(1,:,:) = pattern003
;=========================================================================================================
  fil002    = systemfunc("ls " + ncdir + "CESM2_HR_monthly_heat_flux_convergence_10-50_1degree.nc")
  f002    = addfile(fil002, "r")
  tiw_hfc_PD0     = f002->hfc_PD(:,0:5,:,:)
  tiw_hfc_2xCO20  = f002->hfc_2xCO2(:,0:5,:,:)
  tiw_hfc_4xCO20  = f002->hfc_4xCO2(:,0:5,:,:)
  tiw_hfc_PD      = dim_avg_n_Wrap(tiw_hfc_PD0, 1)
  tiw_hfc_2xCO2   = dim_avg_n_Wrap(tiw_hfc_2xCO20, 1)
  tiw_hfc_4xCO2   = dim_avg_n_Wrap(tiw_hfc_4xCO20, 1)
  printVarSummary(tiw_hfc_PD)

  tiw_hfc_PD_ac1    = clmMonTLL(tiw_hfc_PD)
  tiw_hfc_2xCO2_ac1 = clmMonTLL(tiw_hfc_2xCO2)
  tiw_hfc_4xCO2_ac1 = clmMonTLL(tiw_hfc_4xCO2)
  printVarSummary(tiw_hfc_PD_ac1)

  tiw_hfc_PD_ac      = (tiw_hfc_PD_ac1(8,:,:)   +tiw_hfc_PD_ac1(9,:,:)   +tiw_hfc_PD_ac1(10,:,:)   -tiw_hfc_PD_ac1(2,:,:)   -tiw_hfc_PD_ac1(3,:,:)   -tiw_hfc_PD_ac1(4,:,:))/6.
  tiw_hfc_2xCO2_ac   = (tiw_hfc_2xCO2_ac1(8,:,:)+tiw_hfc_2xCO2_ac1(9,:,:)+tiw_hfc_2xCO2_ac1(10,:,:)-tiw_hfc_2xCO2_ac1(2,:,:)-tiw_hfc_2xCO2_ac1(3,:,:)-tiw_hfc_2xCO2_ac1(4,:,:))/6.
  tiw_hfc_4xCO2_ac   = (tiw_hfc_4xCO2_ac1(8,:,:)+tiw_hfc_4xCO2_ac1(9,:,:)+tiw_hfc_4xCO2_ac1(10,:,:)-tiw_hfc_4xCO2_ac1(2,:,:)-tiw_hfc_4xCO2_ac1(3,:,:)-tiw_hfc_4xCO2_ac1(4,:,:))/6.
  
  copy_VarMeta(tiw_hfc_PD(0,:,:), tiw_hfc_PD_ac)
  copy_VarMeta(tiw_hfc_2xCO2(0,:,:), tiw_hfc_2xCO2_ac)
  copy_VarMeta(tiw_hfc_4xCO2(0,:,:), tiw_hfc_4xCO2_ac)

  pattern0001 = smth9_Wrap(tiw_hfc_PD_ac, 0.50, 0.25, False)
  pattern0002 = smth9_Wrap(tiw_hfc_2xCO2_ac, 0.50, 0.25, False)
  pattern0003 = smth9_Wrap(tiw_hfc_4xCO2_ac, 0.50, 0.25, False)

  pattern_hfc_ac = new((/2, dimsizes(pattern0001&lat), dimsizes(pattern0001&lon)/), typeof(pattern0001))
  pattern_hfc_ac(0,:,:) = pattern0001
  pattern_hfc_ac(1,:,:) = pattern0003

  latS   =  -5.
  latN   =   10. 
  lonL   =  180
  lonR   =  280
  leftstring00  = (/"~F22~g","~F22~h"/) 
  leftstring01  = (/"~F22~i","~F22~j"/) 
  leftstring02  = (/"~F22~k","~F22~l"/) 
  text=(/"PD","4xCO2 -PD"/)

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
  mres@vpYF                  = 0.58
  mres@cnLevelSelectionMode  = "ExplicitLevels"
  mres@cnLevels                := fspan(-0.2,0.2,21)
  
  mres@tiYAxisString= "PD"
  mres@tiYAxisOffsetXF=-0.03
  mres@tiYAxisFontHeightF=0.012
  
  mres@gsnLeftString  = leftstring00(0)
  mres@gsnLeftStringFontHeightF = 0.015
  mres@gsnCenterString = "PD MFC_AC"
  mres@gsnCenterStringFontHeightF = 0.012
  mres@tiMainString =""
  mres@lbLabelBarOn    = True
  mres@tmXBLabels = (/"","","","","",""/)
  mplot0  = gsn_csm_contour_map(wks, pattern_mfc_ac(0,:,:),mres)
  
  mres@vpHeightF            = 0.1
  mres@vpWidthF             = 0.25
  mres@vpXF                  = 0.10
  mres@vpYF                  = 0.43
  
  mres@tiYAxisString= "4xCO2"
  mres@tiYAxisOffsetXF=-0.03
  mres@tiYAxisFontHeightF=0.012
  
  mres@gsnLeftString  = leftstring00(1)
  mres@gsnLeftStringFontHeightF = 0.015
  mres@gsnCenterString = "4xCO2 MFC_AC"
  mres@gsnCenterStringFontHeightF = 0.012
  mres@tiMainString =""
  mres@lbLabelBarOn    = True 
  mres@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W","80~S~o~N~W"/)
  mplot1  = gsn_csm_contour_map(wks, pattern_mfc_ac(1,:,:),mres)
  delete([/mres@tiYAxisString,mres@tiYAxisOffsetXF,mres@tiYAxisFontHeightF/])
  ;---------------------------------------------------------------------------
  ; ----->2 plot hfc
  delete(mres@cnLevels)
  mres@cnLevels                = fspan(-0.5,0.5,21)
  mres@tmYLLabels = (/"","","","","","",""/)
  
  mres@vpHeightF            = 0.1
  mres@vpWidthF             = 0.25
  mres@vpXF                  = 0.40
  mres@vpYF                  = 0.58
  mres@gsnLeftString  = leftstring01(0)
  mres@gsnLeftStringFontHeightF = 0.015
  mres@gsnCenterString = "PD HFC_AC"
  mres@gsnCenterStringFontHeightF = 0.012
  
  mres@lbLabelBarOn    = True 
  mres@tiYAxisString   = ""
  mres@tmXBLabels = (/"","","","","",""/)
  mplot2  = gsn_csm_contour_map(wks,pattern_hfc_ac(0,:,:),mres)
  
  mres@vpHeightF            = 0.1
  mres@vpWidthF             = 0.25
  mres@vpXF                  = 0.40
  mres@vpYF                  = 0.43
  mres@gsnLeftString  = leftstring01(1)
  mres@gsnLeftStringFontHeightF = 0.015
  mres@gsnCenterString = "4xCO2 HFC_AC"
  mres@gsnCenterStringFontHeightF = 0.012
  
  mres@lbLabelBarOn    = True
  mres@tiYAxisString   = ""
  mres@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W","80~S~o~N~W"/)
  mplot3  = gsn_csm_contour_map(wks,pattern_hfc_ac(1,:,:),mres)
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

ddum1=gsn_add_polyline(wks, mplot0, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum2=gsn_add_polyline(wks, mplot1, (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)


w2=360-160
e2=360-100

w02=360-160
e02=360-120

w002=360-150
e002=360-100

s2= -1
n2= 0.5
plres=True
plres@gsLineThicknessF=2
plres@gsLineDashPattern= 2
plres@gsLineColor="red"
plres@gsLineOpacityF=100

eum1=gsn_add_polyline(wks, mplot0, (/w002,e002,e002,w002,w002/), (/s2,s2,n2,n2,s2/), plres)
eum2=gsn_add_polyline(wks, mplot1, (/w02,e02,e02,w02,w02/), (/s2,s2,n2,n2,s2/), plres)

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

plres=True
plres@gsLineThicknessF=2
plres@gsLineDashPattern= 2
plres@gsLineColor="blue"
plres@gsLineOpacityF=100

ddum3=gsn_add_polyline(wks, mplot2, (/w03,e03,e03,w03,w03/), (/s3,s3,n3,n3,s3/), plres)
ddum4=gsn_add_polyline(wks, mplot3, (/w03,e03,e03,w03,w03/), (/s3,s3,n3,n3,s3/), plres)

plres=True
plres@gsLineThicknessF=2
plres@gsLineDashPattern= 2
plres@gsLineColor="red"
plres@gsLineOpacityF=100

eum03=gsn_add_polyline(wks, mplot2, (/w4,e4,e4,w4,w4/), (/s4,s4,n4,n4,s4/), plres)
eum04=gsn_add_polyline(wks, mplot3, (/w4,e4,e4,w4,w4/), (/s4,s4,n4,n4,s4/), plres)

;=====================================================================================================
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

   hfc_PD_ac0      = (hfc_PD(8::12,:,:)   + hfc_PD(9::12,:,:)   + hfc_PD(10::12,:,:)   -hfc_PD(2::12,:,:)   -hfc_PD(3::12,:,:)   -hfc_PD(4::12,:,:))/6.
   hfc_2xCO2_ac0   = (hfc_2xCO2(8::12,:,:)+ hfc_2xCO2(9::12,:,:)+ hfc_2xCO2(10::12,:,:)-hfc_2xCO2(2::12,:,:)-hfc_2xCO2(3::12,:,:)-hfc_2xCO2(4::12,:,:))/6.
   hfc_4xCO2_ac0   = (hfc_4xCO2(8::12,:,:)+ hfc_4xCO2(9::12,:,:)+ hfc_4xCO2(10::12,:,:)-hfc_4xCO2(2::12,:,:)-hfc_4xCO2(3::12,:,:)-hfc_4xCO2(4::12,:,:))/6.
   copy_VarMeta(hfc_PD(0::12,:,:),    hfc_PD_ac0)
   copy_VarMeta(hfc_2xCO2(0::12,:,:), hfc_2xCO2_ac0)
   copy_VarMeta(hfc_4xCO2(0::12,:,:), hfc_4xCO2_ac0)

   mfc_PD_ac0      = (mfc_PD(8::12,:,:)   + mfc_PD(9::12,:,:)   + mfc_PD(10::12,:,:)   -mfc_PD(2::12,:,:)   -mfc_PD(3::12,:,:)   -mfc_PD(4::12,:,:))/6.
   mfc_2xCO2_ac0   = (mfc_2xCO2(8::12,:,:)+ mfc_2xCO2(9::12,:,:)+ mfc_2xCO2(10::12,:,:)-mfc_2xCO2(2::12,:,:)-mfc_2xCO2(3::12,:,:)-mfc_2xCO2(4::12,:,:))/6.
   mfc_4xCO2_ac0   = (mfc_4xCO2(8::12,:,:)+ mfc_4xCO2(9::12,:,:)+ mfc_4xCO2(10::12,:,:)-mfc_4xCO2(2::12,:,:)-mfc_4xCO2(3::12,:,:)-mfc_4xCO2(4::12,:,:))/6.
   copy_VarMeta(mfc_PD(0::12,:,:),    mfc_PD_ac0)
   copy_VarMeta(mfc_2xCO2(0::12,:,:), mfc_2xCO2_ac0)
   copy_VarMeta(mfc_4xCO2(0::12,:,:), mfc_4xCO2_ac0)

   hfc1_PD_ac    = wgt_areaave_Wrap(hfc_PD_ac0(:,{4:6},{210:250}), 1.0, 1.0, 0)
hfc1_2xCO2_ac = wgt_areaave_Wrap(hfc_2xCO2_ac0(:,{4:6},{210:250}), 1.0, 1.0, 0)
hfc1_4xCO2_ac = wgt_areaave_Wrap(hfc_4xCO2_ac0(:,{4:6},{210:250}), 1.0, 1.0, 0)

hfc2_PD_ac    = wgt_areaave_Wrap(hfc_PD_ac0(:,{-1:1},{200:250}), 1.0, 1.0, 0)
hfc2_2xCO2_ac = wgt_areaave_Wrap(hfc_2xCO2_ac0(:,{-1:1},{200:250}), 1.0, 1.0, 0)
hfc2_4xCO2_ac = wgt_areaave_Wrap(hfc_4xCO2_ac0(:,{-1:1},{200:250}), 1.0, 1.0, 0)

mfc1_PD_ac    = wgt_areaave_Wrap(mfc_PD_ac0(:,{1:3},{200:250}), 1.0, 1.0, 0)
mfc1_2xCO2_ac = wgt_areaave_Wrap(mfc_2xCO2_ac0(:,{1:3},{200:250}), 1.0, 1.0, 0)
mfc1_4xCO2_ac = wgt_areaave_Wrap(mfc_4xCO2_ac0(:,{1:3},{200:250}), 1.0, 1.0, 0)

mfc2_PD_ac    = wgt_areaave_Wrap(mfc_PD_ac0(:,{-1:0.5},{200:240}), 1.0, 1.0, 0)
mfc2_2xCO2_ac = wgt_areaave_Wrap(mfc_2xCO2_ac0(:,{-1:0.5},{200:240}), 1.0, 1.0, 0)
mfc2_4xCO2_ac = wgt_areaave_Wrap(mfc_4xCO2_ac0(:,{-1:0.5},{200:240}), 1.0, 1.0, 0)


stat  = 0                         ; 0=mean; 1=variance; 2=std. dev; ...
n     = 360                        ; sub-sample size;  default is n=N
nBoot = 10000                      ; bootstrap replications
N     = 360
opt   = False                     ; use defaults
if (isvar("n") .and. (n.lt.N)) then
    opt   = True
    opt@sample_size = n           ; default N
end if

; seasonal cycle 
BootStrap_hfc1_PD_ac     = bootstrap_stat(hfc1_PD_ac, stat, nBoot, 0, opt)
BootStrap_hfc1_2xCO2_ac  = bootstrap_stat(hfc1_2xCO2_ac, stat, nBoot, 0, opt)
BootStrap_hfc1_4xCO2_ac  = bootstrap_stat(hfc1_4xCO2_ac, stat, nBoot, 0, opt)

BootStrap_hfc2_PD_ac     = bootstrap_stat(hfc2_PD_ac, stat, nBoot, 0, opt)
BootStrap_hfc2_2xCO2_ac  = bootstrap_stat(hfc2_2xCO2_ac, stat, nBoot, 0, opt)
BootStrap_hfc2_4xCO2_ac  = bootstrap_stat(hfc2_4xCO2_ac, stat, nBoot, 0, opt)

BootStrap_mfc1_PD_ac     = bootstrap_stat(mfc1_PD_ac, stat, nBoot, 0, opt)
BootStrap_mfc1_2xCO2_ac  = bootstrap_stat(mfc1_2xCO2_ac, stat, nBoot, 0, opt)
BootStrap_mfc1_4xCO2_ac  = bootstrap_stat(mfc1_4xCO2_ac, stat, nBoot, 0, opt)

BootStrap_mfc2_PD_ac     = bootstrap_stat(mfc2_PD_ac, stat, nBoot, 0, opt)
BootStrap_mfc2_2xCO2_ac  = bootstrap_stat(mfc2_2xCO2_ac, stat, nBoot, 0, opt)
BootStrap_mfc2_4xCO2_ac  = bootstrap_stat(mfc2_4xCO2_ac, stat, nBoot, 0, opt)

;=======================================================================
xBoot_hfc1_PD_ac    = BootStrap_hfc1_PD_ac[1]   ; mean value
xBoot_hfc2_PD_ac    = BootStrap_hfc2_PD_ac[1]   ; mean value
xBoot_mfc1_PD_ac    = BootStrap_mfc1_PD_ac[1]   ; mean value
xBoot_mfc2_PD_ac    = BootStrap_mfc2_PD_ac[1]   ; mean value  
;----------------------------------------------------
xBoot_hfc1_2xCO2_ac    = BootStrap_hfc1_2xCO2_ac[1]   ; mean value
xBoot_hfc2_2xCO2_ac    = BootStrap_hfc2_2xCO2_ac[1]   ; mean value
xBoot_mfc1_2xCO2_ac    = BootStrap_mfc1_2xCO2_ac[1]   ; mean value
xBoot_mfc2_2xCO2_ac    = BootStrap_mfc2_2xCO2_ac[1]   ; mean value  
;----------------------------------------------------   
xBoot_hfc1_4xCO2_ac      = BootStrap_hfc1_4xCO2_ac[1]   ; mean value 
xBoot_hfc2_4xCO2_ac      = BootStrap_hfc2_4xCO2_ac[1]   ; mean value 
xBoot_mfc1_4xCO2_ac      = BootStrap_mfc1_4xCO2_ac[1]   ; mean value 
xBoot_mfc2_4xCO2_ac      = BootStrap_mfc2_4xCO2_ac[1]   ; mean value 
;------------------------------------------------
xBoot_hfc01_PD_ac      = BootStrap_hfc1_PD_ac[0]      
xBoot_hfc02_PD_ac      = BootStrap_hfc2_PD_ac[0]      
xBoot_mfc01_PD_ac      = BootStrap_mfc1_PD_ac[0]      
xBoot_mfc02_PD_ac      = BootStrap_mfc2_PD_ac[0]     
;------------------------------------------------
xBoot_hfc01_2xCO2_ac   = BootStrap_hfc1_2xCO2_ac[0]      
xBoot_hfc02_2xCO2_ac   = BootStrap_hfc2_2xCO2_ac[0]      
xBoot_mfc01_2xCO2_ac   = BootStrap_mfc1_2xCO2_ac[0]      
xBoot_mfc02_2xCO2_ac   = BootStrap_mfc2_2xCO2_ac[0] 
;------------------------------------------------
xBoot_hfc01_4xCO2_ac   = BootStrap_hfc1_4xCO2_ac[0]   
xBoot_hfc02_4xCO2_ac   = BootStrap_hfc2_4xCO2_ac[0]   
xBoot_mfc01_4xCO2_ac   = BootStrap_mfc1_4xCO2_ac[0]   
xBoot_mfc02_4xCO2_ac   = BootStrap_mfc2_4xCO2_ac[0] 
;======================================================
; seasonal cycle
xBootLow_hfc1_PD_ac    = bootstrap_estimate(xBoot_hfc01_PD_ac, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_hfc1_PD_ac    = bootstrap_estimate(xBoot_hfc01_PD_ac, 0.500, False)   ; 50.0% median of bootstrapped estimates  
 xBootHi_hfc1_PD_ac    = bootstrap_estimate(xBoot_hfc01_PD_ac, 0.95, False)   ; 97.5% upper confidence bound 

xBootLow_hfc2_PD_ac    = bootstrap_estimate(xBoot_hfc02_PD_ac, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_hfc2_PD_ac    = bootstrap_estimate(xBoot_hfc02_PD_ac, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_hfc2_PD_ac    = bootstrap_estimate(xBoot_hfc02_PD_ac, 0.95, False)   ; 97.5% upper confidence bound

xBootLow_mfc1_PD_ac    = bootstrap_estimate(xBoot_mfc01_PD_ac, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_mfc1_PD_ac    = bootstrap_estimate(xBoot_mfc01_PD_ac, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_mfc1_PD_ac    = bootstrap_estimate(xBoot_mfc01_PD_ac, 0.95, False)   ; 97.5% upper confidence bound

xBootLow_mfc2_PD_ac    = bootstrap_estimate(xBoot_mfc02_PD_ac, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_mfc2_PD_ac    = bootstrap_estimate(xBoot_mfc02_PD_ac, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_mfc2_PD_ac    = bootstrap_estimate(xBoot_mfc02_PD_ac, 0.95, False)   ; 97.5% upper confidence bound

 ;======================================
 xBootLow_hfc1_2xCO2_ac    = bootstrap_estimate(xBoot_hfc01_2xCO2_ac, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_hfc1_2xCO2_ac    = bootstrap_estimate(xBoot_hfc01_2xCO2_ac, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_hfc1_2xCO2_ac    = bootstrap_estimate(xBoot_hfc01_2xCO2_ac, 0.95, False)   ; 97.5% upper confidence bound

xBootLow_hfc2_2xCO2_ac    = bootstrap_estimate(xBoot_hfc02_2xCO2_ac, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_hfc2_2xCO2_ac    = bootstrap_estimate(xBoot_hfc02_2xCO2_ac, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_hfc2_2xCO2_ac    = bootstrap_estimate(xBoot_hfc02_2xCO2_ac, 0.95, False)   ; 97.5% upper confidence bound

xBootLow_mfc1_2xCO2_ac    = bootstrap_estimate(xBoot_mfc01_2xCO2_ac, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_mfc1_2xCO2_ac    = bootstrap_estimate(xBoot_mfc01_2xCO2_ac, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_mfc1_2xCO2_ac    = bootstrap_estimate(xBoot_mfc01_2xCO2_ac, 0.95, False)   ; 97.5% upper confidence bound

xBootLow_mfc2_2xCO2_ac    = bootstrap_estimate(xBoot_mfc02_2xCO2_ac, 0.05, False)   ;  2.5% lower confidence bound
xBootMed_mfc2_2xCO2_ac    = bootstrap_estimate(xBoot_mfc02_2xCO2_ac, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_mfc2_2xCO2_ac    = bootstrap_estimate(xBoot_mfc02_2xCO2_ac, 0.95, False)   ; 97.5% upper confidence bound
;=======================================
 xBootLow_hfc1_4xCO2_ac    = bootstrap_estimate(xBoot_hfc01_4xCO2_ac, 0.05, False)   ;  2.5% lower confidence bound
 xBootMed_hfc1_4xCO2_ac    = bootstrap_estimate(xBoot_hfc01_4xCO2_ac, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_hfc1_4xCO2_ac    = bootstrap_estimate(xBoot_hfc01_4xCO2_ac, 0.95, False)   ; 97.5% upper confidence bound
 
 xBootLow_hfc2_4xCO2_ac    = bootstrap_estimate(xBoot_hfc02_4xCO2_ac, 0.05, False)   ;  2.5% lower confidence bound
 xBootMed_hfc2_4xCO2_ac    = bootstrap_estimate(xBoot_hfc02_4xCO2_ac, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_hfc2_4xCO2_ac    = bootstrap_estimate(xBoot_hfc02_4xCO2_ac, 0.95, False)   ; 97.5% upper confidence bound
 
 xBootLow_mfc1_4xCO2_ac    = bootstrap_estimate(xBoot_mfc01_4xCO2_ac, 0.05, False)   ;  2.5% lower confidence bound
 xBootMed_mfc1_4xCO2_ac    = bootstrap_estimate(xBoot_mfc01_4xCO2_ac, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_mfc1_4xCO2_ac    = bootstrap_estimate(xBoot_mfc01_4xCO2_ac, 0.95, False)   ; 97.5% upper confidence bound
 
 xBootLow_mfc2_4xCO2_ac    = bootstrap_estimate(xBoot_mfc02_4xCO2_ac, 0.05, False)   ;  2.5% lower confidence bound
 xBootMed_mfc2_4xCO2_ac    = bootstrap_estimate(xBoot_mfc02_4xCO2_ac, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_mfc2_4xCO2_ac    = bootstrap_estimate(xBoot_mfc02_4xCO2_ac, 0.95, False)   ; 97.5% upper confidence bound

  ;======================================
stdarr01 = new((/3/),"double") ; mfc1 ;NEQ
stdarr01(0) =  xBoot_mfc1_PD_ac
stdarr01(1) =  xBoot_mfc1_2xCO2_ac
stdarr01(2) =  xBoot_mfc1_4xCO2_ac

stdarr02 = new((/3/),"double") ; mfc2
stdarr02(0) =  -1 * xBoot_mfc2_PD_ac
stdarr02(1) =  -1 * xBoot_mfc2_2xCO2_ac
stdarr02(2) =  -1 * xBoot_mfc2_4xCO2_ac

stdarr03 = new((/3/),"double") ; hfc1
stdarr03(0) =  -1*xBoot_hfc1_PD_ac
stdarr03(1) =  -1*xBoot_hfc1_2xCO2_ac
stdarr03(2) =  -1*xBoot_hfc1_4xCO2_ac

stdarr04 = new((/3/),"double") ; hfc2
stdarr04(0) =  xBoot_hfc2_PD_ac
stdarr04(1) =  xBoot_hfc2_2xCO2_ac
stdarr04(2) =  xBoot_hfc2_4xCO2_ac


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
 xyres@vpYF    = 0.58
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
 mplot4  = gsn_csm_xy(wks,ispan(1,3,1),stdarr01(:),xyres)   ; advancing the frame
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
 mplot5  = gsn_csm_xy(wks,ispan(1,3,1),stdarr02(:),xyres)   ; advancing the frame
 xyres@vpXF = 0.73
 xyres@vpYF = 0.43
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
 mplot6  = gsn_csm_xy(wks,ispan(1,3,1),stdarr03(:),xyres)   ; advancing the frame
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@tmYROn        = False  ; turn off right tickmarks
 xyres@tmYLOn        = False
 xyres@tmYRBorderOn  = False                       ; turn off right axis
 xyres@tmYLBorderOn  = False
 xyres@vpXF = 0.83
 xyres@vpYF = 0.43
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
 mplot7  = gsn_csm_xy(wks,ispan(1,3,1),stdarr04(:),xyres)   ; advancing the frame


 polyres=True
 polyres@gsLineDashPattern=0
 polyres@gsLineColor="black"
 polyres@gsLineThicknessF=1.0
 polyres@tfPolyDrawOrder = "postDraw"
;
; plot seasonal cycle 
aumm0001     = gsn_add_polyline(wks,mplot4,(/1,1/),(/xBootLow_mfc1_PD_ac,   xBootHi_mfc1_PD_ac/), polyres)
aumm00001    = gsn_add_polyline(wks,mplot4,(/2,2/),(/xBootLow_mfc1_2xCO2_ac,xBootHi_mfc1_2xCO2_ac/), polyres)
aumm000001   = gsn_add_polyline(wks,mplot4,(/3,3/),(/xBootLow_mfc1_4xCO2_ac,xBootHi_mfc1_4xCO2_ac/), polyres)
;
bumm0001    = gsn_add_polyline(wks,mplot5,(/1,1/),(/-1*xBootHi_mfc2_PD_ac,   -1*xBootLow_mfc2_PD_ac/), polyres)
bumm00001   = gsn_add_polyline(wks,mplot5,(/2,2/),(/-1*xBootHi_mfc2_2xCO2_ac,-1*xBootLow_mfc2_2xCO2_ac/), polyres)
bumm000001  = gsn_add_polyline(wks,mplot5,(/3,3/),(/-1*xBootHi_mfc2_4xCO2_ac,-1*xBootLow_mfc2_4xCO2_ac/), polyres)

cumm0001    = gsn_add_polyline(wks,mplot6,(/1,1/),(/-1*xBootHi_hfc1_PD_ac,   -1*xBootLow_hfc1_PD_ac/), polyres)
cumm00001   = gsn_add_polyline(wks,mplot6,(/2,2/),(/-1*xBootHi_hfc1_2xCO2_ac,-1*xBootLow_hfc1_2xCO2_ac/), polyres)
cumm000001  = gsn_add_polyline(wks,mplot6,(/3,3/),(/-1*xBootHi_hfc1_4xCO2_ac,-1*xBootLow_hfc1_4xCO2_ac/), polyres)

dumm0001    = gsn_add_polyline(wks,mplot7,(/1,1/),(/xBootLow_hfc2_PD_ac,   xBootHi_hfc2_PD_ac/), polyres)
dumm00001   = gsn_add_polyline(wks,mplot7,(/2,2/),(/xBootLow_hfc2_2xCO2_ac,xBootHi_hfc2_2xCO2_ac/), polyres)
dumm000001  = gsn_add_polyline(wks,mplot7,(/3,3/),(/xBootLow_hfc2_4xCO2_ac,xBootHi_hfc2_4xCO2_ac/), polyres)
;
wid = 0.24/2
;-------------------------------------------------
; plot seasonal cycle
aum001   = gsn_add_polyline(wks,mplot4,(/1-wid,1+wid/),(/xBootLow_mfc1_PD_ac,xBootLow_mfc1_PD_ac/), polyres)
aum0001  = gsn_add_polyline(wks,mplot4,(/1-wid,1+wid/),(/xBootHi_mfc1_PD_ac,xBootHi_mfc1_PD_ac/), polyres)

bum001   = gsn_add_polyline(wks,mplot4,(/2-wid,2+wid/),(/xBootLow_mfc1_2xCO2_ac,xBootLow_mfc1_2xCO2_ac/), polyres)
bum0001  = gsn_add_polyline(wks,mplot4,(/2-wid,2+wid/),(/xBootHi_mfc1_2xCO2_ac,xBootHi_mfc1_2xCO2_ac/), polyres)

cum001   = gsn_add_polyline(wks,mplot4,(/3-wid,3+wid/),(/xBootLow_mfc1_4xCO2_ac,xBootLow_mfc1_4xCO2_ac/), polyres)
cum0001  = gsn_add_polyline(wks,mplot4,(/3-wid,3+wid/),(/xBootHi_mfc1_4xCO2_ac,xBootHi_mfc1_4xCO2_ac/), polyres)

dum001   = gsn_add_polyline(wks,mplot5,(/1-wid,1+wid/),(/-1*xBootLow_mfc2_PD_ac,-1*xBootLow_mfc2_PD_ac/), polyres)
dum0001  = gsn_add_polyline(wks,mplot5,(/1-wid,1+wid/),(/-1*xBootHi_mfc2_PD_ac, -1*xBootHi_mfc2_PD_ac/), polyres)

eum001   = gsn_add_polyline(wks,mplot5,(/2-wid,2+wid/),(/-1*xBootLow_mfc2_2xCO2_ac,-1*xBootLow_mfc2_2xCO2_ac/), polyres)
eum0001  = gsn_add_polyline(wks,mplot5,(/2-wid,2+wid/),(/-1*xBootHi_mfc2_2xCO2_ac, -1*xBootHi_mfc2_2xCO2_ac/), polyres)

fum001   = gsn_add_polyline(wks,mplot5,(/3-wid,3+wid/),(/-1*xBootLow_mfc2_4xCO2_ac,-1*xBootLow_mfc2_4xCO2_ac/), polyres)
fum0001  = gsn_add_polyline(wks,mplot5,(/3-wid,3+wid/),(/-1*xBootHi_mfc2_4xCO2_ac, -1*xBootHi_mfc2_4xCO2_ac/), polyres)

gum001   = gsn_add_polyline(wks,mplot6,(/1-wid,1+wid/),(/-1*xBootLow_hfc1_PD_ac,-1*xBootLow_hfc1_PD_ac/), polyres)
gum0001  = gsn_add_polyline(wks,mplot6,(/1-wid,1+wid/),(/-1*xBootHi_hfc1_PD_ac, -1*xBootHi_hfc1_PD_ac/), polyres)

hum001   = gsn_add_polyline(wks,mplot6,(/2-wid,2+wid/),(/-1*xBootLow_hfc1_2xCO2_ac,-1*xBootLow_hfc1_2xCO2_ac/), polyres)
hum0001  = gsn_add_polyline(wks,mplot6,(/2-wid,2+wid/),(/-1*xBootHi_hfc1_2xCO2_ac, -1*xBootHi_hfc1_2xCO2_ac/), polyres)

ium001   = gsn_add_polyline(wks,mplot6,(/3-wid,3+wid/),(/-1*xBootLow_hfc1_4xCO2_ac,-1*xBootLow_hfc1_4xCO2_ac/), polyres)
ium0001  = gsn_add_polyline(wks,mplot6,(/3-wid,3+wid/),(/-1*xBootHi_hfc1_4xCO2_ac, -1*xBootHi_hfc1_4xCO2_ac/), polyres)

jum001   = gsn_add_polyline(wks,mplot7,(/1-wid,1+wid/),(/xBootLow_hfc2_PD_ac,xBootLow_hfc2_PD_ac/), polyres)
jum0001  = gsn_add_polyline(wks,mplot7,(/1-wid,1+wid/),(/xBootHi_hfc2_PD_ac, xBootHi_hfc2_PD_ac/), polyres)

kum001   = gsn_add_polyline(wks,mplot7,(/2-wid,2+wid/),(/xBootLow_hfc2_2xCO2_ac,xBootLow_hfc2_2xCO2_ac/), polyres)
kum0001  = gsn_add_polyline(wks,mplot7,(/2-wid,2+wid/),(/xBootHi_hfc2_2xCO2_ac, xBootHi_hfc2_2xCO2_ac/), polyres)

lum001   = gsn_add_polyline(wks,mplot7,(/3-wid,3+wid/),(/xBootLow_hfc2_4xCO2_ac,xBootLow_hfc2_4xCO2_ac/), polyres)
lum0001  = gsn_add_polyline(wks,mplot7,(/3-wid,3+wid/),(/xBootHi_hfc2_4xCO2_ac, xBootHi_hfc2_4xCO2_ac/), polyres)
draw(mplot0)
draw(mplot1)
draw(mplot2)
draw(mplot3)
draw(mplot4)
draw(mplot5)
draw(mplot6)
draw(mplot7)
end 
