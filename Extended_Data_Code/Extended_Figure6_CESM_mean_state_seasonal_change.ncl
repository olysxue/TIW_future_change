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
yrStrt  =  (/71,141,141/)
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

;---->vertical mean 
temp1_xy  = dim_avg_n_Wrap(temp1(:,ilev1:ilev2,:,:),1)
temp2_xy  = dim_avg_n_Wrap(temp2(:,ilev1:ilev2,:,:),1)
temp3_xy  = dim_avg_n_Wrap(temp3(:,ilev1:ilev2,:,:),1)

dlat   = (temp1&lat(2)-temp1&lat(1))*rad
dY     = R*dlat

dTdy1_xy   = center_finite_diff_n(temp1_xy, dY, False, 0, 1)*1000000.
dTdy2_xy   = center_finite_diff_n(temp2_xy, dY, False, 0, 1)*1000000.
dTdy3_xy   = center_finite_diff_n(temp3_xy, dY, False, 0, 1)*1000000.


copy_VarMeta(temp1(:,0,:,:), dTdy1_xy)
copy_VarMeta(temp2(:,0,:,:), dTdy2_xy)
copy_VarMeta(temp3(:,0,:,:), dTdy3_xy)

dTdy1_xy_index  = runave(wgt_areaave_Wrap(dTdy1_xy(:,{2:5},{200:260}), 1.0, 1.0, 0),3,0)
dTdy2_xy_index  = runave(wgt_areaave_Wrap(dTdy2_xy(:,{2:5},{200:260}), 1.0, 1.0, 0),3,0)
dTdy3_xy_index  = runave(wgt_areaave_Wrap(dTdy3_xy(:,{2:5},{200:260}), 1.0, 1.0, 0),3,0)

ntim1  = dimsizes(temp1&time)
ntim2  = dimsizes(temp2&time)
ntim3  = dimsizes(temp3&time)
nyears1 = ntim1/12
nyears2 = ntim2/12
nyears3 = ntim3/12

dTdy1_index_ann  = reshape(dTdy1_xy_index,(/nyears1,12/))
dTdy2_index_ann  = reshape(dTdy2_xy_index,(/nyears2,12/))
dTdy3_index_ann  = reshape(dTdy3_xy_index,(/nyears3,12/))

; annual cycle
dTdy1_index_ave = dim_avg_n_Wrap(dTdy1_index_ann, 0)
dTdy2_index_ave = dim_avg_n_Wrap(dTdy2_index_ann, 0)
dTdy3_index_ave = dim_avg_n_Wrap(dTdy3_index_ann, 0)

; std 
dTdy1_index_std = dim_stddev_n_Wrap(dTdy1_index_ann, 0) 
dTdy2_index_std = dim_stddev_n_Wrap(dTdy2_index_ann, 0) 
dTdy3_index_std = dim_stddev_n_Wrap(dTdy3_index_ann, 0) 

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
;==============================================================================
; calculate domain-avg
temp1_nino3 = wgt_areaave_Wrap(temp1_xy(:,{-2:2},{200:270}), 1.0, 1.0, 0)
temp2_nino3 = wgt_areaave_Wrap(temp2_xy(:,{-2:2},{200:270}), 1.0, 1.0, 0)
temp3_nino3 = wgt_areaave_Wrap(temp3_xy(:,{-2:2},{200:270}), 1.0, 1.0, 0)


temp1_nino3_ann  = reshape(temp1_nino3,(/nyears1,12/))
temp2_nino3_ann  = reshape(temp2_nino3,(/nyears2,12/))
temp3_nino3_ann  = reshape(temp3_nino3,(/nyears3,12/))

; annual cycle
temp1_nino3_ave = dim_avg_n_Wrap(temp1_nino3_ann, 0)
temp2_nino3_ave = dim_avg_n_Wrap(temp2_nino3_ann, 0)
temp3_nino3_ave = dim_avg_n_Wrap(temp3_nino3_ann, 0)

; std 
temp1_nino3_std = dim_stddev_n_Wrap(temp1_nino3_ann, 0) 
temp2_nino3_std = dim_stddev_n_Wrap(temp2_nino3_ann, 0) 
temp3_nino3_std = dim_stddev_n_Wrap(temp3_nino3_ann, 0) 

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

; calculate clm 
clm_temp1 = clmMonTLLL(temp1) ; 12month,lev,lat,lon
clm_temp2 = clmMonTLLL(temp2)
clm_temp3 = clmMonTLLL(temp3)

dlat   = (temp1&lat(2)-temp1&lat(1))*rad
dY     = R*dlat

; month_to_season
clm_temp1_DJF = (clm_temp1(11,:,:,:)+clm_temp1(0,:,:,:)+clm_temp1(1,:,:,:))/3.
clm_temp1_MAM = (clm_temp1(2,:,:,:)+clm_temp1(3,:,:,:)+clm_temp1(4,:,:,:))/3.
clm_temp1_JJA = (clm_temp1(5,:,:,:)+clm_temp1(6,:,:,:)+clm_temp1(7,:,:,:))/3.
clm_temp1_SON = (clm_temp1(8,:,:,:)+clm_temp1(9,:,:,:)+clm_temp1(10,:,:,:))/3.

dTdy1_DJF   = center_finite_diff_n(clm_temp1_DJF, dY, False, 0, 1)*1000000.
dTdy1_MAM   = center_finite_diff_n(clm_temp1_MAM, dY, False, 0, 1)*1000000.
dTdy1_JJA   = center_finite_diff_n(clm_temp1_JJA, dY, False, 0, 1)*1000000.
dTdy1_SON   = center_finite_diff_n(clm_temp1_SON, dY, False, 0, 1)*1000000.

copy_VarMeta(clm_temp1(0,:,:,:), clm_temp1_DJF)
copy_VarMeta(clm_temp1(0,:,:,:), clm_temp1_MAM)
copy_VarMeta(clm_temp1(0,:,:,:), clm_temp1_JJA)
copy_VarMeta(clm_temp1(0,:,:,:), clm_temp1_SON)

copy_VarMeta(clm_temp1(0,:,:,:), dTdy1_DJF)
copy_VarMeta(clm_temp1(0,:,:,:), dTdy1_MAM)
copy_VarMeta(clm_temp1(0,:,:,:), dTdy1_JJA)
copy_VarMeta(clm_temp1(0,:,:,:), dTdy1_SON)

;---->vertical mean 
clm_temp1_DJF_xy  = dim_avg_n_Wrap(clm_temp1_DJF(ilev1:ilev2,:,:),0)
clm_temp1_MAM_xy  = dim_avg_n_Wrap(clm_temp1_MAM(ilev1:ilev2,:,:),0)
clm_temp1_JJA_xy  = dim_avg_n_Wrap(clm_temp1_JJA(ilev1:ilev2,:,:),0)
clm_temp1_SON_xy  = dim_avg_n_Wrap(clm_temp1_SON(ilev1:ilev2,:,:),0)

dTdy1_DJF_xy  = dim_avg_n_Wrap(dTdy1_DJF(ilev1:ilev2,:,:),0)
dTdy1_MAM_xy  = dim_avg_n_Wrap(dTdy1_MAM(ilev1:ilev2,:,:),0)
dTdy1_JJA_xy  = dim_avg_n_Wrap(dTdy1_JJA(ilev1:ilev2,:,:),0)
dTdy1_SON_xy  = dim_avg_n_Wrap(dTdy1_SON(ilev1:ilev2,:,:),0)

;---->longitude mean 
clm_temp1_DJF_yz  = dim_avg_n_Wrap(clm_temp1_DJF(:,:,{lonw:lone}),2)
clm_temp1_MAM_yz  = dim_avg_n_Wrap(clm_temp1_MAM(:,:,{lonw:lone}),2)
clm_temp1_JJA_yz  = dim_avg_n_Wrap(clm_temp1_JJA(:,:,{lonw:lone}),2)
clm_temp1_SON_yz  = dim_avg_n_Wrap(clm_temp1_SON(:,:,{lonw:lone}),2)

dTdy1_DJF_yz  = dim_avg_n_Wrap(dTdy1_DJF(:,:,{lonw:lone}),2)
dTdy1_MAM_yz  = dim_avg_n_Wrap(dTdy1_MAM(:,:,{lonw:lone}),2)
dTdy1_JJA_yz  = dim_avg_n_Wrap(dTdy1_JJA(:,:,{lonw:lone}),2)
dTdy1_SON_yz  = dim_avg_n_Wrap(dTdy1_SON(:,:,{lonw:lone}),2)

;========================================================================
clm_temp2_DJF = (clm_temp2(11,:,:,:)+clm_temp2(0,:,:,:)+clm_temp2(1,:,:,:))/3.
clm_temp2_MAM = (clm_temp2(2,:,:,:) +clm_temp2(3,:,:,:)+clm_temp2(4,:,:,:))/3.
clm_temp2_JJA = (clm_temp2(5,:,:,:) +clm_temp2(6,:,:,:)+clm_temp2(7,:,:,:))/3.
clm_temp2_SON = (clm_temp2(8,:,:,:) +clm_temp2(9,:,:,:)+clm_temp2(10,:,:,:))/3.

copy_VarMeta(clm_temp2(0,:,:,:), clm_temp2_DJF)
copy_VarMeta(clm_temp2(0,:,:,:), clm_temp2_MAM)
copy_VarMeta(clm_temp2(0,:,:,:), clm_temp2_JJA)
copy_VarMeta(clm_temp2(0,:,:,:), clm_temp2_SON)

dTdy2_DJF   = center_finite_diff_n(clm_temp2_DJF, dY, False, 0, 1)*1000000.
dTdy2_MAM   = center_finite_diff_n(clm_temp2_MAM, dY, False, 0, 1)*1000000.
dTdy2_JJA   = center_finite_diff_n(clm_temp2_JJA, dY, False, 0, 1)*1000000.
dTdy2_SON   = center_finite_diff_n(clm_temp2_SON, dY, False, 0, 1)*1000000.

copy_VarMeta(clm_temp2(0,:,:,:), dTdy2_DJF)
copy_VarMeta(clm_temp2(0,:,:,:), dTdy2_MAM)
copy_VarMeta(clm_temp2(0,:,:,:), dTdy2_JJA)
copy_VarMeta(clm_temp2(0,:,:,:), dTdy2_SON)

;---->vertical mean 
clm_temp2_DJF_xy  = dim_avg_n_Wrap(clm_temp2_DJF(ilev1:ilev2,:,:),0)
clm_temp2_MAM_xy  = dim_avg_n_Wrap(clm_temp2_MAM(ilev1:ilev2,:,:),0)
clm_temp2_JJA_xy  = dim_avg_n_Wrap(clm_temp2_JJA(ilev1:ilev2,:,:),0)
clm_temp2_SON_xy  = dim_avg_n_Wrap(clm_temp2_SON(ilev1:ilev2,:,:),0)

dTdy2_DJF_xy  = dim_avg_n_Wrap(dTdy2_DJF(ilev1:ilev2,:,:),0)
dTdy2_MAM_xy  = dim_avg_n_Wrap(dTdy2_MAM(ilev1:ilev2,:,:),0)
dTdy2_JJA_xy  = dim_avg_n_Wrap(dTdy2_JJA(ilev1:ilev2,:,:),0)
dTdy2_SON_xy  = dim_avg_n_Wrap(dTdy2_SON(ilev1:ilev2,:,:),0)

;---->longitude mean 
clm_temp2_DJF_yz  = dim_avg_n_Wrap(clm_temp2_DJF(:,:,{lonw:lone}),2)
clm_temp2_MAM_yz  = dim_avg_n_Wrap(clm_temp2_MAM(:,:,{lonw:lone}),2)
clm_temp2_JJA_yz  = dim_avg_n_Wrap(clm_temp2_JJA(:,:,{lonw:lone}),2)
clm_temp2_SON_yz  = dim_avg_n_Wrap(clm_temp2_SON(:,:,{lonw:lone}),2)

dTdy2_DJF_yz  = dim_avg_n_Wrap(dTdy2_DJF(:,:,{lonw:lone}),2)
dTdy2_MAM_yz  = dim_avg_n_Wrap(dTdy2_MAM(:,:,{lonw:lone}),2)
dTdy2_JJA_yz  = dim_avg_n_Wrap(dTdy2_JJA(:,:,{lonw:lone}),2)
dTdy2_SON_yz  = dim_avg_n_Wrap(dTdy2_SON(:,:,{lonw:lone}),2)
;=======================================================================
clm_temp3_DJF = (clm_temp3(11,:,:,:)+clm_temp3(0,:,:,:)+clm_temp3(1,:,:,:))/3.
clm_temp3_MAM = (clm_temp3(2,:,:,:) +clm_temp3(3,:,:,:)+clm_temp3(4,:,:,:))/3.
clm_temp3_JJA = (clm_temp3(5,:,:,:) +clm_temp3(6,:,:,:)+clm_temp3(7,:,:,:))/3.
clm_temp3_SON = (clm_temp3(8,:,:,:) +clm_temp3(9,:,:,:)+clm_temp3(10,:,:,:))/3.

copy_VarMeta(clm_temp3(0,:,:,:), clm_temp3_DJF)
copy_VarMeta(clm_temp3(0,:,:,:), clm_temp3_MAM)
copy_VarMeta(clm_temp3(0,:,:,:), clm_temp3_JJA)
copy_VarMeta(clm_temp3(0,:,:,:), clm_temp3_SON)

dTdy3_DJF   = center_finite_diff_n(clm_temp3_DJF, dY, False, 0, 1)*1000000.
dTdy3_MAM   = center_finite_diff_n(clm_temp3_MAM, dY, False, 0, 1)*1000000.
dTdy3_JJA   = center_finite_diff_n(clm_temp3_JJA, dY, False, 0, 1)*1000000.
dTdy3_SON   = center_finite_diff_n(clm_temp3_SON, dY, False, 0, 1)*1000000.

copy_VarMeta(clm_temp3(0,:,:,:), dTdy3_DJF)
copy_VarMeta(clm_temp3(0,:,:,:), dTdy3_MAM)
copy_VarMeta(clm_temp3(0,:,:,:), dTdy3_JJA)
copy_VarMeta(clm_temp3(0,:,:,:), dTdy3_SON)

;---->vertical mean 
clm_temp3_DJF_xy  = dim_avg_n_Wrap(clm_temp3_DJF(ilev1:ilev2,:,:),0)
clm_temp3_MAM_xy  = dim_avg_n_Wrap(clm_temp3_MAM(ilev1:ilev2,:,:),0)
clm_temp3_JJA_xy  = dim_avg_n_Wrap(clm_temp3_JJA(ilev1:ilev2,:,:),0)
clm_temp3_SON_xy  = dim_avg_n_Wrap(clm_temp3_SON(ilev1:ilev2,:,:),0)

dTdy3_DJF_xy  = dim_avg_n_Wrap(dTdy3_DJF(ilev1:ilev2,:,:),0)
dTdy3_MAM_xy  = dim_avg_n_Wrap(dTdy3_MAM(ilev1:ilev2,:,:),0)
dTdy3_JJA_xy  = dim_avg_n_Wrap(dTdy3_JJA(ilev1:ilev2,:,:),0)
dTdy3_SON_xy  = dim_avg_n_Wrap(dTdy3_SON(ilev1:ilev2,:,:),0)

;---->longitude mean 
clm_temp3_DJF_yz  = dim_avg_n_Wrap(clm_temp3_DJF(:,:,{lonw:lone}),2)
clm_temp3_MAM_yz  = dim_avg_n_Wrap(clm_temp3_MAM(:,:,{lonw:lone}),2)
clm_temp3_JJA_yz  = dim_avg_n_Wrap(clm_temp3_JJA(:,:,{lonw:lone}),2)
clm_temp3_SON_yz  = dim_avg_n_Wrap(clm_temp3_SON(:,:,{lonw:lone}),2)

dTdy3_DJF_yz  = dim_avg_n_Wrap(dTdy3_DJF(:,:,{lonw:lone}),2)
dTdy3_MAM_yz  = dim_avg_n_Wrap(dTdy3_MAM(:,:,{lonw:lone}),2)
dTdy3_JJA_yz  = dim_avg_n_Wrap(dTdy3_JJA(:,:,{lonw:lone}),2)
dTdy3_SON_yz  = dim_avg_n_Wrap(dTdy3_SON(:,:,{lonw:lone}),2)
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
ncdir2  = "/Volumes/oly2/CESM/monthly/"
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
printVarSummary(uvel1)
ntim1  = dimsizes(uvel1&time)
ntim2  = dimsizes(uvel2&time)
ntim3  = dimsizes(uvel3&time)
nyears1 = ntim1/12
nyears2 = ntim2/12
nyears3 = ntim3/12
uvel1_xy  = dim_avg_n_Wrap(uvel1(:,ilev1:ilev2,:,:),1)
uvel2_xy  = dim_avg_n_Wrap(uvel2(:,ilev1:ilev2,:,:),1)
uvel3_xy  = dim_avg_n_Wrap(uvel3(:,ilev1:ilev2,:,:),1)
uvel01_xy  = dim_avg_n_Wrap(uvel1(:,{6000:15000},:,:),1)
uvel02_xy  = dim_avg_n_Wrap(uvel2(:,{6000:15000},:,:),1)
uvel03_xy  = dim_avg_n_Wrap(uvel3(:,{6000:15000},:,:),1)
printVarSummary(uvel1_xy)
; calculate EUC domain-avg
EUC1 = wgt_areaave_Wrap(uvel01_xy(:,{-1:1},{210:250}), 1.0, 1.0, 0)
EUC2 = wgt_areaave_Wrap(uvel02_xy(:,{-1:1},{210:250}), 1.0, 1.0, 0)
EUC3 = wgt_areaave_Wrap(uvel03_xy(:,{-1:1},{210:250}), 1.0, 1.0, 0)
printVarSummary(EUC1)

EUC1_ann  = reshape(EUC1,(/nyears1,12/))
EUC2_ann  = reshape(EUC2,(/nyears2,12/))
EUC3_ann  = reshape(EUC3,(/nyears3,12/))

; annual cycle
EUC1_ave = dim_avg_n_Wrap(EUC1_ann, 0)
EUC2_ave = dim_avg_n_Wrap(EUC2_ann, 0)
EUC3_ave = dim_avg_n_Wrap(EUC3_ann, 0)

; std 
EUC1_std = dim_stddev_n_Wrap(EUC1_ann, 0) 
EUC2_std = dim_stddev_n_Wrap(EUC2_ann, 0) 
EUC3_std = dim_stddev_n_Wrap(EUC3_ann, 0) 

xp1_EUC_PD = new((/24/), float)
yp1_EUC_PD = new((/24/), float)

do i=0,11
   xp1_EUC_PD(i)          = i+1
   xp1_EUC_PD(2*12-1-i)   = i+1
   yp1_EUC_PD(i)          = EUC1_ave(i) + EUC1_std(i)/2.
   yp1_EUC_PD(2*12-1-i)   = EUC1_ave(i) - EUC1_std(i)/2. 
 end do

 xp1_EUC_4xCO2 = new((/24/), float)
 yp1_EUC_4xCO2 = new((/24/), float)
 
 do i=0,11
    xp1_EUC_4xCO2(i)          = i+1
    xp1_EUC_4xCO2(2*12-1-i)   = i+1
    yp1_EUC_4xCO2(i)          = EUC3_ave(i) + EUC3_std(i)/2.
    yp1_EUC_4xCO2(2*12-1-i)   = EUC3_ave(i) - EUC3_std(i)/2. 
  end do
;---------------------------------------------------------------------------------
; calculate NECC domain-avg
NECC1 = wgt_areaave_Wrap(uvel1_xy(:,{5:7},{210:250}), 1.0, 1.0, 0)
NECC2 = wgt_areaave_Wrap(uvel2_xy(:,{5:7},{210:250}), 1.0, 1.0, 0)
NECC3 = wgt_areaave_Wrap(uvel3_xy(:,{5:7},{210:250}), 1.0, 1.0, 0)

NECC1_ann  = reshape(NECC1,(/nyears1,12/))
NECC2_ann  = reshape(NECC2,(/nyears2,12/))
NECC3_ann  = reshape(NECC3,(/nyears3,12/))

; annual cycle
NECC1_ave = dim_avg_n_Wrap(NECC1_ann, 0)
NECC2_ave = dim_avg_n_Wrap(NECC2_ann, 0)
NECC3_ave = dim_avg_n_Wrap(NECC3_ann, 0)

; std 
NECC1_std = dim_stddev_n_Wrap(NECC1_ann, 0) 
NECC2_std = dim_stddev_n_Wrap(NECC2_ann, 0) 
NECC3_std = dim_stddev_n_Wrap(NECC3_ann, 0) 

xp1_NECC_PD = new((/24/), float)
yp1_NECC_PD = new((/24/), float)

do i=0,11
   xp1_NECC_PD(i)          = i+1
   xp1_NECC_PD(2*12-1-i)   = i+1
   yp1_NECC_PD(i)          = NECC1_ave(i) + NECC1_std(i)/2.
   yp1_NECC_PD(2*12-1-i)   = NECC1_ave(i) - NECC1_std(i)/2. 
 end do

 xp1_NECC_4xCO2 = new((/24/), float)
 yp1_NECC_4xCO2 = new((/24/), float)
 
 do i=0,11
    xp1_NECC_4xCO2(i)          = i+1
    xp1_NECC_4xCO2(2*12-1-i)   = i+1
    yp1_NECC_4xCO2(i)          = NECC3_ave(i) + NECC3_std(i)/2.
    yp1_NECC_4xCO2(2*12-1-i)   = NECC3_ave(i) - NECC3_std(i)/2. 
  end do
;---------------------------------------------------------------------------
; calculate SEC domain-avg
SEC1 = wgt_areaave_Wrap(uvel1_xy(:,{1:3},{210:250}), 1.0, 1.0, 0)
SEC2 = wgt_areaave_Wrap(uvel2_xy(:,{1:3},{210:250}), 1.0, 1.0, 0)
SEC3 = wgt_areaave_Wrap(uvel3_xy(:,{1:3},{210:250}), 1.0, 1.0, 0)

SEC1_ann  = reshape(SEC1,(/nyears1,12/))
SEC2_ann  = reshape(SEC2,(/nyears2,12/))
SEC3_ann  = reshape(SEC3,(/nyears3,12/))

; annual cycle
SEC1_ave = dim_avg_n_Wrap(SEC1_ann, 0)
SEC2_ave = dim_avg_n_Wrap(SEC2_ann, 0)
SEC3_ave = dim_avg_n_Wrap(SEC3_ann, 0)

; std 
SEC1_std = dim_stddev_n_Wrap(SEC1_ann, 0) 
SEC2_std = dim_stddev_n_Wrap(SEC2_ann, 0) 
SEC3_std = dim_stddev_n_Wrap(SEC3_ann, 0) 

xp1_SEC_PD = new((/24/), float)
yp1_SEC_PD = new((/24/), float)

do i=0,11
   xp1_SEC_PD(i)          = i+1
   xp1_SEC_PD(2*12-1-i)   = i+1
   yp1_SEC_PD(i)          = SEC1_ave(i) + SEC1_std(i)/2.
   yp1_SEC_PD(2*12-1-i)   = SEC1_ave(i) - SEC1_std(i)/2. 
 end do

 xp1_SEC_4xCO2 = new((/24/), float)
 yp1_SEC_4xCO2 = new((/24/), float)
 
 do i=0,11
    xp1_SEC_4xCO2(i)          = i+1
    xp1_SEC_4xCO2(2*12-1-i)   = i+1
    yp1_SEC_4xCO2(i)          = SEC3_ave(i) + SEC3_std(i)/2.
    yp1_SEC_4xCO2(2*12-1-i)   = SEC3_ave(i) - SEC3_std(i)/2. 
  end do
;-------------------------------------------------------------------------
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

ntim1  = dimsizes(uvel1&time)
nyears01 = ntim1/12
ntim2  = dimsizes(uvel2&time)
nyears02 = ntim2/12
ntim3  = dimsizes(uvel3&time)
nyears03 = ntim3/12

dUdy1_index_ann  = reshape(dUdy1_xy_index,(/nyears01,12/))
dUdy2_index_ann  = reshape(dUdy2_xy_index,(/nyears02,12/))
dUdy3_index_ann  = reshape(dUdy3_xy_index,(/nyears03,12/))

; annual cycle
dUdy1_index_ave = dim_avg_n_Wrap(dUdy1_index_ann, 0)
dUdy2_index_ave = dim_avg_n_Wrap(dUdy2_index_ann, 0)
dUdy3_index_ave = dim_avg_n_Wrap(dUdy3_index_ann, 0)

; std 
dUdy1_index_std = dim_stddev_n_Wrap(dUdy1_index_ann, 0) 
dUdy2_index_std = dim_stddev_n_Wrap(dUdy2_index_ann, 0) 
dUdy3_index_std = dim_stddev_n_Wrap(dUdy3_index_ann, 0) 

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

; calculate clm 
clm_uvel1 = clmMonTLLL(uvel1) ; 12month,lev,lat,lon
clm_uvel2 = clmMonTLLL(uvel2)
clm_uvel3 = clmMonTLLL(uvel3)

dlat   = (uvel1&lat(2)-uvel1&lat(1))*rad
dY     = R*dlat

; month_to_season
clm_uvel1_DJF = (clm_uvel1(11,:,:,:)+clm_uvel1(0,:,:,:)+clm_uvel1(1,:,:,:))/3.
clm_uvel1_MAM = (clm_uvel1(2,:,:,:) +clm_uvel1(3,:,:,:)+clm_uvel1(4,:,:,:))/3.
clm_uvel1_JJA = (clm_uvel1(5,:,:,:) +clm_uvel1(6,:,:,:)+clm_uvel1(7,:,:,:))/3.
clm_uvel1_SON = (clm_uvel1(8,:,:,:) +clm_uvel1(9,:,:,:)+clm_uvel1(10,:,:,:))/3.

dUdy1_DJF   = center_finite_diff_n(clm_uvel1_DJF, dY, False, 0, 1)*1000000.
dUdy1_MAM   = center_finite_diff_n(clm_uvel1_MAM, dY, False, 0, 1)*1000000.
dUdy1_JJA   = center_finite_diff_n(clm_uvel1_JJA, dY, False, 0, 1)*1000000.
dUdy1_SON   = center_finite_diff_n(clm_uvel1_SON, dY, False, 0, 1)*1000000.

copy_VarMeta(clm_uvel1(0,:,:,:), clm_uvel1_DJF)
copy_VarMeta(clm_uvel1(0,:,:,:), clm_uvel1_MAM)
copy_VarMeta(clm_uvel1(0,:,:,:), clm_uvel1_JJA)
copy_VarMeta(clm_uvel1(0,:,:,:), clm_uvel1_SON)

copy_VarMeta(clm_uvel1(0,:,:,:), dUdy1_DJF)
copy_VarMeta(clm_uvel1(0,:,:,:), dUdy1_MAM)
copy_VarMeta(clm_uvel1(0,:,:,:), dUdy1_JJA)
copy_VarMeta(clm_uvel1(0,:,:,:), dUdy1_SON)

;---->vertical mean 
clm_uvel1_DJF_xy  = dim_avg_n_Wrap(clm_uvel1_DJF(ilev1:ilev2,:,:),0)
clm_uvel1_MAM_xy  = dim_avg_n_Wrap(clm_uvel1_MAM(ilev1:ilev2,:,:),0)
clm_uvel1_JJA_xy  = dim_avg_n_Wrap(clm_uvel1_JJA(ilev1:ilev2,:,:),0)
clm_uvel1_SON_xy  = dim_avg_n_Wrap(clm_uvel1_SON(ilev1:ilev2,:,:),0)

dUdy1_DJF_xy  = dim_avg_n_Wrap(dUdy1_DJF(ilev1:ilev2,:,:),0)
dUdy1_MAM_xy  = dim_avg_n_Wrap(dUdy1_MAM(ilev1:ilev2,:,:),0)
dUdy1_JJA_xy  = dim_avg_n_Wrap(dUdy1_JJA(ilev1:ilev2,:,:),0)
dUdy1_SON_xy  = dim_avg_n_Wrap(dUdy1_SON(ilev1:ilev2,:,:),0)

;---->longitude mean 
clm_uvel1_DJF_yz  = dim_avg_n_Wrap(clm_uvel1_DJF(:,:,{lonw:lone}),2)
clm_uvel1_MAM_yz  = dim_avg_n_Wrap(clm_uvel1_MAM(:,:,{lonw:lone}),2)
clm_uvel1_JJA_yz  = dim_avg_n_Wrap(clm_uvel1_JJA(:,:,{lonw:lone}),2)
clm_uvel1_SON_yz  = dim_avg_n_Wrap(clm_uvel1_SON(:,:,{lonw:lone}),2)

dUdy1_DJF_yz  = dim_avg_n_Wrap(dUdy1_DJF(:,:,{lonw:lone}),2)
dUdy1_MAM_yz  = dim_avg_n_Wrap(dUdy1_MAM(:,:,{lonw:lone}),2)
dUdy1_JJA_yz  = dim_avg_n_Wrap(dUdy1_JJA(:,:,{lonw:lone}),2)
dUdy1_SON_yz  = dim_avg_n_Wrap(dUdy1_SON(:,:,{lonw:lone}),2)

;========================================================================
clm_uvel2_DJF = (clm_uvel2(11,:,:,:)+clm_uvel2(0,:,:,:)+clm_uvel2(1,:,:,:))/3.
clm_uvel2_MAM = (clm_uvel2(2,:,:,:) +clm_uvel2(3,:,:,:)+clm_uvel2(4,:,:,:))/3.
clm_uvel2_JJA = (clm_uvel2(5,:,:,:) +clm_uvel2(6,:,:,:)+clm_uvel2(7,:,:,:))/3.
clm_uvel2_SON = (clm_uvel2(8,:,:,:) +clm_uvel2(9,:,:,:)+clm_uvel2(10,:,:,:))/3.

copy_VarMeta(clm_uvel2(0,:,:,:), clm_uvel2_DJF)
copy_VarMeta(clm_uvel2(0,:,:,:), clm_uvel2_MAM)
copy_VarMeta(clm_uvel2(0,:,:,:), clm_uvel2_JJA)
copy_VarMeta(clm_uvel2(0,:,:,:), clm_uvel2_SON)

dUdy2_DJF   = center_finite_diff_n(clm_uvel2_DJF, dY, False, 0, 1)*1000000.
dUdy2_MAM   = center_finite_diff_n(clm_uvel2_MAM, dY, False, 0, 1)*1000000.
dUdy2_JJA   = center_finite_diff_n(clm_uvel2_JJA, dY, False, 0, 1)*1000000.
dUdy2_SON   = center_finite_diff_n(clm_uvel2_SON, dY, False, 0, 1)*1000000.

copy_VarMeta(clm_uvel2(0,:,:,:), dUdy2_DJF)
copy_VarMeta(clm_uvel2(0,:,:,:), dUdy2_MAM)
copy_VarMeta(clm_uvel2(0,:,:,:), dUdy2_JJA)
copy_VarMeta(clm_uvel2(0,:,:,:), dUdy2_SON)

;---->vertical mean 
clm_uvel2_DJF_xy  = dim_avg_n_Wrap(clm_uvel2_DJF(ilev1:ilev2,:,:),0)
clm_uvel2_MAM_xy  = dim_avg_n_Wrap(clm_uvel2_MAM(ilev1:ilev2,:,:),0)
clm_uvel2_JJA_xy  = dim_avg_n_Wrap(clm_uvel2_JJA(ilev1:ilev2,:,:),0)
clm_uvel2_SON_xy  = dim_avg_n_Wrap(clm_uvel2_SON(ilev1:ilev2,:,:),0)

dUdy2_DJF_xy  = dim_avg_n_Wrap(dUdy2_DJF(ilev1:ilev2,:,:),0)
dUdy2_MAM_xy  = dim_avg_n_Wrap(dUdy2_MAM(ilev1:ilev2,:,:),0)
dUdy2_JJA_xy  = dim_avg_n_Wrap(dUdy2_JJA(ilev1:ilev2,:,:),0)
dUdy2_SON_xy  = dim_avg_n_Wrap(dUdy2_SON(ilev1:ilev2,:,:),0)

;---->longitude mean 
clm_uvel2_DJF_yz  = dim_avg_n_Wrap(clm_uvel2_DJF(:,:,{lonw:lone}),2)
clm_uvel2_MAM_yz  = dim_avg_n_Wrap(clm_uvel2_MAM(:,:,{lonw:lone}),2)
clm_uvel2_JJA_yz  = dim_avg_n_Wrap(clm_uvel2_JJA(:,:,{lonw:lone}),2)
clm_uvel2_SON_yz  = dim_avg_n_Wrap(clm_uvel2_SON(:,:,{lonw:lone}),2)

dUdy2_DJF_yz  = dim_avg_n_Wrap(dUdy2_DJF(:,:,{lonw:lone}),2)
dUdy2_MAM_yz  = dim_avg_n_Wrap(dUdy2_MAM(:,:,{lonw:lone}),2)
dUdy2_JJA_yz  = dim_avg_n_Wrap(dUdy2_JJA(:,:,{lonw:lone}),2)
dUdy2_SON_yz  = dim_avg_n_Wrap(dUdy2_SON(:,:,{lonw:lone}),2)
;=======================================================================
clm_uvel3_DJF = (clm_uvel3(11,:,:,:)+clm_uvel3(0,:,:,:)+clm_uvel3(1,:,:,:))/3.
clm_uvel3_MAM = (clm_uvel3(2,:,:,:) +clm_uvel3(3,:,:,:)+clm_uvel3(4,:,:,:))/3.
clm_uvel3_JJA = (clm_uvel3(5,:,:,:) +clm_uvel3(6,:,:,:)+clm_uvel3(7,:,:,:))/3.
clm_uvel3_SON = (clm_uvel3(8,:,:,:) +clm_uvel3(9,:,:,:)+clm_uvel3(10,:,:,:))/3.

copy_VarMeta(clm_uvel3(0,:,:,:), clm_uvel3_DJF)
copy_VarMeta(clm_uvel3(0,:,:,:), clm_uvel3_MAM)
copy_VarMeta(clm_uvel3(0,:,:,:), clm_uvel3_JJA)
copy_VarMeta(clm_uvel3(0,:,:,:), clm_uvel3_SON)

dUdy3_DJF   = center_finite_diff_n(clm_uvel3_DJF, dY, False, 0, 1)*1000000.
dUdy3_MAM   = center_finite_diff_n(clm_uvel3_MAM, dY, False, 0, 1)*1000000.
dUdy3_JJA   = center_finite_diff_n(clm_uvel3_JJA, dY, False, 0, 1)*1000000.
dUdy3_SON   = center_finite_diff_n(clm_uvel3_SON, dY, False, 0, 1)*1000000.

copy_VarMeta(clm_uvel3(0,:,:,:), dUdy3_DJF)
copy_VarMeta(clm_uvel3(0,:,:,:), dUdy3_MAM)
copy_VarMeta(clm_uvel3(0,:,:,:), dUdy3_JJA)
copy_VarMeta(clm_uvel3(0,:,:,:), dUdy3_SON)

;---->vertical mean 
clm_uvel3_DJF_xy  = dim_avg_n_Wrap(clm_uvel3_DJF(ilev1:ilev2,:,:),0)
clm_uvel3_MAM_xy  = dim_avg_n_Wrap(clm_uvel3_MAM(ilev1:ilev2,:,:),0)
clm_uvel3_JJA_xy  = dim_avg_n_Wrap(clm_uvel3_JJA(ilev1:ilev2,:,:),0)
clm_uvel3_SON_xy  = dim_avg_n_Wrap(clm_uvel3_SON(ilev1:ilev2,:,:),0)

dUdy3_DJF_xy  = dim_avg_n_Wrap(dUdy3_DJF(ilev1:ilev2,:,:),0)
dUdy3_MAM_xy  = dim_avg_n_Wrap(dUdy3_MAM(ilev1:ilev2,:,:),0)
dUdy3_JJA_xy  = dim_avg_n_Wrap(dUdy3_JJA(ilev1:ilev2,:,:),0)
dUdy3_SON_xy  = dim_avg_n_Wrap(dUdy3_SON(ilev1:ilev2,:,:),0)

;---->longitude mean 
clm_uvel3_DJF_yz  = dim_avg_n_Wrap(clm_uvel3_DJF(:,:,{lonw:lone}),2)
clm_uvel3_MAM_yz  = dim_avg_n_Wrap(clm_uvel3_MAM(:,:,{lonw:lone}),2)
clm_uvel3_JJA_yz  = dim_avg_n_Wrap(clm_uvel3_JJA(:,:,{lonw:lone}),2)
clm_uvel3_SON_yz  = dim_avg_n_Wrap(clm_uvel3_SON(:,:,{lonw:lone}),2)

dUdy3_DJF_yz  = dim_avg_n_Wrap(dUdy3_DJF(:,:,{lonw:lone}),2)
dUdy3_MAM_yz  = dim_avg_n_Wrap(dUdy3_MAM(:,:,{lonw:lone}),2)
dUdy3_JJA_yz  = dim_avg_n_Wrap(dUdy3_JJA(:,:,{lonw:lone}),2)
dUdy3_SON_yz  = dim_avg_n_Wrap(dUdy3_SON(:,:,{lonw:lone}),2)

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

ncdir3  = "/Volumes/oly2/CESM/monthly/"
fil001    = systemfunc("ls "+ ncdir3 + "TAUX.mon.0071-0140.20S_20N.120E_80W_1degree.nc")
fil002    = systemfunc("ls "+ ncdir3 + "TAUX.2xCO2.mon.0071-0170.20S_20N.120E_80W_1degree.nc")
fil003    = systemfunc("ls "+ ncdir3 + "TAUX.4xCO2.mon.0071-0170.20S_20N.120E_80W_1degree.nc")
   
f001      = addfile(fil001, "r")
f002      = addfile(fil002, "r")
f003      = addfile(fil003, "r")

YYYY001= cd_calendar(f001->time,-1)/100
iyr001 = ind(YYYY001.ge.yrStrt(0) .and. YYYY001.le.yrLast(0))

YYYY002= cd_calendar(f002->time,-1)/100
iyr002 = ind(YYYY002.ge.yrStrt(1) .and. YYYY002.le.yrLast(1))

YYYY003= cd_calendar(f003->time,-1)/100
iyr003 = ind(YYYY003.ge.yrStrt(2) .and. YYYY003.le.yrLast(2))

taux1  = f001->TAUX(iyr001,:,:)
taux2  = f002->TAUX(iyr002,:,:)
taux3  = f003->TAUX(iyr003,:,:)

; calculate clm 
clm_taux1 = clmMonTLL(taux1) ; 12month,lev,lat,lon
clm_taux2 = clmMonTLL(taux2)
clm_taux3 = clmMonTLL(taux3)

; month_to_season
clm_taux1_DJF = (clm_taux1(11,:,:)+clm_taux1(0,:,:)+clm_taux1(1,:,:))/3.
clm_taux1_MAM = (clm_taux1(2,:,:) +clm_taux1(3,:,:)+clm_taux1(4,:,:))/3.
clm_taux1_JJA = (clm_taux1(5,:,:) +clm_taux1(6,:,:)+clm_taux1(7,:,:))/3.
clm_taux1_SON = (clm_taux1(8,:,:) +clm_taux1(9,:,:)+clm_taux1(10,:,:))/3.

copy_VarMeta(clm_taux1(0,:,:), clm_taux1_DJF)
copy_VarMeta(clm_taux1(0,:,:), clm_taux1_MAM)
copy_VarMeta(clm_taux1(0,:,:), clm_taux1_JJA)
copy_VarMeta(clm_taux1(0,:,:), clm_taux1_SON)

clm_taux2_DJF = (clm_taux2(11,:,:)+clm_taux2(0,:,:)+clm_taux2(1,:,:))/3.
clm_taux2_MAM = (clm_taux2(2,:,:) +clm_taux2(3,:,:)+clm_taux2(4,:,:))/3.
clm_taux2_JJA = (clm_taux2(5,:,:) +clm_taux2(6,:,:)+clm_taux2(7,:,:))/3.
clm_taux2_SON = (clm_taux2(8,:,:) +clm_taux2(9,:,:)+clm_taux2(10,:,:))/3.

copy_VarMeta(clm_taux2(0,:,:), clm_taux2_DJF)
copy_VarMeta(clm_taux2(0,:,:), clm_taux2_MAM)
copy_VarMeta(clm_taux2(0,:,:), clm_taux2_JJA)
copy_VarMeta(clm_taux2(0,:,:), clm_taux2_SON)

clm_taux3_DJF = (clm_taux3(11,:,:)+clm_taux3(0,:,:)+clm_taux3(1,:,:))/3.
clm_taux3_MAM = (clm_taux3(2,:,:) +clm_taux3(3,:,:)+clm_taux3(4,:,:))/3.
clm_taux3_JJA = (clm_taux3(5,:,:) +clm_taux3(6,:,:)+clm_taux3(7,:,:))/3.
clm_taux3_SON = (clm_taux3(8,:,:) +clm_taux3(9,:,:)+clm_taux3(10,:,:))/3.

copy_VarMeta(clm_taux3(0,:,:), clm_taux3_DJF)
copy_VarMeta(clm_taux3(0,:,:), clm_taux3_MAM)
copy_VarMeta(clm_taux3(0,:,:), clm_taux3_JJA)
copy_VarMeta(clm_taux3(0,:,:), clm_taux3_SON)
;========================================================================
yrStrt  =  (/41,141,141/)
yrLast  =  (/70,170,170/)
ilev1   = 0
ilev2   =5
lonw    =210
lone    =250
varname="TAUY"

ncdir4  = "/Volumes/oly2/CESM/monthly/"
fil0001    = systemfunc("ls "+ ncdir4 + "TAUY.mon.0071-0140.20S_20N.120E_80W_1degree.nc")
fil0002    = systemfunc("ls "+ ncdir4 + "TAUY.2xCO2.mon.0071-0170.20S_20N.120E_80W_1degree.nc")
fil0003    = systemfunc("ls "+ ncdir4 + "TAUY.4xCO2.mon.0071-0170.20S_20N.120E_80W_1degree.nc")
   
f0001      = addfile(fil0001, "r")
f0002      = addfile(fil0002, "r")
f0003      = addfile(fil0003, "r")

YYYY0001= cd_calendar(f0001->time,-1)/100
iyr0001 = ind(YYYY0001.ge.yrStrt(0) .and. YYYY0001.le.yrLast(0))

YYYY0002= cd_calendar(f0002->time,-1)/100
iyr0002 = ind(YYYY0002.ge.yrStrt(1) .and. YYYY0002.le.yrLast(1))

YYYY0003= cd_calendar(f0003->time,-1)/100
iyr0003 = ind(YYYY0003.ge.yrStrt(2) .and. YYYY0003.le.yrLast(2))

tauy1  = f0001->TAUY(iyr0001,:,:)
tauy2  = f0002->TAUY(iyr0002,:,:)
tauy3  = f0003->TAUY(iyr0003,:,:)

; calculate clm 
clm_tauy1 = clmMonTLL(tauy1) ; 12month,lev,lat,lon
clm_tauy2 = clmMonTLL(tauy2)
clm_tauy3 = clmMonTLL(tauy3)

; month_to_season
clm_tauy1_DJF = (clm_tauy1(11,:,:)+clm_tauy1(0,:,:)+clm_tauy1(1,:,:))/3.
clm_tauy1_MAM = (clm_tauy1(2,:,:) +clm_tauy1(3,:,:)+clm_tauy1(4,:,:))/3.
clm_tauy1_JJA = (clm_tauy1(5,:,:) +clm_tauy1(6,:,:)+clm_tauy1(7,:,:))/3.
clm_tauy1_SON = (clm_tauy1(8,:,:) +clm_tauy1(9,:,:)+clm_tauy1(10,:,:))/3.

clm_tauy2_DJF = (clm_tauy2(11,:,:)+clm_tauy2(0,:,:)+clm_tauy2(1,:,:))/3.
clm_tauy2_MAM = (clm_tauy2(2,:,:) +clm_tauy2(3,:,:)+clm_tauy2(4,:,:))/3.
clm_tauy2_JJA = (clm_tauy2(5,:,:) +clm_tauy2(6,:,:)+clm_tauy2(7,:,:))/3.
clm_tauy2_SON = (clm_tauy2(8,:,:) +clm_tauy2(9,:,:)+clm_tauy2(10,:,:))/3.

clm_tauy3_DJF = (clm_tauy3(11,:,:)+clm_tauy3(0,:,:)+clm_tauy3(1,:,:))/3.
clm_tauy3_MAM = (clm_tauy3(2,:,:) +clm_tauy3(3,:,:)+clm_tauy3(4,:,:))/3.
clm_tauy3_JJA = (clm_tauy3(5,:,:) +clm_tauy3(6,:,:)+clm_tauy3(7,:,:))/3.
clm_tauy3_SON = (clm_tauy3(8,:,:) +clm_tauy3(9,:,:)+clm_tauy3(10,:,:))/3.

;==========================================================================================
leftstring1  = (/"~F22~a",  "~F22~b",  "~F22~c", "~F22~d"/) 
leftstring2  = (/"~F22~e",  "~F22~f",  "~F22~g", "~F22~h"/) 
leftstring3  = (/"~F22~i",  "~F22~j",  "~F22~k", "~F22~l"/) 
leftstring4  = (/"~F22~m",  "~F22~n",  "~F22~o", "~F22~p"/) 


Centerstring1  = (/"MAM SST & wind","JJA SST & wind","SON SST & wind","DJF SST & wind"/) 
Centerstring2  = (/"MAM uvel",        "JJA uvel",        "SON uvel",        "DJF uvel"/) 
Centerstring3  = (/"MAM SST & wind","JJA SST & wind","SON SST & wind","DJF SST & wind"/) 
Centerstring4  = (/"MAM uvel",        "JJA uvel",        "SON uvel",        "DJF uvel"/) 

text=(/"PD","4xCO2"/)
outdir = "/Users/xueaoyun/Documents/work/TIW_future/revision1/Extend_figure/"
wks = gsn_open_wks("pdf",outdir+"CESM_PD_2xCO2_4xCO2_SST_gradient_zonal_current_seasonal_cycle")
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
res@tmXBLabels            = (/"","","","",""/)
res@tmXBMinorValues       = ispan(120,280,20)
res@tmYLMode              = "Explicit"
res@tmYLValues            = (/-10,-5,0,5,10/)
res@tmYLLabels            = (/"10~S~o~N~S","5~S~o~N~S","0","5~S~o~N~N","10~S~o~N~N"/)
res@tmYLMinorValues       = (/-5,5/)
res@tmXBLabelFontHeightF     = 0.008
res@tmYLLabelFontHeightF     = 0.008
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
res@lbLabelFontHeightF   = 0.008
res@pmLabelBarWidthF     = 0.15
res@pmLabelBarHeightF    = 0.03
res@pmLabelBarOrthogonalPosF = 0.15
res@pmLabelBarParallelPosF   = 0.5
res@gsnRightString       = ""
res@gsnLeftStringFontHeightF = 0.015
res@gsnCenterStringFontHeightF = 0.012

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
vcres@vcRefAnnoOrthogonalPosF = -1.20            ; move ref vector
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
sres@cnLevels             = (/-0.5/)
sres@cnLineThicknessF     = 2
; sres@cnLineLabelDensityF=-0.2
; sres@cnLineLabelInterval  = 1
sres@cnLineLabelsOn       = True
sres@cnLineLabelFontColor = "black"
sres@cnLineLabelFontHeightF = 0.018
sres@cnLineLabelFontThicknessF = 2.0
sres@cnLineLabelPlacementMode = "constant"
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

;=============================================================
; plot
;=============================================================
plot =new(4, graphic)
dumb =plot
vector = plot
;res@cnFillPalette="GMT_wysiwygcont"
;res@cnFillPalette="BlueWhiteOrangeRed"
res@cnFillPalette="NCV_blu_red"
res@gsnSpreadColorStart=15
res@gsnSpreadColorEnd  =245
res@vpHeightF            = 0.1
res@vpWidthF             = 0.2

; res@cnFillColors       = (/20,40,60,80,110,170,190,210,220,230/)
; res@cnLevels             = fspan(-0.6,0.6,31)
;--------------->
;--------------->
res@vpXF           = 0.10
res@vpYF           = 0.95
res@lbLabelBarOn   = False
res@gsnLeftString  = leftstring1(0)
res@gsnCenterString  = Centerstring1(0)

res@cnLevels       = fspan(22,32,21)
res@lbLabelStride  = 4
vcres@vcRefMagnitudeF           = 0.5 
plot(0) = gsn_csm_contour_map_ce(wks,clm_temp1_MAM_xy(:,:),res)
vector(0) = gsn_csm_vector(wks, clm_taux1_MAM(:,:), clm_tauy1_MAM(:,:),vcres)
overlay(plot(0), vector(0))
yspace = 0.03
vcres@vcRefAnnoOn          = False 
res@vpYF           = res@vpYF - res@vpHeightF -yspace
res@gsnLeftString  = leftstring1(1)
res@gsnCenterString  = Centerstring1(1)

plot(1) = gsn_csm_contour_map_ce(wks,clm_temp1_JJA_xy(:,:),res)
vector(1) = gsn_csm_vector(wks, clm_taux1_JJA(:,:), clm_tauy1_JJA(:,:),vcres)
overlay(plot(1), vector(1))

res@vpYF           = res@vpYF - res@vpHeightF -yspace
res@gsnLeftString  = leftstring1(2)
res@gsnCenterString  = Centerstring1(2)

plot(2) = gsn_csm_contour_map_ce(wks,clm_temp1_SON_xy(:,:),res)
vector(2) = gsn_csm_vector(wks, clm_taux1_SON(:,:), clm_tauy1_SON(:,:),vcres)
overlay(plot(2), vector(2))

res@vpYF           = res@vpYF - res@vpHeightF -yspace
res@gsnLeftString  = leftstring1(3)
res@gsnCenterString  = Centerstring1(3)

res@tmXBLabels            = (/"","160~S~o~N~E","160~S~o~N~W","120~S~o~N~W","80~S~o~N~W"/)
res@lbLabelBarOn   =  True 
plot(3) = gsn_csm_contour_map_ce(wks,clm_temp1_DJF_xy(:,:),res)
vector(3) = gsn_csm_vector(wks, clm_taux1_DJF(:,:), clm_tauy1_DJF(:,:),vcres)
overlay(plot(3), vector(3))

draw(plot(0))
draw(plot(1))
draw(plot(2))
draw(plot(3))

txres               = True         ; Text resources
txres@txFontHeightF = 0.015
txres@txFontColor  = "black"
txres@txFontThicknessF = 4.0
txres@txFont    = 22
text1 = (/ "Present-Day", "4xCO2"/)

  gsn_text_ndc(wks, text1(0), 0.25+0.06, 0.99, txres)
  gsn_text_ndc(wks, text1(1), 0.75+0.04, 0.99, txres)
;------------------------------------------------------------
plot1=plot
dumb1=plot
vector1 = vector
res@vpXF           = 0.57
res@vpYF           = 0.95
res@lbLabelBarOn   = False
res@gsnLeftString  = leftstring3(0)
res@gsnCenterString  = Centerstring3(0)
res@cnLevels       = fspan(22,32,21)+4.
vcres@vcRefMagnitudeF           = 0.5 
res@tmXBLabels            = (/"","","","",""/)
vcres@vcRefAnnoOn          = True
plot1(0) = gsn_csm_contour_map_ce(wks,clm_temp3_MAM_xy(:,:),res)
vector1(0) = gsn_csm_vector(wks, clm_taux3_MAM(:,:), clm_tauy3_MAM(:,:),vcres)
overlay(plot1(0), vector1(0))
yspace = 0.03
res@vpYF           = res@vpYF - res@vpHeightF -yspace
res@gsnLeftString  = leftstring3(1)
res@gsnCenterString  = Centerstring3(1)

vcres@vcRefAnnoOn          = False
plot1(1) = gsn_csm_contour_map_ce(wks,clm_temp3_JJA_xy(:,:),res)
vector1(1) = gsn_csm_vector(wks, clm_taux3_JJA(:,:), clm_tauy3_JJA(:,:),vcres)
overlay(plot1(1), vector1(1))

res@vpYF           = res@vpYF - res@vpHeightF -yspace
res@gsnLeftString  = leftstring3(2)
res@gsnCenterString  = Centerstring3(2)

plot1(2) = gsn_csm_contour_map_ce(wks,clm_temp3_SON_xy(:,:),res)
vector1(2) = gsn_csm_vector(wks, clm_taux3_SON(:,:), clm_tauy3_SON(:,:),vcres)
overlay(plot1(2), vector1(2))

res@vpYF           = res@vpYF - res@vpHeightF -yspace
res@gsnLeftString  = leftstring3(3)
res@gsnCenterString  = Centerstring3(3)

res@lbLabelBarOn   = True 
res@tmXBLabels            = (/"","160~S~o~N~E","160~S~o~N~W","120~S~o~N~W","80~S~o~N~W"/)
plot1(3) = gsn_csm_contour_map_ce(wks,clm_temp3_DJF_xy(:,:),res)
vector1(3) = gsn_csm_vector(wks, clm_taux3_DJF(:,:), clm_tauy3_DJF(:,:),vcres)
overlay(plot1(3), vector1(3))

draw(plot1(0))
draw(plot1(1))
draw(plot1(2))
draw(plot1(3))
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
ares@tmXBLabels = (/"","","","",""/)

ares@tmXBLabelFontHeightF  = 0.008
ares@tmYLLabelFontHeightF  = 0.008
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
;ares@cnFillPalette="BlueWhiteOrangeRed"
ares@cnFillPalette="NCV_blu_red"
ares@lbLabelBarOn          = True
ares@lbOrientation      = "Horizontal"
ares@lbBoxEndCapStyle     = "TriangleBothEnds"
ares@pmLabelBarWidthF     = 0.15
ares@pmLabelBarHeightF    = 0.03
ares@pmLabelBarOrthogonalPosF = 0.01
ares@pmLabelBarParallelPosF   = 0.5
ares@lbLabelFontHeightF       = 0.008

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
bres@cnLineLabelInterval  = 1
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
ares@vpWidthF             = 0.18
ares@vpXF                 = 0.33
ares@vpYF                 = 0.95
ares@lbLabelBarOn   = False
ares@gsnLeftString =leftstring2(0)
res@gsnCenterString  = Centerstring2(0)
ares@cnLevels             = fspan(-1.0, 1.0, 21)
bres@cnLevels             = (/-2,-1,0,0.5,1,2/)
plot2(0) = gsn_csm_contour(wks, clm_uvel1_MAM_yz(:,:),ares) 
dumb2(0) = gsn_csm_contour(wks, clm_temp1_MAM_yz(:,:),bres)
overlay(plot2(0), dumb2(0))

ares@vpHeightF            = 0.1
ares@vpWidthF             = 0.18
ares@vpXF                 = 0.33
ares@vpYF                 = ares@vpYF - ares@vpHeightF - yspace
ares@lbLabelBarOn   = False
ares@gsnLeftString =leftstring2(1)
res@gsnCenterString  = Centerstring2(1)

plot2(1) = gsn_csm_contour(wks, clm_uvel1_JJA_yz(:,:),ares) 
dumb2(1) = gsn_csm_contour(wks, clm_temp1_JJA_yz(:,:),bres)
overlay(plot2(1), dumb2(1))

ares@vpYF                 = ares@vpYF - ares@vpHeightF - yspace
ares@lbLabelBarOn   = False
ares@gsnLeftString =leftstring2(2)
res@gsnCenterString  = Centerstring2(2)

plot2(2) = gsn_csm_contour(wks, clm_uvel1_SON_yz(:,:),ares) 
dumb2(2) = gsn_csm_contour(wks, clm_temp1_SON_yz(:,:),bres)
overlay(plot2(2), dumb2(2))

ares@vpYF                 = ares@vpYF - ares@vpHeightF - yspace
ares@lbLabelBarOn   = True 
ares@gsnLeftString =leftstring2(3)
res@gsnCenterString  = Centerstring2(3)

ares@tmXBLabels = (/"10~S~o~N~S","5~S~o~N~S","Eq","5~S~o~N~N","10~S~o~N~N"/)  
plot2(3) = gsn_csm_contour(wks, clm_uvel1_DJF_yz(:,:),ares) 
dumb2(3) = gsn_csm_contour(wks, clm_temp1_DJF_yz(:,:),bres)
overlay(plot2(3), dumb2(3))

w1=1
e1=3

s1= 5
n1= 10000 

plres=True
plres@gsLineThicknessF= 3
plres@gsLineDashPattern= 2
plres@gsLineColor="blue"
plres@gsLineOpacityF=100

ddum1=gsn_add_polyline(wks, plot2(0), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum2=gsn_add_polyline(wks, plot2(1), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum3=gsn_add_polyline(wks, plot2(2), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum4=gsn_add_polyline(wks, plot2(3), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)

w2=5
e2=7

s2= 5
n2= 10000 

plres@gsLineThicknessF= 3
plres@gsLineDashPattern= 2
plres@gsLineColor="red"
plres@gsLineOpacityF=100

ddum01=gsn_add_polyline(wks, plot2(0), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum02=gsn_add_polyline(wks, plot2(1), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum03=gsn_add_polyline(wks, plot2(2), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum04=gsn_add_polyline(wks, plot2(3), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)

draw(plot2(0))
draw(plot2(1))
draw(plot2(2))
draw(plot2(3))

;=============================================================
; plot
;=============================================================
plot3 =new(4, graphic)
dumb3 =new(4, graphic)

; ares@cnFillColors         =(/20,40,60,80,110,170,190,210,220,230/)

ares@vpHeightF            = 0.1
ares@vpWidthF             = 0.18
ares@vpXF                 = 0.80
ares@vpYF                 = 0.95
ares@lbLabelBarOn   = False
ares@gsnLeftString =leftstring4(0)
res@gsnCenterString  = Centerstring4(0)

ares@cnLevels             = fspan(-1.0, 1.0, 21)
bres@cnLevels             = (/-2,-1,0,0.5,1,2/)
ares@tmXBLabels = (/"","","","",""/)
plot3(0) = gsn_csm_contour(wks, clm_uvel3_MAM_yz(:,:),ares) 
dumb3(0) = gsn_csm_contour(wks, clm_temp3_MAM_yz(:,:),bres)
overlay(plot3(0), dumb3(0))

ares@vpHeightF            = 0.1
ares@vpWidthF             = 0.18
ares@vpXF                 = 0.80
ares@vpYF                 = ares@vpYF - ares@vpHeightF - yspace
ares@lbLabelBarOn   = False
ares@gsnLeftString =leftstring4(1)
res@gsnCenterString  = Centerstring4(1)

plot3(1) = gsn_csm_contour(wks, clm_uvel3_JJA_yz(:,:),ares) 
dumb3(1) = gsn_csm_contour(wks, clm_temp3_JJA_yz(:,:),bres)
overlay(plot3(1), dumb3(1))

ares@vpYF                 = ares@vpYF - ares@vpHeightF - yspace
ares@lbLabelBarOn   = False
ares@gsnLeftString =leftstring4(2)
res@gsnCenterString  = Centerstring4(2)

plot3(2) = gsn_csm_contour(wks, clm_uvel3_SON_yz(:,:),ares) 
dumb3(2) = gsn_csm_contour(wks, clm_temp3_SON_yz(:,:),bres)
overlay(plot3(2), dumb3(2))

ares@vpYF                 = ares@vpYF - ares@vpHeightF - yspace
ares@lbLabelBarOn   = True 
ares@gsnLeftString =leftstring4(3)
res@gsnCenterString  = Centerstring4(3)

ares@tmXBLabels = (/"10~S~o~N~S","5~S~o~N~S","Eq","5~S~o~N~N","10~S~o~N~N"/)  
plot3(3) = gsn_csm_contour(wks, clm_uvel3_DJF_yz(:,:),ares) 
dumb3(3) = gsn_csm_contour(wks, clm_temp3_DJF_yz(:,:),bres)
overlay(plot3(3), dumb3(3))

w1=1
e1=3

s1= 5
n1= 10000 

plres=True
plres@gsLineThicknessF= 3
plres@gsLineDashPattern= 2
plres@gsLineColor="blue"
plres@gsLineOpacityF=100

ddum1=gsn_add_polyline(wks, plot3(0), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum2=gsn_add_polyline(wks, plot3(1), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum3=gsn_add_polyline(wks, plot3(2), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum4=gsn_add_polyline(wks, plot3(3), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)

w2=5
e2=7

s2= 5
n2= 10000 

plres@gsLineThicknessF= 3
plres@gsLineDashPattern= 2
plres@gsLineColor="red"
plres@gsLineOpacityF=100

ddum01=gsn_add_polyline(wks, plot3(0), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum02=gsn_add_polyline(wks, plot3(1), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum03=gsn_add_polyline(wks, plot3(2), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum04=gsn_add_polyline(wks, plot3(3), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)

draw(plot3(0))
draw(plot3(1))
draw(plot3(2))
draw(plot3(3))
;===================================================
data1 = new((/2,12/), float)
data1(0,:) = temp1_nino3_ave
data1(1,:) = temp3_nino3_ave

data2 = new((/2,12/), float)
data2(0,:) = EUC1_ave
data2(1,:) = EUC3_ave

data3 = new((/2,12/), float)
data3(0,:) = SEC1_ave
data3(1,:) = SEC3_ave

data4 = new((/2,12/), float)
data4(0,:) = NECC1_ave
data4(1,:) = NECC3_ave

mres               = True 
mres@gsnFrame      = False
mres@gsnDraw       = False
mres@trXMinF       = 0.
mres@trXMaxF       = 13.
mres@trYMinF       = 20.
mres@trYMaxF       = 34.
mres@vpXF          = 0.10                  ; Move slightly to left
mres@vpYF          = 0.38
mres@vpHeightF     = 0.1    
mres@vpWidthF      = 0.18 
mres@gsnYRefLine              = 0. 
mres@gsnYRefLineColors        = (/"black"/)
mres@gsnYRefLineThicknessF    = 2.0
mres@gsnYRefLineDashPattern   = 0 
;mres@gsnXRefLine              = 11.0 
;mres@gsnXRefLineColors        = (/"black"/)
;mres@gsnXRefLineThicknessF    = 2.0
;mres@gsnXRefLineDashPattern   = 1
mres@tfPolyDrawOrder          = "PreDraw" 
mres@tmXBMode                 = "Explicit"                 
mres@tmXBValues               = ispan(0,13,1)
mres@tmXBLabels               = (/"","J","F","M","A","M","J","J","A","S","O","N","D",""/)  
mres@tmYLMode                 = "Explicit"   
mres@tmYLValues               = ispan(20,34,2)
mres@tmYLLabels               = (/"20","22","24","26","28","30","32","34"/)

mres@tmYRMode                 = "Explicit"   
mres@tmYRValues               = ispan(24,34,1)
mres@tmYRLabels               = (/"24","25","26","27","28","29","30","31","32","33","34"/) 
mres@tmYRLabelFontColor       = "red"

mres@tmXBLabelFontHeightF     = 0.001
mres@xyMarkLineModes          = (/"MarkLines","MarkLines"/)  
mres@xyMarkers                = (/16,16/)
mres@xyLineColors             = (/"blue","red"/)
mres@xyDashPatterns           = (/0/)
mres@xyLineThicknesses        = (/2.0,2.0/)
mres@xyMarkerSizes            = (/0.005,0.005/) 
mres@tmXBLabelFontHeightF     = 0.008
mres@tmYLLabelFontHeightF     = 0.008
mres@tmYRLabelFontHeightF     = 0.008
mres@tmXBLabelFont            = 21
mres@tmYLLabelFont            = 21
mres@tmYRLabelFont            = 21

mres@tmXBLabelFontThicknessF  = 2.0
mres@tmYLLabelFontThicknessF  = 2.0 
mres@tmYRLabelFontThicknessF  = 2.0
mres@tmBorderThicknessF       = 2.5
mres@tmXBMajorLengthF         = 0.005
mres@tmYLMajorLengthF         = 0.005
mres@tmXBMinorLengthF         = 0.006
mres@tmYLMinorLengthF         = 0.006 
mres@tmXBMajorThicknessF      = 2.0
mres@tmYLMajorThicknessF      = 2.0
mres@tmXBMinorThicknessF      = 1.5
mres@tmYLMinorThicknessF      = 1.5 
mres@tmXBLabelDeltaF          = -0.2
mres@tmYLLabelDeltaF          = -0.2 
mres@tmYRLabelDeltaF          = -0.2
mres@tiYAxisString            = " "
mres@tiXAxisString            = " "
mres@tiMainString             = " "
mres@tmXBLabelAngleF          = 0
mres@pmLegendDisplayMode    = "Always"            ; turn on legend
mres@pmLegendSide           = "Top"               ; Change location of 
mres@pmLegendZone           = 0        
mres@pmLegendParallelPosF   = -0.20                ; move units right
mres@pmLegendOrthogonalPosF = -0.35                ; move units down
mres@pmLegendWidthF         = 0.06                ; Change width and
mres@pmLegendHeightF        = 0.04                ; height of legend.
mres@lgPerimOn              = False               ; turn off box around
mres@lgLabelFontHeightF     = 0.008              ; label font height
mres@lgLabelFont            = 21                  ; label font height
mres@lgLabelFontThicknessF  = 4                   ; label font height
mres@xyExplicitLegendLabels = (/"PD","4xCO2"/)  
mres@tmXTOn                      = False 
mres@tmYROn                      = False
mres@tmXMajorGrid                = False 
mres@tmXMajorGridLineDashPattern =  2
mres@tmXMajorGridLineColor       = "grey40"
mres@tmXMajorGridThicknessF      = 0.8
mres@tmYMajorGrid                = True
mres@tmYMajorGridLineDashPattern =  2
mres@tmYMajorGridLineColor       = "grey40"
mres@tmYMajorGridThicknessF      = 0.8 
mres@tmYMinorGrid                = False  
mres@tmYMinorGridLineDashPattern =  2
mres@tmYMinorGridLineColor       = "grey40"
mres@tmYMinorGridThicknessF      = 0.3 
mres@tmXMinorGrid                = False 
mres@tmXMinorGridLineDashPattern =  2
mres@tmXMinorGridLineColor       = "grey40"
mres@tmXMinorGridThicknessF      = 0.3 
; mres@tiMainString                ="northward jump of WNPSH"
mres@tiMainFontHeightF =0.02 
mres@gsnLeftString      ="~F22~q"
mres@gsnLeftStringFontHeightF =0.015
mres@gsnCenterString   = "Cold tongue SST"
mres@gsnCenterStringFontHeightF = 0.012
x = ispan(1,12,1)
xspace  = 0.055
plot00 = gsn_csm_xy(wks,x,data1,mres)
delete(mres@tmYLValues)
delete(mres@tmYLLabels)
mres@pmLegendParallelPosF   = -0.20                ; move units right
mres@trYMinF       = 0.
mres@trYMaxF       = 1.
mres@vpXF          = mres@vpXF+mres@vpWidthF+xspace                  ; Move slightly to left
mres@vpYF          = mres@vpYF
mres@vpHeightF     = 0.1    
mres@vpWidthF      = 0.18 
mres@tmYLValues               = (/0,0.2,0.4,0.6,0.8,1.0/)
mres@tmYLLabels               = (/"0","0.2","0.4","0.6","0.8","1"/)  
mres@gsnLeftString      ="~F22~r "
mres@gsnCenterString    = "EUC"
mres@xyExplicitLegendLabels = (/"PD","4xCO2"/)  
plot01 = gsn_csm_xy(wks,x,data2,mres)
delete(mres@tmYLValues)
delete(mres@tmYLLabels)
mres@pmLegendParallelPosF   =-0.20               ; move units right
mres@trYMinF       = -1.
mres@trYMaxF       = 1.
mres@vpXF          = mres@vpXF+mres@vpWidthF+xspace                  ; Move slightly to left
mres@vpYF          = mres@vpYF
mres@vpHeightF     = 0.1    
mres@vpWidthF      = 0.18 
mres@tmYLValues               = (/-1.0,-0.5,0,0.5,1.0/)
mres@tmYLLabels               = (/"-1","-0.5","0","0.5","1"/)  
mres@gsnLeftString      ="~F22~s "
mres@gsnCenterString    = " nSEC"
mres@xyExplicitLegendLabels = (/"PD","4xCO2"/)  
plot02 = gsn_csm_xy(wks,x,data3,mres)

delete(mres@tmYLValues)
delete(mres@tmYLLabels)
mres@pmLegendParallelPosF   = -0.20                ; move units right
mres@trYMinF       = -1.
mres@trYMaxF       = 1.
mres@vpXF          = mres@vpXF+mres@vpWidthF+xspace                  ; Move slightly to left
mres@vpYF          = mres@vpYF
mres@vpHeightF     = 0.1    
mres@vpWidthF      = 0.18
mres@tmYLValues               = (/-1.0,-0.5,0,0.5,1.0/)
mres@tmYLLabels               = (/"-1","-0.5","0","0.5","1"/)  
mres@gsnLeftString      ="~F22~t "
mres@gsnCenterString    = " NECC"
mres@xyExplicitLegendLabels = (/"PD","4xCO2"/)  
plot03 = gsn_csm_xy(wks,x,data4,mres)


 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "royalblue1"                 ; color chosen
 dummy1 = gsn_add_polygon (wks,plot00,xp1_temp_PD,yp1_temp_PD,gsres)

 gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "pink1"                 ; color chosen
 gsres@gsFillOpacityF     = 0.5
 dummy2 = gsn_add_polygon (wks,plot00,xp1_temp_4xCO2,yp1_temp_4xCO2,gsres)

 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "royalblue1"                 ; color chosen
 dummy01 = gsn_add_polygon (wks,plot01,xp1_EUC_PD,yp1_EUC_PD,gsres)

 gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "pink1"                 ; color chosen
 gsres@gsFillOpacityF     = 0.5
 dummy02 = gsn_add_polygon (wks,plot01,xp1_EUC_4xCO2,yp1_EUC_4xCO2,gsres)

 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "royalblue1"                 ; color chosen
 dummy001 = gsn_add_polygon (wks,plot02,xp1_SEC_PD,yp1_SEC_PD,gsres)

 gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "pink1"                 ; color chosen
 gsres@gsFillOpacityF     = 0.5
 dummy002 = gsn_add_polygon (wks,plot02,xp1_SEC_4xCO2,yp1_SEC_4xCO2,gsres)

 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "royalblue1"                 ; color chosen
 dummy0001 = gsn_add_polygon (wks,plot03,xp1_NECC_PD,yp1_NECC_PD,gsres)

 gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "pink1"                 ; color chosen
 gsres@gsFillOpacityF     = 0.5
 dummy0002 = gsn_add_polygon (wks,plot03,xp1_NECC_4xCO2,yp1_NECC_4xCO2,gsres)


draw(plot00)
draw(plot01)
draw(plot02)
draw(plot03)

end 