load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
load "/Users/xueaoyun/Documents/work/libfun/my_utilities.ncl"

begin
;=======================================================================
  datadir1 = "/Volumes/oly2/CESM/energy/mon/"
  fil1    = systemfunc("ls " + datadir1 + "TIW_CKE_01-06.PD.mon.nc") 
  f1      = addfile(fil1, "r")
  CKE_PD  = f1->CKE_PD
  CKE_PD_mon = dim_avg_n_Wrap(CKE_PD(:,0:5,:,:), 1)
  datadir2 = "/Volumes/oly2/CESM/energy/mon/"
  fil2    = systemfunc("ls " + datadir2 + "TIW_CKE_01-06.2xCO2.mon.nc") 
  f2      = addfile(fil2, "r")
  CKE_2xCO2  = f2->CKE_2xCO2
  CKE_2xCO2_mon = dim_avg_n_Wrap(CKE_2xCO2(:,0:5,:,:), 1)
  datadir3 = "/Volumes/oly2/CESM/energy/mon/"
  fil3    = systemfunc("ls " + datadir3 + "TIW_CKE_01-06.4xCO2.mon.nc") 
  f3      = addfile(fil3, "r")
  CKE_4xCO2  =f3->CKE_4xCO2
  CKE_4xCO2_mon = dim_avg_n_Wrap(CKE_4xCO2(:,0:5,:,:), 1)

  CKE_PD_clm    = clmMonTLL(CKE_PD_mon)  
  CKE_2xCO2_clm = clmMonTLL(CKE_2xCO2_mon)
  CKE_4xCO2_clm = clmMonTLL(CKE_4xCO2_mon)

ntim1  = dimsizes(CKE_PD&time)
ntim2  = dimsizes(CKE_2xCO2&time)
ntim3  = dimsizes(CKE_4xCO2&time)
nyears1 = ntim1/12
nyears2 = ntim2/12
nyears3 = ntim3/12
;-------------------------------------------------------------------------
  ; calculate domain-avg BTR NEQ
CKE1 = wgt_areaave_Wrap(CKE_PD_mon(:,{3:5},{200:250}), 1.0, 1.0, 0)
CKE2 = wgt_areaave_Wrap(CKE_2xCO2_mon(:,{3:5},{200:250}), 1.0, 1.0, 0)
CKE3 = wgt_areaave_Wrap(CKE_4xCO2_mon(:,{3:5},{200:250}), 1.0, 1.0, 0)

CKE1_ann  = reshape(CKE1,(/nyears1,12/))
CKE2_ann  = reshape(CKE2,(/nyears2,12/))
CKE3_ann  = reshape(CKE3,(/nyears3,12/))

; annual cycle
CKE1_ave = dim_avg_n_Wrap(CKE1_ann, 0)
CKE2_ave = dim_avg_n_Wrap(CKE2_ann, 0)
CKE3_ave = dim_avg_n_Wrap(CKE3_ann, 0)

; std 
CKE1_std = dim_stddev_n_Wrap(CKE1_ann, 0) 
CKE2_std = dim_stddev_n_Wrap(CKE2_ann, 0) 
CKE3_std = dim_stddev_n_Wrap(CKE3_ann, 0) 

xp1_CKE_PD = new((/24/), double)
yp1_CKE_PD = new((/24/), double)

do i=0,11
   xp1_CKE_PD(i)          = i+1
   xp1_CKE_PD(2*12-1-i)   = i+1
   yp1_CKE_PD(i)          = CKE1_ave(i) + CKE1_std(i)/2.
   yp1_CKE_PD(2*12-1-i)   = CKE1_ave(i) - CKE1_std(i)/2. 
 end do

 xp1_CKE_4xCO2 = new((/24/), double)
 yp1_CKE_4xCO2 = new((/24/), double)
 
 do i=0,11
    xp1_CKE_4xCO2(i)          = i+1
    xp1_CKE_4xCO2(2*12-1-i)   = i+1
    yp1_CKE_4xCO2(i)          = CKE3_ave(i) + CKE3_std(i)/2.
    yp1_CKE_4xCO2(2*12-1-i)   = CKE3_ave(i) - CKE3_std(i)/2. 
  end do
;===================================================================
 ; calculate domain-avg BTR EQ
 CKE01 = wgt_areaave_Wrap(CKE_PD_mon(:,{0:2},{200:250}), 1.0, 1.0, 0)
 CKE02 = wgt_areaave_Wrap(CKE_2xCO2_mon(:,{0:2},{200:250}), 1.0, 1.0, 0)
 CKE03 = wgt_areaave_Wrap(CKE_4xCO2_mon(:,{0:2},{200:250}), 1.0, 1.0, 0)
 
 CKE01_ann  = reshape(CKE01,(/nyears1,12/))
 CKE02_ann  = reshape(CKE02,(/nyears2,12/))
 CKE03_ann  = reshape(CKE03,(/nyears3,12/))
 
 ; annual cycle
 CKE01_ave = dim_avg_n_Wrap(CKE01_ann, 0)
 CKE02_ave = dim_avg_n_Wrap(CKE02_ann, 0)
 CKE03_ave = dim_avg_n_Wrap(CKE03_ann, 0)
 
 ; std 
 CKE01_std = dim_stddev_n_Wrap(CKE01_ann, 0) 
 CKE02_std = dim_stddev_n_Wrap(CKE02_ann, 0) 
 CKE03_std = dim_stddev_n_Wrap(CKE03_ann, 0) 
 
 xp2_CKE_PD = new((/24/), double)
 yp2_CKE_PD = new((/24/), double)
 
 do i=0,11
    xp2_CKE_PD(i)          = i+1
    xp2_CKE_PD(2*12-1-i)   = i+1
    yp2_CKE_PD(i)          = CKE01_ave(i) + CKE01_std(i)/2.
    yp2_CKE_PD(2*12-1-i)   = CKE01_ave(i) - CKE01_std(i)/2. 
  end do
 
  xp2_CKE_4xCO2 = new((/24/), double)
  yp2_CKE_4xCO2 = new((/24/), double)
  
  do i=0,11
     xp2_CKE_4xCO2(i)          = i+1
     xp2_CKE_4xCO2(2*12-1-i)   = i+1
     yp2_CKE_4xCO2(i)          = CKE03_ave(i) + CKE03_std(i)/2.
     yp2_CKE_4xCO2(2*12-1-i)   = CKE03_ave(i) - CKE03_std(i)/2. 
   end do
 ;===================================================================
; month_to_season
clm_CKE1_DJF = (CKE_PD_clm(11,:,:)+CKE_PD_clm(0,:,:)+CKE_PD_clm(1,:,:))/3.
clm_CKE1_MAM = (CKE_PD_clm(2,:,:) +CKE_PD_clm(3,:,:)+CKE_PD_clm(4,:,:))/3.
clm_CKE1_JJA = (CKE_PD_clm(5,:,:) +CKE_PD_clm(6,:,:)+CKE_PD_clm(7,:,:))/3.
clm_CKE1_SON = (CKE_PD_clm(8,:,:) +CKE_PD_clm(9,:,:)+CKE_PD_clm(10,:,:))/3.

;clm_CKE2_DJF = (CKE_2xCO2_clm(11,:,:)+CKE_2xCO2_clm(0,:,:)+CKE_2xCO2_clm(1,:,:))/3.
;clm_CKE2_MAM = (CKE_2xCO2_clm(2,:,:) +CKE_2xCO2_clm(3,:,:)+CKE_2xCO2_clm(4,:,:))/3.
;clm_CKE2_JJA = (CKE_2xCO2_clm(5,:,:) +CKE_2xCO2_clm(6,:,:)+CKE_2xCO2_clm(7,:,:))/3.
;clm_CKE2_SON = (CKE_2xCO2_clm(8,:,:) +CKE_2xCO2_clm(9,:,:)+CKE_2xCO2_clm(10,:,:))/3.

clm_CKE3_DJF = (CKE_4xCO2_clm(11,:,:)+CKE_4xCO2_clm(0,:,:)+CKE_4xCO2_clm(1,:,:))/3.
clm_CKE3_MAM = (CKE_4xCO2_clm(2,:,:) +CKE_4xCO2_clm(3,:,:)+CKE_4xCO2_clm(4,:,:))/3.
clm_CKE3_JJA = (CKE_4xCO2_clm(5,:,:) +CKE_4xCO2_clm(6,:,:)+CKE_4xCO2_clm(7,:,:))/3.
clm_CKE3_SON = (CKE_4xCO2_clm(8,:,:) +CKE_4xCO2_clm(9,:,:)+CKE_4xCO2_clm(10,:,:))/3.
copy_VarCoords(CKE_PD_mon(0,:,:), clm_CKE1_DJF)
copy_VarCoords(CKE_PD_mon(0,:,:), clm_CKE1_MAM)
copy_VarCoords(CKE_PD_mon(0,:,:), clm_CKE1_JJA)
copy_VarCoords(CKE_PD_mon(0,:,:), clm_CKE1_SON)

copy_VarCoords(CKE_4xCO2_mon(0,:,:), clm_CKE3_DJF)
copy_VarCoords(CKE_4xCO2_mon(0,:,:), clm_CKE3_MAM)
copy_VarCoords(CKE_4xCO2_mon(0,:,:), clm_CKE3_JJA)
copy_VarCoords(CKE_4xCO2_mon(0,:,:), clm_CKE3_SON)
;=========================================================================
  datadir1 = "/Volumes/oly2/CESM/energy/mon/"
  fil1    = systemfunc("ls " + datadir1 + "TIW_APE_01-06.PD.mon.nc") 
  f1      = addfile(fil1, "r")
  APE_PD0  =f1->APE_PD(:,:,:,:)
  APE_PD  = APE_PD0 * 10000.
  copy_VarMeta(APE_PD0, APE_PD)
  APE_PD_mon = dim_avg_n_Wrap(APE_PD(:,0:5,:,:), 1)

  datadir2 = " /Volumes/oly2/CESM/energy/mon/"
  fil2    = systemfunc("ls " + datadir2 + "TIW_APE_01-06.2xCO2.mon.nc") 
  f2      = addfile(fil2, "r")
  APE_2xCO20 = f2->APE_2xCO2
  APE_2xCO2  = APE_2xCO20 * 10000.
  copy_VarMeta(APE_2xCO20, APE_2xCO2)
  APE_2xCO2_mon = dim_avg_n_Wrap(APE_2xCO2(:,0:5,:,:), 1)

  datadir3 = "/Volumes/oly2/CESM/energy/mon/"
  fil3    = systemfunc("ls " + datadir3 + "TIW_APE_01-06.4xCO2.mon.nc") 
  f3      = addfile(fil3, "r")
  APE_4xCO20  =f3->APE_4xCO2
  APE_4xCO2  = APE_4xCO20 * 10000.
  copy_VarMeta(APE_4xCO20, APE_4xCO2)
  APE_4xCO2_mon = dim_avg_n_Wrap(APE_4xCO2(:,0:5,:,:), 1)

  ;-------------------------------------------------------------------------
  ; calculate domain-avg BCR NEQ
APE1 = wgt_areaave_Wrap(APE_PD_mon(:,{3:5},{200:250}), 1.0, 1.0, 0)
APE2 = wgt_areaave_Wrap(APE_2xCO2_mon(:,{3:5},{200:250}), 1.0, 1.0, 0)
APE3 = wgt_areaave_Wrap(APE_4xCO2_mon(:,{3:5},{200:250}), 1.0, 1.0, 0)

APE1_ann  = reshape(APE1,(/nyears1,12/))
APE2_ann  = reshape(APE2,(/nyears2,12/))
APE3_ann  = reshape(APE3,(/nyears3,12/))

; annual cycle
APE1_ave = dim_avg_n_Wrap(APE1_ann, 0)
APE2_ave = dim_avg_n_Wrap(APE2_ann, 0)
APE3_ave = dim_avg_n_Wrap(APE3_ann, 0)

; std 
APE1_std = dim_stddev_n_Wrap(APE1_ann, 0) 
APE2_std = dim_stddev_n_Wrap(APE2_ann, 0) 
APE3_std = dim_stddev_n_Wrap(APE3_ann, 0) 

xp1_APE_PD = new((/24/), double)
yp1_APE_PD = new((/24/), double)

do i=0,11
   xp1_APE_PD(i)          = i+1
   xp1_APE_PD(2*12-1-i)   = i+1
   yp1_APE_PD(i)          = APE1_ave(i) + APE1_std(i)/2.
   yp1_APE_PD(2*12-1-i)   = APE1_ave(i) - APE1_std(i)/2. 
 end do

 xp1_APE_4xCO2 = new((/24/), double)
 yp1_APE_4xCO2 = new((/24/), double)
 
 do i=0,11
    xp1_APE_4xCO2(i)          = i+1
    xp1_APE_4xCO2(2*12-1-i)   = i+1
    yp1_APE_4xCO2(i)          = APE3_ave(i) + APE3_std(i)/2.
    yp1_APE_4xCO2(2*12-1-i)   = APE3_ave(i) - APE3_std(i)/2. 
  end do
;===================================================================
 ; calculate domain-avg BCR EQ
 APE01 = wgt_areaave_Wrap(APE_PD_mon(:,{0:2},{200:250}), 1.0, 1.0, 0)
 APE02 = wgt_areaave_Wrap(APE_2xCO2_mon(:,{0:2},{200:250}), 1.0, 1.0, 0)
 APE03 = wgt_areaave_Wrap(APE_4xCO2_mon(:,{0:2},{200:250}), 1.0, 1.0, 0)
 
 APE01_ann  = reshape(APE01,(/nyears1,12/))
 APE02_ann  = reshape(APE02,(/nyears2,12/))
 APE03_ann  = reshape(APE03,(/nyears3,12/))
 
 ; annual cycle
 APE01_ave = dim_avg_n_Wrap(APE01_ann, 0)
 APE02_ave = dim_avg_n_Wrap(APE02_ann, 0)
 APE03_ave = dim_avg_n_Wrap(APE03_ann, 0)
 
 ; std 
 APE01_std = dim_stddev_n_Wrap(APE01_ann, 0) 
 APE02_std = dim_stddev_n_Wrap(APE02_ann, 0) 
 APE03_std = dim_stddev_n_Wrap(APE03_ann, 0) 
 
 xp2_APE_PD = new((/24/), double)
 yp2_APE_PD = new((/24/), double)
 
 do i=0,11
    xp2_APE_PD(i)          = i+1
    xp2_APE_PD(2*12-1-i)   = i+1
    yp2_APE_PD(i)          = APE01_ave(i) + APE01_std(i)/2.
    yp2_APE_PD(2*12-1-i)   = APE01_ave(i) - APE01_std(i)/2. 
  end do
 
  xp2_APE_4xCO2 = new((/24/), double)
  yp2_APE_4xCO2 = new((/24/), double)
  
  do i=0,11
     xp2_APE_4xCO2(i)          = i+1
     xp2_APE_4xCO2(2*12-1-i)   = i+1
     yp2_APE_4xCO2(i)          = APE03_ave(i) + APE03_std(i)/2.
     yp2_APE_4xCO2(2*12-1-i)   = APE03_ave(i) - APE03_std(i)/2. 
   end do
;--------------------------------------------------------------------
  APE_PD_clm    = clmMonTLL(APE_PD_mon)  
  APE_2xCO2_clm = clmMonTLL(APE_2xCO2_mon)
  APE_4xCO2_clm = clmMonTLL(APE_4xCO2_mon)

  ; month_to_season
  clm_APE1_DJF = (APE_PD_clm(11,:,:)+ APE_PD_clm(0,:,:) + APE_PD_clm(1,:,:))/3.
  clm_APE1_MAM = (APE_PD_clm(2,:,:) + APE_PD_clm(3,:,:) + APE_PD_clm(4,:,:))/3.
  clm_APE1_JJA = (APE_PD_clm(5,:,:) + APE_PD_clm(6,:,:) + APE_PD_clm(7,:,:))/3.
  clm_APE1_SON = (APE_PD_clm(8,:,:) + APE_PD_clm(9,:,:) + APE_PD_clm(10,:,:))/3.
  
 ; clm_APE2_DJF = (APE_2xCO2_clm(11,:,:)+APE_2xCO2_clm(0,:,:)+APE_2xCO2_clm(1,:,:))/3.
 ; clm_APE2_MAM = (APE_2xCO2_clm(2,:,:) +APE_2xCO2_clm(3,:,:)+APE_2xCO2_clm(4,:,:))/3.
 ; clm_APE2_JJA = (APE_2xCO2_clm(5,:,:) +APE_2xCO2_clm(6,:,:)+APE_2xCO2_clm(7,:,:))/3.
 ; clm_APE2_SON = (APE_2xCO2_clm(8,:,:) +APE_2xCO2_clm(9,:,:)+APE_2xCO2_clm(10,:,:))/3.
  
  clm_APE3_DJF = (APE_4xCO2_clm(11,:,:)+APE_4xCO2_clm(0,:,:)+APE_4xCO2_clm(1,:,:))/3.
  clm_APE3_MAM = (APE_4xCO2_clm(2,:,:) +APE_4xCO2_clm(3,:,:)+APE_4xCO2_clm(4,:,:))/3.
  clm_APE3_JJA = (APE_4xCO2_clm(5,:,:) +APE_4xCO2_clm(6,:,:)+APE_4xCO2_clm(7,:,:))/3.
  clm_APE3_SON = (APE_4xCO2_clm(8,:,:) +APE_4xCO2_clm(9,:,:)+APE_4xCO2_clm(10,:,:))/3.

copy_VarCoords(APE_PD_mon(0,:,:), clm_APE1_DJF)
copy_VarCoords(APE_PD_mon(0,:,:), clm_APE1_MAM)
copy_VarCoords(APE_PD_mon(0,:,:), clm_APE1_JJA)
copy_VarCoords(APE_PD_mon(0,:,:), clm_APE1_SON)

copy_VarCoords(APE_4xCO2_mon(0,:,:), clm_APE3_DJF)
copy_VarCoords(APE_4xCO2_mon(0,:,:), clm_APE3_MAM)
copy_VarCoords(APE_4xCO2_mon(0,:,:), clm_APE3_JJA)
copy_VarCoords(APE_4xCO2_mon(0,:,:), clm_APE3_SON)
;--------------------------------------------------------------------------------------
leftstring1  = (/"~F22~a",  "~F22~b",  "~F22~c", "~F22~d"/) 
leftstring2  = (/"~F22~e",  "~F22~f",  "~F22~g", "~F22~h"/) 
leftstring3  = (/"~F22~i",  "~F22~j",  "~F22~k", "~F22~l"/) 
leftstring4  = (/"~F22~m",  "~F22~n",  "~F22~o", "~F22~p"/) 


Centerstring1  = (/"PD MAM  ","PD JJA","PD SON","PD DJF"/) 
Centerstring2  = (/"4xCO2 MAM  ","4xCO2 JJA","4xCO2 SON","4xCO2 DJF"/) 
Centerstring3  = (/"PD MAM  ","PD JJA","PD SON","PD DJF"/) 
Centerstring4  = (/"4xCO2 MAM  ","4xCO2 JJA","4xCO2 SON","4xCO2 DJF"/)  
  ;=====================================================================================
pltDir =  "/Users/xueaoyun/Documents/work/TIW_future/revision1/Extend_figure/"
   pltName = "Extended_Figure8_CESM_seasonal_CKE_APE"
   pltPath= pltDir+"/"+pltName
   pltTypeLocal = "pdf"
  wks  = gsn_open_wks(pltTypeLocal, pltPath)  
 ;gsn_define_colormap(wks,"sunshine_9lev") 
 gsn_define_colormap(wks,"BlueWhiteOrangeRed") 
 plot =new(4, graphic)
res               = True                      ; plot mods desired
res@gsnDraw      =  False                   ; do not draw the plot
res@gsnFrame     =  False 
res@gsnSpreadColors  = True                      ; use full colormap
res@gsnSpreadColorStart = 140
res@gsnSpreadColorEnd = 240
res@gsnAddCyclic         = False        ; data not cyclic
;res@cnMissingValFillColor = "black"
;res@gsnScalarContour        = True ;important!
res@cnLevelSelectionMode  = "ExplicitLevels"
res@cnLevels              = fspan(1.,40.,40)
;print(res@cnLevels)
;res@cnFillColors         = (/2,11,30,50,70,90,110,130,150,170,190,230/)
;res@cnFillColors          = (/0,2,3,4,5,6,7,8,9,10,11/)
res@cnFillOn         = True                     ; turn on color
res@cnLinesOn        = False 
res@cnLineLabelsOn   = False                     ; no contour line labels
res@cnInfoLabelOn    = False                     ; no info label

res@mpMinLatF     = -5
res@mpMaxLatF     = 10
res@mpMinLonF     = 180
res@mpMaxLonF     = 280
res@mpShapeMode = "FreeAspect"
res@mpCenterLonF  =230.
res@tmXBMinorOn = False
res@tmYLMinorOn = False
res@tmYLTickSpacingF = 2
res@tmXBTickSpacingF = 5
res@tmBorderThicknessF = 2
res@tmXBMajorThicknessF= 2
res@tmYLMajorThicknessF= 2
res@tmXBLabelFontHeightF = 0.006
res@tmYLLabelFontHeightF = 0.006
res@tmXBMajorLengthF = -0.005
res@tmYLMajorLengthF = -0.005
res@tmXBLabelDeltaF = -0.6
res@tmYLLabelDeltaF = -0.8
;res@trXMinF= 142

res@tmYLMode        = "Explicit"  
res@tmYLValues      = (/-10,-8,-6,-4,-2,0,2,4,6,8,10/)
res@tmYLLabels      = (/"10~S~o~N~S","8~S~o~N~S","6~S~o~N~S","4~S~o~N~S","2~S~o~N~S","0","2~S~o~N~N","4~S~o~N~N","6~S~o~N~N","8~S~o~N~N","10~S~o~N~N"/)

res@tmXBMode        = "Explicit"  
res@tmXBValues            = (/180,200,220,240,260,280/)
;res@tmXBLabels            = (/"180","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W","80~S~o~N~W"/)
res@tmXBLabels            = (/"","","","","",""/)
res@gsnLeftString  =""
res@gsnRightString =""
res@tiYAxisFontHeightF = 0.010
res@tiXAxisFontHeightF = 0.010
res@gsnLeftStringFontHeightF = 0.012

res@lbLabelBarOn =False
res@lbOrientation        = "Horizontal"
res@lbBoxEndCapStyle     = "TriangleBothEnds"
res@lbLabelFontHeightF   = 0.008
res@pmLabelBarWidthF     = 0.18
res@pmLabelBarHeightF    = 0.03
res@pmLabelBarOrthogonalPosF = 0.20
res@pmLabelBarParallelPosF   = 0.5
;res@lbTitleString         = "W/m~S~3~N"
;res@lbTitleFontHeightF = 0.015
res@lbBoxLinesOn    = False
res@lbBoxSeparatorLinesOn= False  
res@vpHeightF            = 0.1      ; Changes the aspect ratio
res@vpWidthF             = 0.18
res@vpXF    = 0.05      ; Make plots wider
res@vpYF    = 0.95
res@tmXBLabelsOn = True
res@gsnLeftString  = leftstring1(0)
res@gsnCenterString  = Centerstring1(0)
plot(0) = gsn_csm_contour_map(wks,clm_CKE1_MAM,res)
;****************************************************
yspace = 0.03
res@vpHeightF            = 0.1      ; Changes the aspect ratio
res@vpWidthF             = 0.18
res@vpXF    = 0.05      ; Make plots wider
res@vpYF    = res@vpYF - res@vpHeightF -yspace
res@gsnLeftString  = leftstring1(1)
res@gsnCenterString  = Centerstring1(1)
plot(1) = gsn_csm_contour_map(wks,clm_CKE1_JJA,res)
;*****************************************************
res@vpHeightF            = 0.1      ; Changes the aspect ratio
res@vpWidthF             = 0.18
res@vpXF    = 0.05      ; Make plots wider
res@vpYF    = res@vpYF - res@vpHeightF -yspace
res@gsnLeftString  = leftstring1(2)
res@gsnCenterString  = Centerstring1(2)
plot(2) = gsn_csm_contour_map(wks,clm_CKE1_SON,res)
;*****************************************************
res@vpHeightF            = 0.1      ; Changes the aspect ratio
res@vpWidthF             = 0.18
res@vpXF    = 0.05      ; Make plots wider
res@vpYF    = res@vpYF - res@vpHeightF -yspace
res@gsnLeftString  = leftstring1(3)
res@gsnCenterString  = Centerstring1(3)
res@tmXBLabels            = (/"180","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W","80~S~o~N~W"/)
res@lbLabelBarOn   =  True 
plot(3) = gsn_csm_contour_map(wks,clm_CKE1_DJF,res)


txres               = True         ; Text resources
txres@txFontHeightF = 0.015
txres@txFontColor  = "black"
txres@txFontThicknessF = 4.0
txres@txFont    = 22
text1 = (/ "BTR", "BCR"/)

  gsn_text_ndc(wks, text1(0), 0.25, 0.99, txres)
  gsn_text_ndc(wks, text1(1), 0.70, 0.99, txres)
;================================================================
plot1=plot
res@vpXF           = 0.27
res@vpYF           = 0.95
res@lbLabelBarOn   = False
res@gsnLeftString  = leftstring2(0)
res@gsnCenterString  = Centerstring2(0)
res@tmXBLabels            = (/"","","","","",""/)
res@tmYLLabelsOn = True
res@tmYLMode        = "Explicit"  
res@tmYLValues      = (/-10,-8,-6,-4,-2,0,2,4,6,8,10/)
res@tmYLLabels      = (/"","","","","","","","","","",""/)
res@lbLabelBarOn =False
plot1(0) = gsn_csm_contour_map(wks,clm_CKE3_MAM,res)
;****************************************************
yspace = 0.03
res@vpHeightF            = 0.1      ; Changes the aspect ratio
res@vpWidthF             = 0.18
res@vpXF    = 0.27                  ; Make plots wider
res@vpYF    = res@vpYF - res@vpHeightF -yspace
res@gsnLeftString  = leftstring2(1)
res@gsnCenterString  = Centerstring2(1)
plot1(1) = gsn_csm_contour_map(wks,clm_CKE3_JJA,res)
;*****************************************************
res@vpHeightF            = 0.1      ; Changes the aspect ratio
res@vpWidthF             = 0.18
res@vpXF    = 0.27      ; Make plots wider
res@vpYF    = res@vpYF - res@vpHeightF -yspace
res@gsnLeftString  = leftstring2(2)
res@gsnCenterString  = Centerstring2(2)
plot1(2) = gsn_csm_contour_map(wks,clm_CKE3_SON,res)
;*****************************************************
res@vpHeightF            = 0.1      ; Changes the aspect ratio
res@vpWidthF             = 0.18
res@vpXF    = 0.27      ; Make plots wider
res@vpYF    = res@vpYF - res@vpHeightF -yspace
res@gsnLeftString  = leftstring2(3)
res@gsnCenterString  = Centerstring2(3)
res@tmXBLabels            = (/"180","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W","80~S~o~N~W"/)
res@lbLabelBarOn   =  True 
plot1(3) = gsn_csm_contour_map(wks,clm_CKE3_DJF,res)
;-------------------------------------------------------------------------  
delete(res@cnLevels) 
plot2=plot
res@vpXF           = 0.50
res@vpYF           = 0.95
res@lbLabelBarOn   = False
res@gsnLeftString  = leftstring3(0)
res@gsnCenterString  = Centerstring3(0)
res@cnLevels              = fspan(1.,30.,30)
res@tmXBLabels            =  (/"","","","","",""/)
res@tmYLLabelsOn = True
res@tmYLMode        = "Explicit"  
res@tmYLValues      = (/-10,-8,-6,-4,-2,0,2,4,6,8,10/)
res@tmYLLabels      = (/"","","","","","","","","","",""/)
res@lbLabelBarOn =False
plot2(0) = gsn_csm_contour_map(wks,clm_APE1_MAM,res)
;****************************************************
yspace = 0.03
res@vpHeightF            = 0.1      ; Changes the aspect ratio
res@vpWidthF             = 0.18
res@vpXF    = 0.50                  ; Make plots wider
res@vpYF    = res@vpYF - res@vpHeightF -yspace
res@gsnLeftString  = leftstring3(1)
res@gsnCenterString  = Centerstring3(1)
plot2(1) = gsn_csm_contour_map(wks,clm_APE1_JJA,res)
;*****************************************************
res@vpHeightF            = 0.1      ; Changes the aspect ratio
res@vpWidthF             = 0.18
res@vpXF    = 0.50      ; Make plots wider
res@vpYF    = res@vpYF - res@vpHeightF -yspace
res@gsnLeftString  = leftstring3(2)
res@gsnCenterString  = Centerstring3(2)
plot2(2) = gsn_csm_contour_map(wks,clm_APE1_SON,res)
;*****************************************************
res@vpHeightF            = 0.1      ; Changes the aspect ratio
res@vpWidthF             = 0.18
res@vpXF    = 0.50      ; Make plots wider
res@vpYF    = res@vpYF - res@vpHeightF -yspace
res@gsnLeftString  = leftstring3(3)
res@gsnCenterString  = Centerstring3(3)
res@tmXBLabels            = (/"180","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W","80~S~o~N~W"/)
res@lbLabelBarOn   =  True 
plot2(3) = gsn_csm_contour_map(wks,clm_APE1_DJF,res)

;==============================================================
plot3=plot
res@vpXF           = 0.72
res@vpYF           = 0.95
res@lbLabelBarOn   = False
res@gsnLeftString  = leftstring4(0)
res@gsnCenterString  = Centerstring4(0)

res@tmXBLabels            =  (/"","","","","",""/)
res@tmYLLabelsOn = True
res@tmYLMode        = "Explicit"  
res@tmYLValues      = (/-10,-8,-6,-4,-2,0,2,4,6,8,10/)
res@tmYLLabels      = (/"","","","","","","","","","",""/)
res@lbLabelBarOn    = False
plot3(0) = gsn_csm_contour_map(wks,clm_APE3_MAM,res)
;****************************************************
yspace = 0.03
res@vpHeightF            = 0.1      ; Changes the aspect ratio
res@vpWidthF             = 0.18
res@vpXF    = 0.72                 ; Make plots wider
res@vpYF    = res@vpYF - res@vpHeightF -yspace
res@gsnLeftString  = leftstring4(1)
res@gsnCenterString  = Centerstring4(1)
plot3(1) = gsn_csm_contour_map(wks,clm_APE3_JJA,res)
;*****************************************************
res@vpHeightF            = 0.1      ; Changes the aspect ratio
res@vpWidthF             = 0.18
res@vpXF    = 0.72             ; Make plots wider
res@vpYF    = res@vpYF - res@vpHeightF -yspace
res@gsnLeftString  = leftstring4(2)
res@gsnCenterString  = Centerstring4(2)
plot3(2) = gsn_csm_contour_map(wks,clm_APE3_SON,res)
;*****************************************************
res@vpHeightF            = 0.1      ; Changes the aspect ratio
res@vpWidthF             = 0.18
res@vpXF    = 0.72      ; Make plots wider
res@vpYF    = res@vpYF - res@vpHeightF -yspace
res@gsnLeftString  = leftstring4(3)
res@gsnCenterString  = Centerstring4(3)
res@tmXBLabels            = (/"180","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W","80~S~o~N~W"/)
res@lbLabelBarOn   =  True 
plot3(3) = gsn_csm_contour_map(wks,clm_APE3_DJF,res)

;==============================================================
  pres                  = True           ; polyline mods desired
  pres@gsLineThicknessF =2.0            ; line thickness   
  pres@gsLineColor      = "black"         ; color of lines
  pres@gsLineDashPattern =1
  pres@tfPolyDrawOrder       = "PostDraw" 
  
  dum1 = gsn_add_polyline(wks,plot(0), (/180,280/) , (/0,0/) ,pres)
  dum2 = gsn_add_polyline(wks,plot(1), (/180,280/) , (/0,0/) ,pres)
  dum3 = gsn_add_polyline(wks,plot(2), (/180,280/) , (/0,0/) ,pres)
  dum4 = gsn_add_polyline(wks,plot(3), (/180,280/) , (/0,0/) ,pres)
  dum01 = gsn_add_polyline(wks,plot1(0), (/180,280/) , (/0,0/) ,pres)
  dum02 = gsn_add_polyline(wks,plot1(1), (/180,280/) , (/0,0/) ,pres)
  dum03 = gsn_add_polyline(wks,plot1(2), (/180,280/) , (/0,0/) ,pres)
  dum04 = gsn_add_polyline(wks,plot1(3), (/180,280/) , (/0,0/) ,pres)
  dum001 = gsn_add_polyline(wks,plot2(0), (/180,280/) , (/0,0/) ,pres)
  dum002 = gsn_add_polyline(wks,plot2(1), (/180,280/) , (/0,0/) ,pres)
  dum003 = gsn_add_polyline(wks,plot2(2), (/180,280/) , (/0,0/) ,pres)
  dum004 = gsn_add_polyline(wks,plot2(3), (/180,280/) , (/0,0/) ,pres)
  dum0001 = gsn_add_polyline(wks,plot3(0), (/180,280/) , (/0,0/) ,pres)
  dum0002 = gsn_add_polyline(wks,plot3(1), (/180,280/) , (/0,0/) ,pres)
  dum0003 = gsn_add_polyline(wks,plot3(2), (/180,280/) , (/0,0/) ,pres)
  dum0004 = gsn_add_polyline(wks,plot3(3), (/180,280/) , (/0,0/) ,pres)

  eum1 = gsn_add_polyline(wks,plot(0), (/180,280/) , (/5,5/) ,pres)
  eum2 = gsn_add_polyline(wks,plot(1), (/180,280/) , (/5,5/) ,pres)
  eum3 = gsn_add_polyline(wks,plot(2), (/180,280/) , (/5,5/) ,pres)
  eum4 = gsn_add_polyline(wks,plot(3), (/180,280/) , (/5,5/) ,pres)
  eum01 = gsn_add_polyline(wks,plot1(0), (/180,280/) , (/5,5/) ,pres)
  eum02 = gsn_add_polyline(wks,plot1(1), (/180,280/) , (/5,5/) ,pres)
  eum03 = gsn_add_polyline(wks,plot1(2), (/180,280/) , (/5,5/) ,pres)
  eum04 = gsn_add_polyline(wks,plot1(3), (/180,280/) , (/5,5/) ,pres)
  eum001 = gsn_add_polyline(wks,plot2(0), (/180,280/) , (/5,5/) ,pres)
  eum002 = gsn_add_polyline(wks,plot2(1), (/180,280/) , (/5,5/) ,pres)
  eum003 = gsn_add_polyline(wks,plot2(2), (/180,280/) , (/5,5/) ,pres)
  eum004 = gsn_add_polyline(wks,plot2(3), (/180,280/) , (/5,5/) ,pres)
  eum0001 = gsn_add_polyline(wks,plot3(0), (/180,280/) , (/5,5/) ,pres)
  eum0002 = gsn_add_polyline(wks,plot3(1), (/180,280/) , (/5,5/) ,pres)
  eum0003 = gsn_add_polyline(wks,plot3(2), (/180,280/) , (/5,5/) ,pres)
  eum0004 = gsn_add_polyline(wks,plot3(3), (/180,280/) , (/5,5/) ,pres)

;txres                       = True
;txres@txPerimOn             = False  
;txres@txBackgroundFillColor = "White"
;txres@txFontHeightF         = 0.013
;txres@txFontThicknessF      = 1.5
;txid1     = gsn_create_text(wks, "~F22~a    PD-BTR", txres)
;txid2     = gsn_create_text(wks, "~F22~b    2xCO2-BTR", txres)
;txid3     = gsn_create_text(wks, "~F22~c    4xCO2-BTR", txres)
;txid03    = gsn_create_text(wks, "10~S~-6~N~Jm~S~-3~N~s", txres)
;txid4     = gsn_create_text(wks, "~F22~d    PD-BCR", txres)
;txid5     = gsn_create_text(wks, "~F22~e    2xCO2-BCR", txres)
;txid6     = gsn_create_text(wks, "~F22~f    4xCO2-BCR", txres)
;txid06    = gsn_create_text(wks,  "10~S~-6~N~Jm~S~-3~N~s", txres)
;
;amres                  = True
;amres@amParallelPosF   = -0.485     ; 0.5 This is the right edge of the plot.
;amres@amOrthogonalPosF = -0.70     ; 0.5 This is the bottom edge of the plot.
;amres@amJust           = "TopLeft"
;
;annoid1   = gsn_add_annotation(plot0, txid1, amres)
;annoid2   = gsn_add_annotation(plot1, txid2, amres)
;annoid3   = gsn_add_annotation(plot2, txid3, amres) 
;annoid01  = gsn_add_annotation(plot3, txid4, amres)
;annoid02  = gsn_add_annotation(plot4, txid5, amres)
;annoid03  = gsn_add_annotation(plot5, txid6, amres) 
;
;txres@txFontHeightF         = 0.015
;amres@amParallelPosF   = 0.51     ; 0.5 This is the right edge of the plot.
;amres@amOrthogonalPosF = -0.72     ; 0.5 This is the bottom edge of the plot.
;txres@txFontHeightF         = 0.008
;annoid03 = gsn_add_annotation(plot2, txid03, amres)  
;annoid03 = gsn_add_annotation(plot5, txid06, amres)  

plres=True
plres@gsLineThicknessF=2.5
; plres@gsFillColor     = -1
plres@gsLineDashPattern= 2
plres@gsLineColor="blue"
plres@gsLineOpacityF=100
w1=360-160
e1=360-110

w2=360-160
e2=360-110


s1=3
n1=5

s2 = 0
n2 = 2

ddum1=gsn_add_polyline(wks, plot(0), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum2=gsn_add_polyline(wks, plot(1), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum3=gsn_add_polyline(wks, plot(2), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum4=gsn_add_polyline(wks, plot(3), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)

ddum01=gsn_add_polyline(wks, plot1(0), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum02=gsn_add_polyline(wks, plot1(1), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum03=gsn_add_polyline(wks, plot1(2), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum04=gsn_add_polyline(wks, plot1(3), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)

ddum001=gsn_add_polyline(wks, plot2(0), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum002=gsn_add_polyline(wks, plot2(1), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum003=gsn_add_polyline(wks, plot2(2), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum004=gsn_add_polyline(wks, plot2(3), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)

ddum0001=gsn_add_polyline(wks, plot3(0), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum0002=gsn_add_polyline(wks, plot3(1), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum0003=gsn_add_polyline(wks, plot3(2), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
ddum0004=gsn_add_polyline(wks, plot3(3), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)

ddum6=gsn_add_polyline(wks, plot(0), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum7=gsn_add_polyline(wks, plot(1), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum8=gsn_add_polyline(wks, plot(2), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum9=gsn_add_polyline(wks, plot(3), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)

ddum06=gsn_add_polyline(wks, plot1(0), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum07=gsn_add_polyline(wks, plot1(1), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum08=gsn_add_polyline(wks, plot1(2), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum09=gsn_add_polyline(wks, plot1(3), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)

ddum006=gsn_add_polyline(wks, plot2(0), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum007=gsn_add_polyline(wks, plot2(1), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum008=gsn_add_polyline(wks, plot2(2), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum009=gsn_add_polyline(wks, plot2(3), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)

ddum0006=gsn_add_polyline(wks, plot3(0), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum0007=gsn_add_polyline(wks, plot3(1), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum0008=gsn_add_polyline(wks, plot3(2), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
ddum0009=gsn_add_polyline(wks, plot3(3), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
draw(plot(0))
draw(plot(1))
draw(plot(2))
draw(plot(3))

draw(plot1(0))
draw(plot1(1))
draw(plot1(2))
draw(plot1(3))

draw(plot2(0))
draw(plot2(1))
draw(plot2(2))
draw(plot2(3))

draw(plot3(0))
draw(plot3(1))
draw(plot3(2))
draw(plot3(3))

;===================================================
data1 = new((/2,12/), double)
data1(0,:) = CKE1_ave
data1(1,:) = CKE3_ave

data2 = new((/2,12/), double)
data2(0,:) = CKE01_ave
data2(1,:) = CKE03_ave

data3 = new((/2,12/), double)
data3(0,:) = APE1_ave
data3(1,:) = APE3_ave

data4 = new((/2,12/), double)
data4(0,:) = APE01_ave
data4(1,:) = APE03_ave

mres               = True 
mres@gsnFrame      = False
mres@gsnDraw       = False
mres@trXMinF       = 0.
mres@trXMaxF       = 13.
mres@trYMinF       = -10.
mres@trYMaxF       = 80.
mres@vpXF          = 0.05                  ; Move slightly to left
mres@vpYF          = 0.38
mres@vpHeightF     = 0.1    
mres@vpWidthF      = 0.18 
mres@gsnYRefLine              = 0. 
mres@gsnYRefLineColors        = (/"black"/)
mres@gsnYRefLineThicknessF    = 1.0
mres@gsnYRefLineDashPattern   = 2 
;mres@gsnXRefLine              = 11.0 
;mres@gsnXRefLineColors        = (/"black"/)
;mres@gsnXRefLineThicknessF    = 2.0
;mres@gsnXRefLineDashPattern   = 1
mres@tfPolyDrawOrder          = "PreDraw" 
mres@tmXBMode                 = "Explicit"                 
mres@tmXBValues               = ispan(0,13,1)
mres@tmXBLabels               = (/"","J","F","M","A","M","J","J","A","S","O","N","D",""/)  

mres@tmYLMode                 = "Explicit"   
mres@tmYLValues               = ispan(-5,40,5)*2
mres@tmYLLabels               = (/"-10","0","10","20","30","40","50","60","70","80"/) 
;mres@tmYLLabelFontColor       = "red"

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
mres@pmLegendOrthogonalPosF =  0.35                ; move units down
mres@pmLegendWidthF         = 0.06                ; Change width and
mres@pmLegendHeightF        = 0.04                ; height of legend.
mres@lgPerimOn              = False               ; turn off box around
mres@lgLabelFontHeightF     = 0.010              ; label font height
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
mres@tiMainFontHeightF =0.012 
mres@gsnLeftString      ="~F22~q"
mres@gsnLeftStringFontHeightF =0.010
mres@gsnCenterString   = "BTR_NEQ"
mres@gsnCenterStringFontHeightF = 0.008
x = ispan(1,12,1)
plot00 = gsn_csm_xy(wks,x,data1,mres)
delete(mres@tmYLValues)
delete(mres@tmYLLabels)
mres@pmLegendParallelPosF   = -0.20                ; move units right
mres@pmLegendOrthogonalPosF =  0.35                ; move units down

mres@trYMinF       = -10.
mres@trYMaxF       = 60.
mres@vpXF          = 0.27                 ; Move slightly to left
mres@vpYF          = mres@vpYF
mres@vpHeightF     = 0.1    
mres@vpWidthF      = 0.18 
mres@tmYLValues               = ispan(-5,40,5)*2
mres@tmYLLabels               = (/"-10","0","10","20","30","40","50","60","70","80"/) 
mres@gsnLeftString      ="~F22~r "
mres@gsnCenterString    = "BTR_EQ"
mres@xyExplicitLegendLabels = (/"PD","4xCO2"/)  
plot01 = gsn_csm_xy(wks,x,data2,mres)
delete(mres@tmYLValues)
delete(mres@tmYLLabels)
mres@pmLegendParallelPosF   =-0.20               ; move units right
mres@pmLegendOrthogonalPosF =  0.35                ; move units down

mres@trYMinF       = -5.
mres@trYMaxF       = 40.
mres@vpXF          = 0.50                 ; Move slightly to left
mres@vpYF          = mres@vpYF
mres@vpHeightF     = 0.1    
mres@vpWidthF      = 0.18 
mres@tmYLValues               = ispan(-5,40,5)
mres@tmYLLabels               = (/"-5","0","5","10","15","20","25","30","35","40"/) 
mres@gsnLeftString      ="~F22~s "
mres@gsnCenterString    = " BCR NEQ"
mres@xyExplicitLegendLabels = (/"PD","4xCO2"/)  
plot02 = gsn_csm_xy(wks,x,data3,mres)

delete(mres@tmYLValues)
delete(mres@tmYLLabels)
mres@pmLegendParallelPosF   = -0.20                ; move units right
mres@pmLegendOrthogonalPosF =  0.35                ; move units down

mres@trYMinF       = -5.
mres@trYMaxF       = 30.
mres@vpXF          = 0.72                 ; Move slightly to left
mres@vpYF          = mres@vpYF
mres@vpHeightF     = 0.1    
mres@vpWidthF      = 0.18
mres@tmYLValues               = ispan(-5,40,5)
mres@tmYLLabels               = (/"-5","0","5","10","15","20","25","30","35","40"/) 
mres@gsnLeftString      ="~F22~t "
mres@gsnCenterString    = " BCR EQ"
mres@xyExplicitLegendLabels = (/"PD","4xCO2"/)  
plot03 = gsn_csm_xy(wks,x,data4,mres)


 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "royalblue1"                 ; color chosen
 dummy1 = gsn_add_polygon (wks,plot00,xp1_CKE_PD,yp1_CKE_PD,gsres)

 gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "royalblue1"               ; color chosen
 gsres@gsFillOpacityF     = 0.5
 dummy2 = gsn_add_polygon(wks,plot01,xp2_CKE_PD,yp2_CKE_PD,gsres)

 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       =    "pink1"                ; color chosen
 dummy01 = gsn_add_polygon (wks,plot00,xp1_CKE_4xCO2,yp1_CKE_4xCO2,gsres) 

 gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "pink1"                 ; color chosen
 gsres@gsFillOpacityF     = 0.5
 dummy02 = gsn_add_polygon (wks,plot01,xp2_CKE_4xCO2,yp2_CKE_4xCO2,gsres)

 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "royalblue1"                 ; color chosen
 dummy001 = gsn_add_polygon (wks,plot02,xp1_APE_PD,yp1_APE_PD,gsres)

 gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "pink1"                 ; color chosen
 gsres@gsFillOpacityF     = 0.5
 dummy002 = gsn_add_polygon (wks,plot02,xp1_APE_4xCO2,yp1_APE_4xCO2,gsres)

 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "royalblue1"                 ; color chosen
 dummy0001 = gsn_add_polygon (wks,plot03,xp2_APE_PD,yp2_APE_PD,gsres)

 gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "pink1"                 ; color chosen
 gsres@gsFillOpacityF     = 0.5
 dummy0002 = gsn_add_polygon (wks,plot03,xp2_APE_4xCO2,yp2_APE_4xCO2,gsres)


draw(plot00)
draw(plot01)
draw(plot02)
draw(plot03)

end 


 
