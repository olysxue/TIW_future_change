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
begin 
;=============================================================
; read data
;=============================================================
varname="tiw_std"
ncdir  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/mean/"
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
sst_PD    = wgt_areaave_Wrap(sst(0,:,{2:5},{200:250}), 1.0, 1.0, 0)
sst_2xCO2 = wgt_areaave_Wrap(sst(1,:,{2:5},{200:250}), 1.0, 1.0, 0)
sst_4xCO2 = wgt_areaave_Wrap(sst(2,:,{2:5},{200:250}), 1.0, 1.0, 0)

ssh_PD    = wgt_areaave_Wrap(ssh(0,:,{3:8},{210:250}), 1.0, 1.0, 0)
ssh_2xCO2 = wgt_areaave_Wrap(ssh(1,:,{3:8},{210:250}), 1.0, 1.0, 0)
ssh_4xCO2 = wgt_areaave_Wrap(ssh(2,:,{3:8},{210:250}), 1.0, 1.0, 0)

va1_PD    = wgt_areaave_Wrap(va(0,:,{3:6},{210:250}), 1.0, 1.0, 0)
va1_2xCO2 = wgt_areaave_Wrap(va(1,:,{3:6},{210:250}), 1.0, 1.0, 0)
va1_4xCO2 = wgt_areaave_Wrap(va(2,:,{3:6},{210:250}), 1.0, 1.0, 0)

va2_PD    = wgt_areaave_Wrap(va(0,:,{-2:1},{220:260}), 1.0, 1.0, 0)
va2_2xCO2 = wgt_areaave_Wrap(va(1,:,{-2:1},{220:250}), 1.0, 1.0, 0)
va2_4xCO2 = wgt_areaave_Wrap(va(2,:,{-2:1},{210:250}), 1.0, 1.0, 0)

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

data1 = new((/2,12/), float)
data1(0,:) = sst_PD_ann_ave
data1(1,:) = sst_4xCO2_ann_ave

data2 = new((/2,12/), float)
data2(0,:) = ssh_PD_ann_ave
data2(1,:) = ssh_4xCO2_ann_ave

data3 = new((/2,12/), float)
data3(0,:) = va1_PD_ann_ave
data3(1,:) = va1_4xCO2_ann_ave

data4 = new((/2,12/), float)
data4(0,:) = va2_PD_ann_ave
data4(1,:) = va2_4xCO2_ann_ave

data5 = new((/2,12/), float)
data5(0,:) = eke_PD_ann_ave
data5(1,:) = eke_4xCO2_ann_ave

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

 BootStrap_eke1_PD     = bootstrap_stat(eke_PD_SON, stat, nBoot, 0, opt)
 BootStrap_eke1_2xCO2  = bootstrap_stat(eke_2xCO2_SON, stat, nBoot, 0, opt)
 BootStrap_eke1_4xCO2  = bootstrap_stat(eke_4xCO2_SON, stat, nBoot, 0, opt)
 ;=====================================================
 xBoot_sst1_PD       = BootStrap_sst1_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sstavg1_PD    = BootStrap_sst1_PD[1]       ; Mean of bootstrapped estimates

 xBoot_ssh1_PD       = BootStrap_ssh1_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sshavg1_PD    = BootStrap_ssh1_PD[1]       ; Mean of bootstrapped estimates

 xBoot_eke1_PD       = BootStrap_eke1_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_ekeavg1_PD    = BootStrap_eke1_PD[1]       ; Mean of bootstrapped estimates
;=====================================================
    xBoot_sst1_2xCO2    = BootStrap_sst1_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sstavg1_2xCO2    = BootStrap_sst1_2xCO2[1]       ; Mean of bootstrapped estimates

    xBoot_ssh1_2xCO2    = BootStrap_ssh1_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sshavg1_2xCO2    = BootStrap_ssh1_2xCO2[1]       ; Mean of bootstrapped estimates

    xBoot_eke1_2xCO2    = BootStrap_eke1_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_ekeavg1_2xCO2    = BootStrap_eke1_2xCO2[1]       ; Mean of bootstrapped estimates
;-----------------------------------------------------
    xBoot_sst1_4xCO2    = BootStrap_sst1_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sstavg1_4xCO2    = BootStrap_sst1_4xCO2[1]       ; Mean of bootstrapped estimates

    xBoot_ssh1_4xCO2    = BootStrap_ssh1_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sshavg1_4xCO2    = BootStrap_ssh1_4xCO2[1]       ; Mean of bootstrapped estimates

    xBoot_eke1_4xCO2    = BootStrap_eke1_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_ekeavg1_4xCO2    = BootStrap_eke1_4xCO2[1]       ; Mean of bootstrapped estimates
;=====================================================
 xBootLow_sst1_PD    = bootstrap_estimate(xBoot_sst1_PD, 0.025, False)   ;  2.5% lower confidence bound
  xBootHi_sst1_PD    = bootstrap_estimate(xBoot_sst1_PD, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_ssh1_PD    = bootstrap_estimate(xBoot_ssh1_PD, 0.025, False)   ;  2.5% lower confidence bound
  xBootHi_ssh1_PD    = bootstrap_estimate(xBoot_ssh1_PD, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_eke1_PD    = bootstrap_estimate(xBoot_eke1_PD, 0.025, False)   ;  2.5% lower confidence bound
  xBootHi_eke1_PD    = bootstrap_estimate(xBoot_eke1_PD, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------
 xBootLow_sst1_2xCO2    = bootstrap_estimate(xBoot_sst1_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
  xBootHi_sst1_2xCO2    = bootstrap_estimate(xBoot_sst1_2xCO2, 0.975, False)   ; 97.5% upper confidence bound
 
 xBootLow_ssh1_2xCO2    = bootstrap_estimate(xBoot_ssh1_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
  xBootHi_ssh1_2xCO2    = bootstrap_estimate(xBoot_ssh1_2xCO2, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_eke1_2xCO2    = bootstrap_estimate(xBoot_eke1_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
  xBootHi_eke1_2xCO2    = bootstrap_estimate(xBoot_eke1_2xCO2, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------
xBootLow_sst1_4xCO2    = bootstrap_estimate(xBoot_sst1_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_sst1_4xCO2    = bootstrap_estimate(xBoot_sst1_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_ssh1_4xCO2    = bootstrap_estimate(xBoot_ssh1_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_ssh1_4xCO2    = bootstrap_estimate(xBoot_ssh1_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_eke1_4xCO2    = bootstrap_estimate(xBoot_eke1_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_eke1_4xCO2    = bootstrap_estimate(xBoot_eke1_4xCO2, 0.975, False)   ; 97.5% upper confidence bound
;---------------------------------------------------------------------------------------------------------
;---------------------------------------------------
BootStrap_sst2_PD     = bootstrap_stat(sst_PD_MAM, stat, nBoot, 0, opt)
BootStrap_sst2_2xCO2  = bootstrap_stat(sst_2xCO2_MAM, stat, nBoot, 0, opt)
BootStrap_sst2_4xCO2  = bootstrap_stat(sst_4xCO2_MAM, stat, nBoot, 0, opt)

BootStrap_ssh2_PD     = bootstrap_stat(ssh_PD_MAM, stat, nBoot, 0, opt)
BootStrap_ssh2_2xCO2  = bootstrap_stat(ssh_2xCO2_MAM, stat, nBoot, 0, opt)
BootStrap_ssh2_4xCO2  = bootstrap_stat(ssh_4xCO2_MAM, stat, nBoot, 0, opt)

BootStrap_eke2_PD     = bootstrap_stat(eke_PD_MAM, stat, nBoot, 0, opt)
BootStrap_eke2_2xCO2  = bootstrap_stat(eke_2xCO2_MAM, stat, nBoot, 0, opt)
BootStrap_eke2_4xCO2  = bootstrap_stat(eke_4xCO2_MAM, stat, nBoot, 0, opt)
;=====================================================
xBoot_sst2_PD       = BootStrap_sst2_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_sstavg2_PD    = BootStrap_sst2_PD[1]       ; Mean of bootstrapped estimates

xBoot_ssh2_PD       = BootStrap_ssh2_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_sshavg2_PD    = BootStrap_ssh2_PD[1]       ; Mean of bootstrapped estimates

xBoot_eke2_PD       = BootStrap_eke2_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_ekeavg2_PD    = BootStrap_eke2_PD[1]       ; Mean of bootstrapped estimates
;=====================================================
   xBoot_sst2_2xCO2    = BootStrap_sst2_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_sstavg2_2xCO2    = BootStrap_sst2_2xCO2[1]       ; Mean of bootstrapped estimates

   xBoot_ssh2_2xCO2    = BootStrap_ssh2_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_sshavg2_2xCO2    = BootStrap_ssh2_2xCO2[1]       ; Mean of bootstrapped estimates

   xBoot_eke2_2xCO2    = BootStrap_eke2_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_ekeavg2_2xCO2    = BootStrap_eke2_2xCO2[1]       ; Mean of bootstrapped estimates
;-----------------------------------------------------
   xBoot_sst2_4xCO2    = BootStrap_sst2_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_sstavg2_4xCO2    = BootStrap_sst2_4xCO2[1]       ; Mean of bootstrapped estimates

   xBoot_ssh2_4xCO2    = BootStrap_ssh2_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_sshavg2_4xCO2    = BootStrap_ssh2_4xCO2[1]       ; Mean of bootstrapped estimates

   xBoot_eke2_4xCO2    = BootStrap_eke2_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_ekeavg2_4xCO2    = BootStrap_eke2_4xCO2[1]       ; Mean of bootstrapped estimates
;=====================================================
xBootLow_sst2_PD    = bootstrap_estimate(xBoot_sst2_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_sst2_PD    = bootstrap_estimate(xBoot_sst2_PD, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_ssh2_PD    = bootstrap_estimate(xBoot_ssh2_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_ssh2_PD    = bootstrap_estimate(xBoot_ssh2_PD, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_eke2_PD    = bootstrap_estimate(xBoot_eke2_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_eke2_PD    = bootstrap_estimate(xBoot_eke2_PD, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------
xBootLow_sst2_2xCO2    = bootstrap_estimate(xBoot_sst2_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_sst2_2xCO2    = bootstrap_estimate(xBoot_sst2_2xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_ssh2_2xCO2    = bootstrap_estimate(xBoot_ssh2_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_ssh2_2xCO2    = bootstrap_estimate(xBoot_ssh2_2xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_eke2_2xCO2    = bootstrap_estimate(xBoot_eke2_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_eke2_2xCO2    = bootstrap_estimate(xBoot_eke2_2xCO2, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------
xBootLow_sst2_4xCO2    = bootstrap_estimate(xBoot_sst2_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_sst2_4xCO2    = bootstrap_estimate(xBoot_sst2_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_ssh2_4xCO2    = bootstrap_estimate(xBoot_ssh2_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_ssh2_4xCO2    = bootstrap_estimate(xBoot_ssh2_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_eke2_4xCO2    = bootstrap_estimate(xBoot_eke2_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_eke2_4xCO2    = bootstrap_estimate(xBoot_eke2_4xCO2, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------------------------------------------------------------
; boreal FALL
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

; boreal spring
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

;===========================================================================================================
outdir = "/Users/xueaoyun/Dropbox/TIW_future/figures_annual/"
wks = gsn_open_wks("pdf",outdir+"plot_CESM_PD_2xCO2_4xCO2_TIW_annual_cycle_box_plot")
wks_type = "eps"
yspace   = 0.05
bres               = True 
bres@gsnFrame      = False
bres@gsnDraw       = False
bres@trXMinF       = 0.
bres@trXMaxF       = 13.
bres@trYMinF       = 0.
bres@trYMaxF       = 0.6
bres@vpXF          = 0.1                  ; Move slightly to left
bres@vpYF          = 0.9
bres@vpHeightF     = 0.1    
bres@vpWidthF      = 0.2 
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
bres@gsnLeftString      ="~F22~a "
bres@gsnCenterString    ="SST"
bres@gsnLeftStringFontHeightF =0.012
bres@gsnCenterStringFontHeightF = 0.012
x = ispan(1,12,1)
plot0 = gsn_csm_xy(wks,x,data1,bres)
delete(bres@tmYLValues)
delete(bres@tmYLLabels)
bres@trYMinF       = 0.
bres@trYMaxF       = 5.
bres@vpXF          = 0.1                  ; Move slightly to left
bres@vpYF          = bres@vpYF-bres@vpHeightF-yspace
bres@vpHeightF     = 0.1    
bres@vpWidthF      = 0.2 
bres@tmYLValues               = (/0,1,2,3,4,5/)
bres@tmYLLabels               = (/"0","1","2","3","4","5"/) 
bres@gsnLeftString      ="~F22~c"
bres@gsnCenterString    ="SSH"
plot1 = gsn_csm_xy(wks,x,data2,bres)
delete(bres@tmYLValues)
delete(bres@tmYLLabels)
bres@trYMinF       = 0.
bres@trYMaxF       = 80.
bres@vpXF          = 0.1                  ; Move slightly to left
bres@vpYF          = bres@vpYF-bres@vpHeightF-yspace
bres@vpHeightF     = 0.1    
bres@vpWidthF      = 0.2 
bres@tmYLValues               = (/0,20,40,60,80/)
bres@tmYLLabels               = (/"0","20","40","60","80"/)
bres@gsnLeftString      ="~F22~e"
bres@gsnCenterString    ="EKE"
plot2 = gsn_csm_xy(wks,x,data5,bres)

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
 dummy1 = gsn_add_polygon (wks,plot0,xp1_sst_PD,yp1_sst_PD,gsres)

 gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "pink1"                 ; color chosen
 gsres@gsFillOpacityF     = 0.5
 dummy2 = gsn_add_polygon (wks,plot0,xp1_sst_4xCO2,yp1_sst_4xCO2,gsres)

 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "royalblue1"                 ; color chosen
 dummy01 = gsn_add_polygon (wks,plot1,xp1_ssh_PD,yp1_ssh_PD,gsres)

 gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "pink1"                 ; color chosen
 gsres@gsFillOpacityF     = 0.5
 dummy02 = gsn_add_polygon (wks,plot1,xp1_ssh_4xCO2,yp1_ssh_4xCO2,gsres)

 gsres                   = True                        ; poly res
 gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "royalblue1"                 ; color chosen
 dummy001 = gsn_add_polygon (wks,plot2,xp1_eke_PD,yp1_eke_PD,gsres)

 gsres                   = True                        ; poly res
 ;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
 gsres@gsFillColor       = "pink1"                 ; color chosen
 gsres@gsFillOpacityF     = 0.5
 dummy002 = gsn_add_polygon (wks,plot2,xp1_eke_4xCO2,yp1_eke_4xCO2,gsres)

txres                       = True
txres@txPerimOn             = False  
txres@txBackgroundFillColor = "transparent"
txres@txFontHeightF         = 0.02
txres@txFontThicknessF      = 4
txres@txFont                = 22
txid1     = gsn_create_text(wks, "CESM-UHR", txres)

amres                  = True
amres@amParallelPosF   =  0.2    ; 0.5 This is the right edge of the plot.
amres@amOrthogonalPosF = -1.2    ; 0.5 This is the bottom edge of the plot.
amres@amJust           = "TopLeft"

annoid1 = gsn_add_annotation(plot0, txid1, amres)

draw(plot0)
draw(plot1)
draw(plot2)
;==================================================================
 xyres                       = True                 ; plot mods desired
 xyres@gsnFrame              = False                ; don't advance frame yet
 xyres@gsnDraw               = False
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.8                  ; bring bars down to zero
 xyres@trXMinF               = 0                    ; adds space on either end
 xyres@trXMaxF               = 4                   ; of the 1st and last bars
 xyres@vpWidthF              = 0.08
 xyres@vpHeightF             = 0.10

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
 xyres@tmYLValues        = (/0,0.1,0.2,0.3,0.4/)*2
 xyres@tmYLLabels        = (/"0","0.2","0.4","0.6","0.8"/)
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

 xyres@vpXF    = 0.35
 yspace1=0.0
 yf=0.9
 xyres@vpYF    = yf
 xyres@gsnXYBarChartColors = (/"royalblue1"/)
 xyres@gsnCenterString       = "MAM"
 xyres@gsnCenterStringOrthogonalPosF = -0.13
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString      ="~F22~b"
 xyres@gsnLeftStringFontHeightF =0.012
 plot3  = gsn_csm_xy(wks,ispan(1,3,1),stdarr1(:),xyres)   ; advancing the frame
 xyres@tmYLOn                = False 
 xyres@tmYLBorderOn          = False 
 xyres@vpXF = xyres@vpWidthF+xyres@vpXF+yspace1
 xyres@gsnXYBarChartColors = (/"firebrick1"/)
 xyres@gsnCenterString       = "SON"
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.8               ; bring bars down to zero
 xyres@tmYLValues        = (/0,0.2,0.4,0.6,0.8/)
 xyres@tmYLLabels        = (/"0","0.2","0.4","0.6","0.8"/)
 xyres@gsnCenterStringOrthogonalPosF = -0.10
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString      = ""
 plot4  = gsn_csm_xy(wks,ispan(1,3,1),stdarr01(:),xyres)   ; advancing the frame
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@tmYLOn                = True 
 xyres@tmYLBorderOn          = True 
 xyres@vpXF    = 0.35
 yspace=0.05
 yf=0.9-xyres@vpHeightF-yspace
 xyres@vpYF    = yf
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 5                ; bring bars down to zero
 xyres@gsnXYBarChartColors = (/"royalblue1"/)
 xyres@tmYLValues        = (/0,1,2,3,4,5/)
 xyres@tmYLLabels        = (/"0","1","2","3","4","5"/)
 xyres@gsnCenterString       = "MAM"
 xyres@gsnCenterStringOrthogonalPosF = -0.13
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString      ="~F22~d "
 xyres@gsnLeftStringFontHeightF =0.012
 plot5  = gsn_csm_xy(wks,ispan(1,3,1),stdarr2(:),xyres)   ; advancing the frame
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@tmYLOn                = False 
 xyres@tmYLBorderOn          = False 
 xyres@vpXF = xyres@vpWidthF+xyres@vpXF+yspace1
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 5               ; bring bars down to zero
 xyres@gsnXYBarChartColors = (/"firebrick1"/)
 xyres@tmYLValues        = (/0,1,2,3,4,5/)
 xyres@tmYLLabels        = (/"0","1","2","3","4","5"/)
 xyres@gsnCenterString       = "SON"
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnCenterStringOrthogonalPosF = -0.10
  xyres@gsnLeftString      =""
 plot6  = gsn_csm_xy(wks,ispan(1,3,1),stdarr02(:),xyres)   ; advancing the frame
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@tmYLOn                = True
 xyres@vpXF    = 0.35
 yspace=0.05
 yf=xyres@vpYF-xyres@vpHeightF-yspace
 xyres@vpYF    = yf
 xyres@tmYLBorderOn  = True 
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 100                ; bring bars down to zero
 xyres@gsnXYBarChartColors = (/"royalblue1"/)
 xyres@tmYLValues        = (/0,10,20,30,40,50/)*2
 xyres@tmYLLabels        = (/"0","20","40","60","80","100"/)
 xyres@gsnCenterString       = "MAM"
 xyres@gsnCenterStringOrthogonalPosF = -0.13
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString      ="~F22~f"
 plot7  = gsn_csm_xy(wks,ispan(1,3,1),stdarr3(:),xyres)   ; advancing the frame

 xyres@vpXF = xyres@vpWidthF+xyres@vpXF+yspace1
 xyres@tmYLOn        = False
 xyres@tmYLBorderOn  = False                        ; turn off right axis
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 100                ; bring bars down to zero
 xyres@gsnXYBarChartColors = (/"firebrick1"/)
 xyres@tmYLValues        = (/0,10,20,30,40,50/)*2
 xyres@tmYLLabels        = (/"0","20","40","60","80","100"/)
 xyres@gsnCenterString       = "SON"
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnCenterStringOrthogonalPosF = -0.10
 xyres@gsnLeftString      =""
 plot8  = gsn_csm_xy(wks,ispan(1,3,1),stdarr03(:),xyres)   ; advancing the frame

  polyres=True
  polyres@gsLineDashPattern=0
  polyres@gsLineColor="black"
  polyres@gsLineThicknessF=1.0
  polyres@tfPolyDrawOrder = "postDraw"

  aumm1    = gsn_add_polyline(wks,plot4,(/1,1/),(/xBootLow_sst1_PD,xBootHi_sst1_PD/), polyres)
  aumm01   = gsn_add_polyline(wks,plot4,(/2,2/),(/xBootLow_sst1_2xCO2,xBootHi_sst1_2xCO2/), polyres)
  aumm001  = gsn_add_polyline(wks,plot4,(/3,3/),(/xBootLow_sst1_4xCO2,xBootHi_sst1_4xCO2/), polyres)

  bumm1    = gsn_add_polyline(wks,plot6,(/1,1/),(/xBootLow_ssh1_PD,xBootHi_ssh1_PD/), polyres)
  bumm01   = gsn_add_polyline(wks,plot6,(/2,2/),(/xBootLow_ssh1_2xCO2,xBootHi_ssh1_2xCO2/), polyres)
  bumm001  = gsn_add_polyline(wks,plot6,(/3,3/),(/xBootLow_ssh1_4xCO2,xBootHi_ssh1_4xCO2/), polyres)

  cumm1    = gsn_add_polyline(wks,plot8,(/1,1/),(/xBootLow_eke1_PD,xBootHi_eke1_PD/), polyres)
  cumm01   = gsn_add_polyline(wks,plot8,(/2,2/),(/xBootLow_eke1_2xCO2,xBootHi_eke1_2xCO2/), polyres)
  cumm001  = gsn_add_polyline(wks,plot8,(/3,3/),(/xBootLow_eke1_4xCO2,xBootHi_eke1_4xCO2/), polyres)

  aumm2    = gsn_add_polyline(wks,plot3,(/1,1/),(/xBootLow_sst2_PD,xBootHi_sst2_PD/), polyres)
  aumm02   = gsn_add_polyline(wks,plot3,(/2,2/),(/xBootLow_sst2_2xCO2,xBootHi_sst2_2xCO2/), polyres)
  aumm002  = gsn_add_polyline(wks,plot3,(/3,3/),(/xBootLow_sst2_4xCO2,xBootHi_sst2_4xCO2/), polyres)

  bumm2    = gsn_add_polyline(wks,plot5,(/1,1/),(/xBootLow_ssh2_PD,xBootHi_ssh2_PD/), polyres)
  bumm02   = gsn_add_polyline(wks,plot5,(/2,2/),(/xBootLow_ssh2_2xCO2,xBootHi_ssh2_2xCO2/), polyres)
  bumm002  = gsn_add_polyline(wks,plot5,(/3,3/),(/xBootLow_ssh2_4xCO2,xBootHi_ssh2_4xCO2/), polyres)

  cumm2    = gsn_add_polyline(wks,plot7,(/1,1/),(/xBootLow_eke2_PD,xBootHi_eke2_PD/), polyres)
  cumm02   = gsn_add_polyline(wks,plot7,(/2,2/),(/xBootLow_eke2_2xCO2,xBootHi_eke2_2xCO2/), polyres)
  cumm002  = gsn_add_polyline(wks,plot7,(/3,3/),(/xBootLow_eke2_4xCO2,xBootHi_eke2_4xCO2/), polyres)
  
  wid = 0.24/2
  dumm1   = gsn_add_polyline(wks,plot4,(/1-wid,1+wid/),(/xBootLow_sst1_PD,xBootLow_sst1_PD/), polyres)
  dumm01  = gsn_add_polyline(wks,plot4,(/1-wid,1+wid/),(/xBootHi_sst1_PD,xBootHi_sst1_PD/), polyres)

  eumm1   = gsn_add_polyline(wks,plot4,(/2-wid,2+wid/),(/xBootLow_sst1_2xCO2,xBootLow_sst1_2xCO2/), polyres)
  eumm01  = gsn_add_polyline(wks,plot4,(/2-wid,2+wid/),(/xBootHi_sst1_2xCO2,xBootHi_sst1_2xCO2/), polyres)

  fumm1   = gsn_add_polyline(wks,plot4,(/3-wid,3+wid/),(/xBootLow_sst1_4xCO2,xBootLow_sst1_4xCO2/), polyres)
  fumm01  = gsn_add_polyline(wks,plot4,(/3-wid,3+wid/),(/xBootHi_sst1_4xCO2,xBootHi_sst1_4xCO2/), polyres)
;----------------------------------------------------------------------------------------------------
  dumm2   = gsn_add_polyline(wks,plot6,(/1-wid,1+wid/),(/xBootLow_ssh1_PD,xBootLow_ssh1_PD/), polyres)
  dumm02  = gsn_add_polyline(wks,plot6,(/1-wid,1+wid/),(/xBootHi_ssh1_PD,xBootHi_ssh1_PD/), polyres)

  eumm2   = gsn_add_polyline(wks,plot6,(/2-wid,2+wid/),(/xBootLow_ssh1_2xCO2,xBootLow_ssh1_2xCO2/), polyres)
  eumm02  = gsn_add_polyline(wks,plot6,(/2-wid,2+wid/),(/xBootHi_ssh1_2xCO2,xBootHi_ssh1_2xCO2/), polyres)

  fumm2   = gsn_add_polyline(wks,plot6,(/3-wid,3+wid/),(/xBootLow_ssh1_4xCO2,xBootLow_ssh1_4xCO2/), polyres)
  fumm02  = gsn_add_polyline(wks,plot6,(/3-wid,3+wid/),(/xBootHi_ssh1_4xCO2,xBootHi_ssh1_4xCO2/), polyres)
;-----------------------------------------------------------------------------------------------------
  dumm3   = gsn_add_polyline(wks,plot8,(/1-wid,1+wid/),(/xBootLow_eke1_PD,xBootLow_eke1_PD/), polyres)
  dumm03  = gsn_add_polyline(wks,plot8,(/1-wid,1+wid/),(/xBootHi_eke1_PD,xBootHi_eke1_PD/), polyres)

  eumm3   = gsn_add_polyline(wks,plot8,(/2-wid,2+wid/),(/xBootLow_eke1_2xCO2,xBootLow_eke1_2xCO2/), polyres)
  eumm03  = gsn_add_polyline(wks,plot8,(/2-wid,2+wid/),(/xBootHi_eke1_2xCO2,xBootHi_eke1_2xCO2/), polyres)

  fumm3   = gsn_add_polyline(wks,plot8,(/3-wid,3+wid/),(/xBootLow_eke1_4xCO2,xBootLow_eke1_4xCO2/), polyres)
  fumm03  = gsn_add_polyline(wks,plot8,(/3-wid,3+wid/),(/xBootHi_eke1_4xCO2,xBootHi_eke1_4xCO2/), polyres)
;------------------------------------------------------------------------------------------------------
  dumm4   = gsn_add_polyline(wks,plot3,(/1-wid,1+wid/),(/xBootLow_sst2_PD,xBootLow_sst2_PD/), polyres)
  dumm04  = gsn_add_polyline(wks,plot3,(/1-wid,1+wid/),(/xBootHi_sst2_PD,xBootHi_sst2_PD/), polyres)
  
  eumm4   = gsn_add_polyline(wks,plot3,(/2-wid,2+wid/),(/xBootLow_sst2_2xCO2,xBootLow_sst2_2xCO2/), polyres)
  eumm04  = gsn_add_polyline(wks,plot3,(/2-wid,2+wid/),(/xBootHi_sst2_2xCO2,xBootHi_sst2_2xCO2/), polyres)
  
  fumm4   = gsn_add_polyline(wks,plot3,(/3-wid,3+wid/),(/xBootLow_sst2_4xCO2,xBootLow_sst2_4xCO2/), polyres)
  fumm04  = gsn_add_polyline(wks,plot3,(/3-wid,3+wid/),(/xBootHi_sst2_4xCO2,xBootHi_sst2_4xCO2/), polyres)
;----------------------------------------------------------------------------------------------------
  dumm5   = gsn_add_polyline(wks,plot5,(/1-wid,1+wid/),(/xBootLow_ssh2_PD,xBootLow_ssh2_PD/), polyres)
  dumm05  = gsn_add_polyline(wks,plot5,(/1-wid,1+wid/),(/xBootHi_ssh2_PD,xBootHi_ssh2_PD/), polyres)
  
  eumm5   = gsn_add_polyline(wks,plot5,(/2-wid,2+wid/),(/xBootLow_ssh2_2xCO2,xBootLow_ssh2_2xCO2/), polyres)
  eumm05  = gsn_add_polyline(wks,plot5,(/2-wid,2+wid/),(/xBootHi_ssh2_2xCO2,xBootHi_ssh2_2xCO2/), polyres)
  
  fumm5   = gsn_add_polyline(wks,plot5,(/3-wid,3+wid/),(/xBootLow_ssh2_4xCO2,xBootLow_ssh2_4xCO2/), polyres)
  fumm05  = gsn_add_polyline(wks,plot5,(/3-wid,3+wid/),(/xBootHi_ssh2_4xCO2,xBootHi_ssh2_4xCO2/), polyres)
;-----------------------------------------------------------------------------------------------------
  dumm6   = gsn_add_polyline(wks,plot7,(/1-wid,1+wid/),(/xBootLow_eke2_PD,xBootLow_eke2_PD/), polyres)
  dumm06  = gsn_add_polyline(wks,plot7,(/1-wid,1+wid/),(/xBootHi_eke2_PD,xBootHi_eke2_PD/), polyres)
  
  eumm6   = gsn_add_polyline(wks,plot7,(/2-wid,2+wid/),(/xBootLow_eke2_2xCO2,xBootLow_eke2_2xCO2/), polyres)
  eumm06  = gsn_add_polyline(wks,plot7,(/2-wid,2+wid/),(/xBootHi_eke2_2xCO2,xBootHi_eke2_2xCO2/), polyres)
  
  fumm6   = gsn_add_polyline(wks,plot7,(/3-wid,3+wid/),(/xBootLow_eke2_4xCO2,xBootLow_eke2_4xCO2/), polyres)
  fumm06  = gsn_add_polyline(wks,plot7,(/3-wid,3+wid/),(/xBootHi_eke2_4xCO2,xBootHi_eke2_4xCO2/), polyres)
;==================================================================================================================
draw(plot3)
draw(plot4)
draw(plot5)
draw(plot6)
draw(plot7)
draw(plot8)
;==========================================read E3SM======================================================
;read E3SM 
varname="tos"
ncdir  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/E3SM/mean/"
scale=1. ;units degC
filelist=(/"tos_day_E3SM-1-0_picontrol_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mon.nc","tos_day_E3SM-1-0_4xCO2_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mon.nc"/)
sst0=readdata1(ncdir, filelist, varname, scale)
sst0@units="degC"
sst1 = sqrt(sst0)
copy_VarMeta(sst0, sst1)
printVarSummary(sst1)
varname="zos"
scale=1. ; units cm to m by scale
filelist=(/"zos_day_E3SM-1-0_picontrol_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mon.nc","zos_day_E3SM-1-0_4xCO2_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mon.nc"/)
ssh0=readdata1(ncdir, filelist, varname, scale)
ssh0@units="meter"
printVarSummary(ssh0)
ssh1 = sqrt(ssh0)*100.
copy_VarMeta(ssh0, ssh1)
scale=1.

ncdir1  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/E3SM/mean/"
filename1   = "TIW_KE_mon_E3SM-1-0_picontrol_gr_00010101-00831231_regrid_1deg_equator_new.nc"
filename2   = "TIW_KE_mon_E3SM-1-0_4xCO2_gr_00010101-00831231_regrid_1deg_equator_new.nc"
fil1        = systemfunc("ls " + ncdir1 + filename1)
fil2        = systemfunc("ls " + ncdir1 + filename2)
f1          =  addfile(fil1, "r")
f2          =  addfile(fil2, "r")
eke01       = f1->KE(:,:,:)
eke02       = f2->KE(:,:,:)

eke1 = new((/2, dimsizes(eke01&time),dimsizes(eke01&lat), dimsizes(eke01&lon)/), typeof(eke01))
eke1(0,:,:,:) = eke01
eke1(1,:,:,:) = eke02
copy_VarCoords(eke01, eke1(0,:,:,:))
copy_VarCoords(eke01, eke1(1,:,:,:))

sst1_PD    = wgt_areaave_Wrap(sst1(0,0:959,{0:5},{200:260}), 1.0, 1.0, 0)
sst1_4xCO2 = wgt_areaave_Wrap(sst1(1,0:959,{0:5},{200:260}), 1.0, 1.0, 0)

ssh1_PD    = wgt_areaave_Wrap(ssh1(0,0:959,{3:8},{210:250}), 1.0, 1.0, 0)
ssh1_4xCO2 = wgt_areaave_Wrap(ssh1(1,0:959,{3:8},{200:240}), 1.0, 1.0, 0)

eke1_PD    = wgt_areaave_Wrap(eke1(0,0:959,{0:5},{210:260}), 1.0, 1.0, 0)
eke1_4xCO2 = wgt_areaave_Wrap(eke1(1,0:959,{0:5},{200:250}), 1.0, 1.0, 0)
printVarSummary(eke1_PD)
ntim  = 960
delete([/sst1, ssh1, eke1/]);
nyears = ntim/12
sst1_PD_ann  = reshape(sst1_PD,(/nyears,12/))
sst1_4xCO2_ann  = reshape(sst1_4xCO2,(/nyears,12/))

ssh1_PD_ann  = reshape(ssh1_PD,(/nyears,12/))
ssh1_4xCO2_ann  = reshape(ssh1_4xCO2,(/nyears,12/))

eke1_PD_ann  = reshape(eke1_PD,(/nyears,12/))
eke1_4xCO2_ann  = reshape(eke1_4xCO2,(/nyears,12/))

; annual cycle
sst1_PD_ann_ave = dim_avg_n_Wrap(sst1_PD_ann, 0)
sst1_4xCO2_ann_ave = dim_avg_n_Wrap(sst1_4xCO2_ann, 0)

ssh1_PD_ann_ave = dim_avg_n_Wrap(ssh1_PD_ann, 0)
ssh1_4xCO2_ann_ave = dim_avg_n_Wrap(ssh1_4xCO2_ann, 0)

eke1_PD_ann_ave = dim_avg_n_Wrap(eke1_PD_ann, 0)
eke1_4xCO2_ann_ave = dim_avg_n_Wrap(eke1_4xCO2_ann, 0)

; std 
sst1_PD_ann_std = dim_stddev_n_Wrap(sst1_PD_ann, 0)
sst1_4xCO2_ann_std = dim_stddev_n_Wrap(sst1_4xCO2_ann, 0)

ssh1_PD_ann_std = dim_stddev_n_Wrap(ssh1_PD_ann, 0)
ssh1_4xCO2_ann_std = dim_stddev_n_Wrap(ssh1_4xCO2_ann, 0)

eke1_PD_ann_std = dim_stddev_n_Wrap(eke1_PD_ann, 0)
eke1_4xCO2_ann_std = dim_stddev_n_Wrap(eke1_4xCO2_ann, 0)

data01 = new((/2,12/), double)
data01(0,:) = sst1_PD_ann_ave
data01(1,:) = sst1_4xCO2_ann_ave

data02 = new((/2,12/), double)
data02(0,:) = ssh1_PD_ann_ave
data02(1,:) = ssh1_4xCO2_ann_ave

data03 = new((/2,12/), double)
data03(0,:) = eke1_PD_ann_ave
data03(1,:) = eke1_4xCO2_ann_ave
;==========================================================
xp2_sst_PD = new((/24/), double)
yp2_sst_PD = new((/24/), double)

do i=0,11
   xp2_sst_PD(i)          = i+1
   xp2_sst_PD(2*12-1-i)   = i+1
   yp2_sst_PD(i)          = sst1_PD_ann_ave(i) + sst1_PD_ann_std(i)/2.
   yp2_sst_PD(2*12-1-i)   = sst1_PD_ann_ave(i) - sst1_PD_ann_std(i)/2. 
 end do

 xp2_sst_4xCO2 = new((/24/), double)
 yp2_sst_4xCO2 = new((/24/), double)
 
 do i=0,11
    xp2_sst_4xCO2(i)          = i+1
    xp2_sst_4xCO2(2*12-1-i)   = i+1
    yp2_sst_4xCO2(i)          = sst1_4xCO2_ann_ave(i) + sst1_4xCO2_ann_std(i)/2.
    yp2_sst_4xCO2(2*12-1-i)   = sst1_4xCO2_ann_ave(i) - sst1_4xCO2_ann_std(i)/2. 
  end do
;=======================================================
xp2_ssh_PD = new((/24/), double)
yp2_ssh_PD = new((/24/), double)

do i=0,11
   xp2_ssh_PD(i)          = i+1
   xp2_ssh_PD(2*12-1-i)   = i+1
   yp2_ssh_PD(i)          = ssh1_PD_ann_ave(i) + ssh1_PD_ann_std(i)/2.
   yp2_ssh_PD(2*12-1-i)   = ssh1_PD_ann_ave(i) - ssh1_PD_ann_std(i)/2. 
 end do

 xp2_ssh_4xCO2 = new((/24/), double)
 yp2_ssh_4xCO2 = new((/24/), double)
 
 do i=0,11
    xp2_ssh_4xCO2(i)          = i+1
    xp2_ssh_4xCO2(2*12-1-i)   = i+1
    yp2_ssh_4xCO2(i)          = ssh1_4xCO2_ann_ave(i) + ssh1_4xCO2_ann_std(i)/2.
    yp2_ssh_4xCO2(2*12-1-i)   = ssh1_4xCO2_ann_ave(i) - ssh1_4xCO2_ann_std(i)/2. 
  end do
;-------------------------------------------------------
xp2_eke_PD = new((/24/), double)
yp2_eke_PD = new((/24/), double)

do i=0,11
   xp2_eke_PD(i)          = i+1
   xp2_eke_PD(2*12-1-i)   = i+1
   yp2_eke_PD(i)          = eke1_PD_ann_ave(i) + eke1_PD_ann_std(i)/2.
   yp2_eke_PD(2*12-1-i)   = eke1_PD_ann_ave(i) - eke1_PD_ann_std(i)/2. 
 end do

 xp2_eke_4xCO2 = new((/24/), double)
 yp2_eke_4xCO2 = new((/24/), double)
 
 do i=0,11
    xp2_eke_4xCO2(i)          = i+1
    xp2_eke_4xCO2(2*12-1-i)   = i+1
    yp2_eke_4xCO2(i)          = eke1_4xCO2_ann_ave(i) + eke1_4xCO2_ann_std(i)/2.
    yp2_eke_4xCO2(2*12-1-i)   = eke1_4xCO2_ann_ave(i) - eke1_4xCO2_ann_std(i)/2. 
  end do
;=========================================================================================

;SON TIWs intensity 
sst1_PD_SON    = (sst1_PD(8::12)    + sst1_PD(9::12)   + sst1_PD(10::12))/3.
sst1_4xCO2_SON = (sst1_4xCO2(8::12) + sst1_4xCO2(9::12)+ sst1_4xCO2(10::12))/3.

ssh1_PD_SON    = (ssh1_PD(8::12)    + ssh1_PD(9::12)   + ssh1_PD(10::12))/3.
ssh1_4xCO2_SON = (ssh1_4xCO2(8::12) + ssh1_4xCO2(9::12)+ ssh1_4xCO2(10::12))/3.

eke1_PD_SON    = (eke1_PD(8::12)    + eke1_PD(9::12)   + eke1_PD(10::12))/3.
eke1_4xCO2_SON = (eke1_4xCO2(8::12) + eke1_4xCO2(9::12)+ eke1_4xCO2(10::12))/3.

; MAM TIWs activity
sst1_PD_MAM    = (sst1_PD(2::12)    + sst1_PD(3::12)   + sst1_PD(4::12))/3.
sst1_4xCO2_MAM = (sst1_4xCO2(2::12) + sst1_4xCO2(3::12)+ sst1_4xCO2(4::12))/3.

ssh1_PD_MAM    = (ssh1_PD(2::12)    + ssh1_PD(3::12)   + ssh1_PD(4::12))/3.
ssh1_4xCO2_MAM = (ssh1_4xCO2(2::12) + ssh1_4xCO2(3::12)+ ssh1_4xCO2(4::12))/3.

eke1_PD_MAM    = (eke1_PD(2::12)    + eke1_PD(3::12)   + eke1_PD(4::12))/3.
eke1_4xCO2_MAM = (eke1_4xCO2(2::12) + eke1_4xCO2(3::12)+ eke1_4xCO2(4::12))/3.
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

 BootStrap_sst01_PD     = bootstrap_stat(sst1_PD_SON, stat, nBoot, 0, opt)
 BootStrap_sst01_4xCO2  = bootstrap_stat(sst1_4xCO2_SON, stat, nBoot, 0, opt)

 BootStrap_ssh01_PD     = bootstrap_stat(ssh1_PD_SON, stat, nBoot, 0, opt)
 BootStrap_ssh01_4xCO2  = bootstrap_stat(ssh1_4xCO2_SON, stat, nBoot, 0, opt)

 BootStrap_eke01_PD     = bootstrap_stat(eke1_PD_SON, stat, nBoot, 0, opt)
 BootStrap_eke01_4xCO2  = bootstrap_stat(eke1_4xCO2_SON, stat, nBoot, 0, opt)
 ;=====================================================
    xBoot_sst01_PD    = BootStrap_sst01_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sstavg01_PD    = BootStrap_sst01_PD[1]       ; Mean of bootstrapped estimates

    xBoot_ssh01_PD    = BootStrap_ssh01_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sshavg01_PD    = BootStrap_ssh01_PD[1]       ; Mean of bootstrapped estimates

    xBoot_eke01_PD    = BootStrap_eke01_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_ekeavg01_PD    = BootStrap_eke01_PD[1]       ; Mean of bootstrapped estimates
;=====================================================
    xBoot_sst01_4xCO2    = BootStrap_sst01_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sstavg01_4xCO2    = BootStrap_sst01_4xCO2[1]       ; Mean of bootstrapped estimates

    xBoot_ssh01_4xCO2    = BootStrap_ssh01_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sshavg01_4xCO2    = BootStrap_ssh01_4xCO2[1]       ; Mean of bootstrapped estimates

    xBoot_eke01_4xCO2    = BootStrap_eke01_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_ekeavg01_4xCO2    = BootStrap_eke01_4xCO2[1]       ; Mean of bootstrapped estimates
;=====================================================
 xBootLow_sst01_PD    = bootstrap_estimate(xBoot_sst01_PD, 0.025, False)   ;  2.5% lower confidence bound
  xBootHi_sst01_PD    = bootstrap_estimate(xBoot_sst01_PD, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_ssh01_PD    = bootstrap_estimate(xBoot_ssh01_PD, 0.025, False)   ;  2.5% lower confidence bound
  xBootHi_ssh01_PD    = bootstrap_estimate(xBoot_ssh01_PD, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_eke01_PD    = bootstrap_estimate(xBoot_eke01_PD, 0.025, False)   ;  2.5% lower confidence bound
  xBootHi_eke01_PD    = bootstrap_estimate(xBoot_eke01_PD, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------
xBootLow_sst01_4xCO2    = bootstrap_estimate(xBoot_sst01_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_sst01_4xCO2    = bootstrap_estimate(xBoot_sst01_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_ssh01_4xCO2    = bootstrap_estimate(xBoot_ssh01_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_ssh01_4xCO2    = bootstrap_estimate(xBoot_ssh01_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_eke01_4xCO2    = bootstrap_estimate(xBoot_eke01_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_eke01_4xCO2    = bootstrap_estimate(xBoot_eke01_4xCO2, 0.975, False)   ; 97.5% upper confidence bound
;---------------------------------------------------------------------------------------------------------
;---------------------------------------------------
BootStrap_sst02_PD     = bootstrap_stat(sst1_PD_MAM, stat, nBoot, 0, opt)
BootStrap_sst02_4xCO2  = bootstrap_stat(sst1_4xCO2_MAM, stat, nBoot, 0, opt)

BootStrap_ssh02_PD     = bootstrap_stat(ssh1_PD_MAM, stat, nBoot, 0, opt)
BootStrap_ssh02_4xCO2  = bootstrap_stat(ssh1_4xCO2_MAM, stat, nBoot, 0, opt)

BootStrap_eke02_PD     = bootstrap_stat(eke1_PD_MAM, stat, nBoot, 0, opt)
BootStrap_eke02_4xCO2  = bootstrap_stat(eke1_4xCO2_MAM, stat, nBoot, 0, opt)
;=====================================================
xBoot_sst02_PD       = BootStrap_sst02_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_sstavg02_PD    = BootStrap_sst02_PD[1]       ; Mean of bootstrapped estimates

xBoot_ssh02_PD       = BootStrap_ssh02_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_sshavg02_PD    = BootStrap_ssh02_PD[1]       ; Mean of bootstrapped estimates

xBoot_eke02_PD       = BootStrap_eke02_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_ekeavg02_PD    = BootStrap_eke02_PD[1]       ; Mean of bootstrapped estimates
;=====================================================
   xBoot_sst02_4xCO2    = BootStrap_sst02_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_sstavg02_4xCO2    = BootStrap_sst02_4xCO2[1]       ; Mean of bootstrapped estimates

   xBoot_ssh02_4xCO2    = BootStrap_ssh02_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_sshavg02_4xCO2    = BootStrap_ssh02_4xCO2[1]       ; Mean of bootstrapped estimates

   xBoot_eke02_4xCO2    = BootStrap_eke02_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_ekeavg02_4xCO2    = BootStrap_eke02_4xCO2[1]       ; Mean of bootstrapped estimates
;=====================================================
xBootLow_sst02_PD    = bootstrap_estimate(xBoot_sst02_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_sst02_PD    = bootstrap_estimate(xBoot_sst02_PD, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_ssh02_PD    = bootstrap_estimate(xBoot_ssh02_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_ssh02_PD    = bootstrap_estimate(xBoot_ssh02_PD, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_eke02_PD    = bootstrap_estimate(xBoot_eke02_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_eke02_PD    = bootstrap_estimate(xBoot_eke02_PD, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------
xBootLow_sst02_4xCO2    = bootstrap_estimate(xBoot_sst02_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_sst02_4xCO2    = bootstrap_estimate(xBoot_sst02_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_ssh02_4xCO2    = bootstrap_estimate(xBoot_ssh02_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_ssh02_4xCO2    = bootstrap_estimate(xBoot_ssh02_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_eke02_4xCO2    = bootstrap_estimate(xBoot_eke02_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
 xBootHi_eke02_4xCO2    = bootstrap_estimate(xBoot_eke02_4xCO2, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------------------------------------------------------------
; boreal FALL
stdarr4 = new((/2/),"double") ;  SST
stdarr4(0) =  xBoot_sstavg01_PD
stdarr4(1) =  xBoot_sstavg01_4xCO2

stdarr5 = new((/2/),"double") ; SSH
stdarr5(0) =  xBoot_sshavg01_PD
stdarr5(1) =  xBoot_sshavg01_4xCO2

stdarr6 = new((/2/),"double") ; EKE 
stdarr6(0) =  xBoot_ekeavg01_PD
stdarr6(1) =  xBoot_ekeavg01_4xCO2 

; boreal spring
stdarr04 = new((/2/),"double") ;  SST
stdarr04(0) =  xBoot_sstavg02_PD
stdarr04(1) =  xBoot_sstavg02_4xCO2

stdarr05 = new((/2/),"double") ; SSH
stdarr05(0) =  xBoot_sshavg02_PD
stdarr05(1) =  xBoot_sshavg02_4xCO2

stdarr06 = new((/2/),"double") ; EKE 
stdarr06(0) =  xBoot_ekeavg02_PD
stdarr06(1) =  xBoot_ekeavg02_4xCO2 

delete(bres@tmYLValues)
delete(bres@tmYLLabels)   
; draw E3SM
bres@trXMinF       = 0.
bres@trXMaxF       = 13.
bres@trYMinF       = 0.
bres@trYMaxF       = 0.8
bres@vpXF          = 0.6                  ; Move slightly to left
bres@vpYF          = 0.9
bres@vpHeightF     = 0.1    
bres@vpWidthF      = 0.2 
bres@gsnLeftString      ="~F22~g"
bres@gsnCenterString    ="SST"
bres@gsnLeftStringFontHeightF =0.012
bres@tmYLMode                 = "Explicit"   
bres@tmYLValues               = (/0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8/)
bres@tmYLLabels               = (/"0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8"/) 
x = ispan(1,12,1)
plot00 = gsn_csm_xy(wks,x,data01,bres)
delete(bres@tmYLValues)
delete(bres@tmYLLabels)
bres@trYMinF       = 1.
bres@trYMaxF       = 5.
bres@vpXF          = 0.6                 ; Move slightly to left
bres@vpYF          = bres@vpYF-bres@vpHeightF-yspace
bres@vpHeightF     = 0.1    
bres@vpWidthF      = 0.2 
bres@tmYLValues               = (/0,1,2,3,4,5/)
bres@tmYLLabels               = (/"0","1","2","3","4","5"/) 
bres@gsnLeftString      ="~F22~i"
bres@gsnCenterString    ="SSH"
plot01 = gsn_csm_xy(wks,x,data02,bres)
delete(bres@tmYLValues)
delete(bres@tmYLLabels)
bres@trYMinF       = 0.
bres@trYMaxF       = 180.
bres@vpXF          = 0.6                  ; Move slightly to left
bres@vpYF          = bres@vpYF-bres@vpHeightF-yspace
bres@vpHeightF     = 0.1    
bres@vpWidthF      = 0.2 
bres@tmYLValues               = (/0,30,60,90,120,150,180/)
bres@tmYLLabels               = (/"0","30","60","90","120","150","180"/)
bres@gsnLeftString      ="~F22~k"
bres@gsnCenterString    ="EKE"
plot02 = gsn_csm_xy(wks,x,data03,bres)
;=======================================================================
gsres                   = True                        ; poly res
gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
gsres@gsFillColor       = "royalblue1"                 ; color chosen
dummy3 = gsn_add_polygon (wks,plot00,xp2_sst_PD,yp2_sst_PD,gsres)

gsres                   = True                        ; poly res
;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
gsres@gsFillColor       = "pink1"                 ; color chosen
gsres@gsFillOpacityF     = 0.5
dummy4 = gsn_add_polygon (wks,plot00,xp2_sst_4xCO2,yp2_sst_4xCO2,gsres)

gsres                   = True                        ; poly res
gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
gsres@gsFillColor       = "royalblue1"                 ; color chosen
dummy03 = gsn_add_polygon (wks,plot01,xp2_ssh_PD,yp2_ssh_PD,gsres)

gsres                   = True                        ; poly res
;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
gsres@gsFillColor       = "pink1"                 ; color chosen
gsres@gsFillOpacityF     = 0.5
dummy04 = gsn_add_polygon (wks,plot01,xp2_ssh_4xCO2,yp2_ssh_4xCO2,gsres)

gsres                   = True                        ; poly res
gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
gsres@gsFillColor       = "royalblue1"                 ; color chosen
dummy003 = gsn_add_polygon (wks,plot02,xp2_eke_PD,yp2_eke_PD,gsres)

gsres                   = True                        ; poly res
;gsres@tfPolyDrawOrder   = "Draw"                   ; draw this first
gsres@gsFillColor       = "pink1"                 ; color chosen
gsres@gsFillOpacityF     = 0.5
dummy004 = gsn_add_polygon (wks,plot02,xp2_eke_4xCO2,yp2_eke_4xCO2,gsres)


txres                       = True
txres@txPerimOn             = False  
txres@txBackgroundFillColor = "transparent"
txres@txFontHeightF         = 0.02
txres@txFontThicknessF      = 4
txres@txFont                = 22
txid2     = gsn_create_text(wks, "E3SMv1", txres)

amres                  = True
amres@amParallelPosF   =  0.25    ; 0.5 This is the right edge of the plot.
amres@amOrthogonalPosF = -1.2    ; 0.5 This is the bottom edge of the plot.
amres@amJust           = "TopLeft"

annoid2 = gsn_add_annotation(plot00, txid2, amres)

draw(plot00)
draw(plot01)
draw(plot02)

delete(xyres@tmYLValues)
delete(xyres@tmYLLabels)
delete(xyres@tmXBValues)
delete(xyres@tmXBLabels)

xyres@vpXF    = 0.85
yspace1=0.0
 yf=0.9
 xyres@trXMinF               = 0                    ; adds space on either end
 xyres@trXMaxF               = 3                   ; of the 1st and last bars
 xyres@vpWidthF              = 0.06
 xyres@vpHeightF             = 0.10

 xyres@gsnXYBarChart         = True                 ; turn on bar chart
 xyres@gsnXYBarChartBarWidth = 0.24                  ; change bar widths
 xyres@tmXBOn                = True                 ; turn off tickmarks at bot
 xyres@tmYLOn                = True
 xyres@tmYLBorderOn          = True
 xyres@tmYROn                = False                 ; turn off tickmarks at bot
 xyres@tmXTOn                = False                 ; turn off tickmarks at bot
 xyres@tmXBMode          = "Explicit"         ; explicit labels
 xyres@tmXBValues        = (/1,2,3/)
 xyres@tmXBLabels        = (/"PD","4CO2",""/)
 xyres@vpYF    = yf
 xyres@gsnXYBarChartColors = (/"royalblue1"/)
 xyres@gsnCenterString       = "MAM"
 xyres@gsnCenterStringOrthogonalPosF = -0.13
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString      ="~F22~h"
 xyres@gsnLeftStringFontHeightF =0.012
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.8               ; bring bars down to zero
 xyres@tmYLValues               = (/0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8/)
 xyres@tmYLLabels               = (/"0","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8"/) 
 plot03  = gsn_csm_xy(wks,ispan(1,2,1),stdarr04(:),xyres)   ; advancing the frame
 xyres@tmYLOn                = False 
 xyres@tmYLBorderOn          = False 
 xyres@vpXF = xyres@vpWidthF+xyres@vpXF+yspace1
 xyres@gsnXYBarChartColors = (/"firebrick1"/)
 xyres@gsnCenterString       = "SON"
 xyres@gsnCenterStringOrthogonalPosF = -0.10
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString      = ""
 plot04  = gsn_csm_xy(wks,ispan(1,2,1),stdarr4(:),xyres)   ; advancing the frame
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@tmYLOn                = True 
 xyres@tmYLBorderOn          = True 
 xyres@vpXF    = 0.85
 yspace=0.05
 yf=0.9-xyres@vpHeightF-yspace
 xyres@vpYF    = yf
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 5                ; bring bars down to zero
 xyres@gsnXYBarChartColors = (/"royalblue1"/)
 xyres@tmYLValues        = (/0,1,2,3,4,5/)
 xyres@tmYLLabels        = (/"0","1","2","3","4","5"/)
 xyres@gsnCenterString       = "MAM"
 xyres@gsnCenterStringOrthogonalPosF = -0.13
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString      ="~F22~j "
 xyres@gsnLeftStringFontHeightF =0.012
 plot05  = gsn_csm_xy(wks,ispan(1,2,1),stdarr05(:),xyres)   ; advancing the frame
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@tmYLOn                = False 
 xyres@tmYLBorderOn          = False 
 xyres@vpXF = xyres@vpWidthF+xyres@vpXF+yspace1
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 5               ; bring bars down to zero
 xyres@gsnXYBarChartColors = (/"firebrick1"/)
 xyres@tmYLValues        = (/0,1,2,3,4,5/)
 xyres@tmYLLabels        = (/"0","1","2","3","4","5"/)
 xyres@gsnCenterString       = "SON"
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnCenterStringOrthogonalPosF = -0.10
 xyres@gsnLeftString      =""
 plot06  = gsn_csm_xy(wks,ispan(1,2,1),stdarr5(:),xyres)   ; advancing the frame
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@tmYLOn                = True
 xyres@vpXF    = 0.85
 yspace=0.05
 yf=xyres@vpYF-xyres@vpHeightF-yspace
 xyres@vpYF    = yf
 xyres@tmYLBorderOn  = True 
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 180                ; bring bars down to zero
 xyres@gsnXYBarChartColors = (/"royalblue1"/)
 xyres@tmYLValues               = (/0,30,60,90,120,150,180/)
 xyres@tmYLLabels               = (/"0","30","60","90","120","150","180"/)
 xyres@gsnCenterString       = "MAM"
 xyres@gsnCenterStringOrthogonalPosF = -0.13
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString      ="~F22~l "
 plot07  = gsn_csm_xy(wks,ispan(1,2,1),stdarr06(:),xyres)   ; advancing the frame
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@vpXF = xyres@vpWidthF+xyres@vpXF+yspace1
 xyres@tmYLOn        = False
 xyres@tmYLBorderOn  = False                        ; turn off right axis
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 180                ; bring bars down to zero
 xyres@gsnXYBarChartColors = (/"firebrick1"/)
 xyres@tmYLValues               = (/0,30,60,90,120,150,180/)
 xyres@tmYLLabels               = (/"0","30","60","90","120","150","180"/)
 xyres@gsnCenterString       = "SON"
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnCenterStringOrthogonalPosF = -0.10
 xyres@gsnLeftString      =""
 plot08  = gsn_csm_xy(wks,ispan(1,2,1),stdarr6(:),xyres)   ; advancing the frame

  polyres=True
  polyres@gsLineDashPattern=0
  polyres@gsLineColor="black"
  polyres@gsLineThicknessF=1.0
  polyres@tfPolyDrawOrder = "postDraw"

  naumm1    = gsn_add_polyline(wks,plot04,(/1,1/),(/xBootLow_sst01_PD,xBootHi_sst01_PD/), polyres)
  naumm001  = gsn_add_polyline(wks,plot04,(/2,2/),(/xBootLow_sst01_4xCO2,xBootHi_sst01_4xCO2/), polyres)

  nbumm1    = gsn_add_polyline(wks,plot06,(/1,1/),(/xBootLow_ssh01_PD,xBootHi_ssh01_PD/), polyres)
  nbumm001  = gsn_add_polyline(wks,plot06,(/2,2/),(/xBootLow_ssh01_4xCO2,xBootHi_ssh01_4xCO2/), polyres)

  ncumm1    = gsn_add_polyline(wks,plot08,(/1,1/),(/xBootLow_eke01_PD,xBootHi_eke01_PD/), polyres)
  ncumm001  = gsn_add_polyline(wks,plot08,(/2,2/),(/xBootLow_eke01_4xCO2,xBootHi_eke01_4xCO2/), polyres)

  naumm2    = gsn_add_polyline(wks,plot03,(/1,1/),(/xBootLow_sst02_PD,xBootHi_sst02_PD/), polyres)
  naumm002  = gsn_add_polyline(wks,plot03,(/2,2/),(/xBootLow_sst02_4xCO2,xBootHi_sst02_4xCO2/), polyres)

  nbumm2    = gsn_add_polyline(wks,plot05,(/1,1/),(/xBootLow_ssh02_PD,xBootHi_ssh02_PD/), polyres)
  nbumm002  = gsn_add_polyline(wks,plot05,(/2,2/),(/xBootLow_ssh02_4xCO2,xBootHi_ssh02_4xCO2/), polyres)

  ncumm2    = gsn_add_polyline(wks,plot07,(/1,1/),(/xBootLow_eke02_PD,xBootHi_eke02_PD/), polyres)
  ncumm002  = gsn_add_polyline(wks,plot07,(/2,2/),(/xBootLow_eke02_4xCO2,xBootHi_eke02_4xCO2/), polyres)
  
  wid = 0.24/2
  ndumm1   = gsn_add_polyline(wks,plot04,(/1-wid,1+wid/),(/xBootLow_sst01_PD,xBootLow_sst01_PD/), polyres)
  ndumm01  = gsn_add_polyline(wks,plot04,(/1-wid,1+wid/),(/xBootHi_sst01_PD,xBootHi_sst01_PD/), polyres)

  nfumm1   = gsn_add_polyline(wks,plot04,(/2-wid,2+wid/),(/xBootLow_sst01_4xCO2,xBootLow_sst01_4xCO2/), polyres)
  nfumm01  = gsn_add_polyline(wks,plot04,(/2-wid,2+wid/),(/xBootHi_sst01_4xCO2,xBootHi_sst01_4xCO2/), polyres)
;----------------------------------------------------------------------------------------------------
  ndumm2   = gsn_add_polyline(wks,plot06,(/1-wid,1+wid/),(/xBootLow_ssh01_PD,xBootLow_ssh01_PD/), polyres)
  ndumm02  = gsn_add_polyline(wks,plot06,(/1-wid,1+wid/),(/xBootHi_ssh01_PD,xBootHi_ssh01_PD/), polyres)

  nfumm2   = gsn_add_polyline(wks,plot06,(/2-wid,2+wid/),(/xBootLow_ssh01_4xCO2,xBootLow_ssh01_4xCO2/), polyres)
  nfumm02  = gsn_add_polyline(wks,plot06,(/2-wid,2+wid/),(/xBootHi_ssh01_4xCO2,xBootHi_ssh01_4xCO2/), polyres)
;-----------------------------------------------------------------------------------------------------
  ndumm3   = gsn_add_polyline(wks,plot08,(/1-wid,1+wid/),(/xBootLow_eke01_PD,xBootLow_eke01_PD/), polyres)
  ndumm03  = gsn_add_polyline(wks,plot08,(/1-wid,1+wid/),(/xBootHi_eke01_PD,xBootHi_eke01_PD/), polyres)

  nfumm3   = gsn_add_polyline(wks,plot08,(/2-wid,2+wid/),(/xBootLow_eke01_4xCO2,xBootLow_eke01_4xCO2/), polyres)
  nfumm03  = gsn_add_polyline(wks,plot08,(/2-wid,2+wid/),(/xBootHi_eke01_4xCO2,xBootHi_eke01_4xCO2/), polyres)
;------------------------------------------------------------------------------------------------------
  ndumm4   = gsn_add_polyline(wks,plot03,(/1-wid,1+wid/),(/xBootLow_sst02_PD,xBootLow_sst02_PD/), polyres)
  ndumm04  = gsn_add_polyline(wks,plot03,(/1-wid,1+wid/),(/xBootHi_sst02_PD,xBootHi_sst02_PD/), polyres)
  
  nfumm4   = gsn_add_polyline(wks,plot03,(/2-wid,2+wid/),(/xBootLow_sst02_4xCO2,xBootLow_sst02_4xCO2/), polyres)
  nfumm04  = gsn_add_polyline(wks,plot03,(/2-wid,2+wid/),(/xBootHi_sst02_4xCO2,xBootHi_sst02_4xCO2/), polyres)
;----------------------------------------------------------------------------------------------------
  ndumm5   = gsn_add_polyline(wks,plot05,(/1-wid,1+wid/),(/xBootLow_ssh02_PD,xBootLow_ssh02_PD/), polyres)
  ndumm05  = gsn_add_polyline(wks,plot05,(/1-wid,1+wid/),(/xBootHi_ssh02_PD,xBootHi_ssh02_PD/), polyres)
  
  nfumm5   = gsn_add_polyline(wks,plot05,(/2-wid,2+wid/),(/xBootLow_ssh02_4xCO2,xBootLow_ssh02_4xCO2/), polyres)
  nfumm05  = gsn_add_polyline(wks,plot05,(/2-wid,2+wid/),(/xBootHi_ssh02_4xCO2,xBootHi_ssh02_4xCO2/), polyres)
;-----------------------------------------------------------------------------------------------------
  ndumm6   = gsn_add_polyline(wks,plot07,(/1-wid,1+wid/),(/xBootLow_eke02_PD,xBootLow_eke02_PD/), polyres)
  ndumm06  = gsn_add_polyline(wks,plot07,(/1-wid,1+wid/),(/xBootHi_eke02_PD,xBootHi_eke02_PD/), polyres)
  
  nfumm6   = gsn_add_polyline(wks,plot07,(/2-wid,2+wid/),(/xBootLow_eke02_4xCO2,xBootLow_eke02_4xCO2/), polyres)
  nfumm06  = gsn_add_polyline(wks,plot07,(/2-wid,2+wid/),(/xBootHi_eke02_4xCO2,xBootHi_eke02_4xCO2/), polyres)
;==================================================================================================================
draw(plot03)
draw(plot04)
draw(plot05)
draw(plot06)
draw(plot07)
draw(plot08)


end 
