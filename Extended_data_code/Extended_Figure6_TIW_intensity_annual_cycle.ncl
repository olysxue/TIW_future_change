 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"   
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"   
 load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"   
 load "/Users/xueaoyun/Documents/work/libfun/my_utilities.ncl"
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


 undef("readdata")
function readdata(ncdir, filelist, varname, scale)
local dims,temp,f
begin
   f= addfile(ncdir + filelist(0), "r")
   temp=f->$varname$(:,0,0)
   dims=getfilevardimsizes(f, varname)
   var=new((/3,dims(0)/12,12/), float)
   
   do i=0,2
   f:= addfile(ncdir + filelist(i), "r")
   var(i,:,:)=reshape(f->$varname$(:,0,0), (/30,12/)) 
   end do 
   out=dim_avg_n_Wrap(var, 1)
   out1 = dim_stddev_n_Wrap(var, 1)
   out=out/scale
   out1 = out/scale 
   printMinMax(out(0,:),0)
   printMinMax(out(1,:),0)
   printMinMax(out(2,:),0)
   return(out)
end
;----------------------------------------------------
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

begin
;====================================================================================
ncdir  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/mean/"
; read SST var
fil1    = systemfunc("ls " + ncdir + "CESM2_HR_SST_std_var_10-50_fft_PD.nc")
fil2    = systemfunc("ls " + ncdir + "CESM2_HR_SST_std_var_10-50_fft_2xCO2.nc")
fil3    = systemfunc("ls " + ncdir + "CESM2_HR_SST_std_var_10-50_fft_4xCO2.nc")
  f1    = addfile(fil1, "r")
  f2    = addfile(fil2, "r")
  f3    = addfile(fil3, "r") 
  tiw_SST_PD     = f1->tiw_std(:,:,:)
  tiw_SST_2xCO2  = f2->tiw_std(:,:,:)
  tiw_SST_4xCO2  = f3->tiw_std(:,:,:)
  printVarSummary(tiw_SST_PD)

  tiw_SST_PD_ac0      = (tiw_SST_PD(8::12,:,:)+tiw_SST_PD(9::12,:,:)+tiw_SST_PD(10::12,:,:)-tiw_SST_PD(2::12,:,:)-tiw_SST_PD(3::12,:,:)-tiw_SST_PD(4::12,:,:))/6.
  tiw_SST_2xCO2_ac0   = (tiw_SST_2xCO2(8::12,:,:)+tiw_SST_2xCO2(9::12,:,:)+tiw_SST_2xCO2(10::12,:,:)-tiw_SST_2xCO2(2::12,:,:)-tiw_SST_2xCO2(3::12,:,:)-tiw_SST_2xCO2(4::12,:,:))/6.
  tiw_SST_4xCO2_ac0   = (tiw_SST_4xCO2(8::12,:,:)+tiw_SST_4xCO2(9::12,:,:)+tiw_SST_4xCO2(10::12,:,:)-tiw_SST_4xCO2(2::12,:,:)-tiw_SST_4xCO2(3::12,:,:)-tiw_SST_4xCO2(4::12,:,:))/6.
  copy_VarMeta(tiw_SST_PD(0::12,:,:), tiw_SST_PD_ac0)
  copy_VarMeta(tiw_SST_2xCO2(0::12,:,:), tiw_SST_2xCO2_ac0)
  copy_VarMeta(tiw_SST_4xCO2(0::12,:,:), tiw_SST_4xCO2_ac0)
  tiw_SST_PD_ac1    = clmMonTLL(tiw_SST_PD)
  tiw_SST_2xCO2_ac1 = clmMonTLL(tiw_SST_2xCO2)
  tiw_SST_4xCO2_ac1 = clmMonTLL(tiw_SST_4xCO2)
  printVarSummary(tiw_SST_PD_ac1)
  tiw_SST_PD_ac      = (tiw_SST_PD_ac1(8,:,:)+tiw_SST_PD_ac1(9,:,:)+tiw_SST_PD_ac1(10,:,:)-tiw_SST_PD_ac1(2,:,:)-tiw_SST_PD_ac1(3,:,:)-tiw_SST_PD_ac1(4,:,:))/6.
  tiw_SST_2xCO2_ac   = (tiw_SST_2xCO2_ac1(8,:,:)+tiw_SST_2xCO2_ac1(9,:,:)+tiw_SST_2xCO2_ac1(10,:,:)-tiw_SST_2xCO2_ac1(2,:,:)-tiw_SST_2xCO2_ac1(3,:,:)-tiw_SST_2xCO2_ac1(4,:,:))/6.
  tiw_SST_4xCO2_ac   = (tiw_SST_4xCO2_ac1(8,:,:)+tiw_SST_4xCO2_ac1(9,:,:)+tiw_SST_4xCO2_ac1(10,:,:)-tiw_SST_4xCO2_ac1(2,:,:)-tiw_SST_4xCO2_ac1(3,:,:)-tiw_SST_4xCO2_ac1(4,:,:))/6.
  
  copy_VarMeta(tiw_SST_PD(0,:,:), tiw_SST_PD_ac)
  copy_VarMeta(tiw_SST_2xCO2(0,:,:), tiw_SST_2xCO2_ac)
  copy_VarMeta(tiw_SST_4xCO2(0,:,:), tiw_SST_4xCO2_ac)
  ;====================================================================================
 ncdir  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/mean/"
 fil01    = systemfunc("ls " + ncdir + "CESM2_TIW_EKE_10-50_fft_surface_PD.nc")
 fil02    = systemfunc("ls " + ncdir + "CESM2_TIW_EKE_10-50_fft_surface_2xCO2.nc")
 fil03    = systemfunc("ls " + ncdir + "CESM2_TIW_EKE_10-50_fft_surface_4xCO2.nc")
   f01    = addfile(fil01, "r")
   f02    = addfile(fil02, "r")
   f03    = addfile(fil03, "r") 
   tiw_EKE0_PD     = f01->EKE(:,:,:)
   tiw_EKE0_2xCO2  = f02->EKE(:,:,:)
   tiw_EKE0_4xCO2  = f03->EKE(:,:,:)
   tiw_EKE_PD         = tiw_EKE0_PD/10000.
   tiw_EKE_2xCO2      = tiw_EKE0_2xCO2/10000.
   tiw_EKE_4xCO2      = tiw_EKE0_4xCO2/10000.
   copy_VarMeta(tiw_EKE0_PD, tiw_EKE_PD)
   copy_VarMeta(tiw_EKE0_2xCO2, tiw_EKE_2xCO2)
   copy_VarMeta(tiw_EKE0_4xCO2, tiw_EKE_4xCO2)
   
   tiw_EKE_PD_ac0      = (tiw_EKE_PD(8::12,:,:)   +tiw_EKE_PD(9::12,:,:)   +tiw_EKE_PD(10::12,:,:)   -tiw_EKE_PD(2::12,:,:)   -tiw_EKE_PD(3::12,:,:)   -tiw_EKE_PD(4::12,:,:))/6.
   tiw_EKE_2xCO2_ac0   = (tiw_EKE_2xCO2(8::12,:,:)+tiw_EKE_2xCO2(9::12,:,:)+tiw_EKE_2xCO2(10::12,:,:)-tiw_EKE_2xCO2(2::12,:,:)-tiw_EKE_2xCO2(3::12,:,:)-tiw_EKE_2xCO2(4::12,:,:))/6.
   tiw_EKE_4xCO2_ac0   = (tiw_EKE_4xCO2(8::12,:,:)+tiw_EKE_4xCO2(9::12,:,:)+tiw_EKE_4xCO2(10::12,:,:)-tiw_EKE_4xCO2(2::12,:,:)-tiw_EKE_4xCO2(3::12,:,:)-tiw_EKE_4xCO2(4::12,:,:))/6.
   copy_VarMeta(tiw_EKE_PD(0::12,:,:),    tiw_EKE_PD_ac0)
   copy_VarMeta(tiw_EKE_2xCO2(0::12,:,:), tiw_EKE_2xCO2_ac0)
   copy_VarMeta(tiw_EKE_4xCO2(0::12,:,:), tiw_EKE_4xCO2_ac0)
 ;===================================================================
   tiw_EKE_PD_ac1    = clmMonTLL(tiw_EKE_PD)
   tiw_EKE_2xCO2_ac1 = clmMonTLL(tiw_EKE_2xCO2)
   tiw_EKE_4xCO2_ac1 = clmMonTLL(tiw_EKE_4xCO2)
 
   tiw_EKE_PD_ac      = (tiw_EKE_PD_ac1(8,:,:)+tiw_EKE_PD_ac1(9,:,:)+tiw_EKE_PD_ac1(10,:,:) -tiw_EKE_PD_ac1(2,:,:)- tiw_EKE_PD_ac1(3,:,:)-tiw_EKE_PD_ac1(4,:,:))/6.
   tiw_EKE_2xCO2_ac   = (tiw_EKE_2xCO2_ac1(8,:,:)+tiw_EKE_2xCO2_ac1(9,:,:)+tiw_EKE_2xCO2_ac1(10,:,:) -tiw_EKE_2xCO2_ac1(2,:,:)- tiw_EKE_2xCO2_ac1(3,:,:)-tiw_EKE_2xCO2_ac1(4,:,:))/6.
   tiw_EKE_4xCO2_ac   = (tiw_EKE_4xCO2_ac1(8,:,:)+tiw_EKE_4xCO2_ac1(9,:,:)+tiw_EKE_4xCO2_ac1(10,:,:) -tiw_EKE_4xCO2_ac1(2,:,:)- tiw_EKE_4xCO2_ac1(3,:,:)-tiw_EKE_4xCO2_ac1(4,:,:))/6.
   copy_VarMeta(tiw_EKE_PD(0,:,:),    tiw_EKE_PD_ac)
   copy_VarMeta(tiw_EKE_2xCO2(0,:,:), tiw_EKE_2xCO2_ac)
   copy_VarMeta(tiw_EKE_4xCO2(0,:,:), tiw_EKE_4xCO2_ac)
  ;===================================================================
  ; read SSH 
  ncdir  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/CESM/mean/"
; read SSH var
fil001    = systemfunc("ls " + ncdir + "CESM2_HR_SSH_std_var_10-50_fft_PD.nc")
fil002    = systemfunc("ls " + ncdir + "CESM2_HR_SSH_std_var_10-50_fft_2xCO2.nc")
fil003    = systemfunc("ls " + ncdir + "CESM2_HR_SSH_std_var_10-50_fft_4xCO2.nc")
  f001    = addfile(fil001, "r")
  f002    = addfile(fil002, "r")
  f003    = addfile(fil003, "r") 
  tiw_SSH_PD     = f001->tiw_std(:,:,:)
  tiw_SSH_2xCO2  = f002->tiw_std(:,:,:)
  tiw_SSH_4xCO2  = f003->tiw_std(:,:,:)
  printVarSummary(tiw_SSH_PD)

  tiw_SSH_PD_ac0      = (tiw_SSH_PD(8::12,:,:)   +tiw_SSH_PD(9::12,:,:)   +tiw_SSH_PD(10::12,:,:)   -tiw_SSH_PD(2::12,:,:)   -tiw_SSH_PD(3::12,:,:)   -tiw_SSH_PD(4::12,:,:))/6.
  tiw_SSH_2xCO2_ac0   = (tiw_SSH_2xCO2(8::12,:,:)+tiw_SSH_2xCO2(9::12,:,:)+tiw_SSH_2xCO2(10::12,:,:)-tiw_SSH_2xCO2(2::12,:,:)-tiw_SSH_2xCO2(3::12,:,:)-tiw_SSH_2xCO2(4::12,:,:))/6.
  tiw_SSH_4xCO2_ac0   = (tiw_SSH_4xCO2(8::12,:,:)+tiw_SSH_4xCO2(9::12,:,:)+tiw_SSH_4xCO2(10::12,:,:)-tiw_SSH_4xCO2(2::12,:,:)-tiw_SSH_4xCO2(3::12,:,:)-tiw_SSH_4xCO2(4::12,:,:))/6.
  copy_VarMeta(tiw_SSH_PD(0::12,:,:),    tiw_SSH_PD_ac0)
  copy_VarMeta(tiw_SSH_2xCO2(0::12,:,:), tiw_SSH_2xCO2_ac0)
  copy_VarMeta(tiw_SSH_4xCO2(0::12,:,:), tiw_SSH_4xCO2_ac0)
;-------------------------------------------------------------
  tiw_SSH_PD_ac1    = clmMonTLL(tiw_SSH_PD)
  tiw_SSH_2xCO2_ac1 = clmMonTLL(tiw_SSH_2xCO2)
  tiw_SSH_4xCO2_ac1 = clmMonTLL(tiw_SSH_4xCO2)
  printVarSummary(tiw_SSH_PD_ac1)
  tiw_SSH_PD_ac      = (tiw_SSH_PD_ac1(8,:,:)   + tiw_SSH_PD_ac1(9,:,:)   + tiw_SSH_PD_ac1(10,:,:)   -tiw_SSH_PD_ac1(2,:,:)   -tiw_SSH_PD_ac1(3,:,:)   -tiw_SSH_PD_ac1(4,:,:))/6.
  tiw_SSH_2xCO2_ac   = (tiw_SSH_2xCO2_ac1(8,:,:)+ tiw_SSH_2xCO2_ac1(9,:,:)+ tiw_SSH_2xCO2_ac1(10,:,:)-tiw_SSH_2xCO2_ac1(2,:,:)-tiw_SSH_2xCO2_ac1(3,:,:)-tiw_SSH_2xCO2_ac1(4,:,:))/6.
  tiw_SSH_4xCO2_ac   = (tiw_SSH_4xCO2_ac1(8,:,:)+ tiw_SSH_4xCO2_ac1(9,:,:)+ tiw_SSH_4xCO2_ac1(10,:,:)-tiw_SSH_4xCO2_ac1(2,:,:)-tiw_SSH_4xCO2_ac1(3,:,:)-tiw_SSH_4xCO2_ac1(4,:,:))/6.
  
  copy_VarMeta(tiw_SSH_PD(0,:,:),    tiw_SSH_PD_ac)
  copy_VarMeta(tiw_SSH_2xCO2(0,:,:), tiw_SSH_2xCO2_ac)
  copy_VarMeta(tiw_SSH_4xCO2(0,:,:), tiw_SSH_4xCO2_ac)
  ;=========================SST==============================
  CESM_pattern1 = smth9_Wrap(tiw_SST_PD_ac, 0.50, 0.25, False)
  CESM_pattern2 = smth9_Wrap(tiw_SST_2xCO2_ac, 0.50, 0.25, False)
  CESM_pattern3 = smth9_Wrap(tiw_SST_4xCO2_ac, 0.50, 0.25, False)
  printMinMax(CESM_pattern1,0)
  printMinMax(CESM_pattern2,0)
  printMinMax(CESM_pattern3,0)

  CESM_pattern_sst = new((/2, dimsizes(CESM_pattern1&lat), dimsizes(CESM_pattern1&lon)/), typeof(CESM_pattern1))
  CESM_pattern_sst(0,:,:) = CESM_pattern1
  CESM_pattern_sst(1,:,:) = CESM_pattern3
;=========================EKE===============================
  CESM_pattern01 = smth9_Wrap(tiw_EKE_PD_ac, 0.50, 0.25, False)
  CESM_pattern02 = smth9_Wrap(tiw_EKE_2xCO2_ac, 0.50, 0.25, False)
  CESM_pattern03 = smth9_Wrap(tiw_EKE_4xCO2_ac, 0.50, 0.25, False)

  CESM_pattern_eke = new((/2, dimsizes(CESM_pattern01&lat), dimsizes(CESM_pattern01&lon)/), typeof(CESM_pattern01))
  CESM_pattern_eke(0,:,:) = CESM_pattern01
  CESM_pattern_eke(1,:,:) = CESM_pattern03
;=========================SSH===============================
  CESM_pattern001 = smth9_Wrap(tiw_SSH_PD_ac, 0.50, 0.25, False)
  CESM_pattern002 = smth9_Wrap(tiw_SSH_2xCO2_ac, 0.50, 0.25, False)
  CESM_pattern003 = smth9_Wrap(tiw_SSH_4xCO2_ac, 0.50, 0.25, False)

  CESM_pattern_ssh = new((/2, dimsizes(CESM_pattern001&lat), dimsizes(CESM_pattern001&lon)/), typeof(CESM_pattern001))
  CESM_pattern_ssh(0,:,:) = CESM_pattern001
  CESM_pattern_ssh(1,:,:) = CESM_pattern003
;========================================read E3SM======================================================
;read E3SM 
varname="tos"
ncdir  = "/Users/xueaoyun/Library/CloudStorage/Dropbox/TIW_future/data/E3SM/mean/"
scale=1. ;units degC
filelist=(/"tos_day_E3SM-1-0_picontrol_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mon.nc","tos_day_E3SM-1-0_4xCO2_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mon.nc"/)
sst0=readdata1(ncdir, filelist, varname, scale)
sst0@units="degC"
sst = sqrt(sst0)
copy_VarMeta(sst0, sst)
printVarSummary(sst)
delete(sst0)
varname="zos"
scale=1. ; units cm to m by scale
filelist=(/"zos_day_E3SM-1-0_picontrol_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mon.nc","zos_day_E3SM-1-0_4xCO2_gr_00010101-00831231_regrid_05deg_equator_10-50fft_amplitude_mon.nc"/)
ssh0=readdata1(ncdir, filelist, varname, scale)
ssh0@units="meter"
printVarSummary(ssh0)
ssh = sqrt(ssh0)*100.
copy_VarMeta(ssh0, ssh)
delete(ssh0)
scale=1.
varname="KE" ;units cm2/s2 to m2/s2 by scale
filelist:=(/"TIW_KE_mon_E3SM-1-0_picontrol_gr_00010101-00831231_regrid_1deg_equator_new.nc","TIW_KE_mon_E3SM-1-0_4xCO2_gr_00010101-00831231_regrid_1deg_equator_new.nc"/)
eke=readdata1(ncdir, filelist, varname, scale)
eke@units="J/m3"

SST_PD = sst(0,:,:,:)
SST_4xCO2 = sst(1,:,:,:)

SSH_PD = ssh(0,:,:,:)
SSH_4xCO2 = ssh(1,:,:,:)

EKE_PD  = eke(0,:,:,:)
EKE_4xCO2 = eke(1,:,:,:)
;====================================================================================================================
SST_PD_ac0      = (SST_PD(8::12,:,:)   + SST_PD(9::12,:,:)   + SST_PD(10::12,:,:)   - SST_PD(2::12,:,:)   - SST_PD(3::12,:,:)   - SST_PD(4::12,:,:))/6.
SST_4xCO2_ac0   = (SST_4xCO2(8::12,:,:)+ SST_4xCO2(9::12,:,:)+ SST_4xCO2(10::12,:,:)- SST_4xCO2(2::12,:,:)- SST_4xCO2(3::12,:,:)- SST_4xCO2(4::12,:,:))/6.
copy_VarMeta(SST_PD(0::12,:,:), SST_PD_ac0)
copy_VarMeta(SST_4xCO2(0::12,:,:), SST_4xCO2_ac0)
 
SST_PD_ac1    = clmMonTLL(SST_PD)
SST_4xCO2_ac1 = clmMonTLL(SST_4xCO2)
printVarSummary(SST_PD_ac1)

SST_PD_ac      = (SST_PD_ac1(8,:,:)   + SST_PD_ac1(9,:,:)   + SST_PD_ac1(10,:,:)   - SST_PD_ac1(2,:,:)   - SST_PD_ac1(3,:,:)   - SST_PD_ac1(4,:,:))/6.
SST_4xCO2_ac   = (SST_4xCO2_ac1(8,:,:)+ SST_4xCO2_ac1(9,:,:)+ SST_4xCO2_ac1(10,:,:)- SST_4xCO2_ac1(2,:,:)- SST_4xCO2_ac1(3,:,:)- SST_4xCO2_ac1(4,:,:))/6.
  
  copy_VarMeta(SST_PD(0,:,:),    SST_PD_ac)
  copy_VarMeta(SST_4xCO2(0,:,:), SST_4xCO2_ac)
;===================================================================================================================
SSH_PD_ac0      = (SSH_PD(8::12,:,:)   + SSH_PD(9::12,:,:)   + SSH_PD(10::12,:,:)   - SSH_PD(2::12,:,:)   - SSH_PD(3::12,:,:)   - SSH_PD(4::12,:,:))/6.
SSH_4xCO2_ac0   = (SSH_4xCO2(8::12,:,:)+ SSH_4xCO2(9::12,:,:)+ SSH_4xCO2(10::12,:,:)- SSH_4xCO2(2::12,:,:)- SSH_4xCO2(3::12,:,:)- SSH_4xCO2(4::12,:,:))/6.
copy_VarMeta(SSH_PD(0::12,:,:),    SSH_PD_ac0)
copy_VarMeta(SSH_4xCO2(0::12,:,:), SSH_4xCO2_ac0)
 
SSH_PD_ac1    = clmMonTLL(SSH_PD)
SSH_4xCO2_ac1 = clmMonTLL(SSH_4xCO2)
printVarSummary(SSH_PD_ac1)

SSH_PD_ac      = (SSH_PD_ac1(8,:,:)   + SSH_PD_ac1(9,:,:)   + SSH_PD_ac1(10,:,:)   - SSH_PD_ac1(2,:,:)   - SSH_PD_ac1(3,:,:)   - SSH_PD_ac1(4,:,:))/6.
SSH_4xCO2_ac   = (SSH_4xCO2_ac1(8,:,:)+ SSH_4xCO2_ac1(9,:,:)+ SSH_4xCO2_ac1(10,:,:)- SSH_4xCO2_ac1(2,:,:)- SSH_4xCO2_ac1(3,:,:)- SSH_4xCO2_ac1(4,:,:))/6.
  
copy_VarMeta(SSH_PD(0,:,:),   SSH_PD_ac)
copy_VarMeta(SSH_4xCO2(0,:,:),SSH_4xCO2_ac)
;==================================================================================================================

EKE_PD_ac0      = (EKE_PD(8::12,:,:)   +EKE_PD(9::12,:,:)   + EKE_PD(10::12,:,:)   - EKE_PD(2::12,:,:)   - EKE_PD(3::12,:,:)   - EKE_PD(4::12,:,:))/6.
EKE_4xCO2_ac0   = (EKE_4xCO2(8::12,:,:)+EKE_4xCO2(9::12,:,:)+ EKE_4xCO2(10::12,:,:)- EKE_4xCO2(2::12,:,:)- EKE_4xCO2(3::12,:,:)- EKE_4xCO2(4::12,:,:))/6.
copy_VarMeta(EKE_PD(0::12,:,:),    EKE_PD_ac0)
copy_VarMeta(EKE_4xCO2(0::12,:,:), EKE_4xCO2_ac0)

EKE_PD_ac1    = clmMonTLL(EKE_PD)
EKE_4xCO2_ac1 = clmMonTLL(EKE_4xCO2)
 
EKE_PD_ac      = (EKE_PD_ac1(8,:,:)   + EKE_PD_ac1(9,:,:)   + EKE_PD_ac1(10,:,:)    - EKE_PD_ac1(2,:,:)   - EKE_PD_ac1(3,:,:)   - EKE_PD_ac1(4,:,:))/6.
EKE_4xCO2_ac   = (EKE_4xCO2_ac1(8,:,:)+ EKE_4xCO2_ac1(9,:,:)+ EKE_4xCO2_ac1(10,:,:) - EKE_4xCO2_ac1(2,:,:)- EKE_4xCO2_ac1(3,:,:)- EKE_4xCO2_ac1(4,:,:))/6.
copy_VarMeta(EKE_PD(0,:,:),    EKE_PD_ac)
copy_VarMeta(EKE_4xCO2(0,:,:), EKE_4xCO2_ac)
;=========================SST==============================
  E3SM_pattern1 = smth9_Wrap(SST_PD_ac, 0.50, 0.25, False)
  E3SM_pattern2 = smth9_Wrap(SST_4xCO2_ac, 0.50, 0.25, False)

  E3SM_pattern_sst = new((/3, dimsizes(E3SM_pattern1&lat), dimsizes(E3SM_pattern1&lon)/), typeof(E3SM_pattern1))
  E3SM_pattern_sst(0,:,:) = E3SM_pattern1
  E3SM_pattern_sst(1,:,:) = E3SM_pattern2
;=========================EKE===============================
  E3SM_pattern01 = smth9_Wrap(EKE_PD_ac, 0.50, 0.25, False)
  E3SM_pattern02 = smth9_Wrap(EKE_4xCO2_ac, 0.50, 0.25, False)

  E3SM_pattern_eke = new((/2, dimsizes(E3SM_pattern01&lat), dimsizes(E3SM_pattern01&lon)/), typeof(E3SM_pattern01))
  E3SM_pattern_eke(0,:,:) = E3SM_pattern01
  E3SM_pattern_eke(1,:,:) = E3SM_pattern02
;=========================SSH===============================
  E3SM_pattern001 = smth9_Wrap(SSH_PD_ac, 0.50, 0.25, False)
  E3SM_pattern002 = smth9_Wrap(SSH_4xCO2_ac, 0.50, 0.25, False)
  
  E3SM_pattern_ssh = new((/3, dimsizes(E3SM_pattern001&lat), dimsizes(E3SM_pattern001&lon)/), typeof(E3SM_pattern001))
  E3SM_pattern_ssh(0,:,:) = E3SM_pattern001
  E3SM_pattern_ssh(1,:,:) = E3SM_pattern002
;------------------------------------------------------------------

wks= gsn_open_wks("pdf","/Users/xueaoyun/Documents/work/TIW_future/revise_figures/plot_TIW_intensity_PD_2xCO2_4xCO2_annual_cycle")
wks_type = "eps"
wks_type@wkWidth  = 0.8
wks_type@wkHeight = 0.2    
;--------------------------------------
; gsn_define_colormap(wks,"posneg_1")
gsn_define_colormap(wks,"BlueWhiteOrangeRed")

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
; res@tmYLValues = (/-10,-8,-6,-4,-2,0,2,4,6,8,10/)
; res@tmYLLabels = (/"10~S~o~N~S","8~S~o~N~S","6~S~o~N~S","4~S~o~N~S","2~S~o~N~S","EQ","2~S~o~N~N","4~S~o~N~N","6~S~o~N~N","8~S~o~N~N","10~S~o~N~N"/)
; res@tmYLMinorValues         = ispan(-10,10,2) 
res@tmYLValues = (/-10,-5,0,5,10/)
res@tmYLLabels = (/"10~S~o~N~S","5~S~o~N~S","EQ","5~S~o~N~N","10~S~o~N~N"/)
res@tmYLMinorValues         = ispan(-10,10,5) 
res@tmXBMode = "Explicit"
res@tmXBValues = (/180,200,220,240,260,280/)
res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W",""/)  
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
; res@lbOrientation        = "Vertical"
res@lbBoxEndCapStyle     = "TriangleBothEnds"
res@pmLabelBarWidthF     = 0.20
; res@pmLabelBarHeightF    = 0.55
res@pmLabelBarOrthogonalPosF = 0.19
; res@pmLabelBarParallelPosF   = 0.6
res@lbLabelFontHeightF       =0.008
res@gsnLeftString        = ""
res@gsnCenterStringFontHeightF  = 0.014

res@gsnSpreadColorEnd  =230

res@vpHeightF            = 0.12
res@vpWidthF             = 0.20

res@lbTitlePosition  = "Bottom"
res@lbTitleDirection = "Across" 
res@lbTitleFontHeightF= 0.010
res@lbTitleOffsetF    = 0.05

plot=new(18,graphic)
dum=plot

; ----->1 plot eke 
panels = (/"~F22~a", "~F22~c", "~F22~e", "~F22~b", "~F22~d", "~F22~f", "~F22~g", "~F22~i","~F22~k","~F22~h","~F22~j", "~F22~l","~F22~n","~F22~o","","~F22~p"/)
text=(/"PD","4xCO2"/)
yspace=0.02
yspace1 = 0.03
yf=0.9
res@tiDeltaF=3
res@gsnCenterStringOrthogonalPosF=0.10
res@cnLevels                = fspan(-40,40,17)
res@vpXF    = 0.08
res@lbTitleString    = "J/m~S~3~N~"
do i=0,1
res@vpYF    = yf-i*yspace-i*res@vpHeightF
res@tiYAxisString=text(i)
res@tiYAxisOffsetXF=-0.02
res@tiYAxisFontHeightF=0.014
if i.eq.0 then
res@gsnCenterString = "CESM EKE"
res@tiMainString=""
else
res@gsnCenterString = ""
res@tiMainString=""
end if

if i.eq.1 then
res@lbLabelBarOn    = True
res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W",""/)  
else
res@lbLabelBarOn    = False
res@tmXBLabels = (/"","","","","",""/)
end if
res@gsnLeftString  = panels(i*3)
res@gsnLeftStringFontHeightF = 0.015
res@gsnLeftStringOrthogonalPosF = 0.01
res@gsnLeftStringParallelPosF  = 0.03
plot(i*3)  = gsn_csm_contour_map(wks,CESM_pattern_eke(i,:,:),res)
end do
delete([/res@tiYAxisString,res@tiYAxisOffsetXF,res@tiYAxisFontHeightF/])

; ----->2 plot sst
res@lbTitleString    = "degC"
delete(res@cnLevels)
res@vpXF    = res@vpWidthF+res@vpXF+yspace1
res@cnLevels             = fspan(-0.3,0.3,31)
res@cnFillPalette         = "WhiteBlueGreenYellowOrangeRed"
; res@cnFillColors         = (/0,12,13,14,15,17,19/)
res@tmYLLabels = (/"","","","",""/)
do i=0,1
res@vpYF    = yf-i*yspace-i*res@vpHeightF
if i.eq.0 then
res@gsnCenterString = "CESM SST"
else
res@gsnCenterString = ""
end if 
if i.eq.1 then
res@lbLabelBarOn    = True
res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W",""/)  
else
res@lbLabelBarOn    = False
res@tmXBLabels = (/"","","","","",""/)
end if
res@gsnLeftString  = panels(1+i*3)
res@gsnLeftStringFontHeightF = 0.015
res@gsnLeftStringOrthogonalPosF = 0.01
res@gsnLeftStringParallelPosF  = 0.03
plot(1+i*3)  = gsn_csm_contour_map(wks,CESM_pattern_sst(i,:,:),res)
end do
; ----->3 plot ssh 
delete(res@cnLevels)
res@vpXF    = res@vpWidthF+res@vpXF+yspace1
res@lbTitleString    = "10~S~-2~N~ m"
res@cnLevels              = fspan(-2,2,21)
res@tmYLLabels = (/"","","","",""/)
do i=0,1
res@vpYF    = yf-i*yspace-i*res@vpHeightF
if i.eq.0 then
res@gsnCenterString = "CESM SSH"
else
res@gsnCenterString = ""
end if
if i.eq.1 then
res@lbLabelBarOn    = True
res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W",""/)  
else
res@lbLabelBarOn    = False
res@tmXBLabels = (/"","","","","",""/)
end if
res@gsnLeftString  = panels(2+i*3)
res@gsnLeftStringFontHeightF = 0.015
res@gsnLeftStringOrthogonalPosF = 0.01
plot(2+i*3)  = gsn_csm_contour_map(wks,CESM_pattern_ssh(i,:,:),res)
end do
;--------------------------------------------------------------------------------------------------------------
;==========================================E3SM================================================================
delete(res@cnLevels)
yf1 = res@vpYF-res@vpHeightF-0.12 
res@tmYLValues = (/-10,-5,0,5,10/)
res@tmYLLabels = (/"10~S~o~N~S","5~S~o~N~S","EQ","5~S~o~N~N","10~S~o~N~N"/)
; ----->1 plot eke
res@vpXF    = 0.08
res@gsnSpreadColorStart= 0
res@gsnSpreadColorEnd  =230
res@cnLevels        =fspan(-50,50,51)
do i=0,1
res@vpYF    = yf1-i*yspace-i*res@vpHeightF
print(res@vpYF)
res@tiYAxisString=text(i)
res@gsnLeftString  = panels(3+i*3+3)
res@tiYAxisOffsetXF=-0.02
res@tiYAxisFontHeightF=0.014
if i.eq.0 then
res@gsnCenterString = "E3SM EKE"
res@tiMainString=""
res@tmXBLabels = (/"","","","","",""/) 
else
res@gsnCenterString = ""
res@tiMainString=""
end if

if i.eq.1 then
res@lbLabelBarOn    = True
res@lbTitleString    = "J/m~S~3~N~"

res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W",""/)  
else
res@lbLabelBarOn    = False
end if

plot(3+i*3+3)  = gsn_csm_contour_map(wks,E3SM_pattern_eke(i,:,:),res)
end do
delete([/res@cnLevels,res@cnFillColors/])
delete([/res@tiYAxisString,res@tiYAxisOffsetXF,res@tiYAxisFontHeightF/])
;============================================================================
; ----->2 plot sst
res@vpXF    = res@vpWidthF+res@vpXF+yspace1
res@cnLevels                = fspan(-0.3,0.3,31)

;res@cnFillColors         =  (/0,144,168,178,188,204,230,246/)
res@lbTitleString    = ""
; res@cnFillColors         = (/0,12,13,14,15,17,19/)
res@tmYLLabels = (/"","","","",""/)
do i=0,1
res@vpYF    = yf1-i*yspace-i*res@vpHeightF
res@gsnLeftString  = panels(4+i*3+3)
if i.eq.0 then
res@gsnCenterString = "E3SM SST"
res@tmXBLabels = (/"","","","","",""/) 
else
res@gsnCenterString = ""
end if
if i.eq.1 then
res@lbLabelBarOn    = True
res@lbTitleString    = "degC"

res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W",""/)  
else
res@lbLabelBarOn    = False
end if
plot(4+i*3+3)  = gsn_csm_contour_map(wks,E3SM_pattern_sst(i,:,:),res)
end do
delete([/res@cnLevels,res@cnFillColors/])
; ----->3 plot ssh
res@vpXF    = res@vpWidthF+res@vpXF+yspace1
res@tmYLLabels = (/"","","","",""/)
res@cnLevels        =fspan(-1,1,21)
do i=0,1
res@vpYF    = yf1-i*yspace-i*res@vpHeightF
res@gsnLeftString  = panels(5+i*3+3)
if i.eq.0 then
res@gsnCenterString = "E3SM SSH"
res@tmXBLabels = (/"","","","","",""/) 
else
res@gsnCenterString = ""
end if
if i.eq.1 then
res@lbLabelBarOn    = True
res@lbTitleString    = "10~S~-2~N~ m"
res@tmXBLabels = (/"","160~S~o~N~W","140~S~o~N~W","120~S~o~N~W","100~S~o~N~W",""/)  
else
res@lbLabelBarOn    = False
res@tmXBLabels = (/"","","","","",""/)
end if
plot(5+i*3+3)  = gsn_csm_contour_map(wks,E3SM_pattern_ssh(i,:,:),res)
end do
delete([/res@cnLevels,res@cnFillColors/])
;=========================================================================
plres=True
plres@gsLineThicknessF=2.5
; plres@gsFillColor     = -1
plres@gsLineDashPattern= 2
plres@gsLineColor="blue"
plres@gsLineOpacityF=100
w1=360-160
e1=360-100

w2=360-160
e2=360-110

w3=360-160
e3=360-110

w4=360-150
e4=360-100
s1=0
n1=5

s2 = 3
n2 = 8

s3 = 3
n3 = 6

dum(0)=gsn_add_polyline(wks, plot(0), (/w2,e2,e2,w2,w2/), (/s1,s1,n1,n1,s1/), plres)
dum(2)=gsn_add_polyline(wks, plot(1), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum(1)=gsn_add_polyline(wks, plot(2), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
dum(3)=gsn_add_polyline(wks, plot(3), (/w2,e2,e2,w2,w2/), (/s1,s1,n1,n1,s1/), plres)
dum(5)=gsn_add_polyline(wks, plot(4), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum(4)=gsn_add_polyline(wks, plot(5), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)


dum(6)=gsn_add_polyline(wks, plot(6), (/w4,e4,e4,w4,w4/), (/s1,s1,n1,n1,s1/), plres)
dum(7)=gsn_add_polyline(wks, plot(8), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
dum(8)=gsn_add_polyline(wks, plot(7), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)
dum(9) =gsn_add_polyline(wks,plot(9), (/w3,e3,e3,w3,w3/), (/s1,s1,n1,n1,s1/), plres)
dum(10)=gsn_add_polyline(wks, plot(11), (/w2,e2,e2,w2,w2/), (/s2,s2,n2,n2,s2/), plres)
dum(11)=gsn_add_polyline(wks, plot(10), (/w1,e1,e1,w1,w1/), (/s1,s1,n1,n1,s1/), plres)


draw(plot(0))
draw(plot(1))
draw(plot(2))
draw(plot(3))
draw(plot(4))
draw(plot(5))
draw(plot(6))
draw(plot(7))
draw(plot(8))
draw(plot(9))
draw(plot(10))
draw(plot(11))
pos=new((/3,2/), float)
do i=0,2
pos(i,:) = get_title_position(plot(i*4))
end do
delete([/sst, ssh, eke/])
;=============================================================
tiw_sst_PD    = wgt_areaave_Wrap(tiw_SST_PD_ac0(:,{0:5},{210:250}), 1.0, 1.0, 0)
tiw_sst_2xCO2 = wgt_areaave_Wrap(tiw_SST_2xCO2_ac0(:,{0:5},{210:250}), 1.0, 1.0, 0)
tiw_sst_4xCO2 = wgt_areaave_Wrap(tiw_SST_4xCO2_ac0(:,{0:5},{210:250}), 1.0, 1.0, 0)

tiw_ssh_PD    = wgt_areaave_Wrap(tiw_SSH_PD_ac0(:,{3:8},{210:260}), 1.0, 1.0, 0)
tiw_ssh_2xCO2 = wgt_areaave_Wrap(tiw_SSH_2xCO2_ac0(:,{3:8},{210:260}), 1.0, 1.0, 0)
tiw_ssh_4xCO2 = wgt_areaave_Wrap(tiw_SSH_4xCO2_ac0(:,{3:8},{210:260}), 1.0, 1.0, 0)

tiw_eke_PD    = wgt_areaave_Wrap(tiw_EKE_PD_ac0(:,{0:5},{210:250}), 1.0, 1.0, 0)
tiw_eke_2xCO2 = wgt_areaave_Wrap(tiw_EKE_2xCO2_ac0(:,{0:5},{210:250}), 1.0, 1.0, 0)
tiw_eke_4xCO2 = wgt_areaave_Wrap(tiw_EKE_4xCO2_ac0(:,{0:5},{210:250}), 1.0, 1.0, 0)

stat  = 0                         ; 0=mean; 1=variance; 2=std. dev; ...
n     = 30                        ; sub-sample size;  default is n=N
nBoot = 10000                      ; bootstrap replications
N     = 30
opt   = False                     ; use defaults
if (isvar("n") .and. (n.lt.N)) then
    opt   = True
    opt@sample_size = n           ; default N
end if

 BootStrap_sst_PD     = bootstrap_stat(tiw_sst_PD, stat, nBoot, 0, opt)
 BootStrap_sst_2xCO2  = bootstrap_stat(tiw_sst_2xCO2, stat, nBoot, 0, opt)
 BootStrap_sst_4xCO2  = bootstrap_stat(tiw_sst_4xCO2, stat, nBoot, 0, opt)

 BootStrap_ssh_PD     = bootstrap_stat(tiw_ssh_PD, stat, nBoot, 0, opt)
 BootStrap_ssh_2xCO2  = bootstrap_stat(tiw_ssh_2xCO2, stat, nBoot, 0, opt)
 BootStrap_ssh_4xCO2  = bootstrap_stat(tiw_ssh_4xCO2, stat, nBoot, 0, opt)

 BootStrap_eke_PD     = bootstrap_stat(tiw_eke_PD, stat, nBoot, 0, opt)
 BootStrap_eke_2xCO2  = bootstrap_stat(tiw_eke_2xCO2, stat, nBoot, 0, opt)
 BootStrap_eke_4xCO2  = bootstrap_stat(tiw_eke_4xCO2, stat, nBoot, 0, opt)

 xBoot_sst_PD       = BootStrap_sst_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sstavg_PD    = BootStrap_sst_PD[1]       ; Mean of bootstrapped estimates
 xBoot_sststd_PD    = BootStrap_sst_PD[2]       ; Std dev of bootstrapped estimates

 xBoot_ssh_PD       = BootStrap_ssh_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sshavg_PD    = BootStrap_ssh_PD[1]       ; Mean of bootstrapped estimates
 xBoot_sshstd_PD    = BootStrap_ssh_PD[2]       ; Std dev of bootstrapped estimates

 xBoot_eke_PD       = BootStrap_eke_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_ekeavg_PD    = BootStrap_eke_PD[1]       ; Mean of bootstrapped estimates
 xBoot_ekestd_PD    = BootStrap_eke_PD[2]       ; Std dev of bootstrapped estimates
;=====================================================
xBoot_sst_2xCO2       = BootStrap_sst_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_sstavg_2xCO2    = BootStrap_sst_2xCO2[1]       ; Mean of bootstrapped estimates
xBoot_sststd_2xCO2    = BootStrap_sst_2xCO2[2]       ; Std dev of bootstrapped estimates

   xBoot_ssh_2xCO2    = BootStrap_ssh_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_sshavg_2xCO2    = BootStrap_ssh_2xCO2[1]       ; Mean of bootstrapped estimates
xBoot_sshstd_2xCO2    = BootStrap_ssh_2xCO2[2]       ; Std dev of bootstrapped estimates

   xBoot_eke_2xCO2    = BootStrap_eke_2xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
xBoot_ekeavg_2xCO2    = BootStrap_eke_2xCO2[1]       ; Mean of bootstrapped estimates
xBoot_ekestd_2xCO2    = BootStrap_eke_2xCO2[2]       ; Std dev of bootstrapped estimates
;-----------------------------------------------------
    xBoot_sst_4xCO2    = BootStrap_sst_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sstavg_4xCO2    = BootStrap_sst_4xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_sststd_4xCO2    = BootStrap_sst_4xCO2[2]       ; Std dev of bootstrapped estimates

    xBoot_ssh_4xCO2    = BootStrap_ssh_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sshavg_4xCO2    = BootStrap_ssh_4xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_sshstd_4xCO2    = BootStrap_ssh_4xCO2[2]       ; Std dev of bootstrapped estimates

    xBoot_eke_4xCO2    = BootStrap_eke_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_ekeavg_4xCO2    = BootStrap_eke_4xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_ekestd_4xCO2    = BootStrap_eke_4xCO2[2]       ; Std dev of bootstrapped estimates
;=====================================================
 xBootLow_sst_PD    = bootstrap_estimate(xBoot_sst_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_sst_PD    = bootstrap_estimate(xBoot_sst_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_sst_PD    = bootstrap_estimate(xBoot_sst_PD, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_ssh_PD    = bootstrap_estimate(xBoot_ssh_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_ssh_PD    = bootstrap_estimate(xBoot_ssh_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_ssh_PD    = bootstrap_estimate(xBoot_ssh_PD, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_eke_PD    = bootstrap_estimate(xBoot_eke_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_eke_PD    = bootstrap_estimate(xBoot_eke_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_eke_PD    = bootstrap_estimate(xBoot_eke_PD, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------
xBootLow_sst_2xCO2    = bootstrap_estimate(xBoot_sst_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_sst_2xCO2    = bootstrap_estimate(xBoot_sst_2xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_sst_2xCO2    = bootstrap_estimate(xBoot_sst_2xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_ssh_2xCO2    = bootstrap_estimate(xBoot_ssh_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_ssh_2xCO2    = bootstrap_estimate(xBoot_ssh_2xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_ssh_2xCO2    = bootstrap_estimate(xBoot_ssh_2xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_eke_2xCO2    = bootstrap_estimate(xBoot_eke_2xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_eke_2xCO2    = bootstrap_estimate(xBoot_eke_2xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_eke_2xCO2    = bootstrap_estimate(xBoot_eke_2xCO2, 0.975, False)   ; 97.5% upper confidence bound
;----------------------------------------------------------------------
xBootLow_sst_4xCO2    = bootstrap_estimate(xBoot_sst_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_sst_4xCO2    = bootstrap_estimate(xBoot_sst_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_sst_4xCO2    = bootstrap_estimate(xBoot_sst_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_ssh_4xCO2    = bootstrap_estimate(xBoot_ssh_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_ssh_4xCO2    = bootstrap_estimate(xBoot_ssh_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_ssh_4xCO2    = bootstrap_estimate(xBoot_ssh_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_eke_4xCO2    = bootstrap_estimate(xBoot_eke_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_eke_4xCO2    = bootstrap_estimate(xBoot_eke_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_eke_4xCO2    = bootstrap_estimate(xBoot_eke_4xCO2, 0.975, False)   ; 97.5% upper confidence bound
;---------------------------------------------------------------------------------------------------------
;stdarr1 = new((/2/),"float") ;  SST
;stdarr1(0) =  xBoot_sstavg_PD
;stdarr1(1) =  xBoot_sstavg_4xCO2
;
;stdarr2 = new((/2/),"float") ; SSH
;stdarr2(0) =  xBoot_sshavg_PD
;stdarr2(1) =  xBoot_sshavg_4xCO2
;
;stdarr3 = new((/2/),"float") ; EKE 
;stdarr3(0) =  xBoot_ekeavg_PD
;stdarr3(1) =  xBoot_ekeavg_4xCO2 
;----------------------------------
stdarr1 = new((/3/),"float") ;  SST
stdarr1(0) =  xBoot_sstavg_PD
stdarr1(1) =  xBoot_sstavg_2xCO2 
stdarr1(2) =  xBoot_sstavg_4xCO2

stdarr2 = new((/3/),"float") ; SSH
stdarr2(0) =  xBoot_sshavg_PD
stdarr2(1) =  xBoot_sshavg_2xCO2 
stdarr2(2) =  xBoot_sshavg_4xCO2

stdarr3 = new((/3/),"float") ; EKE 
stdarr3(0) =  xBoot_ekeavg_PD
stdarr3(1) =  xBoot_ekeavg_2xCO2 
stdarr3(2) =  xBoot_ekeavg_4xCO2 
;=========================================================
delete(filelist)
; E3SM 
sst1_PD    = wgt_areaave_Wrap(SST_PD_ac0(:,{0:5},{200:260}), 1.0, 1.0, 0)
sst1_4xCO2 = wgt_areaave_Wrap(SST_4xCO2_ac0(:,{0:5},{200:260}), 1.0, 1.0, 0)

ssh1_PD    = wgt_areaave_Wrap(SSH_PD_ac0(:,{3:8},{220:260}), 1.0, 1.0, 0)
ssh1_4xCO2 = wgt_areaave_Wrap(SSH_4xCO2_ac0(:,{3:8},{200:240}), 1.0, 1.0, 0)

eke1_PD    = wgt_areaave_Wrap(EKE_PD_ac0(:,{0:5},{210:250}), 1.0, 1.0, 0)
eke1_4xCO2 = wgt_areaave_Wrap(EKE_4xCO2_ac0(:,{0:5},{200:250}), 1.0, 1.0, 0)

stat  = 0                         ; 0=mean; 1=variance; 2=std. dev; ...
n     = 80                        ; sub-sample size;  default is n=N
nBoot = 10000                      ; bootstrap replications
N     = 80
opt   = False                     ; use defaults
if (isvar("n") .and. (n.lt.N)) then
    opt   = True
    opt@sample_size = n           ; default N
end if

 BootStrap_sst1_PD     = bootstrap_stat(sst1_PD, stat, nBoot, 0, opt)
 BootStrap_sst1_4xCO2  = bootstrap_stat(sst1_4xCO2, stat, nBoot, 0, opt)

 BootStrap_ssh1_PD     = bootstrap_stat(ssh1_PD, stat, nBoot, 0, opt)
 BootStrap_ssh1_4xCO2  = bootstrap_stat(ssh1_4xCO2, stat, nBoot, 0, opt)

 BootStrap_eke1_PD     = bootstrap_stat(eke1_PD, stat, nBoot, 0, opt)
 BootStrap_eke1_4xCO2  = bootstrap_stat(eke1_4xCO2, stat, nBoot, 0, opt)

 xBoot_sst1_PD       = BootStrap_sst1_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sstavg1_PD    = BootStrap_sst1_PD[1]       ; Mean of bootstrapped estimates
 xBoot_sststd1_PD    = BootStrap_sst1_PD[2]       ; Std dev of bootstrapped estimates

 xBoot_ssh1_PD       = BootStrap_ssh1_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sshavg1_PD    = BootStrap_ssh1_PD[1]       ; Mean of bootstrapped estimates
 xBoot_sshstd1_PD    = BootStrap_ssh1_PD[2]       ; Std dev of bootstrapped estimates

 xBoot_eke1_PD       = BootStrap_eke1_PD[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_ekeavg1_PD    = BootStrap_eke1_PD[1]       ; Mean of bootstrapped estimates
 xBoot_ekestd1_PD    = BootStrap_eke1_PD[2]       ; Std dev of bootstrapped estimates
;=====================================================
    xBoot_sst1_4xCO2    = BootStrap_sst1_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sstavg1_4xCO2    = BootStrap_sst1_4xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_sststd1_4xCO2    = BootStrap_sst1_4xCO2[2]       ; Std dev of bootstrapped estimates

    xBoot_ssh1_4xCO2    = BootStrap_ssh1_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_sshavg1_4xCO2    = BootStrap_ssh1_4xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_sshstd1_4xCO2    = BootStrap_ssh1_4xCO2[2]       ; Std dev of bootstrapped estimates

    xBoot_eke1_4xCO2    = BootStrap_eke1_4xCO2[0]       ; Bootstrapped values in ascending order (use for histogram)
 xBoot_ekeavg1_4xCO2    = BootStrap_eke1_4xCO2[1]       ; Mean of bootstrapped estimates
 xBoot_ekestd1_4xCO2    = BootStrap_eke1_4xCO2[2]       ; Std dev of bootstrapped estimates
;=====================================================
 xBootLow_sst1_PD    = bootstrap_estimate(xBoot_sst1_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_sst1_PD    = bootstrap_estimate(xBoot_sst1_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_sst1_PD    = bootstrap_estimate(xBoot_sst1_PD, 0.975, False)   ; 97.5% upper confidence bound

 xBootLow_ssh1_PD    = bootstrap_estimate(xBoot_ssh1_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_ssh1_PD    = bootstrap_estimate(xBoot_ssh1_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_ssh1_PD    = bootstrap_estimate(xBoot_ssh1_PD, 0.975, False)   ; 97.5% upper confidence bound
 
 xBootLow_eke1_PD    = bootstrap_estimate(xBoot_eke1_PD, 0.025, False)   ;  2.5% lower confidence bound
 xBootMed_eke1_PD    = bootstrap_estimate(xBoot_eke1_PD, 0.500, False)   ; 50.0% median of bootstrapped estimates
  xBootHi_eke1_PD    = bootstrap_estimate(xBoot_eke1_PD, 0.975, False)   ; 97.5% upper confidence bound
;-----------------------------------------------------
xBootLow_sst1_4xCO2    = bootstrap_estimate(xBoot_sst1_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_sst1_4xCO2    = bootstrap_estimate(xBoot_sst1_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_sst1_4xCO2    = bootstrap_estimate(xBoot_sst1_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_ssh1_4xCO2    = bootstrap_estimate(xBoot_ssh1_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_ssh1_4xCO2    = bootstrap_estimate(xBoot_ssh1_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_ssh1_4xCO2    = bootstrap_estimate(xBoot_ssh1_4xCO2, 0.975, False)   ; 97.5% upper confidence bound

xBootLow_eke1_4xCO2    = bootstrap_estimate(xBoot_eke1_4xCO2, 0.025, False)   ;  2.5% lower confidence bound
xBootMed_eke1_4xCO2    = bootstrap_estimate(xBoot_eke1_4xCO2, 0.500, False)   ; 50.0% median of bootstrapped estimates
 xBootHi_eke1_4xCO2    = bootstrap_estimate(xBoot_eke1_4xCO2, 0.975, False)   ; 97.5% upper confidence bound
;---------------------------------------------------------------------------------------------------------
stdarr01 = new((/2/),typeof(sst1_PD)) ;  SST
stdarr01(0) =  xBoot_sstavg1_PD
stdarr01(1) =  xBoot_sstavg1_4xCO2
printMinMax(stdarr01, 0)
stdarr02 = new((/2/),typeof(ssh1_PD)) ; SSH
stdarr02(0) =  xBoot_sshavg1_PD
stdarr02(1) =  xBoot_sshavg1_4xCO2

stdarr03 = new((/2/),typeof(eke1_PD)) ; EKE 
stdarr03(0) =  xBoot_ekeavg1_PD
stdarr03(1) =  xBoot_ekeavg1_4xCO2 
;===================================================================================
 xyres                       = True                 ; plot mods desired
 xyres@gsnFrame              = False                ; don't advance frame yet
 xyres@gsnDraw               = False
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 50.                 ; bring bars down to zero
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 4.                   ; of the 1st and last bars
 xyres@vpWidthF              = 0.10
 xyres@vpHeightF             = 0.12

 xyres@gsnXYBarChart         = True                 ; turn on bar chart
 xyres@gsnXYBarChartBarWidth = 0.24                  ; change bar widths
 xyres@tmXBOn                = True                 ; turn off tickmarks at bot
 xyres@tmYLOn                = True
 xyres@tmYROn                = False                 ; turn off tickmarks at bot
 xyres@tmXTOn                = False                 ; turn off tickmarks at bot
 xyres@tmXBMode          = "Explicit"         ; explicit labels
 xyres@tmXBValues        = (/1,2,3,4/)
 xyres@tmXBLabels = (/"PD","2CO2","4CO2",""/)
 xyres@tmXBLabelFontHeightF = 0.006
 xyres@tmYLMode          = "Explicit"         ; explicit labels
 xyres@tmYLValues        = (/0,10,20,30,40,50/)
 xyres@tmYLLabels        = (/"0","10","20","30","40","50"/)
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
 xyres@tmXBLabelsOn  = True 
 xyres@vpXF    = 0.08
 yspace=0.03
 yf2= res@vpYF- res@vpWidthF-0.02
 xyres@vpYF    = yf2
 xyres@gsnXYBarChartColors = (/"orange3"/)
 xyres@gsnCenterString       = "CESM_EKE"
 xyres@gsnCenterStringOrthogonalPosF = -0.125
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString  = "~F22~m"
 xyres@gsnLeftStringFontHeightF = 0.015
 xyres@gsnLeftStringOrthogonalPosF = 0.
 xyres@gsnLeftStringParallelPosF = -0.05
 plot(12)  = gsn_csm_xy(wks,ispan(1,3,1),stdarr3(:),xyres)   ; advancing the frame
 delete(xyres@tmXBValues)
 delete(xyres@tmXBLabels)
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 3.                   ; of the 1st and last bars
 xyres@tmXBValues        = (/1,2,3/)
 xyres@tmXBLabels = (/"PD","4CO2","" /)
 xyres@tmYLOn        = False 
 xyres@tmYRBorderOn  = True                        ; turn off right axis
 xyres@tmYLBorderOn  = False
 xyres@vpXF = xyres@vpWidthF+xyres@vpXF
 xyres@gsnLeftString  = ""
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 50.                ; bring bars down to zero
 xyres@gsnCenterString       = "E3SM_EKE"
 xyres@gsnCenterStringOrthogonalPosF = -0.10
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftStringFontHeightF = 0.015
 xyres@gsnLeftStringOrthogonalPosF = 0.
 plot(13)  = gsn_csm_xy(wks,ispan(1,2,1),stdarr03(:),xyres)   ; advancing the frame
 xyres@vpXF = xyres@vpXF+xyres@vpWidthF+yspace
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 delete(xyres@tmXBValues)
 delete(xyres@tmXBLabels)
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 4.                   ; of the 1st and last bars
 xyres@tmXBValues        = (/1,2,3,4/)
 xyres@tmXBLabels = (/"PD","2CO2","4CO2","" /)
 xyres@tmYLOn        = True
 xyres@tmYRBorderOn  = False                        ; turn off right axis
 xyres@tmYLBorderOn  = True 
 xyres@vpWidthF              = 0.10
 xyres@vpHeightF             = 0.12
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.3                 ; bring bars down to zero                   
 xyres@tmYLOn        = True 
 xyres@tmYLMode          = "Explicit"         ; explicit labels
 xyres@tmYLValues        = (/0,0.1,0.2,0.3,0.4,0.5/)
 xyres@tmYLLabels        = (/"0","0.1","0.2","0.3","0.4","0.5"/)
 xyres@gsnXYBarChartColors = (/"cyan3"/)
 xyres@gsnCenterString       = "CESM_SST"
 xyres@gsnCenterStringOrthogonalPosF = -0.125
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString  = "~F22~n"
 xyres@gsnLeftStringFontHeightF = 0.015
 xyres@gsnLeftStringOrthogonalPosF = 0.
 xyres@gsnLeftStringParallelPosF = -0.05
 plot(14)  = gsn_csm_xy(wks,ispan(1,3,1),stdarr1(:),xyres)   ; advancing the frame
 xyres@tmYLOn        = False 
 xyres@tmYRBorderOn  = True                        ; turn off right axis
 xyres@tmYLBorderOn  = False
 delete(xyres@tmXBValues)
 delete(xyres@tmXBLabels)
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 3.                   ; of the 1st and last bars
 xyres@tmXBValues        = (/1,2,3/)
 xyres@tmXBLabels = (/"PD","4CO2","" /)
 xyres@vpXF = xyres@vpWidthF+xyres@vpXF
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 0.3                 ; bring bars down to zero
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 3.                   ; of the 1st and last bars
 xyres@gsnCenterString       = "E3SM_SST"
 xyres@gsnLeftString  = " "
 xyres@gsnCenterStringOrthogonalPosF = -0.10
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftStringFontHeightF = 0.015
 xyres@gsnLeftStringOrthogonalPosF = 0.
 plot(15)  = gsn_csm_xy(wks,ispan(1,2,1),stdarr01(:),xyres)   ; advancing the frame
 delete(xyres@tmXBValues)
 delete(xyres@tmXBLabels)
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 4.                   ; of the 1st and last bars
 xyres@tmXBValues        = (/1,2,3,4/)
 xyres@tmXBLabels = (/"PD","2CO2","4CO2","" /)
 xyres@vpXF = xyres@vpXF+xyres@vpWidthF+yspace
 delete(xyres@tmYLValues)
 delete(xyres@tmYLLabels)
 xyres@tmYLOn        = True
 xyres@tmYRBorderOn  = False                        ; turn off right axis
 xyres@tmYLBorderOn  = True 
 xyres@vpWidthF              = 0.10
 xyres@vpHeightF             = 0.12
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 1.2                 ; bring bars down to zero                   
 xyres@tmYLOn        = True 
 xyres@tmYLMode          = "Explicit"         ; explicit labels
 xyres@tmYLValues        = (/0,0.2,0.4,0.6,0.8,1,1.2/)
 xyres@tmYLLabels        = (/"0","0.2","0.4","0.6","0.8","1.0","1.2"/)
 xyres@gsnXYBarChartColors = (/"pink"/)
 xyres@gsnCenterString       = "CESM_SSH"
 xyres@gsnCenterStringOrthogonalPosF = -0.125
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftString  = "~F22~o"
 xyres@gsnLeftStringFontHeightF = 0.015
 xyres@gsnLeftStringOrthogonalPosF = 0.
 xyres@gsnLeftStringParallelPosF = -0.05
 plot(16)  = gsn_csm_xy(wks,ispan(1,3,1),stdarr2(:),xyres)   ; advancing the frame
 delete(xyres@tmXBValues)
 delete(xyres@tmXBLabels)
 xyres@trXMinF               = 0.                    ; adds space on either end
 xyres@trXMaxF               = 3.                   ; of the 1st and last bars
 xyres@tmXBValues        = (/1,2,3/)
 xyres@tmXBLabels = (/"PD","4CO2","" /)
 xyres@tmYLOn        = False 
 xyres@tmYRBorderOn  = True                        ; turn off right axis
 xyres@tmYLBorderOn  = False
 xyres@vpXF = xyres@vpWidthF+xyres@vpXF
 xyres@trYMinF               = 0                    ; bring bars down to zero
 xyres@trYMaxF               = 1.2                ; bring bars down to zero
 xyres@gsnCenterString       = "E3SM_SSH"
 xyres@gsnLeftString  = ""
 xyres@gsnCenterStringOrthogonalPosF = -0.10
 xyres@gsnCenterStringParallelPosF = 0.475
 xyres@gsnLeftStringFontHeightF = 0.015
 xyres@gsnLeftStringOrthogonalPosF = 0.
 plot(17)  = gsn_csm_xy(wks,ispan(1,2,1),stdarr02(:),xyres)   ; advancing the frame


  polyres=True
  polyres@gsLineDashPattern=0
  polyres@gsLineColor="black"
  polyres@gsLineThicknessF=1.0
  polyres@tfPolyDrawOrder = "postDraw"

  
  eumm1    = gsn_add_polyline(wks,plot(12),(/1,1/),(/xBootLow_eke_PD,xBootHi_eke_PD/), polyres)
  eumm01   = gsn_add_polyline(wks,plot(12),(/2,2/),(/xBootLow_eke_2xCO2,xBootHi_eke_2xCO2/), polyres)
  eumm001  = gsn_add_polyline(wks,plot(12),(/3,3/),(/xBootLow_eke_4xCO2,xBootHi_eke_4xCO2/), polyres)

  fumm1    = gsn_add_polyline(wks,plot(13),(/1,1/),(/xBootLow_eke1_PD,xBootHi_eke1_PD/), polyres)
  fumm001  = gsn_add_polyline(wks,plot(13),(/2,2/),(/xBootLow_eke1_4xCO2,xBootHi_eke1_4xCO2/), polyres)


  aumm1    = gsn_add_polyline(wks,plot(14),(/1,1/),(/xBootLow_sst_PD,xBootHi_sst_PD/), polyres)
  aumm01   = gsn_add_polyline(wks,plot(14),(/2,2/),(/xBootLow_sst_2xCO2,xBootHi_sst_2xCO2/), polyres)
  aumm001  = gsn_add_polyline(wks,plot(14),(/3,3/),(/xBootLow_sst_4xCO2,xBootHi_sst_4xCO2/), polyres)

  bumm1    = gsn_add_polyline(wks,plot(15),(/1,1/),(/xBootLow_sst1_PD,xBootHi_sst1_PD/), polyres)
  bumm001  = gsn_add_polyline(wks,plot(15),(/2,2/),(/xBootLow_sst1_4xCO2,xBootHi_sst1_4xCO2/), polyres)

  cumm1    = gsn_add_polyline(wks,plot(16),(/1,1/),(/xBootLow_ssh_PD,xBootHi_ssh_PD/), polyres)
  cumm01    = gsn_add_polyline(wks,plot(16),(/2,2/),(/xBootLow_ssh_2xCO2,xBootHi_ssh_2xCO2/), polyres)
  cumm001  = gsn_add_polyline(wks,plot(16),(/3,3 /),(/xBootLow_ssh_4xCO2,xBootHi_ssh_4xCO2/), polyres)

  dumm1    = gsn_add_polyline(wks,plot(17),(/1,1/),(/xBootLow_ssh1_PD,xBootHi_ssh1_PD/), polyres)
  dumm001  = gsn_add_polyline(wks,plot(17),(/2,2/),(/xBootLow_ssh1_4xCO2,xBootHi_ssh1_4xCO2/), polyres)


  wid = 0.24/2

  gumm3   = gsn_add_polyline(wks,plot(12),(/1-wid,1+wid/),(/xBootLow_eke_PD,xBootLow_eke_PD/), polyres)
  gumm03  = gsn_add_polyline(wks,plot(12),(/1-wid,1+wid/),(/xBootHi_eke_PD,xBootHi_eke_PD/), polyres)

  iumm3   = gsn_add_polyline(wks,plot(12),(/2-wid,2+wid/),(/xBootLow_eke_2xCO2,xBootLow_eke_2xCO2/), polyres)
  iumm03  = gsn_add_polyline(wks,plot(12),(/2-wid,2+wid/),(/xBootHi_eke_2xCO2,xBootHi_eke_2xCO2/), polyres)

  humm3   = gsn_add_polyline(wks,plot(12),(/3-wid,3+wid/),(/xBootLow_eke_4xCO2,xBootLow_eke_4xCO2/), polyres)
  humm03  = gsn_add_polyline(wks,plot(12),(/3-wid,3+wid/),(/xBootHi_eke_4xCO2,xBootHi_eke_4xCO2/), polyres)

  gumm1   = gsn_add_polyline(wks,plot(14),(/1-wid,1+wid/),(/xBootLow_sst_PD,xBootLow_sst_PD/), polyres)
  gumm01  = gsn_add_polyline(wks,plot(14),(/1-wid,1+wid/),(/xBootHi_sst_PD,xBootHi_sst_PD/), polyres)

  iumm1   = gsn_add_polyline(wks,plot(14),(/2-wid,2+wid/),(/xBootLow_sst_2xCO2,xBootLow_sst_2xCO2/), polyres)
  iumm01  = gsn_add_polyline(wks,plot(14),(/2-wid,2+wid/),(/xBootHi_sst_2xCO2,xBootHi_sst_2xCO2/), polyres)

  humm1   = gsn_add_polyline(wks,plot(14),(/3-wid,3+wid/),(/xBootLow_sst_4xCO2,xBootLow_sst_4xCO2/), polyres)
  humm01  = gsn_add_polyline(wks,plot(14),(/3-wid,3+wid/),(/xBootHi_sst_4xCO2,xBootHi_sst_4xCO2/), polyres)
;----------------------------------------------------------------------------------------------------
  gumm2   = gsn_add_polyline(wks,plot(16),(/1-wid,1+wid/),(/xBootLow_ssh_PD,xBootLow_ssh_PD/), polyres)
  gumm02  = gsn_add_polyline(wks,plot(16),(/1-wid,1+wid/),(/xBootHi_ssh_PD,xBootHi_ssh_PD/), polyres)
  
  iumm2   = gsn_add_polyline(wks,plot(16),(/2-wid,2+wid/),(/xBootLow_ssh_2xCO2,xBootLow_ssh_2xCO2/), polyres)
  iumm02  = gsn_add_polyline(wks,plot(16),(/2-wid,2+wid/),(/xBootHi_ssh_2xCO2,xBootHi_ssh_2xCO2/), polyres)

  humm2   = gsn_add_polyline(wks,plot(16),(/3-wid,3+wid/),(/xBootLow_ssh_4xCO2,xBootLow_ssh_4xCO2/), polyres)
  humm02  = gsn_add_polyline(wks,plot(16),(/3-wid,3+wid/),(/xBootHi_ssh_4xCO2,xBootHi_ssh_4xCO2/), polyres)
;----------------------------------------------------------------------------------------------------
  gumm6   = gsn_add_polyline(wks,plot(13),(/1-wid,1+wid/),(/xBootLow_eke1_PD,xBootLow_eke1_PD/), polyres)
  gumm06  = gsn_add_polyline(wks,plot(13),(/1-wid,1+wid/),(/xBootHi_eke1_PD,xBootHi_eke1_PD/), polyres)
  
  humm6   = gsn_add_polyline(wks,plot(13),(/2-wid,2+wid/),(/xBootLow_eke1_4xCO2,xBootLow_eke1_4xCO2/), polyres)
  humm06  = gsn_add_polyline(wks,plot(13),(/2-wid,2+wid/),(/xBootHi_eke1_4xCO2,xBootHi_eke1_4xCO2/), polyres)
;------------------------------------------------------------------------------------------------------
  gumm4   = gsn_add_polyline(wks,plot(15),(/1-wid,1+wid/),(/xBootLow_sst1_PD,xBootLow_sst1_PD/), polyres)
  gumm04  = gsn_add_polyline(wks,plot(15),(/1-wid,1+wid/),(/xBootHi_sst1_PD,xBootHi_sst1_PD/), polyres)
  
  humm4   = gsn_add_polyline(wks,plot(15),(/2-wid,2+wid/),(/xBootLow_sst1_4xCO2,xBootLow_sst1_4xCO2/), polyres)
  humm04  = gsn_add_polyline(wks,plot(15),(/2-wid,2+wid/),(/xBootHi_sst1_4xCO2,xBootHi_sst1_4xCO2/), polyres)
;----------------------------------------------------------------------------------------------------
  gumm5  = gsn_add_polyline(wks,plot(17),(/1-wid,1+wid/),(/xBootLow_ssh1_PD,xBootLow_ssh1_PD/), polyres)
  gumm05 = gsn_add_polyline(wks,plot(17),(/1-wid,1+wid/),(/xBootHi_ssh1_PD,xBootHi_ssh1_PD/), polyres)
  
  humm5  = gsn_add_polyline(wks,plot(17),(/2-wid,2+wid/),(/xBootLow_ssh1_4xCO2,xBootLow_ssh1_4xCO2/), polyres)
  humm05 = gsn_add_polyline(wks,plot(17),(/2-wid,2+wid/),(/xBootHi_ssh1_4xCO2,xBootHi_ssh1_4xCO2/), polyres)
;----------------------------------------------------------------------------------------------------
;======================================================================================================

  draw(plot(12))
  draw(plot(13))
  draw(plot(14))
  draw(plot(15))
  draw(plot(16))
  draw(plot(17))

end 