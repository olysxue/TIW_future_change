load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

;======================= FUNCTION DEFINITIONS =====================
undef("readdata")
function readdata(var)
local dims1,var,dady
begin
dims1=dimsizes(var)
;print(dims1)
dady=new((/dims1(0),dims1(1),dims1(2)/), typeof(var))

; parameter constant 
PI  = 3.14159265
R   = 6378388.
rad = PI/180.
dens= 1025 ; reference density of seawater

dlat   = (var&lat(2)-var&lat(1))*rad
dY     = R*dlat
dady   = tofloat(center_finite_diff_n(var, dY, False, 0, 1)*3600.*24.*30.42)
copy_VarMeta(var(:,:,:), dady)
;printVarSummary(dady)
return(dady)
end
;========================================================
  begin 
  ;========================================================================
  ; calculate monthly heat flux and momentum flux 
  ;---------------------------Present day----------------------------------
   ncdir1  = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/PD/heat_flux/"
   ncdir01 = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/PD/momentum_flux/"
   fil1    = systemfunc("ls " + ncdir1 + "CESM.PD.daily.heat.flux.10-50.nc")
     f1    = addfile (fil1, "r")
     hf1   = f1->va
  copy_VarMeta(UVEL1, hf1)
   fil01    = systemfunc("ls " + ncdir01 + "CESM.PD.daily.momentum.flux.10-50.nc")
     f01    = addfile (fil01, "r")
     mf1    = f01->ua
  copy_VarMeta(UVEL1, mf1)
  printVarSummary(mf1)
  printVarSummary(hf1)
  ;============================2XCO2=====================================
  ncdir2   = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/2xCO2/heat_flux/"
   ncdir02 = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/2xCO2/momentum_flux/"
   fil2    = systemfunc("ls " + ncdir2 + "CESM.2xCO2.daily.heat.flux.10-50.nc")
     f2    = addfile (fil2, "r")
     hf2   = f2->va
  copy_VarMeta(UVEL2, hf2)
   fil02    = systemfunc("ls " + ncdir02 + "CESM.2xCO2.daily.momentum.flux.10-50.nc")
     f02    = addfile (fil02, "r")
     mf2    = f02->ua
  copy_VarMeta(UVEL2, mf2)
  ;============================4XCO2=====================================
  ncdir3   = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/4xCO2/heat_flux/"
   ncdir03 = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/4xCO2/momentum_flux/"
   fil3    = systemfunc("ls " + ncdir3 + "CESM.4xCO2.daily.heat.flux.10-50.nc")
     f3    = addfile (fil3, "r")
     hf3   = f3->va
  copy_VarMeta(UVEL3, hf3)
   fil03    = systemfunc("ls " + ncdir03 + "CESM.4xCO2.daily.momentum.flux.10-50.nc")
     f03    = addfile (fil03, "r")
     mf3    = f03->ua
  copy_VarMeta(UVEL3, mf3)
  ;========================================================================
  hf1_mon = calculate_monthly_values(hf1, "avg", 0, False)
  mf1_mon = calculate_monthly_values(mf1, "avg", 0, False)
  
  hf2_mon = calculate_monthly_values(hf2, "avg", 0, False)
  mf2_mon = calculate_monthly_values(mf2, "avg", 0, False)
  
  hf3_mon = calculate_monthly_values(hf3, "avg", 0, False)
  mf3_mon = calculate_monthly_values(mf3, "avg", 0, False)
  printVarSummary(hf1_mon)
  ;========================================================================
    fout1=addfile("/glade/scratch/aoyunxue/TIW_data/CESM_HR/CESM_HR_monthly_heat_flux.nc","c")
        fout1->hf_mon_PD    = hf1_mon
        fout1->hf_mon_2xCO2 = hf2_mon
        fout1->hf_mon_4xCO2 = hf3_mon
    fout2=addfile("/glade/scratch/aoyunxue/TIW_data/CESM_HR/CESM_HR_monthly_momentum_flux.nc","c")
        fout2->mf_mon_PD    = mf1_mon
        fout2->mf_mon_2xCO2 = mf2_mon
        fout2->mf_mon_4xCO2 = mf3_mon
  ;========================================================================
  
hfc1 = readdata(hf1_mon)
hfc2 = readdata(hf2_mon)
hfc3 = readdata(hf3_mon)

mfc1 = readdata(mf1_mon)
mfc2 = readdata(mf2_mon)
mfc3 = readdata(mf3_mon)

fout_h = addfile(basedir + "CESM_HR_monthly_heatflux_convergence.nc", "c")
fout_h->hfc_PD    = hfc1
fout_h->hfc_2xCO2 = hfc2
fout_h->hfc_4xCO2 = hfc3

fout_m = addfile(basedir + "CESM_HR_monthly_momentumflux_convergence.nc", "c")
fout_m->mfc_PD    = mfc1
fout_m->mfc_2xCO2 = mfc2
fout_m->mfc_4xCO2 = mfc3

end
