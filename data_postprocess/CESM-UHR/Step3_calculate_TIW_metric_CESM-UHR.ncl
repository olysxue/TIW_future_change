;-------------------------------------------------------------------
; TIW Intensity Diagnostics Script for CESM——UHR HR
; Computes monthly STD, VAR, and EKE from 10–50d band-pass filtered data
; Supports experiments: PD, 2xCO2, 4xCO2
; Applies to UVEL and VVEL averaged over 0–50m, as well as SST and SSH
;-------------------------------------------------------------------
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"

begin
  dens     = 1025.0      ; seawater density
  basedir  = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/"
  exps     = (/"PD", "2xCO2", "4xCO2"/)

  do i = 0, dimsizes(exps)-1
    exp = exps(i)
    print("==== Processing Experiment: " + exp + " ====")

    ;=== Load UVEL and VVEL (0–50m averaged)
    file_u = basedir + exp + "/UVEL/CESM." + exp + ".day.UVEL_0_50m.anom.10-50.fft.nc"
    file_v = basedir + exp + "/VVEL/CESM." + exp + ".day.VVEL_0_50m.anom.10-50.fft.nc"
    fu = addfile(file_u, "r")
    fv = addfile(file_v, "r")
    ua = fu->ua_0_50m
    va = fv->va_0_50m

    ;=== Load SST and SSH
    file_sst = basedir + exp + "/SST/CESM." + exp + ".day.SST.anom.10-50.fft.nc"
    file_ssh = basedir + exp + "/SSH/CESM." + exp + ".day.SSH.anom.10-50.fft.nc"
    fsst = addfile(file_sst, "r")
    fssh = addfile(file_ssh, "r")
    sst = fsst->sst
    ssh = fssh->ssh

    ;=== Compute monthly std/var
    ua_std  = calculate_monthly_values(ua,  "std", 0, False)
    ua_var  = calculate_monthly_values(ua,  "var", 0, False)
    va_std  = calculate_monthly_values(va,  "std", 0, False)
    va_var  = calculate_monthly_values(va,  "var", 0, False)
    sst_std = calculate_monthly_values(sst, "std", 0, False)
    sst_var = calculate_monthly_values(sst, "var", 0, False)
    ssh_std = calculate_monthly_values(ssh, "std", 0, False)
    ssh_var = calculate_monthly_values(ssh, "var", 0, False)

    ;=== Compute EKE = 0.5 * ρ * (u'^2 + v'^2)
    EKE = 0.5 * dens * (ua^2 + va^2)
    EKE@long_name = "Eddy Kinetic Energy (0–50m, band-pass 10–50d)"
    EKE@units     = "J/m3"
    EKE_mon = calculate_monthly_values(EKE, "avg", 0, False)


    ;=== Output: std/var of all fields
    outdir = basedir + exp + "/"

    ;--- Save UVEL/VVEL STD & VAR
    fout1 = addfile(outdir + "CESM_HR_TIW_rms_var_uv_0_50m_10-50_fft_" + exp + ".nc", "c")
    fout1->uvel_std = ua_std
    fout1->uvel_var = ua_var
    fout1->vvel_std = va_std
    fout1->vvel_var = va_var

    ;--- Save SST/SSH STD & VAR
    fout2 = addfile(outdir + "CESM_HR_TIW_rms_var_sst_ssh_10-50_fft_" + exp + ".nc", "c")
    fout2->sst_std  = sst_std
    fout2->sst_var  = sst_var
    fout2->ssh_std  = ssh_std
    fout2->ssh_var  = ssh_var

    ;--- Save EKE
    fout3 = addfile(outdir + "CESM_HR_TIW_EKE_10-50_fft_0_50m_" + exp + ".nc", "c")
    fout3->EKE = EKE_mon

    delete([/ua, va, sst, ssh, ua_std, ua_var, va_std, va_var, sst_std, sst_var, ssh_std, ssh_var, EKE, EKE_mon/])
  end do
end
