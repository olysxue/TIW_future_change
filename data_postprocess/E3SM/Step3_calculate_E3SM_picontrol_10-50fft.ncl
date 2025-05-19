load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin

  ;===================== SST =====================
  dir1 = "/glade/derecho/scratch/aoyunxue/E3SMv2/picontrol/daily/"
  fil0 = systemfunc("ls " + dir1 + "sst_picontrol_0001-0500_remap_equator_new.nc")
  fil1 = systemfunc("ls " + dir1 + "sst_picontrol_0001-0500_remap_equator_anom.nc")
  f0 = addfile(fil0, "r")
  f1 = addfile(fil1, "r")
  sst = f0->sst
  ta  = f1->sst
  copy_VarMeta(sst, ta)
  ta_time = ta(lat|:, lon|:, time|:)

  fca = 1./2.
  fcb = 1./10.
  cf  = ezfftf(ta_time)
  fcf = fspan(0, 0.5, dimsizes(ta_time&time)/2)
  ifcf = ind(fcf.lt.fcb .or. fcf.gt.fca)
  cf(:,:,:,ifcf) = 0.0
  tanom1 = ezfftb(cf, cf@xbar)
  copy_VarMeta(ta_time, tanom1)
  ta1 = tanom1(time|:, lat|:, lon|:)
  delete(cf)
  delete(ta_time)
  delete(tanom1)

  ;===================== SSH =====================
  fil2 = systemfunc("ls " + dir1 + "ssh_picontrol_0001-0500_remap_equator.nc")
  fil3 = systemfunc("ls " + dir1 + "ssh_picontrol_0001-0500_remap_equator_anom.nc")
  f2 = addfile(fil2, "r")
  f3 = addfile(fil3, "r")
  ssh = f2->ssh
  za  = f3->ssh
  copy_VarMeta(ssh, za)
  za_time = za(lat|:, lon|:, time|:)

  cf = ezfftf(za_time)
  fcf = fspan(0, 0.5, dimsizes(za_time&time)/2)
  ifcf = ind(fcf.lt.fcb .or. fcf.gt.fca)
  cf(:,:,:,ifcf) = 0.0
  zanom1 = ezfftb(cf, cf@xbar)
  copy_VarMeta(za_time, zanom1)
  za1 = zanom1(time|:, lat|:, lon|:)
  delete(cf)
  delete(za_time)
  delete(zanom1)

  ;===================== Output =====================
  fout_sst = addfile(dir1 + "sst_picontrol_0001-0500_remap_equator_anom_10-50fft.nc", "c")
  fout_sst->sst = ta1

  fout_ssh = addfile(dir1 + "ssh_picontrol_0001-0500_remap_equator_anom_10-50fft.nc", "c")
  fout_ssh->ssh = za1

end
