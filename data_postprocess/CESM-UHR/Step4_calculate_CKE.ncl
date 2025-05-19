load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

undef ("barotropic_kinetic_energy")
function barotropic_kinetic_energy(ua:numeric, va:numeric,uave:numeric,vave:numeric,lat[*]:numeric, lon[*]:numeric,density,nave)
local cke, nlat, nlon, rad, re, rr, dlat, dlon, dx, tx, dy, ty, nd, nl, dims, uavedx, uavedy, vavedx, vavedy,\
      uu, uv, vv,uuave, uvave, vvave

begin
 if (.not.all(dimsizes(ua).eq.dimsizes(va))) then
     print("barotropic_kinetic_energy: dimension size mismatch")
     print("         ua: "+dimsizes(ua))
     print("         va: "+dimsizes(va))
     exit
 end if
  tx=new(dimsizes(ua),typeof(ua),ua@_FillValue)
  ty=new(dimsizes(ua),typeof(ua),ua@_FillValue)
        nlat = dimsizes(lat)
        nlon = dimsizes(lon)
   dims=dimsizes(ua)
     nd=dimsizes(dims)
  rad    = 4.0*atan(1.0)/180.0
  re     = 6378388.0
  rr     = re*rad
   dlat =  (lat(2)-lat(1))*rad
   dlon =  (lon(2)-lon(1))*rad
   dy   = abs(lat(2)-lat(1))*rr
   
   if (nd .eq. 4)
   uavedy  = center_finite_diff_n(uave,dy,False,0,2)
   vavedy  = center_finite_diff_n(vave,dy,False,0,2)
     do nl=0,nlat-1                      ; loop over each latitude
       dx = re*cos(rad*lat(nl))*dlon
       uavedx = center_finite_diff_n(uave, dx, False, 0, 3)
       vavedx = center_finite_diff_n(vave, dx, False, 0, 3)
     end do
 end if

  if (nd .eq. 3)
   uavedy  = center_finite_diff_n(uave,dy,False,0,1)
   vavedy  = center_finite_diff_n(vave,dy,False,0,1)
     do nl=0,nlat-1                      ; loop over each latitude
       dx = re*cos(rad*lat(nl))*dlon
       uavedx = center_finite_diff_n(uave, dx, False, 0, 2)
       vavedx = center_finite_diff_n(vave, dx, False, 0, 2)
     end do
 end if
;*********************************************************************************
 uu = ua*ua
 uv = ua*va
 vv = va*va
 uuave = runave_n_Wrap(uu,nave,0,0)
  uvave = runave_n_Wrap(uv,nave,0,0)
  vvave = runave_n_Wrap(vv,nave,0,0)

  cke = -1.*density * (uuave * uavedx + uvave * (uavedy + vavedx) + vvave * vavedy)
  cke@long_name = "barotropic kinetic energy"
  copy_VarCoords(ua, cke)
  return(cke)
end

begin

;*******************************************
;  parameter set
;varname = "EKE"
 dens= 1025   ; kg/m3 reference density of seawater
 nave = 31
;*******************************************
; PD 
; read TIW u & v 

 ; datadir1 = "/Users/xueaoyun/Documents/work/TIW_future/data/filter/"
 ; datadir2 = "/Users/xueaoyun/Documents/work/TIW_future/data/filter/"
 ; datadir3 = "/Users/xueaoyun/Documents/work/TIW_future/data/filter/"

 datadir1 = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/PD/UVEL/"
 datadir2 = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/2xCO2/UVEL/"
 datadir3 = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/4xCO2/UVEL/"

   fil1  = systemfunc("ls "+datadir1 + "CESM.PD.day.UVEL_surface.anom.10-50.fft.nc")
   fil2  = systemfunc("ls "+datadir2 + "CESM.2xCO2.day.UVEL_surface.anom.10-50.fft.nc")
   fil3  = systemfunc("ls "+datadir3 + "CESM.4xCO2.day.UVEL_surface.anom.10-50.fft.nc")
 
   f1    = addfile(fil1, "r")
   f2    = addfile(fil2, "r")
   f3    = addfile(fil3, "r")

   ua_PD     = f1->ua(:,{-10:10},{180:280})  
   ua_2xCO2  = f2->ua(:,{-10:10},{180:280})
   ua_4xCO2  = f3->ua(:,{-10:10},{180:280})
   lat = ua_PD&lat 
   lon = ua_PD&lon 
   ;=======================================
    datadir01 = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/PD/VVEL/"
    datadir02 = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/2xCO2/VVEL/"
    datadir03 = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/4xCO2/VVEL/"
   fil01  = systemfunc("ls "+datadir01 + "CESM.PD.day.VVEL_surface.anom.10-50.fft.nc")
   fil02  = systemfunc("ls "+datadir02 + "CESM.2xCO2.day.VVEL_surface.anom.10-50.fft.nc")
   fil03  = systemfunc("ls "+datadir03 + "CESM.4xCO2.day.VVEL_surface.anom.10-50.fft.nc")
 
   f01    = addfile(fil01, "r")
   f02    = addfile(fil02, "r")
   f03    = addfile(fil03, "r")

   va_PD     = f01->va(:,{-10:10},{180:280})  
   va_2xCO2  = f02->va(:,{-10:10},{180:280})
   va_4xCO2  = f03->va(:,{-10:10},{180:280})

;=====================================read ave======================================
 ; datadir01 = "/Users/xueaoyun/Documents/work/TIW_future/data/anom/"
 ; datadir02 = "/Users/xueaoyun/Documents/work/TIW_future/data/anom/"
 ; datadir03 = "/Users/xueaoyun/Documents/work/TIW_future/data/anom/"

 datadir001 = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/PD/UVEL/"
 datadir002 = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/2xCO2/UVEL/"
 datadir003 = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/4xCO2/UVEL/"

   fil001  = systemfunc("ls "+datadir001 + "UVEL.day.0111-0140.10S_10N.180_80W.levidx01_025degree.nc")
   fil002  = systemfunc("ls "+datadir002 + "UVEL.day.0141-0170.10S_10N.180_80W.levidx01_025degree.nc")
   fil003  = systemfunc("ls "+datadir003 + "UVEL.day.0141-0170.10S_10N.180_80W.levidx01_025degree.nc")
 
   f001    = addfile(fil001, "r")
   f002    = addfile(fil002, "r")
   f003    = addfile(fil003, "r")

   u1  = f001->UVEL(:,0,{-10:10},{180:280})  
   u2  = f002->UVEL(:,0,{-10:10},{180:280})
   u3  = f003->UVEL(:,0,{-10:10},{180:280})

   uave_PD = runave_n_Wrap(u1, 61,0,0)
   uave_2xCO2 = runave_n_Wrap(u2,61,0,0)
   uave_4xCO2 = runave_n_Wrap(u3,61,0,0)

;=====================================================================================
   datadir0001 = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/PD/VVEL/"
   datadir0002 = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/2xCO2/VVEL/"
   datadir0003 = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/4xCO2/VVEL/"

   fil0001  = systemfunc("ls "+datadir0001 + "VVEL.day.0111-0140.10S_10N.180_80W.levidx01_025degree.nc")
   fil0002  = systemfunc("ls "+datadir0002 + "VVEL.day.0141-0170.10S_10N.180_80W.levidx01_025degree.nc")
   fil0003  = systemfunc("ls "+datadir0003 + "VVEL.day.0141-0170.10S_10N.180_80W.levidx01_025degree.nc")
 
   f0001    = addfile(fil0001, "r")
   f0002    = addfile(fil0002, "r")
   f0003    = addfile(fil0003, "r")

   v1  = f0001->VVEL(:,0,{-10:10},{180:280})  
   v2  = f0002->VVEL(:,0,{-10:10},{180:280})
   v3  = f0003->VVEL(:,0,{-10:10},{180:280})

   vave_PD    = runave_n_Wrap(v1,61,0,0)
   vave_2xCO2 = runave_n_Wrap(v2,61,0,0)
   vave_4xCO2 = runave_n_Wrap(v3,61,0,0)
;=====================================================================================

cke_PD  = barotropic_kinetic_energy(ua_PD, va_PD, uave_PD, vave_PD ,lat,lon,dens,nave)
cke_2xCO2  = barotropic_kinetic_energy(ua_2xCO2, va_2xCO2, uave_2xCO2, vave_2xCO2 ,lat,lon,dens,nave)
cke_4xCO2  = barotropic_kinetic_energy(ua_4xCO2, va_4xCO2, uave_4xCO2, vave_4xCO2 ,lat,lon,dens,nave)
;=====================================================================================
; output CKE 

 ;  opdir1    = "/Users/xueaoyun/Documents/work/TIW_future/data/energy/"
    opdir1    = "/glade/scratch/aoyunxue/TIW_data/CESM_HR/energy/"

    fname1    = "TIW_CKE_PD.nc"
    system("/bin/rm -f "+ opdir1 + fname1)
    f1      = addfile(opdir1 + fname1,"c")    ; write netCDF file
    f1->CKE = cke_PD 
  
    fname2    = "TIW_CKE_2xCO2.nc"
    system("/bin/rm -f "+ opdir1 + fname2)
    f2      = addfile(opdir1 + fname2,"c")    ; write netCDF file
    f2->CKE = cke_2xCO2 

    fname3    = "TIW_CKE_4xCO2.nc"
    system("/bin/rm -f "+ opdir1 + fname3)
    f3      = addfile(opdir1 + fname3,"c")    ; write netCDF file
    f3->CKE = cke_4xCO2 

end 



