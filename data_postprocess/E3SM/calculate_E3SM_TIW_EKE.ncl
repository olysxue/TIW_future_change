load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
undef ("pentad_to_monthly_value")
function pentad_to_monthly_value (x:numeric, arith:string, nDim[1]:integer, opt[1]:logical)

; olys in hawaii 2019.01.10
; calculate monthly values from pentad data [avg, sum, min, max]
; x: numeric array of 5D or less [eg: time,lev,lat,lon]
; time
; arith: "avg" [also, "ave"], "sum","min","max"     others may be added later
; nDim : scalar integer that specifies the 'time' dimension [generally, 0]
; opt  : only used to eliminate a warning message
;        opt= 0 ... default ... print Warning
;        opt=-1 ... do not print Warning
;
; Sample usage:  x(time,lat,lon) where time are n-hrly or daily values.
;        xMonthVar = pentad_to_monthly_value(x, "var", 0, False)
;        xMonthStd = pentad_to_monthly_value(x, "std", 0, False)
;                  ************
; Note that NCL's 'dimension reduction' which, in general, I view as a feature,
;      introduces some processing issues in the 'corner case' of one value.
;      This requires extra (nuisance) steps and leads to 'code bloat'.  :-(
;                  ************
local dnamx, TNAME, dimx, rankx, ntim, nmon,  xStat
begin
    if (.not.(arith.eq."ave" .or. arith.eq."avg" .or. arith.eq."sum"  \
                             .or. arith.eq."min" .or. arith.eq."max"  \
                             .or. arith.eq."var" .or. arith.eq."std") ) then
        print("calculate_monthly_values: unrecognizezed 'arith' argument="+arith)
        exit
    end if

    dnamx    = getvardims(x)                       ; dimension names
    TNAME    = dnamx(nDim)                         ; typically TNAME="time"
    dimx     = dimsizes(x)                         ; dimensions sizes
    rankx    = dimsizes(dimx)                      ; # of dimensions
    ntim     = dimx(nDim)                          ; size of input time dimension
    nyear    = ntim/73                             ; month of input time
    nmon     = nyear*12

    if (rankx.eq.3) then
    xStat  = new ( (/nmon,dimx(1),dimx(2)/), typeof(x), getVarFillValue(x))

          if (arith.eq."avg" .or. arith.eq."ave") then
                          do i=0,nyear-1
                             xStat(12*i+0,:,:)=dim_avg_n_Wrap(x(73*i+0:73*i+5,:,:),0)
                             xStat(12*i+1,:,:)=dim_avg_n_Wrap(x(73*i+6:73*i+11,:,:),0)
                             xStat(12*i+2,:,:)=dim_avg_n_Wrap(x(73*i+12:73*i+17,:,:),0)
                             xStat(12*i+3,:,:)=dim_avg_n_Wrap(x(73*i+18:73*i+23,:,:),0)
                             xStat(12*i+4,:,:)=dim_avg_n_Wrap(x(73*i+24:73*i+29,:,:),0)
                             xStat(12*i+5,:,:)=dim_avg_n_Wrap(x(73*i+30:73*i+35,:,:),0)
                             xStat(12*i+6,:,:)=dim_avg_n_Wrap(x(73*i+36:73*i+41,:,:),0)
                             xStat(12*i+7,:,:)=dim_avg_n_Wrap(x(73*i+42:73*i+48,:,:),0)
                             xStat(12*i+8,:,:)=dim_avg_n_Wrap(x(73*i+49:73*i+54,:,:),0)
                             xStat(12*i+9,:,:)=dim_avg_n_Wrap(x(73*i+55:73*i+60,:,:),0)
                             xStat(12*i+10,:,:)=dim_avg_n_Wrap(x(73*i+61:73*i+66,:,:),0)
                             xStat(12*i+11,:,:)=dim_avg_n_Wrap(x(73*i+67:73*i+72,:,:),0)
                          end do
          end if
   
          if (arith.eq."sum") then
            do i=0,nyear-1
                    xStat(12*i+0,:,:)=dim_sum_n_Wrap(x(73*i+0:73*i+5,:,:),0)
                    xStat(12*i+1,:,:)=dim_sum_n_Wrap(x(73*i+6:73*i+11,:,:),0)
                    xStat(12*i+2,:,:)=dim_sum_n_Wrap(x(73*i+12:73*i+17,:,:),0)
                    xStat(12*i+3,:,:)=dim_sum_n_Wrap(x(73*i+18:73*i+23,:,:),0)
                    xStat(12*i+4,:,:)=dim_sum_n_Wrap(x(73*i+24:73*i+29,:,:),0)
                    xStat(12*i+5,:,:)=dim_sum_n_Wrap(x(73*i+30:73*i+35,:,:),0)
                    xStat(12*i+6,:,:)=dim_sum_n_Wrap(x(73*i+36:73*i+41,:,:),0)
                    xStat(12*i+7,:,:)=dim_sum_n_Wrap(x(73*i+42:73*i+48,:,:),0)
                    xStat(12*i+8,:,:)=dim_sum_n_Wrap(x(73*i+49:73*i+54,:,:),0)
                    xStat(12*i+9,:,:)=dim_sum_n_Wrap(x(73*i+55:73*i+60,:,:),0)
                    xStat(12*i+10,:,:)=dim_sum_n_Wrap(x(73*i+61:73*i+66,:,:),0)
                    xStat(12*i+11,:,:)=dim_sum_n_Wrap(x(73*i+67:73*i+72,:,:),0)
            end do
          if

                if (arith.eq."var") then
                            do i=0,nyear-1
                                    xStat(12*i+0,:,:)=dim_variance_n_Wrap(x(73*i+0:73*i+5,:,:),0)
                                    xStat(12*i+1,:,:)=dim_variance_n_Wrap(x(73*i+6:73*i+11,:,:),0)
                                    xStat(12*i+2,:,:)=dim_variance_n_Wrap(x(73*i+12:73*i+17,:,:),0)
                                    xStat(12*i+3,:,:)=dim_variance_n_Wrap(x(73*i+18:73*i+23,:,:),0)
                                    xStat(12*i+4,:,:)=dim_variance_n_Wrap(x(73*i+24:73*i+29,:,:),0)
                                    xStat(12*i+5,:,:)=dim_variance_n_Wrap(x(73*i+30:73*i+35,:,:),0)
                                    xStat(12*i+6,:,:)=dim_variance_n_Wrap(x(73*i+36:73*i+41,:,:),0)
                                    xStat(12*i+7,:,:)=dim_variance_n_Wrap(x(73*i+42:73*i+48,:,:),0)
                                    xStat(12*i+8,:,:)=dim_variance_n_Wrap(x(73*i+49:73*i+54,:,:),0)
                                    xStat(12*i+9,:,:)=dim_variance_n_Wrap(x(73*i+55:73*i+60,:,:),0)
                                    xStat(12*i+10,:,:)=dim_variance_n_Wrap(x(73*i+61:73*i+66,:,:),0)
                                    xStat(12*i+11,:,:)=dim_variance_n_Wrap(x(73*i+67:73*i+72,:,:),0)
                            end do
                 end if

                 if (arith.eq."std") then
                            do i=0,nyear-1
                                    xStat(12*i+0,:,:)=dim_stddev_n_Wrap(x(73*i+0:73*i+5,:,:),0)
                                    xStat(12*i+1,:,:)=dim_stddev_n_Wrap(x(73*i+6:73*i+11,:,:),0)
                                    xStat(12*i+2,:,:)=dim_stddev_n_Wrap(x(73*i+12:73*i+17,:,:),0)
                                    xStat(12*i+3,:,:)=dim_stddev_n_Wrap(x(73*i+18:73*i+23,:,:),0)
                                    xStat(12*i+4,:,:)=dim_stddev_n_Wrap(x(73*i+24:73*i+29,:,:),0)
                                    xStat(12*i+5,:,:)=dim_stddev_n_Wrap(x(73*i+30:73*i+35,:,:),0)
                                    xStat(12*i+6,:,:)=dim_stddev_n_Wrap(x(73*i+36:73*i+41,:,:),0)
                                    xStat(12*i+7,:,:)=dim_stddev_n_Wrap(x(73*i+42:73*i+48,:,:),0)
                                    xStat(12*i+8,:,:)=dim_stddev_n_Wrap(x(73*i+49:73*i+54,:,:),0)
                                    xStat(12*i+9,:,:)=dim_stddev_n_Wrap(x(73*i+55:73*i+60,:,:),0)
                                    xStat(12*i+10,:,:)=dim_stddev_n_Wrap(x(73*i+61:73*i+66,:,:),0)
                                    xStat(12*i+11,:,:)=dim_stddev_n_Wrap(x(73*i+67:73*i+72,:,:),0)
                            end do
                 end if
end if 
  ;==========================================================================================
xStat@NCL_tag = "calculate_monthly_values: arith="+arith
return(xStat)

begin

  scenario  = "picontrol"        ; "picontrol" or "4xCO2"
  model     = "E3SM-1-0"
  grid      = "gr"
  start     = "00010101"
  end       = "01001231"
  res       = "1deg"
  region    = "equator"
  method    = "new"
  freq      = "mon"

  fU = addfile("UVEL_pentad_"+scenario+".nc", "r")
  fV = addfile("VVEL_pentad_"+scenario+".nc", "r")
  u  = fU->UVEL   ; [time, lat, lon]
  v  = fV->VVEL

  ; ===== Step 1: Compute monthly total KE =====
  ke = 0.5 * (u^2 + v^2)
  ke_mon = month_avg(ke)
  ke_mon@long_name = "monthly total kinetic energy from pentad data"
  copy_VarMeta(u(0,:,:), ke_mon)

  ; ===== Step 2: Compute monthly mean u/v and LSKE =====
  u_mon = month_avg(u)
  v_mon = month_avg(v)
  lske = 0.5 * (u_mon^2 + v_mon^2)
  lske@long_name = "monthly large-scale kinetic energy"
  copy_VarMeta(u_mon(0,:,:), lske)

  ; ===== Step 3: Subtract to get TIW-EKE =====
  eke = ke_mon - lske
  eke@long_name = "TIW-induced eddy kinetic energy (EKE = TKE - LSKE)"
  copy_VarMeta(lske, eke)

  ; ===== Step 4: Compute mean over time =====
  eke_mean = dim_avg_n_Wrap(eke, 0)
  copy_VarMeta(eke(0,:,:), eke_mean)
  eke_mean@long_name = "time-mean TIW-induced EKE"

  ; ===== Step 5: output =====
  ofile = "TIW_KE_"+freq+"_"+model+"_"+scenario+"_"+grid+"_"+start+"-"+end+ \
          "_regrid_"+res+"_"+region+"_"+method+".nc"
  print("  ➤ Output file: " + ofile)

  ; ===== Step 6: write NetCDF =====
  fout = addfile(ofile, "c")
  fout->TKE        = ke_mon
  fout->LSKE       = lske
  fout->TIW_EKE    = eke
  fout->TIW_EKE_mean = eke_mean


end

