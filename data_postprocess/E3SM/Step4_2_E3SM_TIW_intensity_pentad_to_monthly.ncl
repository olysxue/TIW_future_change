
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
; read E3SMv2 
ncdir1 = "/Users/xueaoyun/Documents/work/TIW_future/data/E3SM/E3SMv1/"
;====================================================================================
fil1 = systemfunc("ls " + ncdir1+"tos_day_E3SM-1-0_0001-0500_remap_equator_10-50fft_1mon_std_1degree.nc")
f1   = addfile(fil1,"r")
sst1 = f1->sst(0:400*73-1,:,:)
printVarSummary(sst1)
;---------------------------------------------------------
fil2 = systemfunc("ls " + ncdir1 + "sst_abrupt_4xCO2_0001-0155_remap_equator_moving_20yr_detrend_10-50fft_1mon_runstd_1degree.nc")
f2   = addfile(fil2,"r")
sst2 = f2->sst(0:140*73-1,:,:)
printVarSummary(sst2)
;===========================================================================
fil01 = systemfunc("ls " + ncdir1+"zos_day_E3SM-1-0_0001-0500_remap_equator_10-50fft_1mon_std_1degree.nc")
f01   = addfile(fil01,"r")
ssh1  = f01->ssh(0:400*73-1,:,:)
printVarSummary(ssh1)
;---------------------------------------------------------
fil02 = systemfunc("ls " + ncdir1 +"ssh_abrupt_4xCO2_0001-0155_remap_equator_moving_20yr_detrend_10-50fft_1mon_runstd_1degree.nc")
f02   = addfile(fil02,"r")
ssh2  = f02->ssh(0:140*73-1,:,:)
printVarSummary(ssh2)
;===========================================================================
fil001 = systemfunc("ls " + ncdir1+"EKE_picontrol_0001-0500_remap_euqator_runave.nc")
f001   = addfile(fil001,"r")
KE1 = f001->KE(0:400*73-1,:,:)
printVarSummary(KE1)
;---------------------------------------------------------
fil002 = systemfunc("ls " + ncdir1 +"EKE_abrupt_4xCO2_0001-0150_remap_equator_runave.nc")
f002   = addfile(fil002,"r")
KE2 = f002->KE(0:140*73-1,:,:)
printVarSummary(KE2)
;============================================================================
sst1_mon = pentad_to_monthly_value(sst1,"avg",0,0)
sst2_mon = pentad_to_monthly_value(sst2,"avg",0,0)
sst1_mon!0 = "time"
sst1_mon!1 = "lat"
sst1_mon!2 = "lon"
sst1_mon&time = ispan(0,400*12-1,1)
sst1_mon&time@units = "months since 0001-01-01 12:00"
sst1_mon&lat  = sst1&LAT
sst1_mon&lon  = sst1&LON 

sst2_mon!0 = "time"
sst2_mon!1 = "lat"
sst2_mon!2 = "lon"
sst2_mon&time = ispan(0,140*12-1,1)
sst2_mon&time@units = "months since 0001-01-01 12:00"
sst2_mon&lat  = sst1&LAT
sst2_mon&lon  = sst1&LON 

ssh1_mon = pentad_to_monthly_value(ssh1,"avg",0,0)
ssh2_mon = pentad_to_monthly_value(ssh2,"avg",0,0)
copy_VarMeta(sst1_mon, ssh1_mon)
copy_VarMeta(sst2_mon, ssh2_mon)

KE1_mon  = pentad_to_monthly_value(KE1,"avg",0,0)
KE2_mon  = pentad_to_monthly_value(KE2,"avg",0,0)
copy_VarMeta(sst1_mon, KE1_mon)
copy_VarMeta(sst2_mon, KE2_mon)

; output files
diro  = "/Users/xueaoyun/Documents/work/TIW_future/data/E3SM/E3SMv1/"                     ; Output directory
filo1 = "tos_mon_E3SM-1-0_0001-0500_remap_equator_10-50fft_new_1mon_std_1degree.nc"             ; Output file
system("/bin/rm -f " + diro + filo1)    ; remove if exists
fout1  = addfile (diro + filo1, "c")  ; open output file
fout1->sst = sst1_mon

filo2 = "sst_mon_abrupt_4xCO2_0001-0155_remap_equator_moving_20yr_detrend_10-50fft_1mon_runstd_1degree.nc"             ; Output file
system("/bin/rm -f " + diro + filo2)    ; remove if exists
fout2  = addfile (diro + filo2, "c")  ; open output file
fout2->sst = sst2_mon

filo3 = "zos_mon_E3SM-1-0_0001-0500_remap_equator_10-50fft_new_1mon_std_1degree.nc"             ; Output file
system("/bin/rm -f " + diro + filo3)    ; remove if exists
fout3  = addfile (diro + filo3, "c")    ; open output file
fout3->ssh = ssh1_mon

filo4 = "ssh_mon_abrupt_4xCO2_0001-0155_remap_equator_moving_20yr_detrend_10-50fft_1mon_runstd_1degree.nc"             ; Output file
system("/bin/rm -f " + diro + filo4)    ; remove if exists
fout4  = addfile (diro + filo4, "c")  ; open output file
fout4->ssh = ssh2_mon

filo5 = "EKE_mon_picontrol_0001-0500_remap_equator_runave.nc"             ; Output file
system("/bin/rm -f " + diro + filo5)    ; remove if exists
fout5  = addfile (diro + filo5, "c")    ; open output file
fout5->KE = KE1_mon

filo6 = "EKE_mon_abrupt_4xCO2_0001-0155_remap_equator_runave.nc"             ; Output file
system("/bin/rm -f " + diro + filo6)    ; remove if exists
fout6  = addfile (diro + filo6, "c")  ; open output file
fout6->KE = KE2_mon
end 