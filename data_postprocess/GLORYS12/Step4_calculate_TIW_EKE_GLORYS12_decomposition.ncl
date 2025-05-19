load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"

begin


  fU = addfile("cmems_uo_180.00W-60.00W_20.00S-20.00N_surface_1993-2023_025degree_anom.nc", "r")
  fV = addfile("cmems_vo_180.00W-60.00W_20.00S-20.00N_surface_1993-2023_025degree_anom.nc", "r")
  u = fU->uo    ; [time, lat, lon]
  v = fV->vo

  ; ==== Step 1: TKE (total KE) ====
  print("  ➤ Calculating total kinetic energy (TKE)...")
  ke = 0.5 * (u^2 + v^2)
  ke_mon = calculate_monthly_values(ke, "avg", 0, False)
  ke_mon@long_name = "monthly total kinetic energy"
  copy_VarMeta(u(0,:,:), ke_mon)

  ; ==== Step 2: LSKE (large-scale KE) ====
  print("  ➤ Calculating large-scale kinetic energy (LSKE)...")
  u_mon = calculate_monthly_values(u, "avg", 0, False)
  v_mon = calculate_monthly_values(v, "avg", 0, False)
  lske = 0.5 * (u_mon^2 + v_mon^2)
  lske@long_name = "monthly large-scale kinetic energy"
  copy_VarMeta(u_mon(0,:,:), lske)

  ; ==== Step 3: TIW-EKE = TKE − LSKE ====
  print("  ➤ Subtracting to get TIW-EKE...")
  eke = ke_mon - lske
  eke@long_name = "TIW-induced eddy kinetic energy"
  copy_VarMeta(lske, eke)

  ; ==== Step 4: Mean over time ====
  print("  ➤ Calculating time-mean TIW-EKE...")
  eke_mean = dim_avg_n_Wrap(eke, 0)
  copy_VarMeta(eke(0,:,:), eke_mean)
  eke_mean@long_name = "mean TIW-induced eddy kinetic energy"

  ; ==== Step 5: Output files ====
  fname_full = "TIW_KE_mon_GLORYS12V1_obs_gr_19930101-20231231_regrid__equator_new.nc"
  fname_mean = "TIW_KE_mon_GLORYS12V1_obs_gr_19930101-20231231_regrid__equator_new_mean.nc"

  print("  ➤ Writing: " + fname_full)
  fout = addfile(fname_full, "c")
  fout->TKE = ke_mon
  fout->LSKE = lske
  fout->TIW_EKE = eke

  print("  ➤ Writing: " + fname_mean)
  fmean = addfile(fname_mean, "c")
  fmean->TIW_EKE_mean = eke_mean

  print("✅ Done. Files saved:")
  print("   - " + fname_full)
  print("   - " + fname_mean)
end
