

c      $Id$

      parameter (NGRIDMAX=1000)

      real*8 alpha0,beta0,
     + massA_grid(NGRIDMAX), alphaA_grid(NGRIDMAX),
     + betaA_grid(NGRIDMAX), kA_grid(NGRIDMAX), 
     + mA,mB,alphaA,betaA,kA,alphaB,betaB,kB,
     + dalphaA_dm,dbetaA_dm,dkA_dm,dalphaB_dm,dbetaB_dm,dkB_dm,
     + dk_dm,dgamma_dm,ddr_dm,ddth_dm,dpbdot_dm,dsi_dm,
     + dk_dm2,dgamma_dm2,ddr_dm2,ddth_dm2,dpbdot_dm2,dsi_dm2

      real*8 pbdot_m, pbdot_d, pbdot_qphi, pbdot_qg, pbdot_old

      real*8 dgamma_dm_factor, dk_dm_factor, dsi_dm_factor, ddr_dm_factor, ddth_dm_factor, dpbdot_m_dm_factor, dpbdot_d1_dm_factor, dpbdot_d2_dm_factor, dpbdot_d3_dm_factor, dpbdot_qphi_dm_factor, dpbdot_qg_dm_factor, dpbdot_dm_factor
      real*8 dgamma_dm2_factor, dk_dm2_factor, dsi_dm2_factor, ddr_dm2_factor, ddth_dm2_factor, dpbdot_m_dm2_factor, dpbdot_d1_dm2_factor, dpbdot_d2_dm2_factor, dpbdot_d3_dm2_factor, dpbdot_qphi_dm2_factor, dpbdot_qg_dm2_factor, dpbdot_dm2_factor

      real*8 dgamma_dm_ap, dk_dm_ap, dsi_dm_ap, ddr_dm_ap, ddth_dm_ap, dpbdot_dm_ap, dgamma_dm2_ap, dk_dm2_ap, dsi_dm2_ap, ddr_dm2_ap, ddth_dm2_ap, dpbdot_dm2_ap

      integer N_grid

      character eosname*16, companion_type*16

      logical ddstg_read, gr_case

      common/ddstg_keys/ gr_case

      common/ddstg_data/ N_grid,ddstg_read,
     + alpha0,beta0,eosname,companion_type,
     + massA_grid, alphaA_grid, betaA_grid, kA_grid, 
     + mA,mB,alphaA,betaA,kA,alphaB,betaB,kB,
     + dalphaA_dm,dbetaA_dm,dkA_dm,dalphaB_dm,dbetaB_dm,dkB_dm,
     + dk_dm,dgamma_dm,ddr_dm,ddth_dm,dpbdot_dm,dsi_dm,
     + dk_dm2,dgamma_dm2,ddr_dm2,ddth_dm2,dpbdot_dm2,dsi_dm2

      common/ddstg_factor_derivatives/ dgamma_dm_factor, dk_dm_factor, dsi_dm_factor, ddr_dm_factor, ddth_dm_factor, dpbdot_m_dm_factor, dpbdot_d1_dm_factor, dpbdot_d2_dm_factor, dpbdot_d3_dm_factor, dpbdot_qphi_dm_factor, dpbdot_qg_dm_factor, dpbdot_dm_factor, dgamma_dm2_factor, dk_dm2_factor, dsi_dm2_factor, ddr_dm2_factor, ddth_dm2_factor, dpbdot_m_dm2_factor, dpbdot_d1_dm2_factor, dpbdot_d2_dm2_factor, dpbdot_d3_dm2_factor, dpbdot_qphi_dm2_factor, dpbdot_qg_dm2_factor, dpbdot_dm2_factor

      common/ddstg_pbdot/ pbdot_m, pbdot_d, pbdot_qphi, pbdot_qg, pbdot_old

      common/ddstg_approx_derivatives/ dgamma_dm_ap, dk_dm_ap, dsi_dm_ap, ddr_dm_ap, ddth_dm_ap, dpbdot_dm_ap, dgamma_dm2_ap, dk_dm2_ap, dsi_dm2_ap, ddr_dm2_ap, ddth_dm2_ap, dpbdot_dm2_ap