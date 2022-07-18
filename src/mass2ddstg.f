c      $Id$
c
      subroutine mass2ddstg(am,am2,
     +       x,ecc,an,arr,ar,xk,si,gamma,pbdot,dr,dth)

c       given system masses m, m2 and keplerian parameters x, ecc, an
c       calculate values of arr, ar, si, gamma, pbdot under GR
      include 'ddstg.h'

      real*8 TWOPI,BARE_SUNMASS,SUNMASS,ARRTOL
      parameter (TWOPI=6.28318530717958648d0,SUNMASS=4.925490947d-6)
      parameter (ARRTOL=1.d-10)

      real*8 am,am2,
     + x,ecc,an,arr,ar,xk,si,gamma,pbdot,dr,dth
      real*8 m, m2, m1, arr0, arrold, XA, XB, GAB, nu, 
     + bare_m2, bare_m1, bare_m,arr_it,
     + gamma_AB, eps, pbdot_d1, pbdot_d2, pbdot_d3,
     + dpbdot_m_dm, dpbdot_d_dm, dpbdot_d1_dm, dpbdot_d2_dm,
     + dpbdot_d3_dm, dpbdot_qphi_dm, dpbdot_qg_dm,
     + dpbdot_m_dm2, dpbdot_d_dm2, dpbdot_d1_dm2, dpbdot_d2_dm2,
     + dpbdot_d3_dm2, dpbdot_qphi_dm2, dpbdot_qg_dm2
      real*8 n, e, k, W_mult

      real*8 dgamma_dalphaA, dgamma_dbetaA, dgamma_dkA, dk_dalphaA, dk_dbetaA, dk_dkA, dsi_dalphaA, dsi_dbetaA, dsi_dkA, ddr_dalphaA, ddr_dbetaA, ddr_dkA, ddth_dalphaA, ddth_dbetaA, ddth_dkA, dpbdot_m_dalphaA, dpbdot_m_dbetaA, dpbdot_m_dkA, dpbdot_d1_dalphaA, dpbdot_d1_dbetaA, dpbdot_d1_dkA, dpbdot_d2_dalphaA, dpbdot_d2_dbetaA, dpbdot_d2_dkA, dpbdot_d3_dalphaA, dpbdot_d3_dbetaA, dpbdot_d3_dkA, dpbdot_qphi_dalphaA, dpbdot_qphi_dbetaA, dpbdot_qphi_dkA, dpbdot_qg_dalphaA, dpbdot_qg_dbetaA, dpbdot_qg_dkA 

      logical exact_der, exact_factor_der


      exact_der = .true.
      exact_factor_der = .true.

      if (ddstg_read .eqv. .false.) then
        write(*,*) "Reading data for DDSTG"
        
        if (alpha0 /= 0.d0 .or. beta0 /= 0.d0) then
          call get_ns_grids()
        endif

        ddstg_read = .true.

        mA = am-am2
        mB = am2
        call interpolate_ddstg_factors()

        write(*,'('' alpha0, beta0:         '', 2d14.5)') alpha0, beta0
        write(*,'('' mA, alphaA, betaA, kA: '', 4d14.5)') mA, alphaA, betaA, kA
        write(*,'('' dalphaA_dm, dbetaA_dm, dkA_dm: '', 4d14.5)') dalphaA_dm, dbetaA_dm, dkA_dm
        write(*,*) "Companion type:           ", companion_type
        write(*,'('' mB, alphaB, betaB, kB: '', 4d14.5)') mB, alphaB, betaB, kB
        write(*,'('' dalphaB_dm, dbetaB_dm, dkB_dm: '', 4d14.5)') dalphaB_dm, dbetaB_dm, dkB_dm

      endif

      if ((mA /= am-am2) .or. (mB /= am2)) then
        mA = am-am2
        mB = am2
        call interpolate_ddstg_factors()
      end if

      BARE_SUNMASS=SUNMASS/(1+alpha0**2)

      m=am*SUNMASS
      m2=am2*SUNMASS
      m1=m-m2

      bare_m=am*BARE_SUNMASS
      bare_m2=am2*BARE_SUNMASS
      bare_m1=bare_m-bare_m2      

      arr0=(m/(an**2)) ** (1.d0/3)
      arr = arr0

      XA = m1/m
      XB = m2/m
      GAB = (1d0+alphaA*alphaB)/(1d0+alpha0**2)
      nu = XA * XB

c     making some aliases
      k = xk
      n = an
      e = ecc

c5	continue
c      arrold = arr
c      arr=arr0*(1+(m1*m2/m**2 - 9)*0.5d0*m/arr) ** (2.d0/3)
c      if (dabs((arr-arrold)/arr).gt.ARRTOL) goto 5
c      arr=arr0*(1+(m1*m2/m**2 - 9)*0.5d0*m/arr) ** (2.d0/3)
c      arr_it = arr 

      gamma_AB = 1 - 2*alphaA*alphaB/(1+alphaA*alphaB)                  !Damour2009 76
      eps = 2*gamma_AB + 1                                              !Damour Tailor1992 near 3.4

c      W_mult = 1d0/(1d0 - 1d0/6d0*(5*eps + 3.0 - 2*nu) * (GAB*m*an) ** (2d0/3.0))

      arr = (GAB*m/an**2) ** (1d0/3) *(1 - 1d0/6*(5*eps + 3 - 2*nu) *   !Damour Tailor1992 3.9 correct
     +     (GAB*m*an) ** (2d0/3)  )

c      write(*,*) "arr st vs gr: ", arr, arr_it
c      arr = arr_it

      ar=arr*m2/m                                                       !correct

      si=x/ar                                                           !the exact value

c      if (si > 1.0) then
c            si = 1.0
c      end if


c     si = an*x/XB * (GAB*m*an) ** (-1.d0/3)                            !Damour1996 5.8 approximation

c      write(*,*) "si st vs gr: ", si, x/(arr_it*m2/m)

C     xk=3*m/(arr*(1-ecc**2))                 
C     IHS based on 060327 changing to non-tw89 defn for omdot.
      xk=3d0/(1-ecc**2)*(GAB*m*an) ** (2.d0/3) *                          !Damour1996 5.3
     2 ((1-1d0/3*alphaA*alphaB)/(1+alphaA*alphaB) - 
     3 (XA*betaB*alphaA**2 + XB*betaA*alphaB**2) /
     4 (6d0 * (1+alphaA*alphaB)**2)  )

      k = xk

c      xk = xk * W_mult
c      write(*,*) "xk st vs gr: ", xk, 3*m/(arr0*(1-ecc**2))
c      xk=3*m/(arr0*(1-ecc**2))         

C     IHS 100420 changing gamma usage as well as per Norbert's advice
C	gamma=ecc*m2*(m1+2*m2)/(an*arr*m)
      gamma=ecc/an*XB/(1+alphaA*alphaB)*(GAB*m*an) ** (2.d0/3) *        !Damour1996 4.9
     + (XB*(1+alphaA*alphaB) + 1d0 + alphaB*kA)

c      gamma = xk * W_mult
c      write(*,*) "gamma st vs gr: ", gamma, ecc*m2*(m1+2*m2)/(an*arr0*m)
c      gamma=ecc*m2*(m1+2*m2)/(an*arr0*m)

      pbdot_m = -1.5d0*twopi/(1+alphaA*alphaB)*nu*(GAB*m*an)**(5.d0/3) *!Damour1992 6.52a
     + ecc**2*(1+ecc**2/4) * (1-ecc**2)**(-3.5d0) *
     + (5d0/3*(alphaA+alphaB) - 2d0/3*(alphaA*XA + alphaB*XB) +
     + (betaA*alphaB + betaB*alphaA)/(1+alphaA*alphaB)) ** 2

      pbdot_d = -twopi/(1+alphaA*alphaB)*nu*(GAB*m*an) *                !Damour1992 6.52b
     + (1+ecc**2/2) * (1-ecc**2)**(-2.5d0)*(alphaA - alphaB) **2 -
     + 2*twopi/(1+alphaA*alphaB)*nu*(GAB*m*an) ** (5.d0/3) * 
     + (1-ecc**2)**(-3.5d0) * (8d0/5 *(1 + (31d0/8)*ecc**2 + 
     +  (19d0/32)*ecc**4) * (alphaA-alphaB) * (alphaA*XA + alphaB*XB) *
     + (XA - XB) + (1 + 3*ecc**2 + (3d0/8)*ecc**4) * (alphaA-alphaB) * 
     + (betaB*alphaA*XA - betaA*alphaA*XB) / (1+alphaA*alphaB))

      pbdot_d1 = -twopi/(1+alphaA*alphaB)*nu*(GAB*m*an) *               
     + (1+ecc**2/2) * (1-ecc**2)**(-2.5d0)*(alphaA - alphaB) **2

      pbdot_d2 = -2*twopi/(1+alphaA*alphaB)*nu*(GAB*m*an) ** (5.d0/3) * 
     + (1-ecc**2)**(-3.5d0) * (8d0/5 *(1 + (31d0/8)*ecc**2 + 
     +  (19d0/32)*ecc**4) * (alphaA-alphaB) * (alphaA*XA + alphaB*XB) *
     + (XA - XB))

      pbdot_d3 = -2*twopi/(1+alphaA*alphaB)*nu*(GAB*m*an) ** (5.d0/3) * 
     + (1-ecc**2)**(-3.5d0) * 
     + ((1 + 3*ecc**2 + (3d0/8)*ecc**4) * (alphaA-alphaB) * 
     + (betaB*alphaA*XA - betaA*alphaA*XB) / (1+alphaA*alphaB))

      pbdot_qphi = -16d0/5*twopi/(1+alphaA*alphaB) * nu *               !Damour1992 6.52c
     + (GAB*m*an)**(5.d0/3) *(1+ (73.d0/24)*ecc**2 + (37.d0/96)*ecc**4)*
     + (1-ecc**2)**(-3.5d0) * 
     + ((alphaA+alphaB) - (alphaA*XA + alphaB*XB)) ** 2

      pbdot_qg = -96d0/5*twopi/(1+alphaA*alphaB) * nu *                 !Damour1992 6.52d
     + (GAB*m*an)**(5.d0/3) *(1+ (73.d0/24)*ecc**2 + (37.d0/96)*ecc**4)*
     + (1-ecc**2)**(-3.5d0)

      pbdot = pbdot_m + pbdot_d + pbdot_qphi + pbdot_qg
      pbdot_old=-(96*twopi/5) * an**(5.d0/3) * (1-ecc**2)**(-3.5d0) * 
     +    (1+ (73.d0/24)*ecc**2 + (37.d0/96)*ecc**4) * 
     +    m1*m2*m**(-1.d0/3)

c      write(*,*) "pbdot st vs gr: ",pbdot, pbdot_old
c      write(*,*) "pbdot_m, pbdot_d: ",pbdot_m, pbdot_d
c      write(*,*) "pbdot_qphi, pbdot_qg: ",pbdot_qphi, pbdot_qg
c      pbdot = pbdot_old

      dr = (GAB*m*an) ** (2.d0/3) / (m**2 * (1+alphaA*alphaB)) * 
     + ((3-alphaA*alphaB)*m1**2 + (6-alphaA*alphaB-kA*alphaB)*m1*m2 + 
     + (2-alphaA*alphaB-kA*alphaB)*m2**2)

      dth = (GAB*m*an) ** (2.d0/3) / (m**2 * (1+alphaA*alphaB)) * 
     + ((3.5d0-0.5*alphaA*alphaB)*m1**2 + (6-alphaA*alphaB-kA*alphaB)*m1*m2 + 
     + (2-alphaA*alphaB-kA*alphaB)*m2**2)

c      write(*,*) "dr st vs gr: ",dr, (3*m1**2 + 6*m1*m2 + 2*m2**2)/(arr*m)
c      write(*,*) "dth st vs gr: ",dth, (3.5d0*m1**2 + 6*m1*m2 + 2*m2**2)/(arr*m)

        
c     Now caclulating derivatives

c     here are approximated expressions for derivatives

      dk_dm = k * 2d0/(3*m)
      dgamma_dm = gamma * (-4d0/(3*m) + 1d0/(m+m2))
      ddr_dm = dr * (-4d0/(3*m) - 6*m/(-3*m**2+m2**2))
      ddth_dm = dth * (-4d0/(3*m) - 2*(7*m-m2)/(-7*m**2+2*m*m2+m2**2))
      dpbdot_dm = pbdot * (2*m+m2)/(3*m**2-3*m*m2)
      dsi_dm = si * 2d0/(3*m)

      dk_dm_ap = k * 2d0/(3*m)
      dgamma_dm_ap = gamma * (-4d0/(3*m) + 1d0/(m+m2))
      ddr_dm_ap = dr * (-4d0/(3*m) - 6*m/(-3*m**2+m2**2))
      ddth_dm_ap = dth * (-4d0/(3*m) - 2*(7*m-m2)/(-7*m**2+2*m*m2+m2**2))
      dpbdot_dm_ap = pbdot * (2*m+m2)/(3*m**2-3*m*m2)
      dsi_dm_ap = si * 2d0/(3*m)

      dpbdot_dm_ap = pbdot_qg * (2*m+m2)/(3*m**2-3*m*m2) + pbdot_d1 * m2/(m**2 - m*m2) + pbdot_d2 * (-7/(3.*m) + 1/m1 + 1/(m - 2*m2)) + pbdot_d3 * (-4/(3.*m) + 1/m1)

      dpbdot_dm = dpbdot_dm_ap

c     here are exact expressions for derivatives

c      dk_dm = k * (-1d0/(3*m) + (-6+alphaA*(2*alphaB*(-2+alphaA*alphaB)+alphaA*alphaB)) / ((-6+alphaA*(2*alphaB*(-2+alphaA*alphaB)+alphaA*alphaB))*m + (alphaB**2*betaA-alphaA**2*betaB)*m2))
c      dgamma_dm = gamma * (-4d0/(3*m) + (1d0+alphaB*kA)/(m*(1+alphaB*kA)+m2*(1+alphaA*alphaB))) 
c      ddr_dm = dr * (-4d0/(3*m) + (2*(-3+alphaA*alphaB)*m + alphaB*(-alphaA+kA)*m2) / ((-3+alphaA*alphaB)*m**2 + alphaB*(-alphaA+kA)*m*m2 + (1+alphaA*alphaB)*m2**2))
c      ddth_dm = dth * (-4d0/(3*m) + 2*((-7+alphaA*alphaB)*m + (1+alphaB+kA)*m2) / ((-7+alphaA*alphaB)*m**2 + 2*(1+alphaB+kA)*m*m2 + (1+alphaA*alphaB)*m2**2))
c      dsi_dm = si * 2d0/(3*m)

c     here are approximated expressions for derivatives

      dk_dm2 = 0d0
      dgamma_dm2 = gamma * (1d0/m2 + 1d0/(m+m2))
      dsi_dm2 = -si/m2
      ddr_dm2 = dr * 2*m2/(-3*m**2+m2**2)
      ddth_dm2 = dth * 2*(m+m2)/(-7*m**2+2*m*m2+m2**2)
      dpbdot_dm2 = pbdot * (1d0/m2 + 1d0/(-m+m2))

      dk_dm2_ap = 0d0
      dgamma_dm2_ap = gamma * (1d0/m2 + 1d0/(m+m2))
      dsi_dm2_ap = -si/m2
      ddr_dm2_ap = dr * 2*m2/(-3*m**2+m2**2)
      ddth_dm2_ap = dth * 2*(m+m2)/(-7*m**2+2*m*m2+m2**2)
      dpbdot_dm2_ap = pbdot * (1d0/m2 + 1d0/(-m+m2))

c     here are exact expressions for derivatives

c      dk_dm2 = k * (alphaB**2*betaA-alphaA**2*betaB) / ((-6+alphaA*(2*alphaB*(-2+alphaA*alphaB)+alphaA*alphaB)*m + (alphaB**2*betaA-alphaA**2*betaB)*m2))
c      dgamma_dm2 = gamma / m2 * (m*(1+alphaB*kA) + 2*m2*(1+alphaA*alphaB)) / (m*(1+alphaB*kA) + m2*(1+alphaA*alphaB))
c      dsi_dm2 = -si/m2
c      ddr_dm2 = dr * (alphaB*(-alphaA+kA)*m + 2*(1+alphaA*alphaB)*m2) / ((-3+alphaA*alphaB)*m**2 + alphaB*(-alphaA+kA)*m*m2 + (1+alphaA*alphaB)*m2**2)
c      ddth_dm2 = dth * (2*((1+alphaB*kA)*m + (1+alphaA*alphaB)*m2)) / ((-7+alphaA*alphaB)*m**2 + 2*(1+alphaB+kA)*m*m2 + (1+alphaA*alphaB)*m2**2)

      if ((exact_der .eqv. .true.) .and. (gr_case .eqv. .false.)) then

c     Caclulating exact derivatives with respect to m

      dk_dm = ((2*(-6 + alphaA*(2*alphaB*(-2 + alphaA*alphaB) + alphaA*betaB))*m - 
     -      alphaB**2*betaA*m2 + alphaA**2*betaB*m2)*
     -    (GAB*m*n)**(2d0/3))/
     -  (6.*(-1 + e**2)*(m + alphaA*alphaB*m)**2)

      dgamma_dm = -(e*m2*(m + alphaB*kA*m + 4*(m2 + alphaA*alphaB*m2))*
     -     (GAB*m*n)**(2d0/3))/(3.*m**3*(n + alphaA*alphaB*n))

      ddr_dm = (((6 - 2*alphaA*alphaB)*m**2 + alphaB*(-alphaA + kA)*m*m2 + 
     -      4*(1 + alphaA*alphaB)*m2**2)*(GAB*m*n)**(2d0/3))/
     -  (3.*(1 + alphaA*alphaB)*m**3)

      ddth_dm = (((7 - alphaA*alphaB)*m**2 + 2*(1 + alphaA*alphaB)*m2**2 + 
     -      m*(m2 + alphaB*kA*m2))*(GAB*m*n)**(2d0/3))/
     -  (3.*(1 + alphaA*alphaB)*m**3)

      dsi_dm = si * 2d0/(3*m)

      dpbdot_m_dm = -(e**2*(4 + e**2)*GAB*m2*((3*alphaA**2*alphaB + 
     - alphaB*(5 + 3*betaA) + alphaA*(3 + 5*alphaB**2 + 3*betaB))*m + 
     -       2*(alphaA - alphaB)*(1 + alphaA*alphaB)*m2)*
     -     (2*(3*alphaA**2*alphaB + alphaB*(5 + 3*betaA) + 
     -          alphaA*(3 + 5*alphaB**2 + 3*betaB))*m**2 + 
     -       (-5*alphaA**2*alphaB + alphaB*(13 + 3*betaA) + 
     -          alphaA*(-5 + 13*alphaB**2 + 3*betaB))*m*m2 + 
     -  14*(alphaA - alphaB)*(1 + alphaA*alphaB)*m2**2)*n*(GAB*m*n)**(2d0/3)*
     -     twopi)/(72.*(1 + alphaA*alphaB)**3*(1 - e**2)**3.5*m**4)
      
      dpbdot_d1_dm = -((alphaA - alphaB)**2*(2 + e**2)*GAB*m2**2*n*twopi)/
     -  (2.*(1 + alphaA*alphaB)*(1 - e**2)**2.5*m**2)
      
      dpbdot_d2_dm =  -((alphaA - alphaB)*(32 + 124*e**2 + 19*e**4)*m2*
     -     (2*alphaA*(m - m2)*(m**2 + 3*m*m2 - 7*m2**2) - 
     -       alphaB*m2*(m**2 - 12*m*m2 + 14*m2**2))*(GAB*m*n)**(5d0/3)*
     -     twopi)/(30.*(1 + alphaA*alphaB)*(1 - e**2)**3.5*m**5)
      
      dpbdot_d3_dm =   -((alphaA - alphaB)*(8 + 3*e**2*(8 + e**2))*m2*
     - (alphaB*betaA*(m - 4*m2)*m2 + 2*alphaA*betaB*(m - m2)*(m + 2*m2))*
     -     (GAB*m*n)**(5d0/3)*twopi)/
     -  (12.*(1 + alphaA*alphaB)**2*(1 - e**2)**3.5*m**4) 
      
      dpbdot_d_dm = dpbdot_d1_dm + dpbdot_d2_dm + dpbdot_d3_dm
      
      dpbdot_qphi_dm = -((96 + 292*e**2 + 37*e**4)*m2*(alphaB*(m - m2) + alphaA*m2)*
     -     (alphaA*m2*(-4*m + 7*m2) + alphaB*(m - m2)*(2*m + 7*m2))*
     -     (GAB*m*n)**(5d0/3)*twopi)/
     -  (90.*(1 + alphaA*alphaB)*(1 - e**2)**3.5*m**5)
      
      dpbdot_qg_dm = -((96 + 292*e**2 + 37*e**4)*m2*(2*m + m2)*(GAB*m*n)**(5d0/3)*
     -     twopi)/(15.*(1 + alphaA*alphaB)*(1 - e**2)**3.5*m**3)
      
      dpbdot_dm = dpbdot_m_dm + dpbdot_d_dm + dpbdot_qphi_dm + 
     + dpbdot_qg_dm

c     Caclulating exact derivatives with respect to m2

      dk_dm2 = ((alphaB**2*betaA - alphaA**2*betaB)*(GAB*m*n)**(2d0/3))/
     -  (2.*(1 + alphaA*alphaB)**2*(-1 + e**2)*m)

      dgamma_dm2 = (e*(m + alphaB*kA*m + 2*(m2 + alphaA*alphaB*m2))*
     -    (GAB*m*n)**(2d0/3))/(m**2*(n + alphaA*alphaB*n))

      dsi_dm2 = -si/m2


      ddr_dm2 = ((alphaB*(alphaA - kA)*m - 2*(1 + alphaA*alphaB)*m2)*
     -    (GAB*m*n)**(2d0/3))/((1 + alphaA*alphaB)*m**2)

      ddth_dm2 = -(((m + alphaB*kA*m + m2 + alphaA*alphaB*m2)*
     -      (GAB*m*n)**(2d0/3))/((1 + alphaA*alphaB)*m**2))

      dpbdot_m_dm2 = -(e**2*(4 + e**2)*((3*alphaA**2*alphaB + alphaB*(5 + 3*betaA) + 
     -          alphaA*(3 + 5*alphaB**2 + 3*betaB))*m + 
     -       2*(alphaA - alphaB)*(1 + alphaA*alphaB)*m2)*
     -     ((3*alphaA**2*alphaB + alphaB*(5 + 3*betaA) + 
     -          alphaA*(3 + 5*alphaB**2 + 3*betaB))*m**2 - 
     -       2*(alphaB*(8 + 8*alphaA*alphaB + 3*betaA) + 3*alphaA*betaB)*m*m2 - 
     -       8*(alphaA - alphaB)*(1 + alphaA*alphaB)*m2**2)*
     -     (GAB*m*n)**(5d0/3)*twopi)/
     -  (24.*(1 + alphaA*alphaB)**3*(1 - e**2)**3.5*m**4)

      dpbdot_d1_dm2 = -((alphaA - alphaB)**2*(2 + e**2)*GAB*(m - 2*m2)*n*twopi)/
     -  (2.*(1 - e**2)**2.5*(m + alphaA*alphaB*m))

      dpbdot_d2_dm2 = -((alphaA - alphaB)*(32 + 124*e**2 + 19*e**4)*
     -     (alphaB*m2*(2*m**2 - 9*m*m2 + 8*m2**2) + 
     -       alphaA*(m - m2)*(m**2 - 7*m*m2 + 8*m2**2))*
     -     (GAB*m*n)**(5d0/3)*twopi)/
     -  (10.*(1 + alphaA*alphaB)*(1 - e**2)**3.5*m**4)

      dpbdot_d3_dm2 = -((alphaA - alphaB)*(8 + 3*e**2*(8 + e**2))*
     -     (alphaA*betaB*(m - 3*m2)*(m - m2) + alphaB*betaA*m2*(-2*m + 3*m2))*
     -     (GAB*m*n)**(5d0/3)*twopi)/
     -  (4.*(1 - e**2)**3.5*m*(m + alphaA*alphaB*m)**2)

      dpbdot_d_dm2 = dpbdot_d1_dm2 + dpbdot_d2_dm2 + dpbdot_d3_dm2

      dpbdot_qphi_dm2 = -((96 + 292*e**2 + 37*e**4)*(alphaB*(m - m2) + alphaA*m2)*
     -     (alphaB*(m - 4*m2)*(m - m2) + alphaA*(3*m - 4*m2)*m2)*
     -     (GAB*m*n)**(5d0/3)*twopi)/
     -  (30.*(1 + alphaA*alphaB)*(1 - e**2)**3.5*m**4)

      dpbdot_qg_dm2 = -((96 + 292*e**2 + 37*e**4)*(m - 2*m2)*(GAB*m*n)**(5d0/3)*
     -     twopi)/(5.*(1 - e**2)**3.5*m*(m + alphaA*alphaB*m))

      dpbdot_dm2 = dpbdot_m_dm2 + dpbdot_d_dm2 + dpbdot_qphi_dm2 + 
     + dpbdot_qg_dm2

      endif  

c      write(*,*) "dk_dm GR/STG:", k * 2d0/(3*m) / dk_dm
c      write(*,*) "dgamma_dm GR/STG:", gamma * (-4d0/(3*m) + 1d0/(m+m2)) / dgamma_dm
c      write(*,*) "dsi_dm GR/STG:", si * 2d0/(3*m) / dsi_dm
c      write(*,*) "ddr_dm GR/STG:", dr * (-4d0/(3*m) - 6*m/(-3*m**2+m2**2)) / ddr_dm
c      write(*,*) "ddth_dm GR/STG:", dth * (-4d0/(3*m) - 2*(7*m-m2)/(-7*m**2+2*m*m2+m2**2)) / ddth_dm
c      write(*,*) "dpbdot_dm GR/STG:", pbdot * (2*m+m2)/(3*m**2-3*m*m2) / dpbdot_dm
c      write(*,*) "dpbdot_dm GR/STG:", dpbdot_dm_ap / dpbdot_dm

c      write(*,*) "dk_dm2 GR/STG:", 0d0 / dk_dm2
c      write(*,*) "dgamma_dm2 GR/STG:", gamma * (1d0/m2 + 1d0/(m+m2)) / dgamma_dm2
c      write(*,*) "dsi_dm2 GR/STG:", -si/m2 / dsi_dm2
c      write(*,*) "ddr_dm2 GR/STG:", dr * 2*m2/(-3*m**2+m2**2) / ddr_dm2
c      write(*,*) "ddth_dm2 GR/STG:", dth * 2*(m+m2)/(-7*m**2+2*m*m2+m2**2) / ddth_dm2
c      write(*,*) "dpbdot_dm2 GR/STG:", pbdot * (1d0/m2 + 1d0/(-m+m2)) / dpbdot_dm2 
c      write(*,*) "dpbdot_dm2 GR/STG:", dpbdot_dm2_ap / dpbdot_dm2    

      if ((exact_factor_der .eqv. .true.) .and. (gr_case .eqv. .false.)  .and. (beta0 /= 0.0d0) .and. (companion_type(1:2).ne.'BH')) then

      dgamma_dalphaA = gamma * (-(alphaB*(m + alphaB*kA*m - 2*(m2 + alphaA*alphaB*m2)))/
     -  (3.*(1 + alphaA*alphaB)*(m + alphaB*kA*m + m2 + alphaA*alphaB*m2)))
      dgamma_dbetaA = 0d0
      dgamma_dkA = gamma * (alphaB*m)/(m + alphaB*kA*m + m2 + alphaA*alphaB*m2)
      dgamma_dm_factor = dalphaA_dm * dgamma_dalphaA + dbetaB_dm * dgamma_dbetaA + dkA_dm * dgamma_dkA
      dgamma_dm2_factor = - dgamma_dm_factor

      dk_dalphaA = k * (2*(3 + alphaA*alphaB)*(2*alphaB*(1 + alphaA*alphaB) + alphaA*betaB)*m - 
     -    2*(2*alphaB**3*betaA + alphaA*(3 + alphaA*alphaB)*betaB)*m2)/
     -  (3.*(1 + alphaA*alphaB)*((-6 + alphaA*(2*alphaB*(-2 + alphaA*alphaB) + alphaA*betaB))*m + 
     -      (alphaB**2*betaA - alphaA**2*betaB)*m2))
      dk_dbetaA = k * (alphaB**2*m2)/((-6 + alphaA*(2*alphaB*(-2 + alphaA*alphaB) + alphaA*betaB))*m + 
     -    (alphaB**2*betaA - alphaA**2*betaB)*m2)
      dk_dkA = 0d0
      dk_dm_factor = dalphaA_dm * dk_dalphaA + dbetaB_dm * dk_dbetaA + dkA_dm * dk_dkA
      dk_dm2_factor = - dk_dm_factor

      dsi_dalphaA = si * (-(alphaB/(3 + 3*alphaA*alphaB)))
      dsi_dbetaA = 0d0
      dsi_dkA = 0d0
      dsi_dm_factor = dalphaA_dm * dsi_dalphaA + dbetaB_dm * dsi_dbetaA + dkA_dm * dsi_dkA
      dsi_dm2_factor = - dsi_dm_factor

      ddr_dalphaA = dr * (alphaB*(2*(3 + alphaA*alphaB)*m**2 - (3 + 2*alphaA*alphaB + alphaB*kA)*m*m2 + 
     -      2*(1 + alphaA*alphaB)*m2**2))/
     -  (3.*(1 + alphaA*alphaB)*((-3 + alphaA*alphaB)*m**2 + alphaB*(-alphaA + kA)*m*m2 + 
     -      (1 + alphaA*alphaB)*m2**2))
      ddr_dbetaA = 0d0
      ddr_dkA = dr * (alphaB*m*m2)/((-3 + alphaA*alphaB)*m**2 + alphaB*(-alphaA + kA)*m*m2 + (1 + alphaA*alphaB)*m2**2)
      ddr_dm_factor = dalphaA_dm * ddr_dalphaA + dbetaB_dm * ddr_dbetaA + dkA_dm * ddr_dkA
      ddr_dm2_factor = - ddr_dm_factor

      ddth_dalphaA = dth * (2*alphaB*((5 + alphaA*alphaB)*m**2 + (1 + alphaA*alphaB)*m2**2 - m*(m2 + alphaB*kA*m2)))/
     -  (3.*(1 + alphaA*alphaB)*((-7 + alphaA*alphaB)*m**2 + (1 + alphaA*alphaB)*m2**2 + 2*m*(m2 + alphaB*kA*m2)))
      ddth_dbetaA = 0d0
      ddth_dkA = dth * (2*alphaB*m*m2)/((-7 + alphaA*alphaB)*m**2 + (1 + alphaA*alphaB)*m2**2 + 2*m*(m2 + alphaB*kA*m2))
      ddth_dm_factor = dalphaA_dm * ddth_dalphaA + dbetaB_dm * ddth_dbetaA + dkA_dm * ddth_dkA
      ddth_dm2_factor = - ddth_dm_factor

      dpbdot_m_dalphaA = pbdot_m * (2*(9*(1 + betaB) + alphaB*(12*alphaA**2*alphaB + 5*alphaA*alphaB**2 + alphaB*(5 - 6*betaA) + 
     -          3*alphaA*(7 + betaB)))*m + 4*(1 + alphaA*alphaB)*(3 + 4*alphaA*alphaB - alphaB**2)*m2)/
     -  (3.*(1 + alphaA*alphaB)*((3*alphaA**2*alphaB + alphaB*(5 + 3*betaA) + alphaA*(3 + 5*alphaB**2 + 3*betaB))*
     -       m + 2*(alphaA - alphaB)*(1 + alphaA*alphaB)*m2))
      dpbdot_m_dbetaA = pbdot_m * (6*alphaB*m)/((3*alphaA**2*alphaB + alphaB*(5 + 3*betaA) + alphaA*(3 + 5*alphaB**2 + 3*betaB))*m + 
     -    2*(alphaA - alphaB)*(1 + alphaA*alphaB)*m2)
      dpbdot_m_dkA = 0.0
      dpbdot_m_dm_factor = dalphaA_dm * dpbdot_m_dalphaA + dbetaB_dm * dpbdot_m_dbetaA + dkA_dm * dpbdot_m_dkA
      dpbdot_m_dm2_factor = - dpbdot_m_dm_factor

      dpbdot_d1_dalphaA = pbdot_d1 * 2/(alphaA - alphaB)
      dpbdot_d1_dbetaA = 0.0
      dpbdot_d1_dkA = 0.0
      dpbdot_d1_dm_factor = dalphaA_dm * dpbdot_d1_dalphaA + dbetaB_dm * dpbdot_d1_dbetaA + dkA_dm * dpbdot_d1_dkA
      dpbdot_d1_dm2_factor = - dpbdot_d1_dm_factor

      dpbdot_d2_dalphaA = pbdot_d2 * ((6 + 8*alphaA*alphaB - 2*alphaB**2)/(1 + alphaA*alphaB) - (3*alphaB*m)/(alphaA*m1 + alphaB*m2))/
     -  (3.*(alphaA - alphaB))
      dpbdot_d2_dbetaA = 0.0
      dpbdot_d2_dkA = 0.0
      dpbdot_d2_dm_factor = dalphaA_dm * dpbdot_d2_dalphaA + dbetaB_dm * dpbdot_d2_dbetaA + dkA_dm * dpbdot_d2_dkA
      dpbdot_d2_dm2_factor = - dpbdot_d2_dm_factor

      dpbdot_d3_dalphaA = pbdot_d3 * (-(3*alphaB*betaB*m1 + 2*alphaA*(-3 + alphaB**2)*betaB*m1 + 2*alphaA*alphaB**2*betaA*m2 + 
     -     alphaB*(3 + alphaB**2)*betaA*m2 + 5*alphaA**2*alphaB*betaB*(-m + m2))/
     -  (3.*(alphaA - alphaB)*(1 + alphaA*alphaB)*(alphaA*betaB*m1 - alphaB*betaA*m2)))
      dpbdot_d3_dbetaA = pbdot_d3 * (alphaB*m2)/(alphaB*betaA*m2 + alphaA*betaB*(-m + m2))
      dpbdot_d3_dkA = 0.0
      dpbdot_d3_dm_factor = dalphaA_dm * dpbdot_d3_dalphaA + dbetaB_dm * dpbdot_d3_dbetaA + dkA_dm * dpbdot_d3_dkA
      dpbdot_d3_dm2_factor = - dpbdot_d3_dm_factor

      dpbdot_d_dm_factor = dpbdot_d1_dm_factor + dpbdot_d2_dm_factor + dpbdot_d3_dm_factor
      dpbdot_d_dm2_factor = dpbdot_d1_dm2_factor + dpbdot_d2_dm2_factor + dpbdot_d3_dm2_factor

      dpbdot_qphi_dalphaA = pbdot_qphi * (2*alphaB**2*m1 + 6*m2 + 8*alphaA*alphaB*m2)/(3.*(1 + alphaA*alphaB)*(alphaB*m1 + alphaA*m2))
      dpbdot_qphi_dbetaA = 0.0
      dpbdot_qphi_dkA = 0.0
      dpbdot_qphi_dm_factor = dalphaA_dm * dpbdot_qphi_dalphaA + dbetaB_dm * dpbdot_qphi_dbetaA + dkA_dm * dpbdot_qphi_dkA
      dpbdot_qphi_dm2_factor = - dpbdot_qphi_dm_factor

      dpbdot_qg_dalphaA = pbdot_qg * (2*alphaB)/(3 + 3*alphaA*alphaB)
      dpbdot_qg_dbetaA = 0.0
      dpbdot_qg_dkA = 0.0
      dpbdot_qg_dm_factor = dalphaA_dm * dpbdot_qg_dalphaA + dbetaB_dm * dpbdot_qg_dbetaA + dkA_dm * dpbdot_qg_dkA
      dpbdot_qg_dm2_factor = - dpbdot_qg_dm_factor

      dpbdot_dm_factor = dpbdot_m_dm_factor + dpbdot_d_dm_factor + dpbdot_qphi_dm_factor + dpbdot_qg_dm_factor
      dpbdot_dm2_factor = dpbdot_m_dm2_factor + dpbdot_d_dm2_factor + dpbdot_qphi_dm2_factor + dpbdot_qg_dm2_factor

      dgamma_dm = dgamma_dm + dgamma_dm_factor
      dgamma_dm2 = dgamma_dm2 + dgamma_dm2_factor

      dk_dm = dk_dm + dk_dm_factor
      dk_dm2 = dk_dm2 + dk_dm2_factor

      dsi_dm = dsi_dm + dsi_dm_factor
      dsi_dm2 = dsi_dm2 + dsi_dm2_factor

      ddr_dm = ddr_dm + ddr_dm_factor
      ddr_dm2 = ddr_dm2 + ddr_dm2_factor

      ddth_dm = ddth_dm + ddth_dm_factor
      ddth_dm2 = ddth_dm2 + ddth_dm2_factor

      dpbdot_dm = dpbdot_dm + dpbdot_dm_factor
      dpbdot_dm2 = dpbdot_dm2 + dpbdot_dm2_factor

      endif

c      dk_dm2 = dk_dm2_ap
c      dgamma_dm2 = dgamma_dm2_ap
c      dsi_dm2 = dsi_dm2_ap
c      ddr_dm2 = ddr_dm2_ap 
c      ddth_dm2 = ddth_dm2_ap
c      dpbdot_dm2 = dpbdot_dm2_ap

c      dk_dm = dk_dm_ap
c      dgamma_dm = dgamma_dm_ap
c      dsi_dm = dsi_dm_ap
c      ddr_dm = ddr_dm_ap 
c      ddth_dm = ddth_dm_ap
c      dpbdot_dm = dpbdot_dm_ap
      
C      write(*,*) "mtot:", am
C      write(*,*) "m2:", am2
C      write(*,*) "alphaA:", alphaA
C      write(*,*) "betaA:", betaA
C      write(*,*) "k:", k
C      write(*,*) "gamma:", gamma
C      write(*,*) "si", si
C      write(*,*) "dr:", dr
C      write(*,*) "dth:", dth
C      write(*,*) "pbdot:", pbdot

C      write(*,*) "dk_dm:", dk_dm * SUNMASS
C      write(*,*) "dgamma_dm:", dgamma_dm * SUNMASS
C      write(*,*) "dsi_dm", dsi_dm * SUNMASS
C      write(*,*) "ddr_dm:", ddr_dm * SUNMASS
C      write(*,*) "ddth_dm:", ddth_dm * SUNMASS
C      write(*,*) "dpbdot_dm:", dpbdot_dm * SUNMASS

C      write(*,*) "dk_dm2:", dk_dm2 * SUNMASS
C      write(*,*) "dgamma_dm2:", dgamma_dm2 * SUNMASS
C      write(*,*) "dsi_dm2:", dsi_dm2 * SUNMASS
C      write(*,*) "ddr_dm2:", ddr_dm2 * SUNMASS
C      write(*,*) "ddth_dm2:", ddth_dm2 * SUNMASS
C      write(*,*) "dpbdot_dm2:", dpbdot_dm2 * SUNMASS
  

        return
      end
