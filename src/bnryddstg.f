c      $Id$
      subroutine bnryddstg(ct,f0,n,torb,fctn)

c  Damour et Deruelle modele pour le chronometrage des temps d'arrive 
c  au barycentre du systeme solaire, a la premiere approximation
c  post-newtonienne de la Relativite Generale.

c  Computes pulsar orbit time, torb, at time of observation t=ct(n)-pepoch.
c  Pulsar proper time is then TP=T+TORB
c  Units are such that c=g=1. Thus masses have units of seconds, with
c  one solar mass = 4.925 490 947d-6 sec.

c  Also computes the binary orbit-related values of fctn: partial
c  derivatives of each arrival time residual with respect to the model
c  parameters.

c  Initial guesses for all parameters must be placed in common/orbit/ by the
c  calling program.  

      implicit real*8 (a-h,o-z)
      include 'dim.h'
      real*8 fctn(NPAP1),k,m,m1,m2,bare_m, bare_m1, bare_m2
      real*8 BARE_SUNMASS,SUNMASS
      parameter (TWOPI=6.28318530717958648d0,SUNMASS=4.925490947d-6)
      parameter (RAD=360.d0/twopi)
      parameter (ARRTOL=1.d-10)
      include 'orbit.h'
      include 'ddstg.h'
      integer loop_counter1, loop_counter2, fake_counter

      BARE_SUNMASS=SUNMASS/(1+alpha0**2)

      tt0=(ct-t0(1))*86400.d0
      an=twopi/pb(1)
      x=a1(1)+xdot*tt0
      ecc=e(1)+edot*tt0
      m=am*SUNMASS
      m2=am2*SUNMASS
      m1=m-m2
      bare_m=am*BARE_SUNMASS
      bare_m2=am2*BARE_SUNMASS
      bare_m1=bare_m-bare_m2

      call mass2ddstg(am,am2,x,ecc,an,arr,ar,xk,si,gamma,pbdot,dr,dth)

c      omdot=360.d0*365.25d0*xk/pb(1)*86400.d0                           !instead of newbin.f
c      write(*,*) omdot

      fake_counter = 0

      k=xk
      frb = 1.d0/pb(1)
      xx = x

c  DD equations 36, 37:
c      dr=(3*m1**2 + 6*m1*m2 + 2*m2**2)/(arr*m)                          !calculated in mass2ddstg.f
      er=ecc*(1+dr)
c      dth=(3.5d0*m1**2 + 6*m1*m2 + 2*m2**2)/(arr*m)                     !calculated in mass2ddstg.f
      eth=ecc*(1+dth)

c  --> Inversion of timing model by iteration: begin of loop
        loop_counter1 = 0
        epsNum = 1.0d-10
        delta  = 0.0d0
 10     delta_old = delta        
        tt  = tt0 - delta
        orbits  = tt*frb - 0.5d0*(pbdot+xpbdot)*(tt*frb)**2
        norbits = orbits
        if(orbits.lt.0.d0) norbits = norbits - 1
        phase = twopi*(orbits - norbits)

        loop_counter1 = loop_counter1 + 1

c  Compute eccentric anomaly u by iterating Keplers equation.
      loop_counter2 = 0
      u=phase+ecc*dsin(phase)*(1+ecc*dcos(phase))
100   du=(phase-(u-ecc*dsin(u)))/(1.d0-ecc*dcos(u))
      u=u+du
      loop_counter2 = loop_counter2 + 1
      if ((dabs(du).gt.1.d-14) .and. (loop_counter2<50)) go to 100

c  DD equations 17b, 17c, 29, and 46 through 52
      su=dsin(u)
      cu=dcos(u)
      onemecu=1.d0-ecc*cu
      cume    = cu - ecc

      cae=(cu-ecc)/onemecu
      sae=dsqrt(1.d0-ecc**2)*su/onemecu
      ae1=datan2(sae,cae)
      if(ae1.lt.0.0) ae1=ae1+twopi
      ae=twopi*orbits + ae1 - phase

          omega=omz(1)/RAD + k*ae 

        sw=dsin(omega)
        cw=dcos(omega)

        psi  = omega + ae1 ! angle w.r.t. ascending node
        spsi = DSIN(psi)
        cpsi = DCOS(psi)

c  Roemer delay (DD)
        alpha = xx*sw
        beta  = xx*DSQRT(1.d0 - eth**2)*cw
        dRoe  = alpha*(cu - er) + beta*su

c  Einstein delay (DD)
        dEin  = gamma * su

c  Shapiro delay

        sqr1me2 = DSQRT(1.d0 - ecc**2)

        brace  = onemecu - si*(sw*cume + sqr1me2*cw*su)

        dlogbr = DLOG(brace)
        dSha   = -2.d0*bare_m2*dlogbr                                   !change m2 to bare_m2

c  Aberration
          a0aligned=an*ar/(twopi*f0*si*sqr1me2)                         !added a0 and b0
          a0=afac*a0aligned
          b0=0.d0

        dAbe = a0*(spsi+ecc*sw)+b0*(cpsi+ecc*cw)

        delta = dRoe + dEin + dSha + dAbe

        diff  = DABS(delta - delta_old)
        if ((diff.gt.epsNum) .and. (loop_counter1 .lt. 50)) goto 10
c  --> Inversion of timing model by iteration: end of loop
      
        torb = -delta

c  Compute eccentric anomaly u by iterating Kepler's equation.

      if (.false.) then

      orbits=tt0/pb(1) - 0.5d0*(pbdot+xpbdot)*(tt0/pb(1))**2
      norbits=orbits
      if(orbits.lt.0.d0) norbits=norbits-1
      phase=twopi*(orbits-norbits)
      u=phase+ecc*dsin(phase)*(1+ecc*dcos(phase))
11    fac=1/(1-ecc*dcos(u))
      du=fac*(phase-(u-ecc*dsin(u)))
      u=u+du
      if(dabs(du).gt.1.d-14) go to 11

C  DD equations 17a, 29:
      ae=2*datan(dsqrt((1+ecc)/(1-ecc))*dtan(0.5d0*u))
      if(ae.lt.0.0) ae=ae+twopi
      if(dabs(ae-phase).gt.3.14d0)  then
        print *,ae,phase,orbits,ct
        stop 'ae error'
        endif
      ae=twopi*orbits + ae - phase
      omega=omz(1)/RAD + (k + xomdot/(an*RAD*365.25d0*86400.d0))*ae

C  DD equations 46 through 52.
      su=dsin(u)
      cu=dcos(u)
      sw=dsin(omega)
      cw=dcos(omega)
      alpha=x*sw
      beta=x*dsqrt(1-eth**2)*cw
      bg=beta+gamma
      dre=alpha*(cu-er) + bg*su
      drep=-alpha*su + bg*cu
      drepp=-alpha*cu - bg*su
      onemecu=1-ecc*cu
      anhat=an/onemecu
c  DD equations 26, 27, 57:
      cume=cu-ecc
      sqr1me2=dsqrt(1-ecc**2)
      brace=onemecu-si*(sw*cume+sqr1me2*cw*su)
      dlogbr=dlog(brace)
      ds=-2*bare_m2*dlogbr
C  These will be different if spin axis not aligned:
      a0aligned=an*ar/(twopi*f0*si*sqr1me2)
      a0=afac*a0aligned
      b0=0.d0
      da=a0*(sin(omega+ae) + ecc*sw) + b0*(cos(omega+ae) + ecc*cw)
c  Now compute d2bar, the orbital time correction in DD equation 42.
      d2bar=dre*(1-anhat*drep+(anhat**2)*(drep**2 + 0.5d0*dre*drepp -
     +    0.5d0*ecc*su*dre*drep/onemecu)) + ds + da
      torb=-d2bar

      endif

c      write(*,*) "torb = ",  torb
c      write(*,*) "u = ",  u
c      write(*,*) "ae = ",  ae
c      write(*,*) "brace = ",  brace
c      write(*,*) "omega = ",  omega
c      write(*,*) "dlogbr = ",  dlogbr
c      write(*,*) "ds = ",  dSha
c      write(*,*) "da = ",  dAbe
c      write(*,*) "tt0 = ",  tt0
c      write(*,*) "phase = ",  phase

      ddSha_domega = -2.d0*bare_m2/brace * si * (-(cu - ecc)*cw + sqr1me2*su*sw)

c  Now we need the partial derivatives. Use DD equations 62a - 62k.
      csigma=x*(-sw*su+sqr1me2*cw*cu)/onemecu
      ce=su*csigma-x*sw-ecc*x*cw*su/sqr1me2
      cx=sw*cume+sqr1me2*cw*su
      comega=x*(cw*cume-sqr1me2*sw*su)
      cgamma=su
      csini=2*bare_m2*(sw*cume+sqr1me2*cw*su)/brace                     !made a change from m2 to bare_m2
      ck=ae*comega
      cdr=-ecc*x*sw
      cdth=-ecc**2*x*cw*su/sqr1me2
      cpbdot=-csigma*an*tt0**2/(2*pb(1))

      ddt_dm2=-2*dlogbr * bare_m2/m2                                    !normally is cr * dr_dm2

      cm2 = ddt_dm2 + ck*dk_dm2 + cgamma*dgamma_dm2 + cdr*ddr_dm2 + 
     + cdth*ddth_dm2 + cpbdot*dpbdot_dm2 + csini*dsi_dm2

      cm = ck*dk_dm + cgamma*dgamma_dm + cdr*ddr_dm + cdth*ddth_dm + 
     + cpbdot*dpbdot_dm + csini*dsi_dm

*      write(*,*) "counter", fake_counter
*      write(*,*) "torb = ", -torb
*      write(*,*) "m2", m2
*      write(*,*) "m", m
*      write(*,*) "ecc = ", ecc
*      write(*,*) "x = ", x
*      write(*,*) "gamma = ",  gamma
*      write(*,*) "sigma = ",  tt0*an
*      write(*,*) "omega = ",  omega
*      write(*,*) "sini = ",  si
*      write(*,*) "k = ",  k
*      write(*,*) "dr = ",  dr
*      write(*,*) "dth = ",  dth
*      write(*,*) "pbdot = ",  pbdot
*
*      write(*,*) "ce = ", ce
*      write(*,*) "cx = ", cx
*      write(*,*) "cgamma = ",  cgamma
*      write(*,*) "csigma = ",  csigma
*      write(*,*) "comega = ",  comega
*      write(*,*) "csini = ",  csini
*      write(*,*) "ck = ",  ck
*      write(*,*) "cdr = ",  cdr
*      write(*,*) "cdth = ",  cdth
*      write(*,*) "cpbdot = ",  cpbdot
*      write(*,*) "ddt_dm2 = ",  ddt_dm2
*      write(*,*) "cm: ", cm
*      write(*,*) "cm2: ", cm2

c      tt0 = tt0 * 1.0000001
c      fake_counter = fake_counter + 1
c      if(fake_counter .le. 1) goto 12

      fctn(9)=cx*f0
      fctn(10)=ce*f0
      fctn(11)=-csigma*an*f0*86400.d0
      fctn(12)=fctn(11)*tt0/(pb(1)*86400.d0)
      fctn(13)=comega*f0
c following used to be fctn(14), omdot, but it is really xomdot.
      fctn(37)=ae*fctn(13)/(an*RAD*365.25d0*86400.d0)  
c following used to be fctn(18), pbdot, but it is really xpbdot
      fctn(38)=tt0*fctn(12)*0.5d-6
      fctn(21)=cm*f0*SUNMASS
      fctn(22)=cm2*f0*SUNMASS
      fctn(24)=tt0*fctn(9)
      fctn(25)=tt0*fctn(10)


      return
      end
