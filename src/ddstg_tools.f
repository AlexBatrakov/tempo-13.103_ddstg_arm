      subroutine get_ns_grids()
      include 'ddstg.h'
      integer :: i, j, N_pc, N_a0, N_b0, i_a0=0, i_b0=0, i_bl
      integer, dimension(1:4) :: bl_rows
      real*8 :: i_a0_real, i_b0_real, x, y
      real*8, dimension(:), allocatable :: a0_grid, b0_grid
      real*8, dimension(:,:), allocatable :: massA_bl, alphaA_bl, 
     + betaA_bl, kA_bl
      character :: path*640, path_eos*640

      call getenv('TEMPO',path)

      path_eos = trim(path)//'/data_ddstg/'//trim(eosname)

      open (unit=99, file=trim(path_eos)//'/dims.dat', status='old')
      read(99, *)
      read(99, *) N_pc, N_a0, N_b0
      close(99)
c      write(*, *) "N_pc, N_a0, N_b0 = ", N_pc, N_a0, N_b0

      N_grid = N_pc

      allocate(a0_grid(N_a0), b0_grid(N_b0))

      open (unit=99, file=trim(path_eos)//'/alpha0.dat', status='old')
      do i = 1, N_a0
        read(99, *) a0_grid(i)
        if (a0_grid(i) <= alpha0 .and. i < N_a0) then
            i_a0 = i
        end if
      end do
      close(99)

      open (unit=99, file=trim(path_eos)//'/beta0.dat', status='old')
      do i = 1, N_b0
        read(99, *) b0_grid(i)
        if (b0_grid(i) <= beta0 .and. i < N_b0) then
            i_b0 = i
        end if
      end do
      close(99)

c      write(*,*) "a0_st, b0_st = ", a0_st, b0_st
c      write(*, *) "i_a0, i_b0 = ", i_a0, i_b0
c      write(*, *) "a0, b0 grid(i) = ", a0_grid(i_a0), b0_grid(i_b0)
c      write(*, *) "a0, b0 grid(i+1) = ",a0_grid(i_a0+1),b0_grid(i_b0+1)

      allocate(massA_bl(N_pc,4), alphaA_bl(N_pc,4), 
     + betaA_bl(N_pc,4), kA_bl(N_pc,4))

      bl_rows = (/i_b0 + N_b0 * (i_a0-1), i_b0 + 1 + N_b0 * (i_a0-1), 
     + i_b0 + N_b0 * i_a0, i_b0 + 1 + N_b0 * i_a0/)
c      write(*, *) "bl_rows = ", bl_rows

      open (unit=99, file=trim(path_eos)//'/massA.dat', status='old')
      open (unit=98, file=trim(path_eos)//'/alphaA.dat', status='old')
      open (unit=97, file=trim(path_eos)//'/betaA.dat', status='old')
      open (unit=96, file=trim(path_eos)//'/kA.dat', status='old')

      do i = 1, N_a0*N_b0
        if (i > bl_rows(4)) then
            exit
        else if (i /= bl_rows(1) .and. i /= bl_rows(2) .and.  
     + i /= bl_rows(3) .and. i /= bl_rows(4)) then
            read(99, "()", advance = "yes")
            read(98, "()", advance = "yes")
            read(97, "()", advance = "yes")
            read(96, "()", advance = "yes")
        else
            !write(*, *) "in bl", i
            if (i == bl_rows(1)) then
                i_bl = 1
            else if (i == bl_rows(2)) then
                i_bl = 2
            else if (i == bl_rows(3)) then
                i_bl = 3
            else if (i == bl_rows(4)) then
                i_bl = 4
            end if
c            write(*,*) "read line ", i
            read(99, *) i_a0_real, i_b0_real, massA_bl(:,i_bl)
            read(98, *) i_a0_real, i_b0_real, alphaA_bl(:,i_bl)
            read(97, *) i_a0_real, i_b0_real, betaA_bl(:,i_bl)
            read(96, *) i_a0_real, i_b0_real, kA_bl(:,i_bl)
        end if
      end do

      close(99)
      close(98)
      close(97)
      close(96)

      x = (beta0 - b0_grid(i_b0)) / (b0_grid(i_b0+1) - b0_grid(i_b0))
      y = (alpha0 - a0_grid(i_a0)) / (a0_grid(i_a0+1) - a0_grid(i_a0))
c      y = (dlog(dabs(a0_st)) - dlog(dabs(a0_grid(i_a0)))) 
c     + / dlog(a0_grid(i_a0+1) / a0_grid(i_a0))
c      write(*, *) "x, y = ", x, y

      do i = 1, N_pc
        massA_grid(i) = massA_bl(i,1)*(1.0-x)*(1.0-y) + 
     + massA_bl(i,2)*x*(1.0-y) + massA_bl(i,3)*(1.0-x)*y +
     + massA_bl(i,4)*x*y
        alphaA_grid(i) = alphaA_bl(i,1)*(1.0-x)*(1.0-y) + 
     + alphaA_bl(i,2)*x*(1.0-y) + alphaA_bl(i,3)*(1.0-x)*y + 
     + alphaA_bl(i,4)*x*y
        betaA_grid(i) = betaA_bl(i,1)*(1.0-x)*(1.0-y) + 
     + betaA_bl(i,2)*x*(1.0-y) + betaA_bl(i,3)*(1.0-x)*y + 
     + betaA_bl(i,4)*x*y
        kA_grid(i) = kA_bl(i,1)*(1.0-x)*(1.0-y) + kA_bl(i,2)*x*(1.0-y) +
     +  kA_bl(i,3)*(1.0-x)*y + kA_bl(i,4)*x*y
      end do


      return
      end subroutine get_ns_grids

      subroutine interpolate_ddstg_factors()
        include 'ddstg.h'
        real*8 SUNMASS
        parameter (SUNMASS=4.925490947d-6)
c        real*8 mA,mB

c        write(*,*) "interpolation ", mA+mB,mB

        if ((alpha0 == 0.0) .and. (beta0 == 0.0)) then
          gr_case = .true.
          alphaA = 0.0
          betaA = 0.0
          kA = 0.0
          alphaB = 0.0
          betaB = 0.0
          kB = 0.0
          dalphaA_dm = 0.0
          dbetaA_dm = 0.0
          dkA_dm = 0.0
          dalphaB_dm = 0.0
          dbetaB_dm = 0.0
          dkB_dm = 0.0
        else
          gr_case = .false.
          call make_polint(massA_grid,alphaA_grid,N_GRID,mA,alphaA,dalphaA_dm)
          call make_polint(massA_grid,betaA_grid,N_GRID,mA,betaA,dbetaA_dm)
          call make_polint(massA_grid,kA_grid,N_GRID,mA,kA,dkA_dm)

          if (companion_type(1:2).eq.'WD') then
            alphaB = alpha0
            betaB = beta0
            kB = 0.0
            dalphaB_dm = 0.0
            dbetaB_dm = 0.0
            dkB_dm = 0.0
          else if (companion_type(1:2).eq.'NS') then
            call make_polint(massA_grid,alphaA_grid,N_GRID,mB,alphaB,dalphaB_dm)
            call make_polint(massA_grid,betaA_grid,N_GRID,mB,betaB,dbetaB_dm)
            call make_polint(massA_grid,kA_grid,N_GRID,mB,kB,dkB_dm)
          else if (companion_type(1:2).eq.'BH') then
            alphaB = 0.0
            betaB = 0.0
            kB = 0.0
            dalphaB_dm = 0.0
            dbetaB_dm = 0.0
            dkB_dm = 0.0
          endif
        endif
        dalphaA_dm = dalphaA_dm / SUNMASS
        dbetaA_dm = dbetaA_dm / SUNMASS
        dkA_dm = dkA_dm / SUNMASS
        dalphaB_dm = dalphaB_dm / SUNMASS
        dbetaB_dm = dbetaB_dm / SUNMASS
        dkB_dm = dkB_dm / SUNMASS
        return
      end subroutine interpolate_ddstg_factors


      subroutine make_polint(xa,ya,n,x,y,dy_dx)
      integer :: n, j, k, m=2
      REAL*8 dy,x,y,dy_dx,xa(n),ya(n)
      call locate(xa,n,x,j)
      k = min(max(j-(m-1)/2,1),n+1-m)
      dy_dx = (ya(k+1) - ya(k)) / (xa(k+1) - xa(k))
      call polint(xa(k),ya(k),m,x,y,dy)
      return
      end subroutine make_polint

      SUBROUTINE polint(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      REAL*8 dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10) !Largest anticipated value of n.
c      Given arrays xa and ya, each of length n, and given a value x, this routine returns a value y, and an error estimate dy. If P (x) is the polynomial of degree N − 1 such that P(xai) = yai,i = 1,...,n, then the returned value y = P(x).
      INTEGER i,m,ns
      REAL*8 den,dif,dift,ho,hp,w,c(NMAX),d(NMAX) 
      ns=1
      dif=dabs(x-xa(1))
      do i=1,n 
            dift=dabs(x-xa(i))
            if (dift.lt.dif) then
                  ns=i
                  dif=dift
            endif
c      Here we find the index ns of the closest table entry,
            c(i)=ya(i) !and initialize the tableau of c’s and d’s.
            d(i)=ya(i) 
      enddo
      y=ya(ns) 
      ns=ns-1
      do m=1,n-1
            do i=1,n-m 
                  ho=xa(i)-x
c      This is the initial approximation to y. For each column of the tableau,
c      we loop over the current c’s and d’s and update them.
                  hp=xa(i+m)-x
                  w=c(i+1)-d(i)
                  den=ho-hp
                  if(den.eq.0.)pause 'failure in polint'
c      This error can occur only if two input xa’s are (to within roundoff) identical. 
                  den=w/den
                  d(i)=hp*den
                  c(i)=ho*den 
            enddo
            if (2*ns.lt.n-m)then 
                  dy=c(ns+1)
            else
                  dy=d(ns)
                  ns=ns-1 
            endif
            y=y+dy 
      enddo
      return 
      END SUBROUTINE polint

      SUBROUTINE locate(xx,n,x,j) 
      INTEGER j,n
      REAL*8 x,xx(n)
c      Given an array xx(1:n), and given a value x, returns a value j such that x is between xx(j) and xx(j+1). xx(1:n) must be monotonic, either increasing or decreasing. j=0 or j=n is returned to indicate that x is out of range.
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then 
            jm=(ju+jl)/2
            if((xx(n).ge.xx(1)).eqv.(x.ge.xx(jm)))then
                  jl=jm 
            else
                  ju=jm 
            endif
      goto 10
      endif 
      if(x.eq.xx(1))then
            j=1
      else if(x.eq.xx(n))then
            j=n-1 
      else
            j=jl 
      endif
      return 
      END SUBROUTINE locate

      SUBROUTINE hunt(xx,n,x,jlo)
      INTEGER jlo,n
      REAL x,xx(n)
c      Given an array xx(1:n), and given a value x, returns a value jlo such that x is between
c      xx(jlo) and xx(jlo+1). xx(1:n) must be monotonic, either increasing or decreasing.
c      jlo=0 or jlo=n is returned to indicate that x is out of range. jlo on input is taken as
c      the initial guess for jlo on output.
      INTEGER inc,jhi,jm
      LOGICAL ascnd
      ascnd=xx(n).ge.xx(1) !True if ascending order of table, false otherwise.
      if(jlo.le.0.or.jlo.gt.n)then !Input guess not useful. Go immediately to bisection.
        jlo=0
        jhi=n+1
        goto 3
      endif
      inc=1
      if(x.ge.xx(jlo).eqv.ascnd)then !Hunt up:
1       jhi=jlo+inc
        if(jhi.gt.n)then! Done hunting, since o end of table.
            jhi=n+1
        else if(x.ge.xx(jhi).eqv.ascnd)then! Not done hunting,
            jlo=jhi
            inc=inc+inc !so double the increment
            goto 1! and try again.
        endif! Done hunting, value bracketed.
      else! Hunt down:
        jhi=jlo
2       jlo=jhi-inc
        if(jlo.lt.1)then! Done hunting, since o end of table.
            jlo=0
        else if(x.lt.xx(jlo).eqv.ascnd)then! Not done hunting,
            jhi=jlo
        inc=inc+inc! so double the increment
        goto 2! and try again.
        endif! Done hunting, value bracketed.
      endif
3     if(jhi-jlo.eq.1)then
        if(x.eq.xx(n))jlo=n-1
        if(x.eq.xx(1))jlo=1
            return
      endif
      jm=(jhi+jlo)/2
      if(x.ge.xx(jm).eqv.ascnd)then
        jlo=jm
      else
        jhi=jm
      endif
      goto 3
      END SUBROUTINE hunt