!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Some Roms routines... 
!!
!!(similar to alex routines with addition of the "Cf2py ..." lines
!!
!!compile with "f2py -DF2PY_REPORT_ON_ARRAY_COPY=1 -c -m romstoolsfort romstoolsfort.F" for python use
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Z interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sigma_to_z_introld (Lm,Mm,N, nz, z_r, z_w, rmask, var,
     &                    z_lev, var_zlv, imin,jmin,kmin, FillValue)
!
! Interpolate field "var" defined in sigma-space to a set of several
! fixed z-levels z_lev(1:nz). 
!

#define CUBIC


      implicit none
      integer Lm,Mm,N, nz,  imin,imax,jmin,jmax, kmin, i,j,k,m

      real*8 z_r(0:Lm+1,0:Mm+1,N),   rmask(0:Lm+1,0:Mm+1),
     &       z_w(0:Lm+1,0:Mm+1,0:N), zlev
      real*4 z_lev(nz), var(imin:Lm+1,jmin:Mm+1,kmin:N),  FillValue,
     &              var_zlv(imin:Lm+1,jmin:Mm+1,nz)

      integer km(0:Lm+1)
      real*8 zz(0:Lm+1,0:N+1), dpth
#ifdef CUBIC
     &     , dz(0:Lm+1,kmin-1:N), FC(0:Lm+1,kmin-1:N), p,q,cff
#endif
      integer numthreads, trd, chunk_size, margin, jstr,jend
C$    integer omp_get_num_threads, omp_get_thread_num

      imax=Lm+1
      jmax=Mm+1

      numthreads=1
C$    numthreads=omp_get_num_threads()
      trd=0
C$    trd=omp_get_thread_num()
      chunk_size=(jmax-jmin + numthreads)/numthreads
      margin=(chunk_size*numthreads -jmax+jmin-1)/2
      jstr=jmin !max(  trd   *chunk_size  -margin, jmin )
      jend=jmax !min( (trd+1)*chunk_size-1-margin, jmax )


Cf2py intent(in) Lm,Mm,N, nz, z_r, z_w, rmask, var, z_lev, imin,jmin,kmin, FillValue
Cf2py intent(out) var_zlv


#ifdef VERBOSE
      write(*,'(2x,2A,3I3,2x,A,3I4)') 'sigma_to_z_intr: imin,jmin,',
     &   'kmin =', imin,jmin,kmin, 'trd,jstr,jend =', trd,jstr,jend
      write(*,*) Lm,Mm,N
#endif


      do j=jstr,jend
        if (kmin.eq.1) then
          if (imin.eq.0 .and. jmin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=z_r(i,j,k)
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=z_w(i,j,0)
              zz(i,N+1)=z_w(i,j,N)
            enddo
          elseif (imin.eq.1 .and. jmin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.5D0*(z_r(i,j,k)+z_r(i-1,j,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.5D0*(z_w(i-1,j,0)+z_w(i,j,0))
              zz(i,N+1)=0.5D0*(z_w(i-1,j,N)+z_w(i,j,N))
            enddo
          elseif (imin.eq.0 .and. jmin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.5*(z_r(i,j,k)+z_r(i,j-1,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.5D0*(z_w(i,j,0)+z_w(i,j-1,0))
              zz(i,N+1)=0.5D0*(z_w(i,j,N)+z_w(i,j-1,N))
            enddo
          elseif (imin.eq.1 .and. jmin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.25D0*( z_r(i,j,k)+z_r(i-1,j,k)
     &                      +z_r(i,j-1,k)+z_r(i-1,j-1,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.25D0*( z_w(i,j,0)+z_w(i-1,j,0)
     &                     +z_w(i,j-1,0)+z_w(i-1,j-1,0))

              zz(i,N+1)=0.25D0*( z_w(i,j,N)+z_w(i-1,j,N)
     &                      +z_w(i,j-1,N)+z_w(i-1,j-1,N))
             enddo
          endif
        else
          if (imin.eq.0 .and. jmin.eq.0) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=z_w(i,j,k)
              enddo
            enddo
          elseif (imin.eq.1 .and. jmin.eq.0) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.5D0*(z_w(i,j,k)+z_w(i-1,j,k))
              enddo
            enddo
          elseif (imin.eq.0 .and. jmin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.5*(z_w(i,j,k)+z_w(i,j-1,k))
              enddo
            enddo
          elseif (imin.eq.1 .and. jmin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.25D0*( z_w(i,j,k)+z_w(i-1,j,k)
     &                      +z_w(i,j-1,k)+z_w(i-1,j-1,k))
              enddo
            enddo
          endif
        endif


#ifdef CUBIC
        do k=kmin,N-1
          do i=imin,imax
            dz(i,k)=zz(i,k+1)-zz(i,k)
            FC(i,k)=var(i,j,k+1)-var(i,j,k)
          enddo
        enddo
        do i=imin,imax
          dz(i,kmin-1)=dz(i,kmin)
          FC(i,kmin-1)=FC(i,kmin)

          dz(i,N)=dz(i,N-1)
          FC(i,N)=FC(i,N-1)
        enddo
        do k=N,kmin,-1   !--> irreversible
          do i=imin,imax
            cff=FC(i,k)*FC(i,k-1)
            if (cff.gt.0.D0) then
              FC(i,k)=cff*(dz(i,k)+dz(i,k-1))/( (FC(i,k)+FC(i,k-1))
     &                                          *dz(i,k)*dz(i,k-1) )
            else
              FC(i,k)=0.D0
            endif
          enddo
        enddo
#endif
        do m=1,nz
          zlev=z_lev(m)
c**       write(*,*) 'm =',m, ' zlev =', zlev
          if (kmin.eq.0) then                      !
            do i=imin,imax                         !
              dpth=zz(i,N)-zz(i,0) 
              if (rmask(i,j).lt.0.5) then
                km(i)=-3     !--> masked out
              elseif (dpth*(zlev-zz(i,N)).gt.0.) then
                km(i)=N+2    !<-- above surface
              elseif (dpth*(zz(i,0)-zlev).gt.0.) then
                km(i)=-2     !<-- below bottom
              else
                km(i)=-1     !--> to search
              endif
            enddo
          else
            do i=imin,imax
              dpth=zz(i,N+1)-zz(i,0)     
              if (rmask(i,j).lt.0.5) then
                km(i)=-3     !--> masked out
              elseif (dpth*(zlev-zz(i,N+1)).gt.0.) then
                km(i)=N+2    !<-- above surface

              elseif (dpth*(zlev-zz(i,N)).gt.0.) then
                km(i)=N      !<-- below surface, but above z_r(N)
              elseif (dpth*(zz(i,0)-zlev).gt.0.) then
                km(i)=-2     !<-- below bottom
              elseif (dpth*(zz(i,1)-zlev).gt.0.) then
                km(i)=0      !<-- above bottom, but below z_r(1)
              else
                km(i)=-1     !--> to search
              endif
            enddo
          endif

          do k=N-1,kmin,-1
            do i=imin,imax
              if (km(i).eq.-1) then
                if((zz(i,k+1)-zlev)*(zlev-zz(i,k)) .ge. 0.) km(i)=k
              endif
            enddo
          enddo

          do i=imin,imax
            if (km(i).eq.-3) then
              var_zlv(i,j,m)=0.             !<-- masked out
            elseif (km(i).eq.-2) then
              var_zlv(i,j,m)=FillValue      !<-- below bottom
            elseif (km(i).eq.N+2) then
              var_zlv(i,j,m)=-FillValue     !<-- above surface
            elseif (km(i).eq.N) then
              var_zlv(i,j,m)=var(i,j,N)     !-> R-point, above z_r(N)
#ifdef CUBIC
     &              +FC(i,N)*(zlev-zz(i,N))
#else
     &              +(zlev-zz(i,N))*(var(i,j,N)-var(i,j,N-1))
     &                                   /(zz(i,N)-zz(i,N-1))
#endif
            elseif (km(i).eq.kmin-1) then   !-> R-point below z_r(1),
              var_zlv(i,j,m)=var(i,j,kmin)  !     but above bottom  
#ifdef CUBIC
     &           -FC(i,kmin)*(zz(i,kmin)-zlev)
#else
     &          -(zz(i,kmin)-zlev)*(var(i,j,kmin+1)-var(i,j,kmin))
     &                                  /(zz(i,kmin+1)-zz(i,kmin))
#endif
            else
              k=km(i)
#ifdef CUBIC
              cff=1.D0/(zz(i,k+1)-zz(i,k))
              p=zlev-zz(i,k)
              q=zz(i,k+1)-zlev

              var_zlv(i,j,m)=cff*( q*var(i,j,k) + p*var(i,j,k+1)
     &             -cff*p*q*( cff*(q-p)*(var(i,j,k+1)-var(i,j,k))
     &                                  +p*FC(i,k+1) -q*FC(i,k) )
     &                                                          )
#else
              var_zlv(i,j,m)=( var(i,j,k)*(zz(i,k+1)-zlev)
     &                        +var(i,j,k+1)*(zlev-zz(i,k))
     &                               )/(zz(i,k+1)-zz(i,k))
#endif
#ifdef VERBOSE
              if ((zz(i,k+1)-zlev)*(zlev-zz(i,k)).lt.0.D0) then
                write(*,'(/1x,2A/)') '### ERROR: sigma_to_z_intr ',
     &                               ':: Search algoritm failure.'
                stop 
              endif
#endif
            endif
          enddo
        enddo ! <-- m
      enddo !<-- j

      return
      end







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute density anomaly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine rho_eos(Lm,Mm,N, T,S, z_r,z_w,rho0, rho1,qp1,rho,bvf)


#define SPLIT_EOS
#define NONLIN_EOS
#define DUKO_2001


!
! Compute density anomaly from T,S via Equation Of State (EOS) for
!-------- ------- ------- ----  for seawater. Following Jackett and 
! McDougall, 1995, physical EOS is assumed to have form
!
!                           rho0 + rho1(T,S)
!           rho(T,S,z) = ------------------------                 (1)
!                         1 - 0.1*|z|/K(T,S,|z|)
!
! where rho1(T,S) is sea-water density perturbation[kg/m^3] at
! standard pressure of 1 Atm (sea surface); |z| is absolute depth,
! i.e. distance from free-surface to the point at which density is
! computed, and
!
!     K(T,S,|z|) = K00 + K01(T,S) + K1(T,S)*|z| + K2(T,S)*|z|^2.  (2)
!
! To reduce errors of pressure-gradient scheme associated with
! nonlinearity of compressibility effects, as well as to reduce
! roundoff errors, the dominant part of density profile,
!
!                           rho0
!                     ----------------                            (3)
!                      1 - 0.1|z|/K00
!
! is removed from from (1). [Since (3) is purely a function of z,
! it does not contribute to pressure gradient.]  This results in
!
!                   rho1 - rho0*[K01+K1*|z|+K2*|z|^2]/[K00-0.1|z|]
!    rho1 + 0.1|z| -----------------------------------------------
!                        K00 + K01 + (K1-0.1)*|z| + K2*|z|^2
!                                                                 (4)
! which is suitable for pressure-gradient calculation.
!
! Optionally, if CPP-switch SPLIT_EOS is defined, term proportional
! to |z| is linearized using smallness 0.1|z|/[K00 + K01] << 1 and
! the resultant EOS has form
!
!              rho(T,S,z) = rho1(T,S) + qp1(T,S)*|z|              (5)
!
! where
!                            rho1(T,S) - rho0*K01(T,S)/K00
!             qp1(T,S)= 0.1 -------------------------------       (6)
!                                   K00 + K01(T,S)
!
! is stored in a special array.
!
! This splitting allows representation of spatial derivatives (and
! also differences) of density as sum of adiabatic derivatives and
! compressible part according to
!
!         d rho     d rho1           d qp1             d |z|
!        ------- = -------- + |z| * -------  +  qp1 * -------     (7)
!         d x,s      d x,s           d x,s             d x,s
!
!                  |<----- adiabatic ----->|   |<- compress ->|
!
! so that constraining of adiabatic derivative for monotonicity is
! equivalent to enforcement of physically stable stratification.
! [This separation and constraining algorithm is subsequently used
! in computation of pressure gradient within prsgrd32ACx-family
! schemes.]
!
! If so prescribed compute the Brunt-Vaisala frequency [1/s^2] at
! horizontal RHO-points and vertical W-points,
!
!                          g    d rho  |
!             bvf^2 = - ------ ------- |                          (8)
!                        rho0    d z   | adiabatic
!
! where density anomaly difference is computed by adiabatically
! rising/lowering the water parcel from RHO point above/below to
! the W-point depth at "z_w".
!
! WARNING: Shared target arrays in the code below: "rho1",
!          "bvf" (if needed), and
!
!          SPLIT_EOS is defined: "qp1"  ["rho" does not exist]
!                   not defined  "rho"  ["qp1" does not exist]
!
!
! Reference:  Jackett, D. R. and T. J. McDougall, 1995, Minimal
!             Adjustment of Hydrostatic Profiles to Achieve Static
!             Stability. J. Atmos. Ocean. Tec., vol. 12, pp. 381-389.
!
! << This equation of state formulation has been derived by Jackett
! and McDougall (1992), unpublished manuscript, CSIRO, Australia. It
! computes in-situ density anomaly as a function of potential
! temperature (Celsius) relative to the surface, salinity (PSU),
! and depth (meters).  It assumes  no  pressure  variation along
! geopotential  surfaces,  that  is,  depth  and  pressure  are
! interchangeable. >>
!                                          John Wilkin, 29 July 92
!


      implicit none
      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k
      real*4 T(0:Lm+1,0:Mm+1,N),  S(0:Lm+1,0:Mm+1,N)
      real*8 z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N),
     &      rho(0:Lm+1,0:Mm+1,N), bvf(0:Lm+1,0:Mm+1,0:N),
     &      rho1(0:Lm+1,0:Mm+1,N), qp1(0:Lm+1,0:Mm+1,N)


      real*8 dpth,cff,cff2, Tt,Ts,sqrtTs, rho0, K0, dr00

      real*8, parameter :: r00=999.842594,   r01=6.793952E-2,
     &          r02=-9.095290E-3, r03=1.001685E-4,  r04=-1.120083E-6,
     &                                              r05=6.536332E-9,
     &          r10=0.824493,     r11=-4.08990E-3,  r12=7.64380E-5,
     &                            r13=-8.24670E-7,  r14=5.38750E-9,
     &          rS0=-5.72466E-3,  rS1=1.02270E-4,   rS2=-1.65460E-6,
     &          r20=4.8314E-4,
     &          K00=19092.56,     K01=209.8925,     K02=-3.041638,
     &                            K03=-1.852732e-3, K04=-1.361629e-5,
     &          K10=104.4077,     K11=-6.500517,    K12=0.1553190,
     &                                              K13=2.326469e-4,
     &          KS0=-5.587545,    KS1=+0.7390729,   KS2=-1.909078e-2

      real*8, parameter :: qp2=0.0000172
      real*8, parameter :: g=9.81

#  ifdef DUKO_2001
      real rho1_0, K0_Duk
#  endif

      integer numthreads, trd, chunk_size, margin, jstr,jend


Cf2py intent(in) Lm,Mm,N, T,S, z_r,z_w, rho0
Cf2py intent(out) rho1,qp1,rho,bvf



!!!!!!!!!!!!!!!!!!



#  ifdef DUKO_2001
      Tt=3.8D0
      Ts=34.5D0
      sqrtTs=sqrt(Ts)

      dr00=r00-1000.D0

      rho1_0=dr00 +Tt*( r01+Tt*( r02+Tt*( r03+Tt*( r04+Tt*r05 ))))
     &                            +Ts*( R10+Tt*( r11+Tt*( r12+Tt*(
     &                                              r13+Tt*r14 )))
     &                   +sqrtTs*( rS0+Tt*( rS1+Tt*rS2 ))+Ts*r20 )

      K0_Duk= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     &       +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     &            +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))
#  endif

      dr00=r00-rho0
        write(*,*) dr00, rho0

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1


      numthreads=1
C$    numthreads=omp_get_num_threads()
      trd=0
C$    trd=omp_get_thread_num()
      chunk_size=(jmax-jmin + numthreads)/numthreads
      margin=(chunk_size*numthreads -jmax+jmin-1)/2
      jstr=max(  trd   *chunk_size  -margin, jmin )
      jend=min( (trd+1)*chunk_size-1-margin, jmax )




      do j=jstr,jend
        do k=1,N
          do i=imin,imax
            Tt=T(i,j,k)                   

            Ts=max(S(i,j,k), 0.)
            sqrtTs=sqrt(Ts)

            rho1(i,j,k)=( dr00 +Tt*( r01+Tt*( r02+Tt*( r03+Tt*(
     &                                           r04+Tt*r05 ))))
     &                         +Ts*( r10+Tt*( r11+Tt*( r12+Tt*(
     &                                            r13+Tt*r14 )))
     &                              +sqrtTs*(rS0+Tt*(
     &                                   rS1+Tt*rS2 ))+Ts*r20 ))

            K0= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     &         +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     &              +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))
 


            qp1(i,j,k)= 0.1D0*(rho0+rho1(i,j,k))*(K0_Duk-K0)
     &                               /((K00+K0)*(K00+K0_Duk))


         rho(i,j,k) = rho1(i,j,k) +  qp1(i,j,k)*(z_w(i,j,N)-z_r(i,j,k))


          enddo
        enddo
 
        cff=g/rho0
        do k=1,N-1
          do i=imin,imax

            dpth=z_w(i,j,N)-0.5*(z_r(i,j,k+1)+z_r(i,j,k))
            cff2=( rho1(i,j,k+1)-rho1(i,j,k)             ! Elementary
     &                        +(qp1(i,j,k+1)-qp1(i,j,k)) ! adiabatic 
     &                            *dpth*(1.-2.*qp2*dpth) ! difference
     &                       )


            bvf(i,j,k)=-cff*cff2 / (z_r(i,j,k+1)-z_r(i,j,k))

          enddo
        enddo

      enddo    ! <-- j

      return
      end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute density anomaly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine rho1_eos(Lm,Mm,N, T,S, z_r,z_w,rho0, rho1)


#define SPLIT_EOS
#define NONLIN_EOS
#define DUKO_2001


!
! Compute density anomaly from T,S via Equation Of State (EOS) for
!-------- ------- ------- ----  for seawater. Following Jackett and 
! McDougall, 1995, physical EOS is assumed to have form
!
!                           rho0 + rho1(T,S)
!           rho(T,S,z) = ------------------------                 (1)
!                         1 - 0.1*|z|/K(T,S,|z|)
!
! where rho1(T,S) is sea-water density perturbation[kg/m^3] at
! standard pressure of 1 Atm (sea surface); |z| is absolute depth,
! i.e. distance from free-surface to the point at which density is
! computed, and
!
!     K(T,S,|z|) = K00 + K01(T,S) + K1(T,S)*|z| + K2(T,S)*|z|^2.  (2)
!
! To reduce errors of pressure-gradient scheme associated with
! nonlinearity of compressibility effects, as well as to reduce
! roundoff errors, the dominant part of density profile,
!
!                           rho0
!                     ----------------                            (3)
!                      1 - 0.1|z|/K00
!
! is removed from from (1). [Since (3) is purely a function of z,
! it does not contribute to pressure gradient.]  This results in
!
!                   rho1 - rho0*[K01+K1*|z|+K2*|z|^2]/[K00-0.1|z|]
!    rho1 + 0.1|z| -----------------------------------------------
!                        K00 + K01 + (K1-0.1)*|z| + K2*|z|^2
!                                                                 (4)
! which is suitable for pressure-gradient calculation.
!
! Optionally, if CPP-switch SPLIT_EOS is defined, term proportional
! to |z| is linearized using smallness 0.1|z|/[K00 + K01] << 1 and
! the resultant EOS has form
!
!              rho(T,S,z) = rho1(T,S) + qp1(T,S)*|z|              (5)
!
! where
!                            rho1(T,S) - rho0*K01(T,S)/K00
!             qp1(T,S)= 0.1 -------------------------------       (6)
!                                   K00 + K01(T,S)
!
! is stored in a special array.
!
! This splitting allows representation of spatial derivatives (and
! also differences) of density as sum of adiabatic derivatives and
! compressible part according to
!
!         d rho     d rho1           d qp1             d |z|
!        ------- = -------- + |z| * -------  +  qp1 * -------     (7)
!         d x,s      d x,s           d x,s             d x,s
!
!                  |<----- adiabatic ----->|   |<- compress ->|
!
! so that constraining of adiabatic derivative for monotonicity is
! equivalent to enforcement of physically stable stratification.
! [This separation and constraining algorithm is subsequently used
! in computation of pressure gradient within prsgrd32ACx-family
! schemes.]
!
! If so prescribed compute the Brunt-Vaisala frequency [1/s^2] at
! horizontal RHO-points and vertical W-points,
!
!                          g    d rho  |
!             bvf^2 = - ------ ------- |                          (8)
!                        rho0    d z   | adiabatic
!
! where density anomaly difference is computed by adiabatically
! rising/lowering the water parcel from RHO point above/below to
! the W-point depth at "z_w".
!
! WARNING: Shared target arrays in the code below: "rho1",
!          "bvf" (if needed), and
!
!          SPLIT_EOS is defined: "qp1"  ["rho" does not exist]
!                   not defined  "rho"  ["qp1" does not exist]
!
!
! Reference:  Jackett, D. R. and T. J. McDougall, 1995, Minimal
!             Adjustment of Hydrostatic Profiles to Achieve Static
!             Stability. J. Atmos. Ocean. Tec., vol. 12, pp. 381-389.
!
! << This equation of state formulation has been derived by Jackett
! and McDougall (1992), unpublished manuscript, CSIRO, Australia. It
! computes in-situ density anomaly as a function of potential
! temperature (Celsius) relative to the surface, salinity (PSU),
! and depth (meters).  It assumes  no  pressure  variation along
! geopotential  surfaces,  that  is,  depth  and  pressure  are
! interchangeable. >>
!                                          John Wilkin, 29 July 92
!


      implicit none
      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k
      real*4 T(0:Lm+1,0:Mm+1,N),  S(0:Lm+1,0:Mm+1,N)
      real*8 z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N),
     &      rho1(0:Lm+1,0:Mm+1,N)


      real*8 dpth,cff,cff2, Tt,Ts,sqrtTs, rho0, K0, dr00

      real*8, parameter :: r00=999.842594,   r01=6.793952E-2,
     &          r02=-9.095290E-3, r03=1.001685E-4,  r04=-1.120083E-6,
     &                                              r05=6.536332E-9,
     &          r10=0.824493,     r11=-4.08990E-3,  r12=7.64380E-5,
     &                            r13=-8.24670E-7,  r14=5.38750E-9,
     &          rS0=-5.72466E-3,  rS1=1.02270E-4,   rS2=-1.65460E-6,
     &          r20=4.8314E-4,
     &          K00=19092.56,     K01=209.8925,     K02=-3.041638,
     &                            K03=-1.852732e-3, K04=-1.361629e-5,
     &          K10=104.4077,     K11=-6.500517,    K12=0.1553190,
     &                                              K13=2.326469e-4,
     &          KS0=-5.587545,    KS1=+0.7390729,   KS2=-1.909078e-2

      real*8, parameter :: qp2=0.0000172
      real*8, parameter :: g=9.81

#  ifdef DUKO_2001
      real rho1_0, K0_Duk
#  endif

      integer numthreads, trd, chunk_size, margin, jstr,jend


Cf2py intent(in) Lm,Mm,N, T,S, z_r,z_w, rho0
Cf2py intent(out) rho1



!!!!!!!!!!!!!!!!!!



#  ifdef DUKO_2001
      Tt=3.8D0
      Ts=34.5D0
      sqrtTs=sqrt(Ts)

      dr00=r00-1000.D0

      rho1_0=dr00 +Tt*( r01+Tt*( r02+Tt*( r03+Tt*( r04+Tt*r05 ))))
     &                            +Ts*( R10+Tt*( r11+Tt*( r12+Tt*(
     &                                              r13+Tt*r14 )))
     &                   +sqrtTs*( rS0+Tt*( rS1+Tt*rS2 ))+Ts*r20 )

      K0_Duk= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     &       +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     &            +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))
#  endif

      dr00=r00-rho0
        write(*,*) dr00, rho0

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1


      numthreads=1
C$    numthreads=omp_get_num_threads()
      trd=0
C$    trd=omp_get_thread_num()
      chunk_size=(jmax-jmin + numthreads)/numthreads
      margin=(chunk_size*numthreads -jmax+jmin-1)/2
      jstr=max(  trd   *chunk_size  -margin, jmin )
      jend=min( (trd+1)*chunk_size-1-margin, jmax )




      do j=jstr,jend
        do k=1,N
          do i=imin,imax
            Tt=T(i,j,k)                   

            Ts=max(S(i,j,k), 0.)
            sqrtTs=sqrt(Ts)

            rho1(i,j,k)=( dr00 +Tt*( r01+Tt*( r02+Tt*( r03+Tt*(
     &                                           r04+Tt*r05 ))))
     &                         +Ts*( r10+Tt*( r11+Tt*( r12+Tt*(
     &                                            r13+Tt*r14 )))
     &                              +sqrtTs*(rS0+Tt*(
     &                                   rS1+Tt*rS2 ))+Ts*r20 ))


          enddo
        enddo
 


      enddo    ! <-- j

      return
      end




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute neutral density gradients
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine rho_grad(Lm,Mm,N, T,S, z_r,z_w,rho0,pm,pn,
     & rho1,qp1,drdz,drdx,drdy)


#define SPLIT_EOS
#define NONLIN_EOS
#define DUKO_2001



      implicit none
      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k
      real*4 T(0:Lm+1,0:Mm+1,N),  S(0:Lm+1,0:Mm+1,N),
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1)
      real*8 z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N),
     &      rho1(0:Lm+1,0:Mm+1,N), drdz(0:Lm+1,0:Mm+1,0:N),
     &      drdx(1:Lm+1,0:Mm+1,N), drdy(0:Lm+1,1:Mm+1,N),
     &      qp1(0:Lm+1,0:Mm+1,N)

      real*8 dpth,cff,cff2, Tt,Ts,sqrtTs, rho0, K0, dr00

      real*8, parameter :: r00=999.842594,   r01=6.793952E-2,
     &          r02=-9.095290E-3, r03=1.001685E-4,  r04=-1.120083E-6,
     &                                              r05=6.536332E-9,
     &          r10=0.824493,     r11=-4.08990E-3,  r12=7.64380E-5,
     &                            r13=-8.24670E-7,  r14=5.38750E-9,
     &          rS0=-5.72466E-3,  rS1=1.02270E-4,   rS2=-1.65460E-6,
     &          r20=4.8314E-4,
     &          K00=19092.56,     K01=209.8925,     K02=-3.041638,
     &                            K03=-1.852732e-3, K04=-1.361629e-5,
     &          K10=104.4077,     K11=-6.500517,    K12=0.1553190,
     &                                              K13=2.326469e-4,
     &          KS0=-5.587545,    KS1=+0.7390729,   KS2=-1.909078e-2

      real*8, parameter :: qp2=0.0000172
      real*8, parameter :: g=9.81

#  ifdef DUKO_2001
      real rho1_0, K0_Duk
#  endif

      integer numthreads, trd, chunk_size, margin, jstr,jend


Cf2py intent(in) Lm,Mm,N, T,S, z_r,z_w, rho0,pm,pn
Cf2py intent(out) rho1,qp1,drdz,drdx,drdy



!!!!!!!!!!!!!!!!!!



#  ifdef DUKO_2001
      Tt=3.8D0
      Ts=34.5D0
      sqrtTs=sqrt(Ts)

      dr00=r00-1000.D0

      rho1_0=dr00 +Tt*( r01+Tt*( r02+Tt*( r03+Tt*( r04+Tt*r05 ))))
     &                            +Ts*( R10+Tt*( r11+Tt*( r12+Tt*(
     &                                              r13+Tt*r14 )))
     &                   +sqrtTs*( rS0+Tt*( rS1+Tt*rS2 ))+Ts*r20 )

      K0_Duk= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     &       +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     &            +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))
#  endif

      dr00=r00-rho0


      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1

      do j=jmin,jmax
        do k=1,N
          do i=imin,imax
            Tt=T(i,j,k)                   

            Ts=max(S(i,j,k), 0.)
            sqrtTs=sqrt(Ts)

            rho1(i,j,k)=( dr00 +Tt*( r01+Tt*( r02+Tt*( r03+Tt*(
     &                                           r04+Tt*r05 ))))
     &                         +Ts*( r10+Tt*( r11+Tt*( r12+Tt*(
     &                                            r13+Tt*r14 )))
     &                              +sqrtTs*(rS0+Tt*(
     &                                   rS1+Tt*rS2 ))+Ts*r20 ))

            K0= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     &         +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     &              +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))
 


            qp1(i,j,k)= 0.1D0*(rho0+rho1(i,j,k))*(K0_Duk-K0)
     &                               /((K00+K0)*(K00+K0_Duk))



          enddo
        enddo
 
        cff=g/rho0

        do k=1,N-1
          do i=imin,imax

            dpth=z_w(i,j,N)-0.5*(z_r(i,j,k+1)+z_r(i,j,k))
            cff2=( rho1(i,j,k+1)-rho1(i,j,k)             ! Elementary
     &                        +(qp1(i,j,k+1)-qp1(i,j,k)) ! adiabatic 
     &                            *dpth*(1.-2.*qp2*dpth) ! difference
     &                       )


            drdz(i,j,k)=-cff*cff2 / (z_r(i,j,k+1)-z_r(i,j,k))

          enddo
        enddo


      enddo    ! <-- j


!---------------------------------------------------------------------------------------
      do k=N,1,-1

      write(*,*) k
        do j=jmin,jmax
          do i=imin+1,imax


            dpth=0.5*( z_w(i,j,N)+z_w(i-1,j,N)
     &                -z_r(i,j,k)-z_r(i-1,j,k))

            drdx(i,j,k)=-cff*( rho1(i,j,k)-rho1(i-1,j,k) ! Elementary
     &                +(qp1(i,j,k)-qp1(i-1,j,k))         ! adiabatic
     &                     *dpth*(1.-qp2*dpth) )         ! difference
     &                    *0.5*(pm(i,j)+pm(i-1,j))


          enddo
        enddo

      write(*,*) 'drdy',k
!---------------------------------------------------------------------------------------
        do j=jmin+1,jmax
          do i=imin,imax


            dpth=0.5*( z_w(i,j,N)+z_w(i,j-1,N)
     &                -z_r(i,j,k)-z_r(i,j-1,k))
            
            drdy(i,j,k)=-cff*( rho1(i,j,k)-rho1(i,j-1,k)   ! Elementary
     &                +(qp1(i,j,k)-qp1(i,j-1,k))         ! adiabatic
     &                     *dpth*(1.-qp2*dpth) )  
     &                    *0.5*(pn(i,j)+pn(i,j-1)) 




          enddo
        enddo

      enddo
!---------------------------------------------------------------------------------------

      return
      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute density anomaly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      subroutine rho_eos_alex (Lm,Mm,N, T,S, z_r,z_w, rp_r,rp_w,rho,bvf)


#define SPLIT_EOS
#define NONLIN_EOS
#define DUKO_2001



      implicit none
      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k
      real*4 T(0:Lm+1,0:Mm+1,N),  S(0:Lm+1,0:Mm+1,N)
      real*8 z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N),
     &      rp_r(0:Lm+1,0:Mm+1,N), rp_w(0:Lm+1,0:Mm+1,0:N),
     &      rho(0:Lm+1,0:Mm+1,N), bvf(0:Lm+1,0:Mm+1,0:N)

      real*8 rho1(0:Lm+1,N), qp1(0:Lm+1,N)

      real*8, parameter ::   r00=999.842594D0,     r01=6.793952D-2,
     &      r02=-9.095290D-3,   r03=1.001685D-4,    r04=-1.120083D-6,
     &                                              r05=6.536332D-9,
     &      r10=0.824493D0,     r11=-4.08990D-3,    r12=7.64380D-5,
     &                          r13=-8.24670D-7,    r14=5.38750D-9,
     &      rS0=-5.72466D-3,    rS1=1.02270D-4,     rS2=-1.65460D-6,
     &      r20=4.8314D-4,

     &      K00=19092.56D0,     K01=209.8925D0,     K02=-3.041638D0,
     &                          K03=-1.852732D-3,   K04=-1.361629D-5,
     &      K10=104.4077D0,     K11=-6.500517D0,    K12=0.1553190D0,
     &                                              K13=2.326469D-4,
     &      KS0=-5.587545D0,    KS1=+0.7390729D0,   KS2=-1.909078D-2


      real*8, parameter :: qp2=0.0000172
      real*8 dpth,cff, Tt,Ts,sqrtTs, rho0, K0, dr00

      integer numthreads, trd, chunk_size, margin, jstr,jend
C$    integer omp_get_num_threads, omp_get_thread_num


Cf2py intent(in) Lm,Mm,N, T,S, z_r,z_w
Cf2py intent(out) rp_r,rp_w,rho,bvf

#  ifdef DUKO_2001
      real K0_Duk, rho1_0

      Tt=3.8D0
      Ts=34.5D0
      sqrtTs=sqrt(Ts)

      cff=r00-1000.D0
      rho1_0=cff +Tt*( r01+Tt*( r02+Tt*( r03+Tt*( r04+Tt*r05 ))))
     &                            +Ts*( R10+Tt*( r11+Tt*( r12+Tt*(
     &                                              r13+Tt*r14 )))
     &                   +sqrtTs*( rS0+Tt*( rS1+Tt*rS2 ))+Ts*r20 )

      K0_Duk= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     &       +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     &            +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))

c**   write(*,'(2(4x,A,F16.8))') 'rho1_0=',rho1_0, 'K0_Duk=',K0_Duk
      rho0=rho1_0 
#  endif

      dr00=r00-rho0

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1


      numthreads=1
C$    numthreads=omp_get_num_threads()
      trd=0
C$    trd=omp_get_thread_num()
      chunk_size=(jmax-jmin + numthreads)/numthreads
      margin=(chunk_size*numthreads -jmax+jmin-1)/2
      jstr=max(  trd   *chunk_size  -margin, jmin )
      jend=min( (trd+1)*chunk_size-1-margin, jmax )

c**    write(*,*) 'trd,jstr,jend =', trd,jstr,jend



      do j=jstr,jend
        do k=1,N
          do i=imin,imax
            Tt=T(i,j,k)
            Ts=max(S(i,j,k),0.)
            sqrtTs=sqrt(Ts)

            rho1(i,k)=( dr00 +Tt*( r01+Tt*( r02+Tt*( r03+Tt*(
     &                                         r04+Tt*r05 ))))

     &                       +Ts*( r10+Tt*( r11+Tt*( r12+Tt*(
     &                                          r13+Tt*r14 )))
     &                            +sqrtTs*(rS0+Tt*(
     &                                 rS1+Tt*rS2 ))+Ts*r20 ))
#  ifdef MASKING
     &                                            *rmask(i,j)
#  endif


            K0= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     &         +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     &              +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))

            qp1(i,k)= 0.1D0*(rho0+rho1(i,k))*(K0_Duk-K0)
     &                          /((K00+K0)*(K00+K0_Duk))
#   ifdef MASKING
     &                   *rmask(i,j)
#   endif

            rho(i,j,k) = rho1(i,k) +  qp1(i,k)*(z_w(i,j,N)-z_r(i,j,k))

          enddo
        enddo

        cff=1.D0
        do k=1,N-1
          do i=imin,imax
            dpth=z_w(i,j,N)-0.5*(z_r(i,j,k+1)+z_r(i,j,k))

            bvf(i,j,k)=-cff*( rho1(i,k+1)-rho1(i,k)   ! Elementary
     &                      +(qp1(i,k+1)-qp1(i,k))  ! adiabatic
     &                      *dpth*(1.-2.*qp2*dpth)  ! difference
     &                 )/(z_r(i,j,k+1)-z_r(i,j,k))
          enddo
        enddo
        do i=imin,imax
          bvf(i,j,N)=bvf(i,j,N-1)
          bvf(i,j,0)=bvf(i,j,1)
          rp_w(i,j,N)=rho1(i,N) + (rho1(i,N)-rho1(i,N-1))
     &                           *(z_w(i,j,N)-z_r(i,j,N))
     &                         /(z_r(i,j,N)-z_r(i,j,N-1))

          rp_r(i,j,N)=rp_w(i,j,N) + bvf(i,j,N)*(z_w(i,j,N)-z_r(i,j,N))
        enddo

        do k=N-1,1,-1
          do i=imin,imax
            rp_r(i,j,k)=rp_r(i,j,k+1) +bvf(i,j,k)*( z_r(i,j,k+1)
     &                                             -z_r(i,j,k))

            rp_w(i,j,k)=( rp_r(i,j,k)*(z_r(i,j,k+1)-z_w(i,j,k))
     &                   +rp_r(i,j,k+1)*(z_w(i,j,k)-z_r(i,j,k))
     &                              )/(z_r(i,j,k+1)-z_r(i,j,k))
          enddo
        enddo

        do i=imin,imax
          rp_w(i,j,0)=rp_r(i,j,1) +bvf(i,j,1)*(z_r(i,j,1)-z_w(i,j,0))
        enddo
      enddo ! <-- j

      return
      end








!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute density anomaly
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      subroutine test_args (Lm,Mm,N, T,S, rp_r)


      implicit none
      integer Lm,Mm,N
      real*4 T(0:Lm+1,0:Mm+1,N),  S(0:Lm+1,0:Mm+1,N)
      real*8 rp_r(0:Lm+1,0:Mm+1,N)


Cf2py intent(in) Lm,Mm,N,T
Cf2py intent(out) rp_r
Cf2py intent(inout) S

      write(*,*) T(1,1,1)

      end







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Z interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sigma_to_z_intr (Lm,Mm,N, nz, z_r, z_w, rmask, var,
     &                    z_lev, var_zlv, imin,jmin,kmin, FillValue)
!
! Interpolate field "var" defined in sigma-space to 3-D z_lev. 
!


#define INTERP_BELOW
#define CUBIC



      implicit none
      integer Lm,Mm,N, nz,  imin,imax,jmin,jmax, kmin, i,j,k,m

      real*8 z_r(0:Lm+1,0:Mm+1,N),   rmask(0:Lm+1,0:Mm+1),
     &       z_w(0:Lm+1,0:Mm+1,0:N), z_lev(0:Lm+1,0:Mm+1,nz)
      real*4 var(imin:Lm+1,jmin:Mm+1,kmin:N),  FillValue,
     &              var_zlv(imin:Lm+1,jmin:Mm+1,nz)

      integer km(0:Lm+1)
      real*8 zz(0:Lm+1,0:N+1), dpth
#ifdef CUBIC
     &     , dz(0:Lm+1,kmin-1:N), FC(0:Lm+1,kmin-1:N), p,q,cff
#endif
      integer numthreads, trd, chunk_size, margin, jstr,jend
C$    integer omp_get_num_threads, omp_get_thread_num

      imax=Lm+1
      jmax=Mm+1

      numthreads=1
C$    numthreads=omp_get_num_threads()
      trd=0
C$    trd=omp_get_thread_num()
      chunk_size=(jmax-jmin + numthreads)/numthreads
      margin=(chunk_size*numthreads -jmax+jmin-1)/2
      jstr=jmin !max(  trd   *chunk_size  -margin, jmin )
      jend=jmax !min( (trd+1)*chunk_size-1-margin, jmax )


Cf2py intent(in) Lm,Mm,N, nz, z_r, z_w, rmask, var, z_lev, imin,jmin,kmin, FillValue
Cf2py intent(out) var_zlv


#ifdef VERBOSE
      write(*,'(2x,2A,3I3,2x,A,3I4)') 'sigma_to_z_intr: imin,jmin,',
     &   'kmin =', imin,jmin,kmin, 'trd,jstr,jend =', trd,jstr,jend
      write(*,*) Lm,Mm,N
#endif


      do j=jstr,jend
        if (kmin.eq.1) then
          if (imin.eq.0 .and. jmin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=z_r(i,j,k)
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=z_w(i,j,0)
              zz(i,N+1)=z_w(i,j,N)
            enddo
          elseif (imin.eq.1 .and. jmin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.5D0*(z_r(i,j,k)+z_r(i-1,j,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.5D0*(z_w(i-1,j,0)+z_w(i,j,0))
              zz(i,N+1)=0.5D0*(z_w(i-1,j,N)+z_w(i,j,N))
            enddo
          elseif (imin.eq.0 .and. jmin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.5*(z_r(i,j,k)+z_r(i,j-1,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.5D0*(z_w(i,j,0)+z_w(i,j-1,0))
              zz(i,N+1)=0.5D0*(z_w(i,j,N)+z_w(i,j-1,N))
            enddo
          elseif (imin.eq.1 .and. jmin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.25D0*( z_r(i,j,k)+z_r(i-1,j,k)
     &                      +z_r(i,j-1,k)+z_r(i-1,j-1,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.25D0*( z_w(i,j,0)+z_w(i-1,j,0)
     &                     +z_w(i,j-1,0)+z_w(i-1,j-1,0))

              zz(i,N+1)=0.25D0*( z_w(i,j,N)+z_w(i-1,j,N)
     &                      +z_w(i,j-1,N)+z_w(i-1,j-1,N))
             enddo
          endif
        else
          if (imin.eq.0 .and. jmin.eq.0) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=z_w(i,j,k)
              enddo
            enddo
          elseif (imin.eq.1 .and. jmin.eq.0) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.5D0*(z_w(i,j,k)+z_w(i-1,j,k))
              enddo
            enddo
          elseif (imin.eq.0 .and. jmin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.5*(z_w(i,j,k)+z_w(i,j-1,k))
              enddo
            enddo
          elseif (imin.eq.1 .and. jmin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.25D0*( z_w(i,j,k)+z_w(i-1,j,k)
     &                      +z_w(i,j-1,k)+z_w(i-1,j-1,k))
              enddo
            enddo
          endif
        endif
#ifdef CUBIC
        do k=kmin,N-1
          do i=imin,imax
            dz(i,k)=zz(i,k+1)-zz(i,k)
            FC(i,k)=var(i,j,k+1)-var(i,j,k)
          enddo
        enddo
        do i=imin,imax
          dz(i,kmin-1)=dz(i,kmin)
          FC(i,kmin-1)=FC(i,kmin)

          dz(i,N)=dz(i,N-1)
          FC(i,N)=FC(i,N-1)
        enddo
        do k=N,kmin,-1   !--> irreversible
          do i=imin,imax
            cff=FC(i,k)*FC(i,k-1)
            if (cff.gt.0.D0) then
              FC(i,k)=cff*(dz(i,k)+dz(i,k-1))/( (FC(i,k)+FC(i,k-1))
     &                                          *dz(i,k)*dz(i,k-1) )
            else
              FC(i,k)=0.D0
            endif
          enddo
        enddo
#endif
        do m=1,nz
c**       write(*,*) 'm =',m, ' zlev =', zlev
          if (kmin.eq.0) then                      !
            do i=imin,imax                         !
              dpth=zz(i,N)-zz(i,0) 
              if (rmask(i,j).lt.0.5) then
                km(i)=-3     !--> masked out
              elseif (dpth*(z_lev(i,j,m)-zz(i,N)).gt.0.) then
                km(i)=N+2    !<-- above surface
              elseif (dpth*(zz(i,0)-z_lev(i,j,m)).gt.0.) then
                km(i)=-2     !<-- below bottom
              else
                km(i)=-1     !--> to search
              endif
            enddo
          else
            do i=imin,imax
              dpth=zz(i,N+1)-zz(i,0)     
              if (rmask(i,j).lt.0.5) then
                km(i)=-3     !--> masked out
              elseif (dpth*(z_lev(i,j,m)-zz(i,N+1)).gt.0.) then
                km(i)=N+2    !<-- above surface

              elseif (dpth*(z_lev(i,j,m)-zz(i,N)).gt.0.) then
                km(i)=N      !<-- below surface, but above z_r(N)
              elseif (dpth*(zz(i,0)-z_lev(i,j,m)).gt.0.) then
                km(i)=-2     !<-- below bottom
              elseif (dpth*(zz(i,1)-z_lev(i,j,m)).gt.0.) then
                km(i)=0      !<-- above bottom, but below z_r(1)
              else
                km(i)=-1     !--> to search
              endif
            enddo
          endif
          do k=N-1,kmin,-1
            do i=imin,imax
              if (km(i).eq.-1) then
                if((zz(i,k+1)-z_lev(i,j,m))*(z_lev(i,j,m)-zz(i,k))
     & .ge. 0.) km(i)=k
              endif
            enddo
          enddo

          do i=imin,imax
            if (km(i).eq.-3) then
              var_zlv(i,j,m)=0.             !<-- masked out
            elseif (km(i).eq.-2) then
#ifdef INTERP_BELOW
              var_zlv(i,j,m)=var(i,j,kmin)  !     but above bottom  
#ifdef CUBIC
     &           -FC(i,kmin)*(zz(i,kmin)-z_lev(i,j,m))
#else
     &          -(zz(i,kmin)-z_lev(i,j,m))*(var(i,j,kmin+1)-var(i,j,kmin))
     &                                  /(zz(i,kmin+1)-zz(i,kmin))
#endif
#else
              var_zlv(i,j,m)=FillValue      !<-- below bottom
#endif
            elseif (km(i).eq.N+2) then
              var_zlv(i,j,m)=-FillValue     !<-- above surface
            elseif (km(i).eq.N) then
              var_zlv(i,j,m)=var(i,j,N)     !-> R-point, above z_r(N)
#ifdef CUBIC
     &              +FC(i,N)*(z_lev(i,j,m)-zz(i,N))
#else
     &              +(z_lev(i,j,m)-zz(i,N))*(var(i,j,N)-var(i,j,N-1))
     &                                   /(zz(i,N)-zz(i,N-1))
#endif
            elseif (km(i).eq.kmin-1) then   !-> R-point below z_r(1),
              var_zlv(i,j,m)=var(i,j,kmin)  !     but above bottom  
#ifdef CUBIC
     &           -FC(i,kmin)*(zz(i,kmin)-z_lev(i,j,m))
#else
     &          -(zz(i,kmin)-z_lev(i,j,m))*(var(i,j,kmin+1)-var(i,j,kmin))
     &                                  /(zz(i,kmin+1)-zz(i,kmin))
#endif
            else
              k=km(i)
#ifdef CUBIC
              cff=1.D0/(zz(i,k+1)-zz(i,k))
              p=z_lev(i,j,m)-zz(i,k)
              q=zz(i,k+1)-z_lev(i,j,m)

              var_zlv(i,j,m)=cff*( q*var(i,j,k) + p*var(i,j,k+1)
     &             -cff*p*q*( cff*(q-p)*(var(i,j,k+1)-var(i,j,k))
     &                                  +p*FC(i,k+1) -q*FC(i,k) )
     &                                                          )
#else
              var_zlv(i,j,m)=( var(i,j,k)*(zz(i,k+1)-z_lev(i,j,m))
     &                        +var(i,j,k+1)*(z_lev(i,j,m)-zz(i,k))
     &                               )/(zz(i,k+1)-zz(i,k))
#endif
#ifdef VERBOSE
              if ((zz(i,k+1)-z_lev(i,j,m))*(z_lev(i,j,m)-zz(i,k))
     & .lt.0.D0) then
                write(*,'(/1x,2A/)') '### ERROR: sigma_to_z_intr ',
     &                               ':: Search algoritm failure.'
                stop 
              endif
#endif
            endif
          enddo
        enddo ! <-- m
      enddo !<-- j

      return
      end







