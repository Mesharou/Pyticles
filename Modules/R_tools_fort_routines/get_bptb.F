!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute Bottom Pressure torque J(Pb,H)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






      subroutine get_bptb (Lm,Mm,N,T,S, z_r,z_w,rho0,pm,pn
     &  ,bpt)    



      implicit none

      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k,
     &      istr,iend,jstr,jend,istrU,jstrV



      real*8 T(0:Lm+1,0:Mm+1,N),  S(0:Lm+1,0:Mm+1,N),
     &      rho1(0:Lm+1,0:Mm+1,N), qp1(0:Lm+1,0:Mm+1,N),
     &      z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N),
     &      Hz(0:Lm+1,0:Mm+1,0:N),
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      Tt,Ts,sqrtTs, rho0,  K0, dr00,
     &      cff ,cff1, cff2, cfr, HalfGRho, GRho,
     &      var1, var2,var3, var4

      real*8 bpt(1:Lm+1,1:Mm+1)


      real*8 P(0:Lm+1,0:Mm+1,N),
     &      ru(1:Lm+1,0:Mm+1), rv(0:Lm+1,1:Mm+1),
     &      rho(0:Lm+1,0:Mm+1,N), dpth,
     &      dR(0:Lm+1,N), dZ(0:Lm+1,N),
     &      FC(0:Lm+2,0:Mm+2), dZx(0:Lm+1,0:Mm+1),
     &      rx(0:Lm+2,0:Mm+2), dRx(0:Lm+1,0:Mm+1)

      real*8, parameter :: OneFifth=0.2, OneTwelfth=1./12., epsil=0.

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
     &          KS0=-5.587545,    KS1=+0.7390729,   KS2=-1.909078e-2,
     &          qp2=0.0000172

#  ifdef DUKO_2001
      real rho1_0, K0_Duk
#  endif
# include "scalars.h"



Cf2py intent(in) Lm,Mm,N, T,S,z_r,z_w,rho0,pm,pn
Cf2py intent(out) bpt




!
! A non-conservative Density-Jacobian scheme using cubic polynomial
! fits for rho and z_r as functions of nondimensianal coordinates xi,
! eta, and s (basically their respective fortran indices). The cubic
! polynomials are constructed by specifying first derivatives of
! interpolated fields on co-located (non-staggered) grid. These
! derivatives are computed using harmonic (rather that algebraic)
! averaging of elementary differences, which guarantees monotonicity
! of the resultant interpolant.
!
! In the code below, if CPP-switch SPLIT_EOS is defined, the Equation
! of State (EOS) is assumed to have form
!
!       rho(T,S,z) = rho1(T,S) + qp1(T,S)*dpth*[1.-qp2*dpth]
!
! where rho1 is potential density at 1 atm and qp1 is compressibility
! coefficient, which does not depend on z, and dpth=zeta-z, and qp2
! is just a constant. In this case
!
!   d rho    d rho1   d qp1                                    d z
!  ------- = ------ + ----- *dpth*[..] - qp1*[1.-2.*qp2*dpth]*------
!   d s,x     d s,x   d s,x                                    d s,x
!
!           |<--- adiabatic part --->|  |<--- compressible part --->|
!
! where the first two terms constitute "adiabatic derivative" of
! density, which is subject to harmonic averaging, while the last
! term is added in later. This approach quarantees that density
! profile reconstructed by cubic polynomial maintains its positive
! statification in physical sense as long as discrete values of
! density are positively stratified.
!
! This scheme retains exact antisymmetry J(rho,z_r)=-J(z_r,rho)
! [with the exception of harmonic averaging algorithm in the case
! when CPP-switch SPLIT_EOS is defined, see above]. If parameter
! OneFifth (see above) is set to zero, the scheme becomes identical
! to standard Jacobian.
!
! NOTE: This routine is an alternative form of prsgrd32 and it
!       produces results identical to that if its prototype.
!


        istr=0
        istrU=1
        iend=Lm+1
        jstr=0
        jstrV=1
        jend=Mm+1

        imin=istrU
        imax=iend
        jmin=jstrV
        jmax=jend


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! COMPUTE DENSITY 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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


!---------------------------------------------------------------------------------------
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


!---------------------------------------------------------------------------------------

            Hz(i,j,k)=z_w(i,j,k)-z_w(i,j,k-1)




          enddo
        enddo
      enddo    ! <-- j



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!
! Preliminary step (same for XI- and ETA-components:
!------------ ---- ----- --- --- --- ---------------
!
      GRho=g/rho0
      HalfGRho=0.5*GRho
 
      do j=jstrV-1,jend
        do k=1,N-1
          do i=istrU-1,iend
            dZ(i,k)=z_r(i,j,k+1)-z_r(i,j,k)
# ifdef SPLIT_EOS
            dpth=z_w(i,j,N)-0.5*(z_r(i,j,k+1)+z_r(i,j,k))

            dR(i,k)=rho1(i,j,k+1)-rho1(i,j,k)            ! Elementary
     &              +(qp1(i,j,k+1)-qp1(i,j,k))           ! adiabatic
     &                     *dpth*(1.-qp2*dpth)           ! difference
# else
            dR(i,k)=rho(i,j,k+1)-rho(i,j,k)
# endif
          enddo
        enddo
        do i=istrU-1,iend
          dR(i,N)=dR(i,N-1)
          dR(i,0)=dR(i,1)
          dZ(i,N)=dZ(i,N-1)
          dZ(i,0)=dZ(i,1)
        enddo
        do k=N,1,-1               !--> irreversible
          do i=istrU-1,iend
            cff=2.*dZ(i,k)*dZ(i,k-1)
            dZ(i,k)=cff/(dZ(i,k)+dZ(i,k-1))
 
            cfr=2.*dR(i,k)*dR(i,k-1)
            if (cfr.gt.epsil) then
              dR(i,k)=cfr/(dR(i,k)+dR(i,k-1))
            else
              dR(i,k)=0.
            endif
# ifdef SPLIT_EOS
            dpth=z_w(i,j,N)-z_r(i,j,k)
            dR(i,k)=dR(i,k)  -qp1(i,j,k)*dZ(i,k)*(1.-2.*qp2*dpth)
            rho(i,j,k)=rho1(i,j,k) +qp1(i,j,k)*dpth*(1.-qp2*dpth)
# endif
          enddo
        enddo
        do i=istrU-1,iend
          P(i,j,N)=g*z_w(i,j,N) + GRho*( rho(i,j,N)
     &       +0.5*(rho(i,j,N)-rho(i,j,N-1))*(z_w(i,j,N)-z_r(i,j,N))
     &          /(z_r(i,j,N)-z_r(i,j,N-1)) )*(z_w(i,j,N)-z_r(i,j,N))
        enddo
        do k=N-1,1,-1
          do i=istrU-1,iend
            P(i,j,k)=P(i,j,k+1)+HalfGRho*( (rho(i,j,k+1)+rho(i,j,k))
     &                                     *(z_r(i,j,k+1)-z_r(i,j,k))
 
     &     -OneFifth*( (dR(i,k+1)-dR(i,k))*( z_r(i,j,k+1)-z_r(i,j,k)
     &                              -OneTwelfth*(dZ(i,k+1)+dZ(i,k)) )
 
     &                -(dZ(i,k+1)-dZ(i,k))*( rho(i,j,k+1)-rho(i,j,k)
     &                              -OneTwelfth*(dR(i,k+1)+dR(i,k)) )
     &                                                             ))
          enddo
        enddo
      enddo   !<-- j




        ru = 0.
        rv = 0.



!
! Compute XI-component of pressure gradient term:
!-------- ------------ -- -------- -------- -----
!
      do k=N,1,-1
        do j=jstr,jend
          do i=imin,imax
            FC(i,j)=(z_r(i,j,k)-z_r(i-1,j,k))
# ifdef SPLIT_EOS
            dpth=0.5*( z_w(i,j,N)+z_w(i-1,j,N)
     &                -z_r(i,j,k)-z_r(i-1,j,k))

            rx(i,j)=( rho1(i,j,k)-rho1(i-1,j,k)          ! Elementary
     &                +(qp1(i,j,k)-qp1(i-1,j,k))         ! adiabatic
     &                     *dpth*(1.-qp2*dpth) )         ! difference
# else
            rx(i,j)=(rho(i,j,k)-rho(i-1,j,k))
# endif

          enddo
        enddo
 
# ifndef EW_PERIODIC
        if (WESTERN_EDGE) then         ! Extrapolate elementary
          do j=jstr,jend               ! differences near physical
            FC(imin-1,j)=FC(imin,j)    ! boundaries to compencate.
            rx(imin-1,j)=rx(imin,j)    ! for reduced loop ranges.
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=jstr,jend
            FC(imax+1,j)=FC(imax,j)
            rx(imax+1,j)=rx(imax,j)
          enddo
        endif
# endif
 
        do j=jstr,jend
          do i=istrU-1,iend
            cff=2.*FC(i,j)*FC(i+1,j)
            if (cff.gt.epsil) then
              dZx(i,j)=cff/(FC(i,j)+FC(i+1,j))
            else
              dZx(i,j)=0.
            endif
 
            cfr=2.*rx(i,j)*rx(i+1,j)
            if (cfr.gt.epsil) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i+1,j))
            else
              dRx(i,j)=0.
            endif
# ifdef SPLIT_EOS
            dRx(i,j)=dRx(i,j) -qp1(i,j,k)*dZx(i,j)
     &         *(1.-2.*qp2*(z_w(i,j,N)-z_r(i,j,k)))
# endif
          enddo               !--> discard FC, rx

          do i=istrU,iend
            ru(i,j)=ru(i,j) + 
     &    0.5*(Hz(i,j,k)+Hz(i-1,j,k))*2./(pn(i,j)+pn(i-1,j))*(
     &                              P(i-1,j,k)-P(i,j,k)-HalfGRho*(

     &            (rho(i,j,k)+rho(i-1,j,k))*(z_r(i,j,k)-z_r(i-1,j,k))
 
     &   -OneFifth*( (dRx(i,j)-dRx(i-1,j))*( z_r(i,j,k)-z_r(i-1,j,k)
     &                            -OneTwelfth*(dZx(i,j)+dZx(i-1,j)) )
 
     &              -(dZx(i,j)-dZx(i-1,j))*( rho(i,j,k)-rho(i-1,j,k)
     &                            -OneTwelfth*(dRx(i,j)+dRx(i-1,j)) )
     &                                                            )))
          enddo
        enddo
!
! ETA-component of pressure gradient term:
!-------------- -- -------- -------- -----
!
        do j=jmin,jmax
          do i=istr,iend
            FC(i,j)=(z_r(i,j,k)-z_r(i,j-1,k))      
# ifdef SPLIT_EOS
            dpth=0.5*( z_w(i,j,N)+z_w(i,j-1,N)
     &                -z_r(i,j,k)-z_r(i,j-1,k))
            
            rx(i,j)=( rho1(i,j,k)-rho1(i,j-1,k)          ! Elementary
     &                +(qp1(i,j,k)-qp1(i,j-1,k))         ! adiabatic
     &                     *dpth*(1.-qp2*dpth) )         ! difference
# else
            rx(i,j)=(rho(i,j,k)-rho(i,j-1,k))
# endif
          enddo
        enddo
 
# ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=istr,iend
            FC(i,jmin-1)=FC(i,jmin)
            rx(i,jmin-1)=rx(i,jmin)

          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=istr,iend
            FC(i,jmax+1)=FC(i,jmax)
            rx(i,jmax+1)=rx(i,jmax)
          enddo
        endif
# endif
 
        do j=jstrV-1,jend
          do i=istr,iend
            cff=2.*FC(i,j)*FC(i,j+1)
            if (cff.gt.epsil) then
c**         if ((FC(i,j).gt.0. .and. FC(i,j+1).gt.0.) .or.
c**  &          (FC(i,j).lt.0. .and. FC(i,j+1).lt.0.)) then
              dZx(i,j)=cff/(FC(i,j)+FC(i,j+1))
            else
              dZx(i,j)=0.
            endif
 
            cfr=2.*rx(i,j)*rx(i,j+1)
            if (cfr.gt.epsil) then
c**         if ((rx(i,j).gt.0. .and. rx(i,j+1).gt.0.) .or.
c**  &          (rx(i,j).lt.0. .and. rx(i,j+1).lt.0.)) then
              dRx(i,j)=cfr/(rx(i,j)+rx(i,j+1))
            else
              dRx(i,j)=0.
            endif
# ifdef SPLIT_EOS
            dRx(i,j)=dRx(i,j) -qp1(i,j,k)*dZx(i,j)
     &         *(1.-2.*qp2*(z_w(i,j,N)-z_r(i,j,k)))
# endif
          enddo               !--> discard FC, rx
 
          if (j.ge.jstrV) then
            do i=istr,iend

            rv(i,j)=rv(i,j) + 
     &            0.5*(Hz(i,j,k)+Hz(i,j-1,k))*2./(pm(i,j)+pm(i,j-1))*(
     &                             P(i,j-1,k)-P(i,j,k) -HalfGRho*(

     &            (rho(i,j,k)+rho(i,j-1,k))*(z_r(i,j,k)-z_r(i,j-1,k))
 
     &   -OneFifth*( (dRx(i,j)-dRx(i,j-1))*( z_r(i,j,k)-z_r(i,j-1,k)
     &                            -OneTwelfth*(dZx(i,j)+dZx(i,j-1)) )
 
     &              -(dZx(i,j)-dZx(i,j-1))*( rho(i,j,k)-rho(i,j-1,k)
     &                            -OneTwelfth*(dRx(i,j)+dRx(i,j-1)) )
     &                                                            )))

            enddo
          endif
        enddo
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute rotational
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! 
! 
!        do j=jmin+1,jmax-1
!           do i=imin+1,imax-1
! 
! 
!             ru(i,j) = ru(i,j)*0.5*(pn(i,j)+pn(i-1,j))
!      &                       *0.5*(pm(i,j)+pm(i-1,j))
! 
! 
!             rv(i,j) = rv(i,j)*0.5*(pm(i,j)+pm(i,j-1))
!      &                       *0.5*(pn(i,j)+pn(i,j-1))
! 
! 
! 
!            cff1 = 0.25*(pm(i,j) + pm(i-1,j) + pm(i,j-1) + pm(i-1,j-1))
!            cff2 = 0.25*(pn(i,j) + pn(i-1,j) + pn(i,j-1) + pn(i-1,j-1))
! 
!             bpt(i,j) = (rv(i,j) - rv(i-1,j)) * cff1
!      &                 - (ru(i,j) - ru(i,j-1)) * cff2
!         enddo
!       enddo
! 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NEW COMPUTATION OF ROT (added 14/09/13)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=jmin+1,jmax-1
        do i=imin+1,imax-1

            cff=0.5*(pn(i,j)+pn(i-1,j))
            ru(i,j)=ru(i,j)*cff

            cff=0.5*(pm(i,j)+pm(i,j-1))
            rv(i,j)=rv(i,j)*cff


           cff1 = 0.25*(pm(i,j) + pm(i-1,j) + pm(i,j-1) + pm(i-1,j-1))
     &            * 0.25*(pn(i,j) + pn(i-1,j) + pn(i,j-1) + pn(i-1,j-1))


             bpt(i,j) = (rv(i,j) - rv(i-1,j)) * cff1
     &                  - (ru(i,j) - ru(i,j-1)) * cff1

        enddo
       enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






      return
      end

 
