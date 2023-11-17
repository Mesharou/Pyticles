





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
      do j=jstrR,jendR
        do k=1,N                                       ! NONLINEAR
          do i=istrR,iendR                             !  EQUATION
            Tt=T(i,j,k)                     !  OF STATE
#  ifdef SALINITY
            Ts=max(S(i,j,k), 0.)
            sqrtTs=sqrt(Ts)
#  else
            Ts=0.0
            sqrtTs=0.0
#  endif
            rho1(i,j,k)=( dr00 +Tt*( r01+Tt*( r02+Tt*( r03+Tt*(
     &                                           r04+Tt*r05 ))))
     &                         +Ts*( r10+Tt*( r11+Tt*( r12+Tt*(
     &                                            r13+Tt*r14 )))
     &                              +sqrtTs*(rS0+Tt*(
     &                                   rS1+Tt*rS2 ))+Ts*r20 ))
#  ifdef MASKING
     &                                              *rmask(i,j)
#  endif
     
          

     
            K0= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     &         +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     &              +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))
 
#  ifdef SPLIT_EOS
#   ifdef DUKO_2001
            qp1(i,j,k)= 0.1D0*(rho0+rho1(i,j,k))*(K0_Duk-K0)
     &                               /((K00+K0)*(K00+K0_Duk))
#   else
            qp1(i,j,k)=( 0.000649D0 +0.1D0*( K00*rho1(i,j,k)
     &                             -rho0*K0)/(K00*(K00+K0)) )
#   endif
#   ifdef MASKING
     &                                            *rmask(i,j)
#   endif
#  else
            K1=B00+Tt*(B01+Tt*(B02+Tt*B03)) +Ts*( B10+Tt*( B11
     &                                    +Tt*B12 )+sqrtTs*BS1 )

            K2=E00+Tt*(E01+Tt*E02) +Ts*(E10+Tt*(E11+Tt*E12))


c**         rho(i,j,k)=rho1(i,j,k) + 0.1*(z_w(i,j,N)-z_r(i,j,k))
c**  &                  *(K00*rho1(i,j,k)-rho0*K0)/(K00*(K00+K0))
 
            dpth=z_w(i,j,N)-z_r(i,j,k)
            cff=K00-0.1*dpth
            cff1=K0+dpth*(K1+K2*dpth)
            rho(i,j,k)=( rho1(i,j,k)*cff*(K00+cff1)
     &                           -0.1*dpth*rho0*cff1
     &                            )/(cff*(cff+cff1))
 
c**  &                    +0.1*dpth*rho0/cff
 
 
c**            rho(i,j,k)=(rho0+rho1(i,j,k))
c**  &            /( 1.-0.1*dpth/(K00+K0+dpth*(K1+K2*dpth)))
c**  &                                                 -rho0
 
c**         rho(i,j,k)=rho1(i,j,k)
 
 
#   ifdef MASKING
            rho(i,j,k)=rho(i,j,k)*rmask(i,j)
#   endif
#  endif

            

            
          enddo

        enddo
      enddo
   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

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

!

# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then     ! Restrict extended ranges one
        imin=istrU               ! point inward near the physical
      else                       ! boundary. Note that this version
        imin=istrU-1             ! of code is suitable for MPI
      endif                      ! configuration too.
      if (EASTERN_EDGE) then
        imax=iend
      else
        imax=iend+1
      endif
# else
      imin=istr-1
      imax=iend+1
# endif
# ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          jmin=jstrV
        else
          jmin=jstrV-1
        endif
        if (NORTHERN_EDGE) then
          jmax=jend
        else
          jmax=jend+1
        endif
# else
        jmin=jstr-1
        jmax=jend+1
# endif


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

!
! Compute XI-component of pressure gradient term:
!-------- ------------ -- -------- -------- -----
!
      do k=N,1,-1
        do j=jstr,jend
          do i=imin,imax
              

              
              
            FCP(i,j)=(z_r(i,j,k)-z_r(i-1,j,k))
# ifdef MASKING
     &                              *umask(i,j)
# endif   
# ifdef SPLIT_EOS
            dpth=0.5*( z_w(i,j,N)+z_w(i-1,j,N)
     &                -z_r(i,j,k)-z_r(i-1,j,k))

            rx(i,j)=( rho1(i,j,k)-rho1(i-1,j,k)          ! Elementary
     &                +(qp1(i,j,k)-qp1(i-1,j,k))         ! adiabatic
     &                     *dpth*(1.-qp2*dpth) )         ! difference
# else
            rx(i,j)=(rho(i,j,k)-rho(i-1,j,k))
# endif
# ifdef MASKING
     &                              *umask(i,j)
# endif
     

         
          enddo
        enddo
 
# ifndef EW_PERIODIC
        if (WESTERN_EDGE) then         ! Extrapolate elementary
          do j=jstr,jend               ! differences near physical
            FCP(imin-1,j)=FCP(imin,j)    ! boundaries to compencate.
            rx(imin-1,j)=rx(imin,j)    ! for reduced loop ranges.
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=jstr,jend
            FCP(imax+1,j)=FCP(imax,j)
            rx(imax+1,j)=rx(imax,j)
          enddo
        endif
# endif
 
        do j=jstr,jend
          do i=istrU-1,iend
            cff=2.*FCP(i,j)*FCP(i+1,j)
            if (cff.gt.epsil) then
              dZx(i,j)=cff/(FCP(i,j)+FCP(i+1,j))
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
          enddo               !--> discard FCP, rx

          do i=istrU,iend
            ru(i,j,k)=0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)*(
     &                              P(i-1,j,k)-P(i,j,k)-HalfGRho*(

     &            (rho(i,j,k)+rho(i-1,j,k))*(z_r(i,j,k)-z_r(i-1,j,k))
 
     &   -OneFifth*( (dRx(i,j)-dRx(i-1,j))*( z_r(i,j,k)-z_r(i-1,j,k)
     &                            -OneTwelfth*(dZx(i,j)+dZx(i-1,j)) )
 
     &              -(dZx(i,j)-dZx(i-1,j))*( rho(i,j,k)-rho(i-1,j,k)
     &                            -OneTwelfth*(dRx(i,j)+dRx(i-1,j)) )
     &                                                            )))
                  

              MPrsgrd(i,j,k,1) = ru(i,j,k)
        

          enddo
        enddo





! ETA-component of pressure gradient term:
!-------------- -- -------- -------- -----
!
        do j=jmin,jmax
          do i=istr,iend
            FCP(i,j)=(z_r(i,j,k)-z_r(i,j-1,k))
# ifdef MASKING
     &                              *vmask(i,j)
# endif  
# ifdef SPLIT_EOS
            dpth=0.5*( z_w(i,j,N)+z_w(i,j-1,N)
     &                -z_r(i,j,k)-z_r(i,j-1,k))
            
            rx(i,j)=( rho1(i,j,k)-rho1(i,j-1,k)          ! Elementary
     &                +(qp1(i,j,k)-qp1(i,j-1,k))         ! adiabatic
     &                     *dpth*(1.-qp2*dpth) )         ! difference
# else
            rx(i,j)=(rho(i,j,k)-rho(i,j-1,k))
# endif
# ifdef MASKING
     &                              *vmask(i,j)
# endif  

          enddo
        enddo
 
# ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=istr,iend
            FCP(i,jmin-1)=FCP(i,jmin)
            rx(i,jmin-1)=rx(i,jmin)

          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=istr,iend
            FCP(i,jmax+1)=FCP(i,jmax)
            rx(i,jmax+1)=rx(i,jmax)

          enddo
        endif
# endif
 
        do j=jstrV-1,jend
          do i=istr,iend
            cff=2.*FCP(i,j)*FCP(i,j+1)
            if (cff.gt.epsil) then
c**         if ((FCP(i,j).gt.0. .and. FCP(i,j+1).gt.0.) .or.
c**  &          (FCP(i,j).lt.0. .and. FCP(i,j+1).lt.0.)) then
              dZx(i,j)=cff/(FCP(i,j)+FCP(i,j+1))
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
          enddo               !--> discard FCP, rx
 
          if (j.ge.jstrV) then
            do i=istr,iend
              rv(i,j,k)=0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)*(
     &                             P(i,j-1,k)-P(i,j,k) -HalfGRho*(

     &            (rho(i,j,k)+rho(i,j-1,k))*(z_r(i,j,k)-z_r(i,j-1,k))
 
     &   -OneFifth*( (dRx(i,j)-dRx(i,j-1))*( z_r(i,j,k)-z_r(i,j-1,k)
     &                            -OneTwelfth*(dZx(i,j)+dZx(i,j-1)) )
 
     &              -(dZx(i,j)-dZx(i,j-1))*( rho(i,j,k)-rho(i,j-1,k)
     &                            -OneTwelfth*(dRx(i,j)+dRx(i,j-1)) )
     &                                                            )))

              MPrsgrd(i,j,k,2) = rv(i,j,k)

            enddo
          endif
        enddo
      enddo


 
