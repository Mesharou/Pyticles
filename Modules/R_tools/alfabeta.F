

      subroutine alfabeta (Lm,Mm,t,s,rho0,alpha,beta)
!
! Compute thermal expansion and saline contraction coefficients
! as functions of potential temperature, salinity from a polynomial
! expression (Jackett & McDougall, 1992). The coefficients are
! evaluated at the surface.
!
!  alpha(Ts,Tt,0)=-d(rho1(Ts,Tt,0))/d(Tt) / rho0
!  beta(Ts,Tt,0) = d(rho1(Ts,Tt,0))/d(Ts) / rho0
!
!  Adapted from original "rati" and "beta" routines.
!
      implicit none

      integer Lm,Mm, i,j

      integer imin,imax,jmin,jmax

      real      t(0:Lm+1,0:Mm+1), s(0:Lm+1,0:Mm+1),
     &          alpha(0:Lm+1,0:Mm+1), beta(0:Lm+1,0:Mm+1),
     &          rho0
# include "scalars.h"

      real Q01, Q02, Q03, Q04, Q05, U00, U01, U02, U03, U04,
     &                                   V00, V01, V02, W00
      parameter(                  Q01=6.793952E-2, Q02=-9.095290E-3,
     &          Q03=+1.001685E-4, Q04=-1.120083E-6, Q05=+6.536332E-9,
     &          U00=+0.824493   , U01=-4.08990E-3 , U02=+7.64380E-5 ,
     &          U03=-8.24670E-7 , U04=+5.38750E-9 , V00=-5.72466E-3 ,
     &          V01=+1.02270E-4 , V02=-1.65460E-6 , W00=+4.8314E-4  )
      real      sqrtTs, cff, Tt, Ts


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Cf2py intent(in) Lm,Mm,t,s, rho0
Cf2py intent(out) alpha,beta

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        imin=0
        imax=Lm+1
        jmin=0
        jmax=Mm+1

      cff=1./rho0
      do j=jmin,jmax
        do i=imin,imax
          Tt=t(i,j)
#  ifdef SALINITY
          Ts=s(i,j)
          sqrtTs=sqrt(Ts)
#  else
          Ts=0.
          sqrtTs=0.
#  endif


! !             rho1=(dr00+Tt*( Q01+Tt*( Q02+Tt*( Q03+Tt*(
! !      &                                           Q04+Tt*Q05 ))))
! !      &                         +Ts*( U00+Tt*( U01+Tt*( U02+Tt*(
! !      &                                            U03+Tt*U04 )))
! !      &                              +sqrtTs*(V00+Tt*(
! !      &                                   V01+Tt*V02 ))+Ts*W00 ))

          alpha(i,j)=-cff*( Q01+Tt*( 2.*Q02+Tt*( 3.*Q03+Tt*(
     &                                  4.*Q04 +Tt*5.*Q05 )))
     &                         +Ts*( U01+Tt*( 2.*U02+Tt*(
     &                                     3.*U03 +Tt*4.*U04 ))
     &                                +sqrtTs*( V01+Tt*2.*V02))
     &                                                       )
 
          beta(i,j)= cff*( U00+Tt*(U01+Tt*(U02+Tt*(U03+Tt*U04)))
     &            +1.5*(V00+Tt*(V01+Tt*V02))*sqrtTs+2.*W00*Ts )
        enddo
      enddo

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      return
      end