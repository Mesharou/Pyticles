
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute bottom  vorticity stretching (r'$- f \vec{u_b} . \vec{\nabla} h $')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






      subroutine get_fwb(Lm,Mm,N,u,v, z_r,z_w,pm,pn,f
     &  ,fwb)    


      implicit none

      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k,
     &      istr,iend,jstr,jend,istrU,jstrV


      real*8 u(1:Lm+1,0:Mm+1,N),  v(0:Lm+1,1:Mm+1,N),
     &      FlxU(1:Lm+1,0:Mm+1,N),  FlxV(0:Lm+1,1:Mm+1,N),
     &      z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N),
     &      Hz(0:Lm+1,0:Mm+1,N),  f(0:Lm+1,0:Mm+1),
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      dn_u(0:Lm+1,0:Mm+1), dm_v(0:Lm+1,0:Mm+1),
     &      var1, var2,var3, var4

      real*8 Wxi(1:Lm+1,0:Mm+1),Weta(0:Lm+1,1:Mm+1)

      real*8 fwb(0:Lm+1,0:Mm+1)

# include "scalars.h"




Cf2py intent(in) Lm,Mm,N, u,v,z_r,z_w,pm,pn,f
Cf2py intent(out) fwb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!
! Compute "omega" vertical velocity by means of integration of mass
! divergence of mass fluxes from bottom up. In this computation,
! unlike that in omega.F, there is (1) immediate multiplication by
! pm*pn so that the result has meaning of velocity, rather than
! finite volume mass flux through vertical facet of tracer grid box;
! and (2, also unlike omega.F) no subtraction of vertical velocity
! of moving grid-box interface (the effect of "breething" of vertical
! grid system due to evolving free surface) is made now.
! Consequently, Wrk(:,N).ne.0, unlike its counterpart W(:,:,N).eqv.0
! in omega.F. Once omega vertical velocity is computed, interpolate
! it to vertical RHO-points.
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!
! Compute and add contributions due to (quasi-)horizontal motions
! along S=const surfaces by multiplying horizontal velocity
! components by slops S-coordinate surfaces:
!

        do j=jmin,jmax
          do i=imin+1,imax
            Wxi(i,j)=u(i,j,1)*(z_r(i,j,1)-z_r(i-1,j,1))
     &                                  *(pm(i,j)+pm(i-1,j))
          enddo
        enddo
        do j=jmin+1,jmax
          do i=imin,imax
            Weta(i,j)=v(i,j,1)*(z_r(i,j,1)-z_r(i,j-1,1))
     &                                   *(pn(i,j)+pn(i,j-1))
          enddo
        enddo
        do j=jmin+1,jmax-1
          do i=imin+1,imax-1
            fwb(i,j)=f(i,j)*0.25*(Wxi(i,j)+Wxi(i+1,j)
     &                                +Weta(i,j)+Weta(i,j+1))
          enddo
        enddo



!         do j=jmin,jmax
!           do i=imin+1,imax
!             Wxi(i,j)=u(i,j,N)*(z_r(i,j,N)-z_r(i-1,j,N))
!      &                                  *(pm(i,j)+pm(i-1,j))
!           enddo
!         enddo
!         do j=jmin+1,jmax
!           do i=imin,imax
!             Weta(i,j)=v(i,j,N)*(z_r(i,j,N)-z_r(i,j-1,N))
!      &                                   *(pn(i,j)+pn(i,j-1))
!           enddo
!         enddo
!         do j=jmin+1,jmax-1
!           do i=imin+1,imax-1
!             fws(i,j)=-f(i,j)*0.25*(Wxi(i,j)+Wxi(i+1,j)
!      &                                +Weta(i,j)+Weta(i,j+1))
!           enddo
!         enddo






      return
      end








