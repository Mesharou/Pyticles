
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute difference between:
!
! -f wb (where wb is vertical velocity)
! and
! bottom  vorticity stretching (r'$- f \vec{u_b} . \vec{\nabla} h $')
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






      subroutine get_fwdivub(Lm,Mm,N,u,v, z_r,z_w,pm,pn,f
     &  ,fdivu)    


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

      real*8 Wrk(0:Lm+1,0:N),Wvlc(0:Lm+1,0:Mm+1,N),
     &      Wxi(1:Lm+1,0:Mm+1),Weta(0:Lm+1,1:Mm+1)

      real*8 fdivu(0:Lm+1,0:Mm+1)

# include "scalars.h"




Cf2py intent(in) Lm,Mm,N, u,v,z_r,z_w,pm,pn,f
Cf2py intent(out) fdivu

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






      do j=jmin,jmax
        do i=imin,imax
          do k=1,N,+1
           Hz(i,j,k) = z_w(i,j,k) - z_w(i,j,k-1)
          enddo
        enddo
      enddo




      do j=jmin,jmax
        do i=imin+1,imax
            dn_u(i,j) = 2./(pn(i,j)+pn(i-1,j))
            do k=1,N,+1
              FlxU(i,j,k) = 0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)
     &      * u(i,j,k) 
            enddo
          enddo
      enddo



      do j=jmin+1,jmax
        do i=imin,imax
            dm_v(i,j) = 2./(pm(i,j)+pm(i,j-1))
            do k=1,N,+1
              FlxV(i,j,k) = 0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)
     &      * v(i,j,k) 
            enddo
          enddo
      enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      do j=jmin,jmax-1

        do i=imin,imax
          Wrk(i,0)=0.
        enddo

        do k=1,N,+1
          do i=imin,imax-1
            Wrk(i,k)=Wrk(i,k-1)-pm(i,j)*pn(i,j)*(
     &                      FlxU(i+1,j,k)-FlxU(i,j,k)
     &                     +FlxV(i,j+1,k)-FlxV(i,j,k))
          enddo
        enddo


        do i=imin+1,imax
          fdivu(i,j)= f(i,j)*( -0.125*Wrk(i,2) +0.75*Wrk(i,1)
     &                                 +0.375*Wrk(i,0) )
        enddo
      enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



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
            fdivu(i,j)=fdivu(i,j)+f(i,j)*0.25*(Wxi(i,j)+Wxi(i+1,j)
     &                                +Weta(i,j)+Weta(i,j+1))
          enddo
        enddo









      return
      end








