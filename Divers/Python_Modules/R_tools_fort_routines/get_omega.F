
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute vertical velocity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






      subroutine get_omega(Lm,Mm,N,u,v, z_r,z_w,pm,pn
     &  ,W)    


      implicit none

      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k,
     &      istr,iend,jstr,jend,istrU,jstrV


      real*8 u(1:Lm+1,0:Mm+1,N),  v(0:Lm+1,1:Mm+1,N),
     &      FlxU(1:Lm+1,0:Mm+1,N),  FlxV(0:Lm+1,1:Mm+1,N),
     &      z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N),
     &      Hz(0:Lm+1,0:Mm+1,N),
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      dn_u(0:Lm+1,0:Mm+1), dm_v(0:Lm+1,0:Mm+1),
     &      var1, var2,var3, var4


      real*8 wrk(0:Lm+1),W(0:Lm+1,0:Mm+1,0:N)

# include "scalars.h"




Cf2py intent(in) Lm,Mm,N, u,v,z_r,z_w,pm,pn
Cf2py intent(out) W

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1


      istr=1
      iend=Lm
      jstr=1
      jend=Mm

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

      do j=jstr,jend
        do i=istr,iend
          W(i,j,0)=0.
        enddo

        do k=1,N,+1        !--> recursive
          do i=istr,iend
            W(i,j,k)=W(i,j,k-1) -FlxU(i+1,j,k)+FlxU(i,j,k)
     &                          -FlxV(i,j+1,k)+FlxV(i,j,k)
          enddo
        enddo

        do i=istr,iend
          wrk(i)=W(i,j,N)/(z_w(i,j,N)-z_w(i,j,0))
          W(i,j,N)=0.
        enddo

        do k=N-1,1,-1
          do i=istr,iend
            W(i,j,k)=W(i,j,k)-wrk(i)*(z_w(i,j,k)-z_w(i,j,0))
          enddo
        enddo
      enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      do j=jstr,jend
        do i=istr,iend
           do k=1,N,+1 

             W(i,j,k)=W(i,j,k)*pm(i,j)*pn(i,j)

          enddo
        enddo
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                       ! Set lateral
        do k=0,N                                   ! boundary 
          do j=jstr,jend                           ! conditions
            W(istr-1,j,k)=W(istr,j,k)
            W(iend+1,j,k)=W(iend,j,k)
          enddo
        enddo
 
        do k=0,N
          do i=istr,iend
            W(i,jstr-1,k)=W(i,jstr,k)
            W(i,jend+1,k)=W(i,jend,k)
          enddo
        enddo

        do k=0,N
          W(istr-1,jstr-1,k)=W(istr,jstr,k)
        enddo

        do k=0,N
          W(istr-1,jend+1,k)=W(istr,jend,k)
        enddo

        do k=0,N
          W(iend+1, jstr-1,k)=W(iend,jstr,k)
        enddo

        do k=0,N
          W(iend+1,jend+1,k)=W(iend,jend,k)
        enddo

      return
      end








