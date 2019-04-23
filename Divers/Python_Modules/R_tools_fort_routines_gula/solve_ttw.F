
      subroutine solve_ttw(Lm,Mm,N,bx,by,Av,sustr,svstr,f,pm,pn,
     &  z_r,z_w,u,v,ug,vg)    


      implicit none

      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k,
     &      istr,iend,jstr,jend,istrU,jstrV


      real*8 u(1:Lm+1,0:Mm+1,N),  v(0:Lm+1,1:Mm+1,N),
     &      ug(1:Lm+1,0:Mm+1,N),  vg(0:Lm+1,1:Mm+1,N),
     &      z_r(0:Lm+1,0:Mm+1,N),z_w(0:Lm+1,0:Mm+1,0:N),
     &      Av(0:Lm+1,0:Mm+1,0:N),
     &      Hz(0:Lm+1,0:Mm+1,N),f(0:Lm+1,0:Mm+1),
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      sustr(0:Lm+1,0:Mm+1), svstr(0:Lm+1,0:Mm+1),
     &      bx(0:Lm+1,0:Mm+1,N), by(0:Lm+1,0:Mm+1,N)
!      &      AA(0:Lm+1,0:Mm+1,0:N),rho(0:Lm+1,0:Mm+1,N)
     
      real  b11,b12,b21,b22, det, TauX,TauY
      real, dimension(0:N) :: FC, c11,c12,c21,c22
      real, dimension(N) ::  d1,d2

      real cff, cff0,cff1,cff2,cff3, xx

# include "scalars.h"


Cf2py intent(in) Lm,Mm,N,bx,by,Av,sustr,svstr,f,pm,pn,z_r,z_w
Cf2py intent(out) u,v,ug,vg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1

      istr=1
      iend=Lm
      jstr=1
      jend=Mm

      istrU=istr+1
      jstrV=jstr+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





      do j=jmin,jmax
        do i=imin,imax
          do k=1,N,+1
           Hz(i,j,k) = z_w(i,j,k) - z_w(i,j,k-1)
          enddo
        enddo
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   
#   define I_EXT_RANGE istrU-1,iend+1
#   define J_EXT_RANGE jstrV-1,jend+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!       do k=1,N,+1             ! Put rho on Psi-grid
!           do j=jstr,jend+1
!             do i=istr,iend+1
!               AA(i,j,k)=0.25*( rho(i,j,k) +rho(i-1,j,k)
!      &                           +rho(i,j-1,k) +rho(i-1,j-1,k)
!      &                                                     )
!             enddo                   
!           enddo
!       enddo !<-- k


!       do k=1,N,+1
!         if (k == 1) then            ! The following segment converts
!           do j=J_EXT_RANGE          ! "AA" from density anomaly to
!             do i=I_EXT_RANGE        ! hydrostatic pressure anomaly.
!               AA(i,j,0)=0.          ! <-- initialize vertical
!             enddo                   !           integration
!           enddo
!         else
!           cff=0.125
! !          cff=0.5
!           do j=J_EXT_RANGE
!             do i=I_EXT_RANGE
!               xx=AA(i,j,0)  !<-- save (k-1)th value
! 
!               AA(i,j,0)=AA(i,j,0) +cff*(AA(i,j,k)+AA(i,j,k-1))*(
!      &            z_r(i,j,k)+z_r(i-1,j,k)+z_r(i,j-1,k)+z_r(i-1,j-1,k)
!      &           -z_r(i,j,k-1)-z_r(i-1,j,k-1)
!      &                               -z_r(i,j-1,k-1)-z_r(i-1,j-1,k-1)
!      &                                                              )
! 
! !              AA(i,j,0)=AA(i,j,0) +cff*(AA(i,j,k)+AA(i,j,k-1))*(
! !     &            z_r(i,j,k)
! !     &           -z_r(i,j,k-1)
! !     &                                                              )
! 
!               AA(i,j,k-1)=xx
! 
!             enddo
!           enddo
!         endif
!       enddo !<-- k
! 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
!       do j=J_EXT_RANGE           ! After the above steps "AA(:,:,k)"
!         do i=I_EXT_RANGE         ! becomes baroclinic pressure anomaly
!           AA(i,j,N)=AA(i,j,0)    ! normalized as P/(g*rho0) and defined
!           do k=1,N
!             AA(i,j,k) = AA(i,j,k) - AA(i,j,N) + z_w(i,j,N)
!           enddo
!         enddo                    ! at horizontal  PSI- and vertical RHO
!       enddo  



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=jstr,jend           ! Compute horizontal pressure gradient
        do i=istr,iend         ! components and place them into d1,d2


!           cff=0.5*g
!           do k=1,N
!             d1(k)=cff*pm(i,j)*Hz(i,j,k)*( AA(i,j  ,k)-AA(i+1,j  ,k)
!      &                                   +AA(i,j+1,k)-AA(i+1,j+1,k))
! 
!             d2(k)=cff*pn(i,j)*Hz(i,j,k)*( AA(i,j  ,k)+AA(i+1,j  ,k)
!      &                                   -AA(i,j+1,k)-AA(i+1,j+1,k))
!           enddo

         do k=1,N
           d1(k)= -Hz(i,j,k)*bx(i,j  ,k)
           d2(k)= -Hz(i,j,k)*by(i,j  ,k)
         enddo

          TauX=sustr(i,j)
          TauY=svstr(i,j)

          do k=1,N-1
            FC(k)=2.*Av(i,j,k)/(Hz(i,j,k)+Hz(i,j,k+1))
          enddo
          FC(0)=0.
          FC(N)=0.

! Block-tri-diagonal problem for u,v-velocity components which are in
! simultaneous geostrophic balance with the density field initialized
! above and in viscous-Coriolis (Ekman) balance with each other.

          b11=FC(1)                     ! Free-slip bottom boundary
          b22=FC(1)                     ! condition: everything is the
          b12=-f(i,j)*Hz(i,j,1)             ! same as it would be in k-loop
          b21= f(i,j)*Hz(i,j,1)             ! for k=1 with all (k-1)-terms
                                        ! vanishing,
          det=1./(b11*b22-b12*b21)      ! FC(k-1) --> FC(0) --> 0

          c11(1)= det*b22*FC(1) ; c12(1)=-det*b12*FC(1)
          c21(1)=-det*b21*FC(1) ; c22(1)= det*b11*FC(1)

!          d1(1)= det*( b22*d1(1) -b12*d2(1))
!          d2(1)= det*(-b21*d1(1) +b11*d2(1))

          cff1 = det*( b22*d1(1) -b12*d2(1))
          d2(1)= det*(-b21*d1(1) +b11*d2(1))
          d1(1) = cff1

          do k=2,N-1,+1                          !--> forward sweep
            b11=FC(k)    +FC(k-1)*(1.-c11(k-1))  ! of block-Gaussian
            b22=FC(k)    +FC(k-1)*(1.-c22(k-1))  ! elimination
            b12=-f(i,j)*Hz(i,j,k) -FC(k-1)*c12(k-1)
            b21= f(i,j)*Hz(i,j,k) -FC(k-1)*c21(k-1)

            det=1./(b11*b22-b12*b21)

            c11(k)= det*b22*FC(k) ; c12(k)=-det*b12*FC(k)
            c21(k)=-det*b21*FC(k) ; c22(k)= det*b11*FC(k)

            cff1=d1(k)+FC(k-1)*d1(k-1)
            cff2=d2(k)+FC(k-1)*d2(k-1)
            d1(k)= det*( b22*cff1 -b12*cff2)
            d2(k)= det*(-b21*cff1 +b11*cff2)
          enddo

          b11=          FC(N-1)*(1.-c11(N-1))    ! surface boundary
          b22=          FC(N-1)*(1.-c22(N-1))    ! condition: apply
          b12=-f(i,j)*Hz(i,j,N) -FC(N-1)*c12(N-1)    ! wind stress, i.e.,
          b21= f(i,j)*Hz(i,j,N) -FC(N-1)*c21(N-1)    ! replace what would
                                                 ! be upper interface
          det=1./(b11*b22-b12*b21)               ! implicit viscous
                                                 ! flux
          cff1=d1(N)+TauX +FC(N-1)*d1(N-1)       ! FC(k)*(u(k+1)-u(k))
          cff2=d2(N)+TauY +FC(N-1)*d2(N-1)       ! at k=N with TauX,y.
          d1(N)= det*( b22*cff1 -b12*cff2)
          d2(N)= det*(-b21*cff1 +b11*cff2)

          do k=N-1,1,-1                          !--> backsubstitution
            d1(k)=d1(k) +c11(k)*d1(k+1)+c12(k)*d2(k+1)
            d2(k)=d2(k) +c21(k)*d1(k+1)+c22(k)*d2(k+1)
          enddo

! After this moment d1,d2 contain horizontal velocity components u,v
! computed at common location on C-grid (at RHO-points).  Interpolate
! them to their native placements.

          if (i < iend .and. j >= jstr) then     ! Note: this code
            do k=1,N                             ! segment relies on
              u(i+1,j,k)=0.5*d1(k)             ! the fact that
            enddo                                ! i,j-indices are
          endif                                  ! increasing in their
          if (i >= istrU .and. j >= jstr) then   ! respective loops
            do k=1,N                             ! enclosing this part,
              u(i,j,k)=u(i,j,k)+0.5*d1(k)    ! hence for each
            enddo                                ! individual value of
          endif                                  ! "i" in computation
                                                 ! of "u" the upper
          if (i >= istr .and. j < jend) then     ! "if" statement is
            do k=1,N                             ! executed first.
              v(i,j+1,k)=0.5*d2(k)             ! Same applies to "j"
            enddo                                ! and "v".
          endif
          if (i >= istr .and. j >= jstrV) then
            do k=1,N
              v(i,j,k)=v(i,j,k)+0.5*d2(k)
            enddo
          endif
        enddo  !<-- i
      enddo  !<-- j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      do k=1,N                 ! Alternative to VTWB: compute velocity
                              ! components from geostrophic balance at
      do j=jstr,jend          ! their native C-grid locations.
        do i=istrU,iend
          ug(i,j,k)=-by(i,j,k)/f(i,j)
        enddo
      enddo
      do j=jstrV,jend
        do i=istr,iend
          vg(i,j,k)=bx(i,j,k)/f(i,j)
        enddo
      enddo
      enddo

     
!       do k=1,N                 ! Alternative to VTWB: compute velocity
!        cff=0.5*g               ! components from geostrophic balance at
!        do j=jstr,jend          ! their native C-grid locations.
!         do i=istrU,iend
!        ug(i,j,k)=-cff*(pn(i,j)+pn(i-1,j))*(AA(i,j+1,k)-AA(i,j,k))/f(i,j)
!         enddo
!        enddo
!        do j=jstrV,jend
!         do i=istr,iend
!        vg(i,j,k)=+cff*(pm(i,j)+pm(i,j-1))*(AA(i+1,j,k)-AA(i,j,k))/f(i,j)
!         enddo
!        enddo
!       enddo


      return
      end




