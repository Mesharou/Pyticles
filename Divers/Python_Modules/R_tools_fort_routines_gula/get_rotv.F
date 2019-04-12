
!======================================================================
!
! Compute rotational in a flux consistent form
!
! 
!
! - updated 16/08/19
!======================================================================

      subroutine get_rotv(Lm,Mm, u,v,pm,pn,rot)


      implicit none
      integer Lm,Mm, imin,imax,jmin,jmax, i,j
      real*8 rot(1:Lm+1,1:Mm+1),
     &      u(1:Lm+1,0:Mm+1), v(0:Lm+1,1:Mm+1),
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      dvdx, dudy, cff


Cf2py intent(in) Lm,Mm,u,v,pm,pn
Cf2py intent(out) rot

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1


!!!!!!!!!!!!!!!!!!!!!

      do j=jmin,jmax
        do i=imin+1,imax

         u(i,j) = u(i,j) * 2./ (pm(i,j) + pm(i-1,j))

        enddo !<- i
      enddo !<- j

!!!!!!!!!!!!!!!!!!!!!

      do j=jmin+1,jmax
        do i=imin,imax

         v(i,j) = v(i,j) * 2./ (pn(i,j) + pn(i,j-1))

        enddo !<- i
      enddo !<- j


        do j=jmin+1,jmax
          do i=imin+1,imax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

          cff = 0.25*(pm(i,j) + pm(i-1,j) + pm(i,j-1) + pm(i-1,j-1))
     &        * 0.25*(pn(i,j) + pn(i-1,j) + pn(i,j-1) + pn(i-1,j-1))

            dvdx = (v(i,j) - v(i-1,j)) * cff
            dudy = (u(i,j) - u(i,j-1)) * cff

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            rot(i,j) = dvdx - dudy

        enddo !<- i
      enddo !<- j



      return
      end







 
 
