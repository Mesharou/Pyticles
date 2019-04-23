
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute ROT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_rot(Lm,Mm, u,v,pm,pn,rot)


      implicit none
      integer Lm,Mm, imin,imax,jmin,jmax, i,j
      real*8 rot(1:Lm+1,1:Mm+1),
     &      u(1:Lm+1,0:Mm+1), v(0:Lm+1,1:Mm+1),
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      dvdx, dudy


Cf2py intent(in) Lm,Mm,u,v,pm,pn
Cf2py intent(out) rot

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1


        do j=jmin+1,jmax
          do i=imin+1,imax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            dvdx = (v(i,j) - v(i-1,j)) * 0.25 * 
     &           (pm(i,j)+pm(i-1,j)+pm(i,j-1)+pm(i-1,j-1))



            dudy = (u(i,j) - u(i,j-1) )* 0.25 * 
     &           (pn(i,j)+pn(i-1,j)+pn(i,j-1)+pn(i-1,j-1))



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            rot(i,j) = dvdx - dudy

            !write(*,*) i,j,dvdx, dudy, rot(i,j)


        enddo !<- i
      enddo !<- j



      return
      end







 
 
