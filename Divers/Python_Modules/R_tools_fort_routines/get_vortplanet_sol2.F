
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute planetary vorticity balance term
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! Has not been modifed yet!!!


      subroutine get_vortplanet_sol2(Lm,Mm, u,v,H,pm,pn,f,rot)


      implicit none
      integer Lm,Mm, imin,imax,jmin,jmax, i,j,k
      real*8 rot(1:Lm+1,1:Mm+1), H(0:Lm+1,0:Mm+1),
     &      u(1:Lm+1,0:Mm+1), v(0:Lm+1,1:Mm+1),
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      f(0:Lm+1,0:Mm+1),
     &      cff1,cff2


Cf2py intent(in) Lm,Mm,u,v,H,pm,pn,f
Cf2py intent(out) rot

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




      do j=jmin+1,jmax
        do i=imin+1,imax

            cff1 = 0.25*(H(i,j)+ H(i-1,j)) * u(i,j) * 
     &          (f(i,j) - f(i-1,j)) * 0.5 * (pm(i,j)+ pm(i-1,j))
     &         + 0.25*(H(i,j-1)+ H(i-1,j-1)) * u(i,j-1) * 
     &          (f(i,j-1) - f(i-1,j-1)) * 0.5 * (pm(i,j-1)+ pm(i-1,j-1))


            cff2 = 0.25*(H(i,j)+ H(i,j-1)) * v(i,j) * 
     &          (f(i,j) - f(i,j-1)) * 0.5 * (pn(i,j)+ pn(i,j-1))
     &         + 0.25*(H(i-1,j)+ H(i-1,j-1)) * v(i-1,j) * 
     &          (f(i-1,j) - f(i-1,j-1)) * 0.5 * (pn(i-1,j)+ pn(i-1,j-1))


            rot(i,j) = cff1 + cff2

        enddo !<- i
      enddo !<- j






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      return
      end







 
 
