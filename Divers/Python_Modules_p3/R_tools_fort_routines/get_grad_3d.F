
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute horizontal GRAD amplitude for a 3d var function
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_grad_3d(Lm,Mm,N, psi,pm,pn,grad)


      implicit none
      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k
      real*8 grad(1:Lm+1,1:Mm+1,N),
     &      psi(0:Lm+1,0:Mm+1,N), 
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      dvdx, dudy


Cf2py intent(in) Lm,Mm,N,psi,pm,pn
Cf2py intent(out) grad

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1

        do k=1,N
          do j=jmin+1,jmax
            do i=imin+1,imax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            dvdx = 0.5*((psi(i,j,k) - psi(i-1,j,k)) * 0.5 * 
     &           (pm(i,j)+pm(i-1,j))
     &           + (psi(i,j-1,k)-psi(i-1,j-1,k)) * 0.5 * 
     &           (pm(i,j-1)+pm(i-1,j-1)))

            dudy = 0.5*((psi(i,j,k) - psi(i,j-1,k) )* 0.5 * 
     &           (pn(i,j)+pn(i,j-1))
     &           + (psi(i-1,j,k) - psi(i-1,j-1,k) )* 0.5 * 
     &           (pn(i-1,j)+pn(i-1,j-1)))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            grad(i,j,k) = sqrt(dvdx**2 + dudy**2)

            !write(*,*) i,j,dvdx, dudy, rot(i,j)


          enddo !<- i
        enddo !<- j
      enddo !<- k


      return
      end







 
 
