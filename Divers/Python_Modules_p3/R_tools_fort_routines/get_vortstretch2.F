
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute planetary stretching term from the vorticity balance
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_vortstretch2(Lm,Mm, u,v,H,pm,pn,f,rot)


      implicit none
      integer Lm,Mm, imin,imax,jmin,jmax, i,j,k
      real*8 rot(1:Lm+1,1:Mm+1),H(0:Lm+1,0:Mm+1),
     &      u(1:Lm+1,0:Mm+1), v(0:Lm+1,1:Mm+1),
     &      Hu(0:Lm+2), Hv(0:Mm+2),
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      f(0:Lm+1,0:Mm+1),
     &      dudx(0:Lm+1,0:Mm+1), dvdy(0:Lm+1,1:Mm+1),
     &      cff1,cff2,cff


Cf2py intent(in) Lm,Mm,u,v,H,pm,pn,f
Cf2py intent(out) rot

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=jmin,jmax
        do i=imin+1,imax

            !cff = 0.5*(pm(i,j)+pm(i-1,j))*0.5*(pn(i,j)+pn(i-1,j))
            cff = 0.5*(pn(i,j)+pn(i-1,j))
            Hu(i) = 0.5*(H(i,j)+ H(i-1,j)) * u(i,j) / cff

        enddo !<- i
            Hu(imin) = Hu(imin+1)
            Hu(imax+1) = Hu(imax)

        do i=imin,imax
            dudx(i,j) = f(i,j) *(Hu(i+1) - Hu(i))* pm(i,j)*pn(i,j)
        enddo !<- i

      enddo !<- j

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i=imin,imax

        do j=jmin+1,jmax

            !cff = 0.5*(pm(i,j)+pm(i,j-1))*0.5*(pn(i,j)+pn(i,j-1))
            cff = 0.5*(pm(i,j)+pm(i,j-1))
            Hv(j) = 0.5*(H(i,j)+ H(i,j-1)) * v(i,j) / cff
        enddo !<- j
            Hv(jmin) = Hv(jmin+1)
            Hv(jmax+1) = Hv(jmax)

        do j=jmin,jmax
            dvdy(i,j) = f(i,j) *(Hv(j+1) - Hv(j))* pn(i,j)*pm(i,j)
        enddo !<- j

      enddo !<- i

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do i=imin+1,imax
        do j=jmin+1,jmax

         cff =0.25*(pn(i,j) + pn(i-1,j) + pn(i,j-1) + pn(i-1,j-1))
     &      * 0.25*(pm(i,j) + pm(i-1,j) + pm(i,j-1) + pm(i-1,j-1))

              rot(i,j) = 0.25*(dudx(i,j) + dudx(i-1,j) 
     &              + dudx(i,j-1) + dudx(i-1,j-1)
     &              + dvdy(i,j) + dvdy(i-1,j)
     &              + dvdy(i,j-1) + dvdy(i-1,j-1))*cff

        enddo !<- i
      enddo !<- j






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
           
      return
      end







 
 
