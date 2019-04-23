!-------------------------------------------------------------------------
!
!      compute pv div flux = [ absvrt(3d) . b] on a psi grid and a rho-level
!
!-------------------------------------------------------------------------
!
!
!
!
!
!-------------------------------------------------------------------------

      subroutine get_pvdiv_flux(Lm,Mm,N,b,
     & u,v, z_r,z_w,pm,pn,f,iz,Jpvdiv)

      implicit none
      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k,iz
      real*8 b(0:Lm+1,0:Mm+1),
     &      u(1:Lm+1,0:Mm+1,N), v(0:Lm+1,1:Mm+1,N),
     &      work(1:Lm+1,1:Mm+1,N),
     &      uz(1:Lm+1,0:Mm+1), vz(0:Lm+1,1:Mm+1),
     &      Jpvdiv(1:Lm+1,1:Mm+1),f(0:Lm+1,0:Mm+1),
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N),
     &      dpth,cff,cff2,cff3,
     &      var1, var2,var3, var4,
     &      absvrt(1:Lm+1,1:Mm+1),
     &      dx(1:Lm+1,1:Mm+1),dy(1:Lm+1,1:Mm+1),
     &      dz(1:Lm+1,1:Mm+1)

      real*8, parameter :: g=9.81

Cf2py intent(in) Lm,Mm,N,b, u,v,z_r,z_w,pm,pn,f,iz
Cf2py intent(out) Jpvdiv

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1

!-------------------------------------------------------------------------
! COMPUTE VERTICAL VORTICITY
!-------------------------------------------------------------------------

       do i=imin+1,imax
        do j=jmin+1,jmax

            dpth=0.25*(z_r(i,j,iz+1)+z_r(i-1,j,iz+1)
     &           +     z_r(i-1,j-1,iz+1)+z_r(i,j-1,iz+1))

            CALL interp_1d(N,v(i,j,:)
     &                ,0.5*(z_r(i,j,:)+z_r(i,j-1,:))
     &                ,0.5*(z_w(i,j,:)+z_w(i,j-1,:))
     &                ,dpth,var1,1,0)

            CALL interp_1d(N,v(i-1,j,:)
     &                ,0.5*(z_r(i-1,j,:)+z_r(i-1,j-1,:))
     &                ,0.5*(z_w(i-1,j,:)+z_w(i-1,j-1,:))
     &                ,dpth,var2,1,0)

            CALL interp_1d(N,u(i,j,:)
     &                ,0.5*(z_r(i-1,j,:)+z_r(i,j,:))
     &                ,0.5*(z_w(i-1,j,:)+z_w(i,j,:))
     &                ,dpth,var3,1,0)

            CALL interp_1d(N,u(i,j-1,:)
     &                ,0.5*(z_r(i,j-1,:)+z_r(i-1,j-1,:))
     &                ,0.5*(z_w(i,j-1,:)+z_w(i-1,j-1,:))
     &                ,dpth,var4,1,0)

          cff = 0.25*(f(i,j) + f(i-1,j) + f(i,j-1) + f(i-1,j-1))
          cff2 = 0.25*(pm(i,j) + pm(i-1,j) + pm(i,j-1) + pm(i-1,j-1))
          cff3 = 0.25*(pn(i,j) + pn(i-1,j) + pn(i,j-1) + pn(i-1,j-1))

            absvrt(i,j)= cff 
     &                + (var1-var2) * cff2
     &                - (var3-var4) * cff3

         enddo
       enddo

!-------------------------------------------------------------------------
! COMPUTE VORTICITY
!-------------------------------------------------------------------------

       do i=imin+1,imax
        do j=jmin+1,jmax

            dpth=0.25*(z_r(i,j,iz+1)+z_r(i-1,j,iz+1)
     &           +     z_r(i-1,j-1,iz+1)+z_r(i,j-1,iz+1))

         do k=1,N-1
         work(i,j,k) = 2.*(v(i,j,k+1)+v(i-1,j,k+1)-v(i,j,k)-v(i-1,j,k))
     &          / ( z_r(i,j,k+1)+z_r(i-1,j,k+1)-z_r(i,j,k)-z_r(i-1,j,k)
     &   + z_r(i,j-1,k+1)+z_r(i-1,j-1,k+1)-z_r(i,j-1,k)-z_r(i-1,j-1,k))
         enddo

            CALL interp_1d(N,work(i,j,:)
     &                ,0.25*(z_w(i,j,1:N-1)+z_w(i-1,j-1,1:N-1)
     &                      +z_w(i-1,j,1:N-1)+z_w(i,j-1,1:N-1))
     &                ,0.25*(z_r(i,j,1:N)+z_r(i-1,j-1,1:N)
     &                      +z_r(i-1,j,1:N)+z_r(i,j-1,1:N))
     &                ,dpth,vz(i,j),1,0)

         do k=1,N-1
         work(i,j,k) = 2.*(u(i,j,k+1)+u(i,j-1,k+1)-u(i,j,k)-u(i,j-1,k))
     &          / ( z_r(i,j,k+1)+z_r(i-1,j,k+1)-z_r(i,j,k)-z_r(i-1,j,k)
     &   + z_r(i,j-1,k+1)+z_r(i-1,j-1,k+1)-z_r(i,j-1,k)-z_r(i-1,j-1,k))
         enddo

            CALL interp_1d(N,work(i,j,:)
     &                ,0.25*(z_w(i,j,1:N-1)+z_w(i-1,j-1,1:N-1)
     &                      +z_w(i-1,j,1:N-1)+z_w(i,j-1,1:N-1))
     &                ,0.25*(z_r(i,j,1:N)+z_r(i-1,j-1,1:N)
     &                      +z_r(i-1,j,1:N)+z_r(i,j-1,1:N))
     &                ,dpth,uz(i,j),1,0)

         enddo
       enddo

!-------------------------------------------------------------------------
! HEIGHT GRADIENT
!-------------------------------------------------------------------------

        do j=jmin+1,jmax
          do i=imin+1,imax

            dx(i,j) = 0.5* (z_r(i,j,iz+1) + z_r(i,j-1,iz+1)
     &               - z_r(i-1,j,iz+1)- z_r(i-1,j-1,iz+1))
     &      * 0.25*(pm(i,j)+pm(i,j-1)+ pm(i-1,j)+pm(i-1,j-1))

            dy(i,j)  = 0.5* (z_r(i,j,iz+1) + z_r(i-1,j,iz+1)
     &               - z_r(i,j-1,iz+1)- z_r(i-1,j-1,iz+1))
     &      * 0.25*(pn(i,j)+pn(i,j-1)+ pn(i-1,j)+pn(i-1,j-1))
 
            cff = 1./sqrt(1 + dx(i,j)**2 + dy(i,j)**2)

            dx(i,j) = - cff * dx(i,j)
            dy(i,j) = - cff * dy(i,j) 
            dz(i,j) = cff

          enddo
        enddo

!-------------------------------------------------------------------------
! COMPUTE VORTICITY
!-------------------------------------------------------------------------

        do j=jmin+1,jmax
          do i=imin+1,imax

            Jpvdiv(i,j) = 0.25 * 
     &           (b(i,j)+b(i-1,j)+b(i,j-1)+b(i-1,j-1))
     &           * ( absvrt(i,j) * dz(i,j)
     &               +  uz(i,j)  * dy(i,j)
     &               -  vz(i,j)  * dx(i,j) )


          enddo
        enddo

      return
      end







 
 
