
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute PV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_Jdiab_fromdiags(Lm,Mm,N, Tdiab, Sdiab, alpha, beta,
     & u,v, z_r,z_w,rho0,pm,pn,f,Jdiab)

      implicit none
      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k
      real*8 Tdiab(0:Lm+1,0:Mm+1),  Sdiab(0:Lm+1,0:Mm+1),
     &      alpha(0:Lm+1,0:Mm+1), beta(0:Lm+1,0:Mm+1),
     &      u(1:Lm+1,0:Mm+1,N), v(0:Lm+1,1:Mm+1,N),
     &      Jdiab(1:Lm+1,1:Mm+1),f(0:Lm+1,0:Mm+1),
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      hbls(0:Lm+1,0:Mm+1),
     &      z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N),
     &      dpth,cff,cff2, Tt,Ts,sqrtTs, rho0, K0, dr00,
     &      var1, var2,var3, var4,cff3,
     &      absvrt(1:Lm+1,1:Mm+1),absvrt0(1:Lm+1,1:Mm+1)

      real*8, parameter :: g=9.81


Cf2py intent(in) Lm,Mm,N,Tdiab, Sdiab, alpha, beta, u,v,z_r,z_w,rho0,pm,pn,f
Cf2py intent(out) Jdiab


      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! COMPUTE VORTICITY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!---------------------------------------------------------------------------------------

       do i=imin+1,imax
        do j=jmin+1,jmax

            dpth=0.25*(z_r(i,j,1)+z_r(i-1,j,1)
     &           +     z_r(i-1,j-1,1)+z_r(i,j-1,1))

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


c            absvrt0(i,j)= cff 
c     &                + (v(i,j,N)-v(i-1,j,N)) * cff2
c     &                - (u(i,j,N)-u(i,j-1,N)) * cff3

c            write(*,*) i,j,dpth,z_r(i,j,1),absvrt(i,j), absvrt0(i,j)

         if (var1.gt.1e3) then
            CALL interp_1d(N,v(i,j,:)
     &                ,0.5*(z_r(i,j,:)+z_r(i,j-1,:))
     &                ,0.5*(z_w(i,j,:)+z_w(i,j-1,:))
     &                ,dpth,var1,1,1)

         write(*,*)
         write(*,*) var1,dpth,N
         write(*,*) v(i,j,:)
         write(*,*) 0.5*(z_r(i,j,:)+z_r(i,j-1,:))
         write(*,*) 0.5*(z_w(i,j,:)+z_w(i,j-1,:))
         write(*,*)
         write(*,*) 
         write(*,*)
         endif
         enddo
       enddo


!---------------------------------------------------------------------------------------

 



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! COMPUTE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        do j=jmin,jmax
          do i=imin,imax
                Tdiab(i,j) = Tdiab(i,j)*alpha(i,j)
                Sdiab(i,j) = Sdiab(i,j)*beta(i,j)
          enddo
        enddo


        do j=jmin+1,jmax
          do i=imin+1,imax


            cff = g/(0.25*(z_w(i,j,1) + z_w(i-1,j,1)
     &                 + z_w(i,j-1,1) + z_w(i-1,j-1,1)
     &                 - z_w(i,j,0) - z_w(i-1,j,0)
     &                 - z_w(i,j-1,0) - z_w(i-1,j-1,0)))



            Jdiab(i,j) = cff * absvrt(i,j) * 0.25 *  (
     &           (Tdiab(i,j)+Tdiab(i-1,j)+Tdiab(i,j-1)+Tdiab(i-1,j-1))
     &         - (Sdiab(i,j)+Sdiab(i-1,j)+Sdiab(i,j-1)+Sdiab(i-1,j-1))
     &           )

c            if (Jdiab(i,j).gt.1e3) then
c             write(*,*) Jdiab(i,j),cff,absvrt(i,j),Tdiab(i,j),Sdiab(i,j)
c            endif
          enddo
        enddo


      return
      end







 
 
