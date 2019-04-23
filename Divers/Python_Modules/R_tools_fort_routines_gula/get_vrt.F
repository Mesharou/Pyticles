
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute PV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_absvrt(Lm,Mm,N,u,v, z_r,z_w,pm,
     & pn,f,absvrt, absvrt0,dpth,var1,var2,var3,var4)


      implicit none
      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k
      real*8 stflx(0:Lm+1,0:Mm+1),  ssflx(0:Lm+1,0:Mm+1),
     &      u(1:Lm+1,0:Mm+1,N), v(0:Lm+1,1:Mm+1,N),
     &      J1(1:Lm+1,1:Mm+1),f(0:Lm+1,0:Mm+1),
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      hbls(0:Lm+1,0:Mm+1),cff3,
     &      z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N),
     &      dpth(1:Lm+1,1:Mm+1),cff,cff2, Tt,Ts,sqrtTs, rho0, K0, dr00,
     &      var1(1:Lm+1,1:Mm+1), var2(1:Lm+1,1:Mm+1),
     &      var3(1:Lm+1,1:Mm+1), var4(1:Lm+1,1:Mm+1),
     &      absvrt(1:Lm+1,1:Mm+1),absvrt0(1:Lm+1,1:Mm+1)



Cf2py intent(in) Lm,Mm,N,  u,v,z_r,z_w,pm,pn,f
Cf2py intent(out) absvrt, absvrt0, dpth,var1,var2,var3,var4


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

            dpth(i,j)=0.25*(z_r(i,j,N)+z_r(i-1,j,N)
     &           +     z_r(i-1,j-1,N)+z_r(i,j-1,N))

            CALL interp_1d(N,v(i,j,:)
     &                ,0.5*(z_r(i,j,:)+z_r(i,j-1,:))
     &                ,0.5*(z_w(i,j,:)+z_w(i,j-1,:))
     &                ,dpth(i,j),var1(i,j),1,0)

            CALL interp_1d(N,v(i-1,j,:)
     &                ,0.5*(z_r(i-1,j,:)+z_r(i-1,j-1,:))
     &                ,0.5*(z_w(i-1,j,:)+z_w(i-1,j-1,:))
     &                ,dpth(i,j),var2(i,j),1,0)


            CALL interp_1d(N,u(i,j,:)
     &                ,0.5*(z_r(i-1,j,:)+z_r(i,j,:))
     &                ,0.5*(z_w(i-1,j,:)+z_w(i,j,:))
     &                ,dpth(i,j),var3(i,j),1,0)

            CALL interp_1d(N,u(i,j-1,:)
     &                ,0.5*(z_r(i,j-1,:)+z_r(i-1,j-1,:))
     &                ,0.5*(z_w(i,j-1,:)+z_w(i-1,j-1,:))
     &                ,dpth(i,j),var4(i,j),1,0)

            cff = 0.25*(f(i,j) + f(i-1,j) + f(i,j-1) + f(i-1,j-1))
            cff2 = 0.25*(pm(i,j) + pm(i-1,j) + pm(i,j-1) + pm(i-1,j-1))
            cff3 = 0.25*(pn(i,j) + pn(i-1,j) + pn(i,j-1) + pn(i-1,j-1))

            absvrt(i,j)= cff 
     &                + (var1(i,j)-var2(i,j)) * cff2
     &                - (var3(i,j)-var4(i,j)) * cff3


            absvrt0(i,j)= cff 
     &                + (v(i,j,N)-v(i-1,j,N)) * cff2
     &                - (u(i,j,N)-u(i,j-1,N)) * cff3

          write(*,*) i,j,dpth(i,j),z_r(i,j,N),absvrt(i,j), absvrt0(i,j)

          write(*,*) i,j,v(i,j,:)
          write(*,*) i,j,0.5*(z_r(i,j,:)+z_r(i,j-1,:))
          write(*,*) i,j,0.5*(z_w(i,j,:)+z_w(i,j-1,:))
          write(*,*) i,j,var1(i,j)

         enddo

       enddo


!---------------------------------------------------------------------------------------

 



      return
      end







 
 
