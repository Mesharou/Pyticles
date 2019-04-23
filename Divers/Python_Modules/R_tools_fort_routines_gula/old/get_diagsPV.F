
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute PV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_diagsPV(Lm,Mm,N, T,S, u,v, z_r,z_w,rho0,pm,pn,f,
     & pv)


      implicit none
      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k
      real*8 T(0:Lm+1,0:Mm+1,N),  S(0:Lm+1,0:Mm+1,N),
     &      u(1:Lm+1,0:Mm+1,N), v(0:Lm+1,1:Mm+1,N),
     &      pv(1:Lm+1,1:Mm+1,0:N),
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      f(0:Lm+1,0:Mm+1),
     &      z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N),
     &      rho1(0:Lm+1,0:Mm+1,N), drdz(0:Lm+1,0:Mm+1,0:N),
     &      drdx(1:Lm+1,0:Mm+1,N), drdy(0:Lm+1,1:Mm+1,N),
     &      dudz(1:Lm+1,0:Mm+1,0:N), dvdz(0:Lm+1,1:Mm+1,0:N),
     &      absvrt(1:Lm+1,1:Mm+1,N), 
     &      qp1(0:Lm+1,0:Mm+1,N),
     &      dpth,cff,cff2, Tt,Ts,sqrtTs, rho0, K0, dr00,
     &      dvdx, dudy,
     &      cffi(0:Lm+1), cffj(0:Mm+1),
     &      var1, var2

      real*8, parameter :: r00=999.842594,   r01=6.793952E-2,
     &          r02=-9.095290E-3, r03=1.001685E-4,  r04=-1.120083E-6,
     &                                              r05=6.536332E-9,
     &          r10=0.824493,     r11=-4.08990E-3,  r12=7.64380E-5,
     &                            r13=-8.24670E-7,  r14=5.38750E-9,
     &          rS0=-5.72466E-3,  rS1=1.02270E-4,   rS2=-1.65460E-6,
     &          r20=4.8314E-4,
     &          K00=19092.56,     K01=209.8925,     K02=-3.041638,
     &                            K03=-1.852732e-3, K04=-1.361629e-5,
     &          K10=104.4077,     K11=-6.500517,    K12=0.1553190,
     &                                              K13=2.326469e-4,
     &          KS0=-5.587545,    KS1=+0.7390729,   KS2=-1.909078e-2,
     &          qp2=0.0000172, g=9.81

#  ifdef DUKO_2001
      real rho1_0, K0_Duk
#  endif

      integer numthreads, trd, chunk_size, margin, jstr,jend


Cf2py intent(in) Lm,Mm,N, T,S, u,v,z_r,z_w,rho0,pm,pn,f
Cf2py intent(out) pv


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! COMPUTE NEUTRAL DENSITY GRADIENTS 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#  ifdef DUKO_2001
      Tt=3.8D0
      Ts=34.5D0
      sqrtTs=sqrt(Ts)

      dr00=r00-1000.D0

      rho1_0=dr00 +Tt*( r01+Tt*( r02+Tt*( r03+Tt*( r04+Tt*r05 ))))
     &                            +Ts*( R10+Tt*( r11+Tt*( r12+Tt*(
     &                                              r13+Tt*r14 )))
     &                   +sqrtTs*( rS0+Tt*( rS1+Tt*rS2 ))+Ts*r20 )

      K0_Duk= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     &       +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     &            +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))
#  endif

      dr00=r00-rho0


      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1



      do j=jmin,jmax



!---------------------------------------------------------------------------------------
        do k=1,N
          do i=imin,imax
            Tt=T(i,j,k)                   

            Ts=max(S(i,j,k), 0.)
            sqrtTs=sqrt(Ts)

            rho1(i,j,k)=( dr00 +Tt*( r01+Tt*( r02+Tt*( r03+Tt*(
     &                                           r04+Tt*r05 ))))
     &                         +Ts*( r10+Tt*( r11+Tt*( r12+Tt*(
     &                                            r13+Tt*r14 )))
     &                              +sqrtTs*(rS0+Tt*(
     &                                   rS1+Tt*rS2 ))+Ts*r20 ))

            K0= Tt*( K01+Tt*( K02+Tt*( K03+Tt*K04 )))
     &         +Ts*( K10+Tt*( K11+Tt*( K12+Tt*K13 ))
     &              +sqrtTs*( KS0+Tt*( KS1+Tt*KS2 )))
 


            qp1(i,j,k)= 0.1D0*(rho0+rho1(i,j,k))*(K0_Duk-K0)
     &                               /((K00+K0)*(K00+K0_Duk))



          enddo
        enddo

 !---------------------------------------------------------------------------------------


        cff=g/rho0

        do k=1,N-1
          do i=imin,imax

            dpth=z_w(i,j,N)-0.5*(z_r(i,j,k+1)+z_r(i,j,k))

            cff2=( rho1(i,j,k+1)-rho1(i,j,k)             ! Elementary
     &                        +(qp1(i,j,k+1)-qp1(i,j,k)) ! adiabatic 
     &                            *dpth*(1.-2.*qp2*dpth) ! difference
     &                       )


            drdz(i,j,k)=-cff*cff2 / (z_r(i,j,k+1)-z_r(i,j,k))

          enddo
        enddo
      enddo    ! <-- j




!---------------------------------------------------------------------------------------
      do k=N,1,-1
       do i=imin+1,imax

        do j=jmin,jmax
            dpth=0.5*( z_w(i,j,N)+z_w(i-1,j,N)
     &                -z_r(i,j,k)-z_r(i-1,j,k))

            cffj(j)=-cff*( rho1(i,j,k)-rho1(i-1,j,k) ! Elementary
     &                +(qp1(i,j,k)-qp1(i-1,j,k))         ! adiabatic
     &                     *dpth*(1.-qp2*dpth) )         ! difference
     &                    *0.5*(pm(i,j)+pm(i-1,j))
         enddo


         do j=jmin+1,jmax
            drdx(i,j,k)= 0.5 * (cffj(j) + cffj(j-1))
         enddo

       enddo


!---------------------------------------------------------------------------------------

        do j=jmin+1,jmax
          do i=imin,imax
            dpth=0.5*( z_w(i,j,N)+z_w(i,j-1,N)
     &                -z_r(i,j,k)-z_r(i,j-1,k))
            
            cffi(i)=-cff*( rho1(i,j,k)-rho1(i,j-1,k)   ! Elementary
     &                +(qp1(i,j,k)-qp1(i,j-1,k))         ! adiabatic
     &                     *dpth*(1.-qp2*dpth) )  
     &                    *0.5*(pn(i,j)+pn(i,j-1)) 
          enddo



          do i=imin+1,imax
            drdy(i,j,k)= 0.5 * (cffi(i) + cffi(i-1))
          enddo



        enddo
      enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! COMPUTE VORTICITY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do k=N,1,-1
       do j=jmin+1,jmax
        do i=imin+1,imax


            dpth=-0.5*( z_w(i,j,N)+z_w(i-1,j,N)
     &                -z_w(i,j,k)-z_w(i-1,j,k))

!            write(*,*) k,z_w(i,j,N),z_w(i,j,k) ,dpth

            CALL interp_1d(N,v(i,j,:),z_r(i,j,:),z_w(i,j,:)
     &                ,dpth,var1,1,1)
            CALL interp_1d(N,v(i-1,j,:),z_r(i-1,j,:),z_w(i-1,j,:)
     &                ,dpth,var2,1,1)
    


            dvdx = (var1 - var2) *0.5*(pm(i,j)+pm(i-1,j))




            
            dpth = -0.5*( z_w(i,j,N)+z_w(i,j-1,N)
     &                -z_w(i,j,k)-z_w(i,j-1,k))

            CALL interp_1d(N,u(i,j,:),z_r(i,j,:),z_w(i,j,:)
     &                ,dpth,var1,1,0)
            CALL interp_1d(N,u(i,j-1,:),z_r(i,j-1,:),z_w(i,j-1,:)
     &                ,dpth,var2,1,0)
            
            dudy = (var1 - var2)*(pn(i,j)+pn(i,j-1))




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          absvrt(i,j,k) = 0.25*(f(i,j)+f(i-1,j)+f(i,j-1)+f(i-1,j-1))
     &                  +dvdx-dudy

            if (k.eq.1) then
                write(*,*) k,dpth,z_w(i,j,k),z_w(i,j-1,k),var1,var2
     &                   ,absvrt(i,j,k)
            end if

        enddo
       enddo
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Vertical velocity gradients
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      do k=1,N-1
       do j=jmin+1,jmax
        do i=imin+1,imax



        ! compute dudz on psi - w grid

        dudz(i,j,k) = 0.5*((u(i,j,k+1)-u(i,j,k)) / 
     &       (z_r(i,j,k+1)+z_r(i-1,j,k+1)-z_r(i,j,k)-z_r(i-1,j,k)) * 2.
     &              + (u(i,j-1,k+1)-u(i,j-1,k)) / 
     &          (z_r(i,j-1,k+1)+z_r(i-1,j-1,k+1)
     &          -z_r(i,j-1,k)-z_r(i-1,j-1,k)) * 2.)


        ! compute dvdz on psi - w grid

        dvdz(i,j,k) = 0.5*((v(i,j,k+1)-v(i,j,k)) /
     &       (z_r(i,j,k+1)+z_r(i,j-1,k+1)-z_r(i,j,k)-z_r(i,j-1,k)) * 2.
     &              + (v(i-1,j,k+1)-v(i-1,j,k)) /
     &          (z_r(i-1,j,k+1)+z_r(i-1,j-1,k+1)
     &          -z_r(i-1,j,k)-z_r(i-1,j-1,k)) * 2.)



        enddo
       enddo
      enddo







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! COMPUTE NEUTRAL DENSITY GRADIENTS 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      do k=1,N-1
        do j=jmin+1,jmax
          do i=imin+1,imax


            pv(i,j,k) = absvrt(i,j,k) *
     &       0.25*(drdz(i,j,k)+drdz(i-1,j,k)+drdz(i,j-1,k)+drdz(i-1,j-1,k))
     &          - dvdz(i,j,k) * 0.5*(drdx(i,j,k) + drdx(i,j,k+1))
     &          + dudz(i,j,k) * 0.5*(drdy(i,j,k) + drdy(i,j,k+1))

          enddo
        enddo
      enddo


      !write(*,*) 'z_w(:,:,0)',z_w(:,:,0)






      return
      end







 
 
