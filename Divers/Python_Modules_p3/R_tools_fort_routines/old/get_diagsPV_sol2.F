
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute PV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_diagsPV_sol2(Lm,Mm,N, T,S, u,v, z_r,z_w,rho0,pm,pn,
     & f,pv)


      implicit none
      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k

      real*4 T(0:Lm+1,0:Mm+1,N),  S(0:Lm+1,0:Mm+1,N),
     &      u(1:Lm+1,0:Mm+1,N), v(0:Lm+1,1:Mm+1,N),
     &      z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N),
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      f(0:Lm+1,0:Mm+1),
     &      pv(1:Lm+1,1:Mm+1,0:N),
     &      rho1(0:Lm+1,0:Mm+1,N), qp1(0:Lm+1,0:Mm+1,N),
     &      drdx(1:Lm+1,0:Mm+1,N), drdy(0:Lm+1,1:Mm+1,N),
     &      drdz(0:Lm+1,0:Mm+1,0:N),
     &      dudz, dvdz, absvrt,
     &      dpth,cff,cff1,cff2, Tt,Ts,sqrtTs, rho0, K0, dr00,
     &      dvdx, dudy, zlev,
     &      cffi(0:Lm+1), cffj(0:Mm+1),
     &      var1, var2, var3, var4

      real*4, parameter :: r00=999.842594,   r01=6.793952E-2,
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

      real*4, parameter :: verbose=0

      real*4, parameter :: FillValue = -9999.

Cf2py intent(in) Lm,Mm,N, T,S, u,v,z_r,z_w,rho0,pm,pn,f
Cf2py intent(out) pv



      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (verbose.eq.1) write(*,*) 'COMPUTE NEUTRAL DENSITY GRADIENTS'
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (verbose.eq.1) write(*,*) 'density'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        cff=g/rho0



      do j=jmin,jmax
!---------------------------------------------------------------------------------------
          do i=imin,imax

            do k=1,N

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





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        cff=g/rho0
        do i=imin,imax
          do k=1,N-1

            dpth=z_w(i,j,N)-0.5*(z_r(i,j,k+1)+z_r(i,j,k))
 

            cff2=( rho1(i,j,k+1)-rho1(i,j,k)             ! Elementary
     &                        +(qp1(i,j,k+1)-qp1(i,j,k)) ! adiabatic 
     &                            *dpth*(1.-2.*qp2*dpth) ! difference
     &                       )


            drdz(i,j,k)=-cff*cff2 / (z_r(i,j,k+1)-z_r(i,j,k))

          enddo
        enddo
      enddo    ! <-- j




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Buoyancy gradients
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do k=N-1,1,-1
        do i=imin+1,imax
          do j=jmin,jmax




            zlev = 0.5*(z_w(i,j,k)+z_w(i-1,j,k))
            dpth = 0.5*(z_w(i,j,N)+z_w(i-1,j,N)) - zlev

         !!write(*,*) 'all',rho1(i,j,k),z_r(i,j,k),z_w(i,j,k)

#ifdef PV_CUBIC
            CALL interp_1d(N,rho1(i,j,:),z_r(i,j,:),z_w(i,j,:)
     &                ,zlev,var1,1,0)
            CALL interp_1d(N,rho1(i-1,j,:),z_r(i-1,j,:),z_w(i-1,j,:)
     &                ,zlev,var2,1,0)
            CALL interp_1d(N,qp1(i,j,:),z_r(i,j,:),z_w(i,j,:)
     &                ,zlev,var3,1,0)
            CALL interp_1d(N,qp1(i-1,j,:),z_r(i-1,j,:),z_w(i-1,j,:)
     &                ,zlev,var4,1,0)
#else
            var1 = ((z_r(i,j,k+1) - zlev)*rho1(i,j,k) + 
     &             (zlev - z_r(i,j,k))*rho1(i,j,k+1) )
     &             /(z_r(i,j,k+1) - z_r(i,j,k))

            var2 = ((z_r(i-1,j,k+1) - zlev)*rho1(i-1,j,k) + 
     &             (zlev - z_r(i-1,j,k))*rho1(i-1,j,k+1) )
     &             /(z_r(i-1,j,k+1) - z_r(i-1,j,k))

            var3= ((z_r(i,j,k+1) - zlev)*qp1(i,j,k) + 
     &             (zlev - z_r(i,j,k))*qp1(i,j,k+1) )
     &             /(z_r(i,j,k+1) - z_r(i,j,k))

            var4 = ((z_r(i-1,j,k+1) - zlev)*qp1(i-1,j,k) + 
     &             (zlev - z_r(i-1,j,k))*qp1(i-1,j,k+1) )
     &             /(z_r(i-1,j,k+1) - z_r(i-1,j,k))

#endif

            cffj(j)=-cff*( var1 - var2 ! Elementary
     &                +(var3-var4)         ! adiabatic
     &                     *dpth*(1.-qp2*dpth) )         ! difference
     &                    *0.5*(pm(i,j)+pm(i-1,j))



          enddo



         do j=jmin+1,jmax

            drdx(i,j,k)= 0.5 * (cffj(j) + cffj(j-1))
            !!write(*,*) 'all3',dpth,cffj(j),cffj(j-1), drdx(i,j,k)

         enddo

       enddo










!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=jmin+1,jmax
        do i=imin,imax


            zlev = 0.5*(z_w(i,j,k)+z_w(i,j-1,k))
            dpth = 0.5*(z_w(i,j,N)+z_w(i,j-1,N)) - zlev


#ifdef PV_CUBIC
            CALL interp_1d(N,rho1(i,j,:),z_r(i,j,:),z_w(i,j,:)
     &                ,zlev,var1,1,0)
            CALL interp_1d(N,rho1(i,j-1,:),z_r(i,j-1,:),z_w(i,j-1,:)
     &                ,zlev,var2,1,0)

            CALL interp_1d(N,qp1(i,j,:),z_r(i,j,:),z_w(i,j,:)
     &                ,zlev,var3,1,0)
            CALL interp_1d(N,qp1(i,j-1,:),z_r(i,j-1,:),z_w(i,j-1,:)
     &                ,zlev,var4,1,0)
#else
            var1 = ((z_r(i,j,k+1) - zlev)*rho1(i,j,k) + 
     &             (zlev - z_r(i,j,k))*rho1(i,j,k+1) )
     &             /(z_r(i,j,k+1) - z_r(i,j,k))

            var2 = ((z_r(i,j-1,k+1) - zlev)*rho1(i,j-1,k) + 
     &             (zlev - z_r(i,j-1,k))*rho1(i,j-1,k+1) )
     &             /(z_r(i,j-1,k+1) - z_r(i,j-1,k))

            var3= ((z_r(i,j,k+1) - zlev)*qp1(i,j,k) + 
     &             (zlev - z_r(i,j,k))*qp1(i,j,k+1) )
     &             /(z_r(i,j,k+1) - z_r(i,j,k))

            var4 = ((z_r(i,j-1,k+1) - zlev)*qp1(i,j-1,k) + 
     &             (zlev - z_r(i,j-1,k))*qp1(i,j-1,k+1) )
     &             /(z_r(i,j-1,k+1) - z_r(i,j-1,k))

#endif



            cffi(i)=-cff*( var1 - var2 ! Elementary
     &                +(var3-var4)         ! adiabatic
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


        do j=jmin+1,jmax
          do i=imin+1,imax


          if ((z_r(i,j,N).eq.0).or.(z_r(i-1,j,N).eq.0)
     &    .or.(z_r(i,j-1,N).eq.0) .or.(z_r(i-1,j-1,N).eq.0)) then

            pv(i,j,:) = FillValue

          else


            do k=N-1,1,-1

            drdz(i,j,k)=0.25*(drdz(i,j,k) + drdz(i-1,j,k)
     &                + drdz(i,j,k)+ drdz(i,j-1,k))






            zlev=0.25*(z_w(i,j,k)+z_w(i-1,j,k)
     &                +z_w(i,j,k)+z_w(i,j-1,k))

!            !!write(*,*) k,z_w(i,j,N),z_w(i,j,k) ,dpth

#ifdef PV_CUBIC
            CALL interp_1d(N,v(i,j,:)
     &                ,0.5*(z_r(i,j,:)+z_r(i,j-1,:))
     &                ,0.5*(z_w(i,j,:)+z_w(i,j-1,:))
     &                ,zlev,var1,1,0)
            CALL interp_1d(N,v(i-1,j,:)
     &                ,0.5*(z_r(i-1,j,:)+z_w(i-1,j-1,:))
     &                ,0.5*(z_w(i-1,j,:)+z_w(i-1,j-1,:))
     &                ,zlev,var2,1,0)
#else


            cff1 = 0.5*(z_r(i,j,k+1)+z_r(i,j-1,k+1)) - zlev
            cff2 = zlev - 0.5*(z_r(i,j,k)+z_r(i,j-1,k))

          var1 = ( cff1*v(i,j,k) + cff2 * v(i,j,k+1) )
     &             /(cff1 + cff2)      

            cff1 = 0.5*(z_r(i-1,j,k+1)+z_r(i-1,j-1,k+1))- zlev
            cff2 = zlev -0.5*(z_r(i-1,j,k)+z_r(i-1,j-1,k))

          var2 = ( cff1*v(i-1,j,k) + cff2 * v(i-1,j,k+1) )
     &             /(cff1 + cff2)    


#endif
  

            dvdx = (var1 - var2) * 0.25 * 
     &           (pm(i,j)+pm(i-1,j)+pm(i,j-1)+pm(i-1,j-1))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


            zlev = 0.5*(z_w(i,j,k)+z_w(i,j-1,k))



#ifdef PV_CUBIC
            CALL interp_1d(N,u(i,j,:),0.5*(z_r(i,j,:)+z_r(i-1,j,:))
     &                ,0.5*(z_w(i,j,:)+z_w(i-1,j,:))
     &                ,zlev,var1,1,0)
            CALL interp_1d(N,u(i,j-1,:),0.5*(z_r(i,j-1,:)+z_r(i,j-1,:))
     &                ,0.5*(z_w(i,j-1,:)+z_w(i-1,j-1,:))
     &                ,zlev,var2,1,0)
#else
            cff1 = 0.5*(z_r(i,j,k+1)+z_r(i-1,j,k+1))- zlev
            cff2 = zlev -0.5*(z_r(i,j,k)+z_r(i-1,j,k))

          var1 = ( cff1*u(i,j,k) + cff2 * u(i,j,k+1) )
     &             /(cff1 + cff2)      

            cff1 = 0.5*(z_r(i,j-1,k+1)+z_r(i-1,j-1,k+1))- zlev
            cff2 = zlev -0.5*(z_r(i,j-1,k)+z_r(i-1,j-1,k))

          var2 = ( cff1*u(i,j-1,k) + cff2 * u(i,j-1,k+1) )
     &             /(cff1 + cff2)    
#endif


            dudy = (var1 - var2)* 0.25 * 
     &           (pn(i,j)+pn(i-1,j)+pn(i,j-1)+pn(i-1,j-1))



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          absvrt = 0.25*(f(i,j)+f(i-1,j)+f(i,j-1)+f(i-1,j-1))
     &                  +dvdx-dudy


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Vertical velocity gradients
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        ! compute dudz on psi - w grid

        dudz = 0.5*((u(i,j,k+1)-u(i,j,k)) / 
     &       (z_r(i,j,k+1)+z_r(i-1,j,k+1)-z_r(i,j,k)-z_r(i-1,j,k)) * 2.
     &              + (u(i,j-1,k+1)-u(i,j-1,k)) / 
     &          (z_r(i,j-1,k+1)+z_r(i-1,j-1,k+1)
     &          -z_r(i,j-1,k)-z_r(i-1,j-1,k)) * 2.)


        ! compute dvdz on psi - w grid

        dvdz = 0.5*((v(i,j,k+1)-v(i,j,k)) /
     &       (z_r(i,j,k+1)+z_r(i,j-1,k+1)-z_r(i,j,k)-z_r(i,j-1,k)) * 2.
     &              + (v(i-1,j,k+1)-v(i-1,j,k)) /
     &          (z_r(i-1,j,k+1)+z_r(i-1,j-1,k+1)
     &          -z_r(i-1,j,k)-z_r(i-1,j-1,k)) * 2.)




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! COMPUTE NEUTRAL DENSITY GRADIENTS 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            pv(i,j,k) = absvrt * drdz(i,j,k)
     &          - dvdz * drdx(i,j,k) 
     &          + dudz * drdy(i,j,k)


            enddo !<-- k




          endif !<-- z_r.eq.0 (land)
        enddo !<- i
      enddo !<- j





      return
      end







 
 
