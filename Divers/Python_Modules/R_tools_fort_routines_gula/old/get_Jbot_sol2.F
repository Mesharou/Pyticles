
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute PV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_Jbot_sol2(Lm,Mm,N, T,S, u,v, z_r,z_w,rho0,pm,pn,
     & hbbls,rdrg,Jbot)


      implicit none
      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k
      real*8 T(0:Lm+1,0:Mm+1,N),  S(0:Lm+1,0:Mm+1,N),
     &      u(1:Lm+1,0:Mm+1), v(0:Lm+1,1:Mm+1),
     &      Jbot(1:Lm+1,1:Mm+1),
     &      dx(1:Lm+1,1:Mm+1),
     &      dy(1:Lm+1,1:Mm+1),
     &      dz(1:Lm+1,1:Mm+1),
     &      rd(0:Lm+1,0:Mm+1),
     &      ubot(1:Lm+1,1:Mm+1),vbot(1:Lm+1,1:Mm+1),
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      hbbls(0:Lm+1,0:Mm+1),
     &      z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N),
     &      Hz(0:Lm+1,0:Mm+1),
     &      rho1(0:Lm+1,0:Mm+1,N), drdz(0:Lm+1,0:Mm+1),
     &      drdx(1:Lm+1,0:Mm+1), drdy(0:Lm+1,1:Mm+1),
     &      qp1(0:Lm+1,0:Mm+1,N),
     &      dpth,cff,Tt,Ts,sqrtTs, rho0,  K0, dr00,
     &      cff2,cff3,zlev,
     &      dvdx, dudy,
     &      cffi(0:Lm+1), cffj(0:Mm+1),
     &      var1, var2,var3, var4


      real*8 Zob, rdrg

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
     &          qp2=0.0000172

#  ifdef DUKO_2001
      real rho1_0, K0_Duk
#  endif

# include "scalars.h"

Cf2py intent(in) Lm,Mm,N, T,S, u,v,z_r,z_w,rho0,pm,pn,hbbls,rdrg
Cf2py intent(out) Jbot









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
      enddo    ! <-- j



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! COMPUTE BUOYANCY GRADIENTS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        cff=g/rho0



!---------------------------------------------------------------------------------------




      do j=jmin,jmax
          do i=imin,imax

            dpth=z_w(i,j,N)-z_r(i,j,1)


            CALL interp_1d(N,rho1(i,j,:),z_r(i,j,:),z_w(i,j,:)
     &                ,z_w(i,j,0),var1,1,0)
            CALL interp_1d(N,rho1(i,j,:),z_r(i,j,:),z_w(i,j,:)
     &                ,z_w(i,j,1),var2,1,0)

            CALL interp_1d(N,qp1(i,j,:),z_r(i,j,:),z_w(i,j,:)
     &                ,z_w(i,j,0),var3,1,0)
            CALL interp_1d(N,qp1(i,j,:),z_r(i,j,:),z_w(i,j,:)
     &                ,z_w(i,j,1),var4,1,0)

            cff2=( var2-var1             ! Elementary
     &                        +(var4-var3) ! adiabatic 
     &                            *dpth*(1.-2.*qp2*dpth) ! difference
     &                       )


            drdz(i,j) =- cff*cff2 / (z_w(i,j,1)-z_w(i,j,0))



          enddo
      enddo    ! <-- j






!---------------------------------------------------------------------------------------

       do i=imin+1,imax
        do j=jmin,jmax

            zlev=0.5*(z_r(i,j,1)+z_r(i-1,j,1))

            dpth=0.5*(z_w(i,j,N)+z_w(i-1,j,N)) - zlev


            !if ((z_r(i,j,N).ne.0).and.(z_r(i-1,j,N).ne.0)) then
            !if (dpth.gt.z_w(i,j,N)) then
             !write(*,*) 'calling', dpth, z_r(i,j,N), z_w(i,j,N)
             !write(*,*) z_w(i,j,:)
             !write(*,*) z_w(i-1,j,:)
            !endif
            !endif

            !if (dpth.gt.z_w(i-1,j,N)) then
            ! write(*,*) 'calling 2',  dpth, z_r(i,j,N), z_w(i,j,N)
            !endif


            !if ((z_r(i,j,N).ne.0).and.(z_r(i-1,j,N).ne.0)) then


            CALL interp_1d(N,rho1(i,j,:),z_r(i,j,:),z_w(i,j,:)
     &                ,zlev,var1,1,0)
            CALL interp_1d(N,rho1(i-1,j,:),z_r(i-1,j,:),z_w(i-1,j,:)
     &                ,zlev,var2,1,0)


            CALL interp_1d(N,qp1(i,j,:),z_r(i,j,:),z_w(i,j,:)
     &                ,zlev,var3,1,0)
            CALL interp_1d(N,qp1(i-1,j,:),z_r(i-1,j,:),z_w(i-1,j,:)
     &                ,zlev,var4,1,0)

            cffj(j)=-cff*( var1 - var2 ! Elementary
     &                +(var3-var4)         ! adiabatic
     &                     *dpth*(1.-qp2*dpth) )         ! difference
     &                    *0.5*(pm(i,j)+pm(i-1,j))

            !else

            !cffj(j)= 9999.

            !endif

         enddo



         do j=jmin+1,jmax

            drdx(i,j)= 0.5 * (cffj(j) + cffj(j-1))
            !!write(*,*) 'all3',dpth,cffj(j),cffj(j-1), drdx(i,j,k)

         enddo

         !!write(*,*) 'all2',i,j,drdx(i,j)


       enddo


!---------------------------------------------------------------------------------------

        do j=jmin+1,jmax

          do i=imin,imax

            zlev=0.5*(z_r(i,j,1)+z_r(i,j-1,1))
            dpth=0.5*(z_w(i,j,N)+z_w(i-1,j,N)) - zlev

            !write(*,*) 'all2',z_w(i,j,N),z_r(i,j,1),i,j

            CALL interp_1d(N,rho1(i,j,:),z_r(i,j,:),z_w(i,j,:)
     &                ,zlev,var1,1,0)
            CALL interp_1d(N,rho1(i,j-1,:),z_r(i,j-1,:),z_w(i,j-1,:)
     &                ,zlev,var2,1,0)

            CALL interp_1d(N,qp1(i,j,:),z_r(i,j,:),z_w(i,j,:)
     &                ,zlev,var3,1,0)
            CALL interp_1d(N,qp1(i,j-1,:),z_r(i,j-1,:),z_w(i,j-1,:)
     &                ,zlev,var4,1,0)




            cffi(i)=-cff*( var1 - var2 ! Elementary
     &                +(var3-var4)         ! adiabatic
     &                     *dpth*(1.-qp2*dpth) )  
     &                    *0.5*(pn(i,j)+pn(i,j-1)) 


         !!!write(*,*) 'all2',i,j,dpth,cffi(i),var1,var2

          enddo



          do i=imin+1,imax

            drdy(i,j)= 0.5 * (cffi(i) + cffi(i-1))


          enddo



        enddo




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bottom Drag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       Zob=0.01


       do j=jmin+1,jmax-1
         do i=imin+1,imax-1

            Hz(i,j) = z_w(i,j,1) - z_w(i,j,0)

# ifdef LINEAR_DRAG_ONLY
            rd(i,j)=rdrg
# else
            cff=sqrt( 0.333333333333*(
     &              u(i,j)**2 +u(i+1,j)**2
     &                    +u(i,j)*u(i+1,j)
     &              +v(i,j)**2+v(i,j+1)**2
     &                    +v(i,j)*v(i,j+1)
     &                                               ))
            rd(i,j)=rdrg + cff*(vonKar/log(Hz(i,j)/Zob))**2
# endif



          enddo
        enddo

       do j=jmin+1,jmax
            rd(imax ,j)=rd(imax-1 ,j)
            rd(imin ,j)=rd(imin+1 ,j)
        enddo

       do i=imin+1,imax
            rd(i ,jmax)=rd(i ,jmax-1)
            rd(i ,jmin)=rd(i ,jmin+1)
        enddo







!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! HEIGHT GRADIENT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        do j=jmin+1,jmax
          do i=imin+1,imax

            dx(i,j) = 0.5* (z_r(i,j,1) + z_r(i,j-1,1)
     &               - z_r(i-1,j,1)- z_r(i-1,j-1,1))
     &      * 0.25*(pm(i,j)+pm(i,j-1)+ pm(i+1,j)+pm(i-1,j-1))

            dy(i,j)  = 0.5* (z_r(i,j,1) + z_r(i-1,j,1)
     &               - z_r(i,j-1,1)- z_r(i-1,j-1,1))
     &      * 0.25*(pn(i,j)+pn(i,j-1)+ pn(i+1,j)+pn(i-1,j-1))
 

            cff = 1./sqrt(1 + dx(i,j)**2 + dy(i,j)**2)

            dx(i,j) = cff * dx(i,j)
            dy(i,j) = cff * dy(i,j) 
            dz(i,j) = cff

            !if (cff.ne.1) then
            !    write(*,*) i,j,dx(i,j), dy(i,j), dz(i,j)
            !endif

          enddo
        enddo




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! COMPUTE 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



        do j=jmin+1,jmax
          do i=imin+1,imax


            cff2 = 0.25 * (drdz(i-1,j) + drdz(i,j) 
     &               + drdz(i-1,j-1) + drdz(i,j-1))
            cff3 =0.25 * (rd(i,j)+rd(i-1,j)
     &               + rd(i,j-1)+rd(i-1,j-1))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            cff = 1./(0.25*(hbbls(i,j) + hbbls(i-1,j)
     &                 + hbbls(i,j-1) + hbbls(i-1,j-1)))

            ubot(i,j) = cff3*0.5*(u(i,j-1) + u(i,j)) 

            vbot(i,j) = cff3*0.5*(v(i-1,j)+v(i,j) )

            Jbot(i,j) = cff * ((drdx(i,j) * vbot(i,j)
     &              - drdy(i,j) * ubot(i,j)) * dz(i,j)
     &              - (cff2 * vbot(i,j)) * dx(i,j)
     &              + (cff2 * ubot(i,j)) * dy(i,j)  
     &               )



          enddo
        enddo






      return
      end







 
 
