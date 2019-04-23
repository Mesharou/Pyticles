


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute Advective part of the barotropic vorticity balance equation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine get_kediss_2d (Lm,Mm,u,v, pm,pn
     &  ,kediss)    

      implicit none

      integer Lm,Mm, imin,imax,jmin,jmax, i,j,k,
     &      istr,iend,jstr,jend,istrU,jstrV


      real*8 u(1:Lm+1,0:Mm+1),  v(0:Lm+1,1:Mm+1),
     &      FlxU(1:Lm+1,0:Mm+1),  FlxV(0:Lm+1,1:Mm+1),
     &      Hz(0:Lm+1,0:Mm+1), gamma,
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      dn_u(0:Lm+1,0:Mm+1), dm_v(0:Lm+1,0:Mm+1),
     &      var1, var2,var3, var4, cff, cff1, cff2



      real*8 wrk1(0:Lm+1,0:Mm+1), wrk2(0:Lm+1,0:Mm+1),
     &      UFx(0:Lm+1,0:Mm+1), UFe(0:Lm+1,0:Mm+1),
     &      VFx(0:Lm+1,0:Mm+1), VFe(0:Lm+1,0:Mm+1)

      real*8 advu(1:Lm+1,0:Mm+1), advv(0:Lm+1,1:Mm+1)

      real*8 kediss(0:Lm+1,0:Mm+1)

# include "scalars.h"

      parameter (gamma=0.25)  

Cf2py intent(in) Lm,Mm,N, u,v,pm,pn
Cf2py intent(out) kediss

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      do j=jmin,jmax
        do i=imin,imax
           Hz(i,j) = 1.
        enddo
      enddo



      do j=jmin,jmax
        do i=imin+1,imax
            dn_u(i,j) = 2./(pn(i,j)+pn(i-1,j))
              FlxU(i,j) = 0.5*(Hz(i,j)+Hz(i-1,j))*dn_u(i,j)
     &      * u(i,j) 
          enddo
      enddo



      do j=jmin+1,jmax
        do i=imin,imax
            dm_v(i,j) = 2./(pm(i,j)+pm(i,j-1))
              FlxV(i,j) = 0.5*(Hz(i,j)+Hz(i,j-1))*dm_v(i,j)
     &      * v(i,j) 
          enddo
      enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!
! Add in horizontal advection of momentum: Compute diagonal [UFx,VFe]
! and off-diagonal [UFe,VFx] components of tensor of momentum flux
! due to horizontal advection; after that add divergence of these
! terms to r.h.s.
!
#  define uxx wrk1
#  define Huxx wrk2



        do j=jmin,jmax
         do i=imin+1,imax-2

            uxx(i,j)=u(i-1,j)-2.*u(i,j)
     &                                +u(i+1,j)
            Huxx(i,j)=FlxU(i-1,j)-2.*FlxU(i,j)
     &                                  +FlxU(i+1,j)
          enddo
        enddo





        do j=jmin,jmax
         do i=imin+2,imax-1
            cff=FlxU(i,j)+FlxU(i+1,j)-0.125*( Huxx(i  ,j)
     &                                           +Huxx(i+1,j))
         UFx(i,j)=0.25*( cff*(u(i,j)+u(i+1,j))
     &                          -gamma*( max(cff,0.)*uxx(i  ,j)
     &                                  +min(cff,0.)*uxx(i+1,j)
     &                                                      ))

        UFx(i,j) = UFx(i,j) - 0.25*( u(i,j)+u(i+1,j)
     &                         -0.125*(uxx(i,j)+uxx(i+1,j))
     &                  )*( FlxU(i,j)+FlxU(i+1,j)
     &                      -0.125*(Huxx(i,j)+Huxx(i+1,j)))

          enddo
        enddo



#  undef Huxx
#  undef uxx




#  define vee wrk1
#  define Hvee wrk2

        do j=jmin+1,jmax-2
         do i=imin,imax

            vee(i,j)=v(i,j-1)-2.*v(i,j)+v(i,j+1)

            Hvee(i,j)=FlxV(i,j-1)-2.*FlxV(i,j)+FlxV(i,j+1)

          enddo
        enddo


        do j=jmin+2,jmax-1
         do i=imin,imax


            cff=FlxV(i,j)+FlxV(i,j+1)-0.125*( Hvee(i,j  )
     &                                           +Hvee(i,j+1))

            VFe(i,j)=0.25*( cff*(v(i,j)+v(i,j+1))
     &                          -gamma*( max(cff,0.)*vee(i,j  )
     &                                  +min(cff,0.)*vee(i,j+1)
     &                                                      )) 


        VFe(i,j)=VFe(i,j) - 0.25*( v(i,j)+v(i,j+1)
     &                        -0.125*(vee(i,j)+vee(i,j+1))
     &                  )*( FlxV(i,j)+FlxV(i,j+1)
     &                      -0.125*(Hvee(i,j)+Hvee(i,j+1)))


          enddo
        enddo
#  undef Hvee
#  undef vee





#  define uee wrk1
#  define Hvxx wrk2

       do j=jmin+1,jmax-1
         do i=imin,imax
            uee(i,j)=u(i,j-1)-2.*u(i,j)
     &                                 +u(i,j+1)
          enddo
        enddo

        do j=jmin,jmax
         do i=imin+1,imax-1
           Hvxx(i,j)=FlxV(i-1,j)-2.*FlxV(i,j)+FlxV(i+1,j)
          enddo
        enddo

        do j=jmin+1,jmax-1
         do i=imin+1,imax-1

            cff=FlxV(i,j)+FlxV(i-1,j)-0.125*( Hvxx(i  ,j)
     &                                           +Hvxx(i-1,j))
            UFe(i,j)=0.25*( cff*(u(i,j)+u(i,j-1))
     &                          -gamma*( max(cff,0.)*uee(i,j-1)
     &                                  +min(cff,0.)*uee(i,j  )
     &                                                      ))

            UFe(i,j)=UFe(i,j) - 0.25*( u(i,j)+u(i,j-1)
     &                        -0.125*(uee(i,j)+uee(i,j-1))
     &                  )*( FlxV(i,j)+FlxV(i-1,j)
     &                     -0.125*(Hvxx(i,j)+Hvxx(i-1,j)))

          enddo
        enddo
#  undef Hvxx
#  undef uee


#  define vxx wrk1
#  define Huee wrk2

        do j=jmin,jmax
         do i=imin+1,imax-1

            vxx(i,j)=v(i-1,j)-2.*v(i,j)
     &                                 +v(i+1,j)
          enddo
        enddo

        do j=jmin+1,jmax-1
         do i=imin,imax

           Huee(i,j)=FlxU(i,j-1)-2.*FlxU(i,j)+FlxU(i,j+1)
          enddo
        enddo

        do j=jmin+1,jmax-1
         do i=imin+1,imax-1


            cff=FlxU(i,j)+FlxU(i,j-1)-0.125*( Huee(i,j  )
     &                                           +Huee(i,j-1))
            VFx(i,j)=0.25*( cff*(v(i,j)+v(i-1,j))
     &                          -gamma*( max(cff,0.)*vxx(i-1,j)
     &                                  +min(cff,0.)*vxx(i  ,j)
     &                                                      ))

            VFx(i,j)=VFx(i,j) - 0.25*( v(i,j)+v(i-1,j)
     &                        -0.125*(vxx(i,j)+vxx(i-1,j))
     &                  )*( FlxU(i,j)+FlxU(i,j-1)
     &                     -0.125*(Huee(i,j)+Huee(i,j-1)))

          enddo
        enddo
#  undef Huee
#  undef vxx










      do j=jmin+2,jmax-1
        do i=imin+2,imax-1

            advu(i,j) = (-UFx(i,j)+UFx(i-1,j)
     &              -UFe(i,j+1)+UFe(i,j))



          enddo
        enddo

      do j=jmin+2,jmax-1
        do i=imin+2,imax-1


            advv(i,j) = (-VFx(i+1,j)+VFx(i,j)
     &              -VFe(i,j)+VFe(i,j-1))


          enddo
        enddo





! Divide all diagnostic terms by (pm*pn).
! There after the unit of these terms are :
!  s-2


      do j=jmin+2,jmax-1
        do i=imin+2,imax-1


            cff=0.25*(pm(i,j)+pm(i-1,j))
     &             *(pn(i,j)+pn(i-1,j))


              advu(i,j)=advu(i,j)*cff
     &             * u(i,j)



            cff=0.25*(pm(i,j)+pm(i,j-1))
     &             *(pn(i,j)+pn(i,j-1))



              advv(i,j)=advv(i,j)*cff
     &             * v(i,j)


          enddo
        enddo



!!!!!!!!!!!!!!!!!


      do j=jmin+2,jmax-1
        do i=imin+2,imax-1


              kediss(i,j) = 0.5*( advu(i-1,j) + advu(i,j)
     &                            + advv(i,j-1) + advv(i,j) )


          enddo
        enddo


!!!!!!!!!!!!!!!!!!

      return
      end











