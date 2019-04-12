


# undef UPSTREAM


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute Advective part of the barotropic vorticity balance equation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine get_adv_4th (Lm,Mm,N,u,v, z_r,z_w,pm,pn
     &  ,adv)    

      implicit none

      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k,
     &      istr,iend,jstr,jend,istrU,jstrV


      real*8 u(1:Lm+1,0:Mm+1,N),  v(0:Lm+1,1:Mm+1,N),
     &      FlxU(1:Lm+1,0:Mm+1,N),  FlxV(0:Lm+1,1:Mm+1,N),
     &      z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N),
     &      Hz(0:Lm+1,0:Mm+1,N), gamma,
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      dn_u(0:Lm+1,0:Mm+1), dm_v(0:Lm+1,0:Mm+1),
     &      var1, var2,var3, var4, cff, cff1, cff2


      real*8 wrkXadv(0:Lm+1,0:Mm+1,2),wrkYadv(0:Lm+1,0:Mm+1,2)

      real*8 wrk1(0:Lm+1,0:Mm+1), wrk2(0:Lm+1,0:Mm+1),
     &      UFx(0:Lm+1,0:Mm+1), UFe(0:Lm+1,0:Mm+1),
     &      VFx(0:Lm+1,0:Mm+1), VFe(0:Lm+1,0:Mm+1)

      real*8 adv(1:Lm+1,1:Mm+1), vrtXadv(1:Lm+1,1:Mm+1), 
     &      vrtYadv(1:Lm+1,1:Mm+1)

# include "scalars.h"

      parameter (gamma=0.25)  

Cf2py intent(in) Lm,Mm,N, u,v,z_r,z_w,pm,pn
Cf2py intent(out) adv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      do j=jmin,jmax
        do i=imin,imax
          do k=1,N,+1
           Hz(i,j,k) = z_w(i,j,k) - z_w(i,j,k-1)
          enddo
        enddo
      enddo



      do j=jmin,jmax
        do i=imin+1,imax
            dn_u(i,j) = 2./(pn(i,j)+pn(i-1,j))
            do k=1,N,+1
              FlxU(i,j,k) = 0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)
     &      * u(i,j,k) 
            enddo
          enddo
      enddo



      do j=jmin+1,jmax
        do i=imin,imax
            dm_v(i,j) = 2./(pm(i,j)+pm(i,j-1))
            do k=1,N,+1
              FlxV(i,j,k) = 0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)
     &      * v(i,j,k) 
            enddo
          enddo
      enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





      do k=1,N

!
! Add in horizontal advection of momentum: Compute diagonal [UFx,VFe]
! and off-diagonal [UFe,VFx] components of tensor of momentum flux
! due to horizontal advection; after that add divergence of these
! terms to r.h.s.
!
#  define uxx wrk1
#  define Huxx wrk2



        do j=jmin,jmax
         do i=imin+2,imax-1

            uxx(i,j)=u(i-1,j,k)-2.*u(i,j,k)
     &                                +u(i+1,j,k)
            Huxx(i,j)=FlxU(i-1,j,k)-2.*FlxU(i,j,k)
     &                                  +FlxU(i+1,j,k)
          enddo
        enddo





        do j=jmin,jmax
         do i=imin+2,imax-1
#  ifdef UPSTREAM
            cff=FlxU(i,j,k)+FlxU(i+1,j,k)-0.125*( Huxx(i  ,j)
     &                                           +Huxx(i+1,j))
            UFx(i,j)=0.25*( cff*(u(i,j,k)+u(i+1,j,k))
     &                          -gamma*( max(cff,0.)*uxx(i  ,j)
     &                                  +min(cff,0.)*uxx(i+1,j)
     &                                                      ))
#  else
            UFx(i,j)=0.25*( u(i,j,k)+u(i+1,j,k)
     &                         -0.125*(uxx(i,j)+uxx(i+1,j))
     &                  )*( FlxU(i,j,k)+FlxU(i+1,j,k)
     &                      -0.125*(Huxx(i,j)+Huxx(i+1,j)))

#  endif
          enddo
        enddo



#  undef Huxx
#  undef uxx




#  define vee wrk1
#  define Hvee wrk2

        do j=jmin+2,jmax-1
         do i=imin,imax

            vee(i,j)=v(i,j-1,k)-2.*v(i,j,k)+v(i,j+1,k)

            Hvee(i,j)=FlxV(i,j-1,k)-2.*FlxV(i,j,k)+FlxV(i,j+1,k)

          enddo
        enddo


        do j=jmin+2,jmax-1
         do i=imin,imax

#  ifdef UPSTREAM
            cff=FlxV(i,j,k)+FlxV(i,j+1,k)-0.125*( Hvee(i,j  )
     &                                           +Hvee(i,j+1))
            VFe(i,j)=0.25*( cff*(v(i,j,k)+v(i,j+1,k))
     &                          -gamma*( max(cff,0.)*vee(i,j  )
     &                                  +min(cff,0.)*vee(i,j+1)
     &                                                      )) 

#  else
            VFe(i,j)=0.25*( v(i,j,k)+v(i,j+1,k)
     &                        -0.125*(vee(i,j)+vee(i,j+1))
     &                  )*( FlxV(i,j,k)+FlxV(i,j+1,k)
     &                      -0.125*(Hvee(i,j)+Hvee(i,j+1)))
#  endif

          enddo
        enddo
#  undef Hvee
#  undef vee





#  define uee wrk1
#  define Hvxx wrk2

       do j=jmin+1,jmax-1
         do i=imin,imax
            uee(i,j)=u(i,j-1,k)-2.*u(i,j,k)
     &                                 +u(i,j+1,k)
          enddo
        enddo

        do j=jmin,jmax
         do i=imin+1,imax-1
           Hvxx(i,j)=FlxV(i-1,j,k)-2.*FlxV(i,j,k)+FlxV(i+1,j,k)
          enddo
        enddo

        do j=jmin+1,jmax-1
         do i=imin+1,imax-1
#  ifdef UPSTREAM
            cff=FlxV(i,j,k)+FlxV(i-1,j,k)-0.125*( Hvxx(i  ,j)
     &                                           +Hvxx(i-1,j))
            UFe(i,j)=0.25*( cff*(u(i,j,k)+u(i,j-1,k))
     &                          -gamma*( max(cff,0.)*uee(i,j-1)
     &                                  +min(cff,0.)*uee(i,j  )
     &                                                      ))
#  else
            UFe(i,j)=0.25*( u(i,j,k)+u(i,j-1,k)
     &                        -0.125*(uee(i,j)+uee(i,j-1))
     &                  )*( FlxV(i,j,k)+FlxV(i-1,j,k)
     &                     -0.125*(Hvxx(i,j)+Hvxx(i-1,j)))
#  endif
          enddo
        enddo
#  undef Hvxx
#  undef uee


#  define vxx wrk1
#  define Huee wrk2

        do j=jmin,jmax
         do i=imin+1,imax-1

            vxx(i,j)=v(i-1,j,k)-2.*v(i,j,k)
     &                                 +v(i+1,j,k)
          enddo
        enddo

        do j=jmin+1,jmax-1
         do i=imin,imax

           Huee(i,j)=FlxU(i,j-1,k)-2.*FlxU(i,j,k)+FlxU(i,j+1,k)
          enddo
        enddo

        do j=jmin+1,jmax-1
         do i=imin+1,imax-1

#  ifdef UPSTREAM
            cff=FlxU(i,j,k)+FlxU(i,j-1,k)-0.125*( Huee(i,j  )
     &                                           +Huee(i,j-1))
            VFx(i,j)=0.25*( cff*(v(i,j,k)+v(i-1,j,k))
     &                          -gamma*( max(cff,0.)*vxx(i-1,j)
     &                                  +min(cff,0.)*vxx(i  ,j)
     &                                                      ))
#  else
            VFx(i,j)=0.25*( v(i,j,k)+v(i-1,j,k)
     &                        -0.125*(vxx(i,j)+vxx(i-1,j))
     &                  )*( FlxU(i,j,k)+FlxU(i,j-1,k)
     &                     -0.125*(Huee(i,j)+Huee(i,j-1)))
#  endif
          enddo
        enddo
#  undef Huee
#  undef vxx


      do j=jmin+2,jmax-1
        do i=imin+2,imax-1

              if (k.eq.1) then
                wrkXadv(i,j,1) = -UFx(i,j)+UFx(i-1,j)
                wrkYadv(i,j,1) = -UFe(i,j+1)+UFe(i,j)
              else
                wrkXadv(i,j,1) = wrkXadv(i,j,1) - UFx(i,j)+UFx(i-1,j)
                wrkYadv(i,j,1) = wrkYadv(i,j,1) - UFe(i,j+1)+UFe(i,j)
              endif


          enddo
        enddo

      do j=jmin+2,jmax-1
        do i=imin+2,imax-1

              if (k.eq.1) then
                wrkXadv(i,j,2) = -VFx(i+1,j)+VFx(i,j)
                wrkYadv(i,j,2) = -VFe(i,j)+VFe(i,j-1)
              else
                wrkXadv(i,j,2) = wrkXadv(i,j,2) -VFx(i+1,j)+VFx(i,j)
                wrkYadv(i,j,2) = wrkYadv(i,j,2) -VFe(i,j)+VFe(i,j-1)
              endif

          enddo
        enddo
      enddo


! 
! ! Divide all diagnostic terms by (pm*pn).
! ! There after the unit of these terms are :
! !  s-2
! 
! 
!       do j=jmin+2,jmax-1
!         do i=imin+2,imax-1
! 
!             cff=0.25*(pm(i,j)+pm(i-1,j))
!      &             *(pn(i,j)+pn(i-1,j))
! 
!             
!             wrkXadv(i,j,1)=wrkXadv(i,j,1)*cff
!             wrkYadv(i,j,1)=wrkYadv(i,j,1)*cff
! 
! 
!             cff=0.25*(pm(i,j)+pm(i,j-1))
!      &             *(pn(i,j)+pn(i,j-1))
! 
! 
!             wrkXadv(i,j,2)=wrkXadv(i,j,2)*cff
!             wrkYadv(i,j,2)=wrkYadv(i,j,2)*cff
! 
!           enddo
!         enddo
! 
! 
! 
!       do j=jmin+2,jmax-2
!         do i=imin+2,imax-2
! 
!            cff1 = 0.25*(pm(i,j) + pm(i-1,j) + pm(i,j-1) + pm(i-1,j-1))
!            cff2 = 0.25*(pn(i,j) + pn(i-1,j) + pn(i,j-1) + pn(i-1,j-1))
! 
!             vrtXadv(i,j) = (wrkXadv(i,j,2) - wrkXadv(i-1,j,2)) * cff1
!      &                     - (wrkXadv(i,j,1) - wrkXadv(i,j-1,1)) * cff2
! 
! 
!             vrtYadv(i,j) = (wrkYadv(i,j,2) - wrkYadv(i-1,j,2)) * cff1
!      &                     - (wrkYadv(i,j,1) - wrkYadv(i,j-1,1)) * cff2
! 
! 
!             adv(i,j) = vrtXadv(i,j) + vrtYadv(i,j)
! 
! 
!           enddo
!          enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! NEW COMPUTATION OF ROT (added 14/09/13)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=jmin+2,jmax-1
        do i=imin+2,imax-1

            cff=0.5*(pn(i,j)+pn(i-1,j))
        
            wrkXadv(i,j,1)=wrkXadv(i,j,1)*cff
            wrkYadv(i,j,1)=wrkYadv(i,j,1)*cff

            cff=0.5*(pm(i,j)+pm(i,j-1))

            wrkXadv(i,j,2)=wrkXadv(i,j,2)*cff
            wrkYadv(i,j,2)=wrkYadv(i,j,2)*cff

          enddo
        enddo



      do j=jmin+2,jmax-2
        do i=imin+2,imax-2

           cff1 = 0.25*(pm(i,j) + pm(i-1,j) + pm(i,j-1) + pm(i-1,j-1))
     &          * 0.25*(pn(i,j) + pn(i-1,j) + pn(i,j-1) + pn(i-1,j-1))

            vrtXadv(i,j) = (wrkXadv(i,j,2) - wrkXadv(i-1,j,2)) * cff1
     &                     - (wrkXadv(i,j,1) - wrkXadv(i,j-1,1)) * cff1


            vrtYadv(i,j) = (wrkYadv(i,j,2) - wrkYadv(i-1,j,2)) * cff1
     &                     - (wrkYadv(i,j,1) - wrkYadv(i,j-1,1)) * cff1


            adv(i,j) = vrtXadv(i,j) + vrtYadv(i,j)


          enddo
         enddo




      return
      end











