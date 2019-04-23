


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! compute Advective part of the barotropic vorticity balance equation
!!
!! Boundary conditions addded on 05/09/14
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine get_kediss (Lm,Mm,N,u,v, z_r,z_w,pm,pn, rmask
     &  ,kediss)    

      implicit none

      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k


      real*8 u(1:Lm+1,0:Mm+1,N),  v(0:Lm+1,1:Mm+1,N),
     &      FlxU(1:Lm+1,0:Mm+1,N),  FlxV(0:Lm+1,1:Mm+1,N),
     &      z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N),
     &      rmask(0:Lm+1,0:Mm+1),
     &      Hz(0:Lm+1,0:Mm+1,N), gamma,
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      dn_u(0:Lm+1,0:Mm+1), dm_v(0:Lm+1,0:Mm+1),
     &      var1, var2,var3, var4, cff, cff1, cff2



      real*8 wrk1(0:Lm+1,0:Mm+1), wrk2(0:Lm+1,0:Mm+1),
     &      UFx(0:Lm+1,0:Mm+1), UFe(0:Lm+1,0:Mm+1),
     &      VFx(0:Lm+1,0:Mm+1), VFe(0:Lm+1,0:Mm+1),
     &      umask(1:Lm+1,0:Mm+1),  vmask(0:Lm+1,1:Mm+1)

      real*8 advu(1:Lm+1,0:Mm+1,N), advv(0:Lm+1,1:Mm+1,N)

      real*8 kediss(0:Lm+1,0:Mm+1,N)

# include "scalars.h"

      parameter (gamma=0.25)  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!# include "compute_tile_bounds.h"
# include "compute_auxiliary_bounds.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Cf2py intent(in) Lm,Mm,N, u,v,z_r,z_w,pm,pn, rmask
Cf2py intent(out) kediss


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Compute grid variables and fluxes


      do j=jstrR,jendR
        do i=istrR,iendR
          do k=1,N,+1
           Hz(i,j,k) = z_w(i,j,k) - z_w(i,j,k-1)
          enddo
        enddo
      enddo


      do j=jstrR,jendR
        do i=istr,iendR
            dn_u(i,j) = 2./(pn(i,j)+pn(i-1,j))
            umask(i,j) = rmask(i,j)*rmask(i-1,j)
            do k=1,N,+1
              u(i,j,k) = u(i,j,k) * umask(i,j)
              FlxU(i,j,k) = 0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)
     &      * u(i,j,k) 
            enddo
          enddo
      enddo


      do j=jstr,jendR
        do i=istrR,iendR
            dm_v(i,j) = 2./(pm(i,j)+pm(i,j-1))
            vmask(i,j) =  rmask(i,j)*rmask(i,j-1)
            do k=1,N,+1
              v(i,j,k) = v(i,j,k) * vmask(i,j)
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
#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then        ! Sort out bounding indices of
          imin=istrU                  ! extended ranges: note that in
        else                          ! the vicinity of physical
          imin=istrU-1                ! boundaries values at the
        endif                         ! extremal points of stencil
        if (EASTERN_EDGE) then        ! are not available, so an
          imax=iend                   ! extrapolation rule needs to
        else                          ! be applied. Also note that
          imax=iend+1                 ! for this purpose periodic
        endif                         ! ghost points and MPI margins 
#  else
        imin=istrU-1                  ! are not considered as
        imax=iend+1                   ! physical boundaries.
#  endif

        do j=jstr,jend
          do i=imin,imax
            uxx(i,j)=u(i-1,j,k)-2.*u(i,j,k)
     &                                +u(i+1,j,k)
            Huxx(i,j)=FlxU(i-1,j,k)-2.*FlxU(i,j,k)
     &                                  +FlxU(i+1,j,k)
          enddo
        enddo
#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=jstr,jend
            uxx(istrU-1,j) =uxx(istrU,j)
            Huxx(istrU-1,j)=Huxx(istrU,j)
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=jstr,jend
            uxx(iend+1,j) =uxx(iend,j)
            Huxx(iend+1,j)=Huxx(iend,j)
          enddo
        endif
#  endif



        do j=jstr,jend
          do i=istrU-1,iend

         cff=FlxU(i,j,k)+FlxU(i+1,j,k)-0.125*( Huxx(i  ,j)
     &                                           +Huxx(i+1,j))
         UFx(i,j)=0.25*( cff*(u(i,j,k)+u(i+1,j,k))
     &                          -gamma*( max(cff,0.)*uxx(i  ,j)
     &                                  +min(cff,0.)*uxx(i+1,j)
     &                                                      ))

         UFx(i,j) = UFx(i,j) - 0.25*( u(i,j,k)+u(i+1,j,k)
     &                         -0.125*(uxx(i,j)+uxx(i+1,j))
     &                  )*( FlxU(i,j,k)+FlxU(i+1,j,k)
     &                      -0.125*(Huxx(i,j)+Huxx(i+1,j)))

          enddo
        enddo
#  undef Huxx
#  undef uxx




#  define vee wrk1
#  define Hvee wrk2
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          jmin=jstrV
        else
          jmin=jstrV-1
        endif
        if (NORTHERN_EDGE) then
          jmax=jend
        else
          jmax=jend+1 
        endif
#  else
        jmin=jstrV-1
        jmax=jend+1
#  endif
        do j=jmin,jmax
          do i=istr,iend
            vee(i,j)=v(i,j-1,k)-2.*v(i,j,k)+v(i,j+1,k)
            Hvee(i,j)=FlxV(i,j-1,k)-2.*FlxV(i,j,k)+FlxV(i,j+1,k)
          enddo
        enddo
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=istr,iend
            vee(i,jstrV-1)=vee(i,jstrV)
            Hvee(i,jstrV-1)=Hvee(i,jstrV)
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=istr,iend
            vee(i,jend+1)=vee(i,jend)
            Hvee(i,jend+1)=Hvee(i,jend)
          enddo
        endif
#  endif



        do j=jstrV-1,jend
          do i=istr,iend

            cff=FlxV(i,j,k)+FlxV(i,j+1,k)-0.125*( Hvee(i,j  )
     &                                           +Hvee(i,j+1))

            VFe(i,j)=0.25*( cff*(v(i,j,k)+v(i,j+1,k))
     &                          -gamma*( max(cff,0.)*vee(i,j  )
     &                                  +min(cff,0.)*vee(i,j+1)
     &                                                      )) 

            VFe(i,j)=VFe(i,j) - 0.25*( v(i,j,k)+v(i,j+1,k)
     &                        -0.125*(vee(i,j)+vee(i,j+1))
     &                  )*( FlxV(i,j,k)+FlxV(i,j+1,k)
     &                      -0.125*(Hvee(i,j)+Hvee(i,j+1)))

          enddo
        enddo
#  undef Hvee
#  undef vee





#  define uee wrk1
#  define Hvxx wrk2
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          jmin=jstr
        else
          jmin=jstr-1
        endif
        if (NORTHERN_EDGE) then
          jmax=jend
        else
          jmax=jend+1
        endif
#  else
        jmin=jstr-1
        jmax=jend+1
#  endif
        do j=jmin,jmax
          do i=istrU,iend
            uee(i,j)=u(i,j-1,k)-2.*u(i,j,k)
     &                                 +u(i,j+1,k)
          enddo
        enddo
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=istrU,iend
            uee(i,jstr-1)=uee(i,jstr)
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=istrU,iend
            uee(i,jend+1)=uee(i,jend)
          enddo
        endif
#  endif
        do j=jstr,jend+1
          do i=istrU-1,iend
           Hvxx(i,j)=FlxV(i-1,j,k)-2.*FlxV(i,j,k)+FlxV(i+1,j,k)
          enddo
        enddo

        do j=jstr,jend+1
          do i=istrU,iend

            cff=FlxV(i,j,k)+FlxV(i-1,j,k)-0.125*( Hvxx(i  ,j)
     &                                           +Hvxx(i-1,j))
            UFe(i,j)=0.25*( cff*(u(i,j,k)+u(i,j-1,k))
     &                          -gamma*( max(cff,0.)*uee(i,j-1)
     &                                  +min(cff,0.)*uee(i,j  )
     &                                                      ))

            UFe(i,j)=UFe(i,j) - 0.25*( u(i,j,k)+u(i,j-1,k)
     &                        -0.125*(uee(i,j)+uee(i,j-1))
     &                  )*( FlxV(i,j,k)+FlxV(i-1,j,k)
     &                     -0.125*(Hvxx(i,j)+Hvxx(i-1,j)))

          enddo
        enddo
#  undef Hvxx
#  undef uee


#  define vxx wrk1
#  define Huee wrk2
#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          imin=istr
        else
          imin=istr-1
        endif
        if (EASTERN_EDGE) then
          imax=iend
        else
          imax=iend+1
        endif
#  else
        imin=istr-1
        imax=iend+1
#  endif
        do j=jstrV,jend
          do i=imin,imax
            vxx(i,j)=v(i-1,j,k)-2.*v(i,j,k)
     &                                 +v(i+1,j,k)
          enddo
        enddo
#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=jstrV,jend
            vxx(istr-1,j)=vxx(istr,j)
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=jstrV,jend
            vxx(iend+1,j)=vxx(iend,j)
          enddo
        endif
#  endif
        do j=jstrV-1,jend
          do i=istr,iend+1
           Huee(i,j)=FlxU(i,j-1,k)-2.*FlxU(i,j,k)+FlxU(i,j+1,k)
          enddo
        enddo
        do j=jstrV,jend
          do i=istr,iend+1
            cff=FlxU(i,j,k)+FlxU(i,j-1,k)-0.125*( Huee(i,j  )
     &                                           +Huee(i,j-1))
            VFx(i,j)=0.25*( cff*(v(i,j,k)+v(i-1,j,k))
     &                          -gamma*( max(cff,0.)*vxx(i-1,j)
     &                                  +min(cff,0.)*vxx(i  ,j)
     &                                                      ))

            VFx(i,j)=VFx(i,j) - 0.25*( v(i,j,k)+v(i-1,j,k)
     &                        -0.125*(vxx(i,j)+vxx(i-1,j))
     &                  )*( FlxU(i,j,k)+FlxU(i,j-1,k)
     &                     -0.125*(Huee(i,j)+Huee(i,j-1)))

          enddo
        enddo
#  undef Huee
#  undef vxx



        do j=jstr,jend
          do i=istrU,iend
            advu(i,j,k) = (-UFx(i,j)+UFx(i-1,j)
     &              -UFe(i,j+1)+UFe(i,j))
          enddo
        enddo

        do j=jstrV,jend
          do i=istr,iend
            advv(i,j,k) = (-VFx(i+1,j)+VFx(i,j)
     &              -VFe(i,j)+VFe(i,j-1))
          enddo
        enddo


      enddo !k




! Divide all diagnostic terms by the cell volume Hz/(pm*pn).
! There after the unit of diag terms are :
! (unit of velocity) * s-1  =  m * s-2


        do j=jstr,jend
          do i=istrU,iend
            cff=0.5*(pm(i,j)+pm(i-1,j))
     &             *(pn(i,j)+pn(i-1,j))


            do k=1,N  
       
              advu(i,j,k)=advu(i,j,k)*cff
     &             * u(i,j,k)
     &             /(Hz(i,j,k)+Hz(i-1,j,k))
            enddo
           enddo
         enddo


        do j=jstrV,jend
          do i=istr,iend

            cff=0.5*(pm(i,j)+pm(i,j-1))
     &             *(pn(i,j)+pn(i,j-1))

            do k=1,N  

              advv(i,j,k)=advv(i,j,k)*cff
     &             * v(i,j,k)
     &             /(Hz(i,j,k)+Hz(i,j-1,k))

            enddo
          enddo
        enddo



!!!!!!!!!!!!!!!!!

        do k=1,N  
          do i=istr,iend
            kediss(i,jstr,k)=0.
            kediss(i,jend,k)=0.
          enddo
        enddo


        do k=1,N  
          do j=jstr,jend
            kediss(istr,j,k)=0.
            kediss(iend,j,k)=0.
          enddo
        enddo

       write(*,*) istr,jstr,kediss(istr,jstr,1)
       write(*,*) istrU,jstrV,kediss(istr,jstr,1)


        do j=jstrV,jend-1
          do i=istrU,iend-1
            do k=1,N  

              kediss(i,j,k) = 0.5*( advu(i,j,k) + advu(i+1,j,k)
     &                            + advv(i,j,k) + advv(i,j+1,k) )
! #  ifdef MASKING
!      &                          *rmask(i,j)
! #  endif
            enddo
          enddo
        enddo

         write(*,*) istr,jstr,kediss(istr,jstr,1)

!!!!!!!!!!!!!!!!!!


      return
      end











