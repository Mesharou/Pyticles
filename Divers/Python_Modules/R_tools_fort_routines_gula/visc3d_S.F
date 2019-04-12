


# define DIAGNOSTICS_UV



      subroutine visc3d_S  (Lm,Mm,N,u,v,z_r,z_w
     &  ,pm,pn,f,dt,rmask,umask, vmask
     &  ,visc2, v_sponge,coord,coordmax
     &  ,MHmix)

!
! Compute horizontal (along geopotential surfaces) viscous terms as
! divergence of symmetric stress tensor.
!
! Compute harmonic mixing of momentum, rotated along geopotentials,
! from the horizontal  divergence  of the  stress  tensor.
! A transverse  isotropy is assumed so the stress tensor is splitted
! into vertical and horizontal subtensors.
!
!  Reference:
!
! [1] Stelling, G. S., and J. A. Th. M. van Kester, 1994: On the
!       approximation of horizontal gradients in sigma-coordinates
!       for bathymetry with steep bottom slopes. Int. J. Num. Meth.
!       in Fluids, v. 18, pp. 915-935.
!
! [2] Wajsowicz, R.C, 1993: A consistent formulation of the
!       anisotropic stress tensor for use in models of the
!       large-scale ocean circulation, JCP, 105, 333-338.
!
! [3] Sadourny, R. and K. Maynard, 1997: Formulations of lateral
!       diffusion in geophysical fluid dynamics models, In
!       "Numerical Methods of Atmospheric and Oceanic Modelling".
!       Lin, Laprise, and Ritchie, Eds., NRC Research Press,
!       547-556.
!
! [4] Griffies, S.M. and R.W. Hallberg, 2000: Biharmonic friction
!       with a Smagorinsky-like viscosity for use in large-scale
!       eddy-permitting ocean models, Monthly Weather Rev.,v. 128,
!       No. 8, pp. 2935-2946.
!
      implicit none



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      integer Lm,Mm,LLm,MMm,N, imin,imax,jmin,jmax, i,j,k,
     &      dt,
     &      itrc, k1,k2,isp,ibnd

        !INPUT
      real*8 u(1:Lm+1,0:Mm+1,N),  v(0:Lm+1,1:Mm+1,N)
     &      ,z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N)
     &      ,pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1)
     &      ,rmask(0:Lm+1,0:Mm+1), pmask(1:Lm+1,1:Mm+1)
     &      ,umask(1:Lm+1,0:Mm+1),  vmask(0:Lm+1,1:Mm+1)
     &      ,f(0:Lm+1,0:Mm+1)

      integer coord(4), coordmax(4)

      ! OUTPUTS
      real*8  MHmix(0:Lm+1,0:Mm+1,N,2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real visc2, v_sponge

      real UFe(0:Lm+1,0:Mm+1),
     &     VFe(0:Lm+1,0:Mm+1),     cff,cff1,cff2,
     &     UFx(0:Lm+1,0:Mm+1),
     &     VFx(0:Lm+1,0:Mm+1),
     &  dm_p(0:Lm+1,0:Mm+1), dn_p(0:Lm+1,0:Mm+1),
     &  dm_r(0:Lm+1,0:Mm+1), dn_r(0:Lm+1,0:Mm+1),
     &     wrk(0:Lm+1,0:Mm+1)


      real*8 Hz(0:Lm+1,0:Mm+1,N),
     &      dn_u(0:Lm+1,0:Mm+1), dm_v(0:Lm+1,0:Mm+1),
     &      dm_u(0:Lm+1,0:Mm+1), dn_v(0:Lm+1,0:Mm+1),
     &      visc2_r(0:Lm+1,0:Mm+1), visc2_p(0:Lm+1,0:Mm+1)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# include "compute_auxiliary_bounds.h"



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


Cf2py intent(in) Lm,Mm,N,u,v,z_r,z_w,pm,pn,f,dt,rmask,umask, vmask,visc2, v_sponge,coord,coordmax
Cf2py intent(out) MHmix


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! define visc coefficient (from set_nudgcof.F )
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      LLm = coordmax(4)-2
      MMm = coordmax(2)-2


       isp=min((LLm+1)/12,(MMm+1)/12)
      ! isp=(LLm+1)/12    !old version of the code


       do j=max(-1,jstr-1),jend
         do i=max(-1,istr-1),iend
          ibnd=isp

#   ifdef OBC_WEST
          ibnd=min(ibnd,i+coord(3))
#   endif
#   ifdef OBC_EAST
          ibnd=min(ibnd,LLm+1-i-coord(3))
#   endif
#   ifdef OBC_SOUTH
          ibnd=min(ibnd,j+coord(1))
#   endif
#   ifdef OBC_NORTH
          ibnd=min(ibnd,MMm+1-j-coord(1))
#   endif

          wrk(i,j)=float(isp-ibnd)/float(isp)
        enddo
       enddo


       do j=jstr-1,jend+1
        do i=istr-1,iend+1
#  ifdef SPONGE
#   ifdef UV_VIS2
          visc2_r(i,j)=visc2 + v_sponge*wrk(i,j)
#   endif
#  endif
        enddo
       enddo

       do j=jstr,jend+1
        do i=istr,iend+1
#  ifdef SPONGE
#   ifdef UV_VIS2
          visc2_p(i,j)=visc2 + 0.25*v_sponge*( wrk(i,j)
     &                  +wrk(i-1,j)+wrk(i,j-1)+wrk(i-1,j-1))
#   endif
#  endif
        enddo
       enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! define grid variables
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      do j=jstrR,jendR
        do i=istrR,iendR
          do k=1,N,+1
           Hz(i,j,k) = z_w(i,j,k) - z_w(i,j,k-1)
          enddo
              dm_r(i,j)=1./pm(i,j)
              dn_r(i,j)=1./pn(i,j)
        enddo
      enddo

      do j=jstrR,jendR
        do i=istr,iendR
            dn_u(i,j) = 2./(pn(i,j)+pn(i-1,j))
            dm_u(i,j) = 2./(pm(i,j)+pm(i-1,j))
          enddo
      enddo



      do j=jstr,jendR
        do i=istrR,iendR
            dm_v(i,j) = 2./(pm(i,j)+pm(i,j-1))
            dn_v(i,j) = 2./(pn(i,j)+pn(i,j-1))
          enddo
      enddo



! Compute n/m and m/n at horizontal PSI-points.
! Set mask according to slipperness parameter gamma.
!
      do j=jstrV-1,jend
        do i=istrU-1,iend


          dm_p(i,j)=4./(pm(i,j)+pm(i,j-1)+pm(i-1,j)+pm(i-1,j-1))
          dn_p(i,j)=4./(pn(i,j)+pn(i,j-1)+pn(i-1,j)+pn(i-1,j-1))

#ifdef MASKING
c          pmask(i,j)=rmask(i,j)*rmask(i-1,j)*rmask(i,j-1)
c     &                                      *rmask(i-1,j-1)
C***      if (gamma2.lt.0.) pmask(i,j)=2.-pmask(i,j)
!
! Set no-slip boundary conditions on land-mask boundaries
! regardless of supplied value of gamma2.
!

          cff1=1.  !<-- computation of off-diagonal nonlinear terms
          cff2=2.

          if (rmask(i-1,j  ).gt.0.5 .and. rmask(i,j  ).gt.0.5 .and.
     &        rmask(i-1,j-1).gt.0.5 .and. rmask(i,j-1).gt.0.5) then
            pmask(i,j)=1.

          elseif(rmask(i-1,j  ).lt.0.5 .and.rmask(i,j  ).gt.0.5 .and.
     &           rmask(i-1,j-1).gt.0.5 .and.rmask(i,j-1).gt.0.5) then
            pmask(i,j)=cff1
          elseif(rmask(i-1,j  ).gt.0.5 .and.rmask(i,j  ).lt.0.5 .and.
     &           rmask(i-1,j-1).gt.0.5 .and.rmask(i,j-1).gt.0.5) then
            pmask(i,j)=cff1
          elseif(rmask(i-1,j  ).gt.0.5 .and.rmask(i,j  ).gt.0.5 .and.
     &           rmask(i-1,j-1).lt.0.5 .and.rmask(i,j-1).gt.0.5) then
            pmask(i,j)=cff1
          elseif(rmask(i-1,j  ).gt.0.5 .and.rmask(i,j  ).gt.0.5 .and.
     &           rmask(i-1,j-1).gt.0.5 .and.rmask(i,j-1).lt.0.5) then
            pmask(i,j)=cff1


          elseif(rmask(i-1,j  ).gt.0.5 .and.rmask(i,j  ).lt.0.5 .and.
     &           rmask(i-1,j-1).gt.0.5 .and.rmask(i,j-1).lt.0.5) then
            pmask(i,j)=cff2
          elseif(rmask(i-1,j  ).lt.0.5 .and.rmask(i,j  ).gt.0.5 .and.
     &           rmask(i-1,j-1).lt.0.5 .and.rmask(i,j-1).gt.0.5) then
            pmask(i,j)=cff2
          elseif(rmask(i-1,j  ).gt.0.5 .and.rmask(i,j  ).gt.0.5 .and.
     &           rmask(i-1,j-1).lt.0.5 .and.rmask(i,j-1).lt.0.5) then
            pmask(i,j)=cff2
          elseif(rmask(i-1,j  ).lt.0.5 .and.rmask(i,j  ).lt.0.5 .and.
     &           rmask(i-1,j-1).gt.0.5 .and.rmask(i,j-1).gt.0.5) then
            pmask(i,j)=cff2

          else
            pmask(i,j)=0.
          endif
#endif

        enddo
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!  Compute horizontal harmonic viscosity along constant S-surfaces.
!
!  Compute flux-components of the horizontal divergence of the stress
!  tensor (m5/s2) in XI- and ETA-directions.
!
      do k=1,N


        do j=jstrV-1,jend
          do i=istrU-1,iend
            cff=0.5*Hz(i,j,k)*visc2_r(i,j)*(
 
     &      dn_r(i,j)*pm(i,j)*( (pn(i  ,j)+pn(i+1,j))*u(i+1,j,k)
     &                         -(pn(i-1,j)+pn(i  ,j))*u(i  ,j,k)
     &                                                              )

     &     -dm_r(i,j)*pn(i,j)*( (pm(i,j  )+pm(i,j+1))*v(i,j+1,k)
     &                         -(pm(i,j-1)+pm(i,j  ))*v(i,j  ,k)
     &                                                             ))
            UFx(i,j)=  cff*dn_r(i,j)*dn_r(i,j)
            VFe(i,j)= -cff*dm_r(i,j)*dm_r(i,j)

          enddo
        enddo



        do j=jstr,jend+1
          do i=istr,iend+1
            cff=0.125*(Hz(i-1,j,k)+Hz(i,j,k)+Hz(i-1,j-1,k)
     &                         +Hz(i,j-1,k))*visc2_p(i,j)*(
 
     &      0.25*(pm(i-1,j)+pm(i,j)+pm(i-1,j-1)+pm(i,j-1))*dn_p(i,j)
     &                     *( (pn(i  ,j-1)+pn(i  ,j))*v(i  ,j,k)
     &                       -(pn(i-1,j-1)+pn(i-1,j))*v(i-1,j,k)
     &                                                              )

     &     +0.25*(pn(i-1,j)+pn(i,j)+pn(i-1,j-1)+pn(i,j-1))*dm_p(i,j)
     &                     *( (pm(i-1,j  )+pm(i,j  ))*u(i,j  ,k)
     &                       -(pm(i-1,j-1)+pm(i,j-1))*u(i,j-1,k)
     &                                                             ))
# ifdef MASKING
     &                                                   *pmask(i,j)
# endif
            UFe(i,j)= cff*dm_p(i,j)*dm_p(i,j)
            VFx(i,j)= cff*dn_p(i,j)*dn_p(i,j)
          enddo
        enddo




! 
! 
!         if (k.eq.1) then
!         write(*,*) 'u(1:5,82,0,nstp)',u(1:5,82,1)
!         write(*,*) 'v(1:5,82,0,nstp)',v(1:5,82,1)
!         write(*,*) 'u(1:5,81,0,nstp)',u(1:5,81,1)
!         write(*,*) 'v(1:5,81,0,nstp)',v(1:5,81,1)
!         write(*,*) 'UFx(i,j)',UFx(1:5,82)
!         write(*,*) 'VFe(i,j)',VFe(1:5,82)
!         write(*,*) 'UFe(i,j)',UFe(1:5,82)
!         write(*,*) 'VFx(i,j)',VFx(1:5,82)
!         write(*,*) 'dm_p(i,j)*dm_p(i,j)',dm_p(5,82),dn_p(5,82)
!         write(*,*) 'pmask(i,j)',pmask(5,82)
!         write(*,*) 'pm',pm(4:5,81:82)
!         write(*,*) 'pn',pn(4:5,81:82)
!         write(*,*) 'visc2_p',visc2_p(5,82)
!         write(*,*) 'Hz',Hz(4,82,1)+Hz(5,82,1)+Hz(4,81,1)+Hz(5,81,1)
!         endif
! 
! 



!
! Apply viscous terms. Note that at this stage arrays u,v(...,3-nstp)
! contain Hz*U and Hz*V with units of [m2/s].   Also compute vertical
! integral of viscous terms and add it into coupling terms for the
! barotropic mode
!

          do j=jstr,jend
            do i=istrU,iend
              cff=0.125*(pm(i-1,j)+pm(i,j))*(pn(i-1,j) +pn(i,j))
     &               *( (pn(i-1,j)+pn(i,j))*(UFx(i,j)-UFx(i-1,j))
     &                 +(pm(i-1,j)+pm(i,j))*(UFe(i,j+1)-UFe(i,j))
     &                                                          )
# ifdef DIAGNOSTICS_UV
              MHmix(i,j,k,1) = cff *dm_u(i,j)*dn_u(i,j)
# endif
            enddo
          enddo


        do j=jstrV,jend
          do i=istr,iend
            cff=0.125*(pm(i,j)+pm(i,j-1))*(pn(i,j) +pn(i,j-1))
     &             *( (pn(i,j-1)+pn(i,j))*(VFx(i+1,j)-VFx(i,j))
     &               +(pm(i,j-1)+pm(i,j))*(VFe(i,j)-VFe(i,j-1))
     &                                                        )

# ifdef DIAGNOSTICS_UV
              MHmix(i,j,k,2) = cff*dm_v(i,j)*dn_v(i,j)
# endif


            enddo
          enddo


      enddo
      return
      end

 
