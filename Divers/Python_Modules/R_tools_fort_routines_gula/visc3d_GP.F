


# define MASKING
# define DIAGNOSTICS_UV



      subroutine visc3d_GP  (Lm,Mm,N,u,v,z_r,z_w
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
     &     VFe(0:Lm+1,0:Mm+1),     cff,
     &     UFx(0:Lm+1,0:Mm+1),     cff1,
     &     VFx(0:Lm+1,0:Mm+1),     cff2,
     &     UFs(0:Lm+1,0:Mm+1,2),   cff3,
     &     VFs(0:Lm+1,0:Mm+1,2),   cff4,
     &   dmUde(0:Lm+1,0:Mm+1,2),   cff5,
     &   dmVde(0:Lm+1,0:Mm+1,2),   cff6,
     &   dnUdx(0:Lm+1,0:Mm+1,2),   cff7,
     &   dnVdx(0:Lm+1,0:Mm+1,2),   cff8,
     &    dUdz(0:Lm+1,0:Mm+1,2),
     &    dVdz(0:Lm+1,0:Mm+1,2),   dmUdz,
     &  dZde_p(0:Lm+1,0:Mm+1,2),   dnUdz,
     &  dZde_r(0:Lm+1,0:Mm+1,2),   dmVdz,
     &  dZdx_p(0:Lm+1,0:Mm+1,2),   dnVdz,
     &  dZdx_r(0:Lm+1,0:Mm+1,2),
     &  dm_p(0:Lm+1,0:Mm+1), dn_p(0:Lm+1,0:Mm+1),
     &  dm_r(0:Lm+1,0:Mm+1), dn_r(0:Lm+1,0:Mm+1),
     &  wrk(0:Lm+1,0:Mm+1)


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


      !isp=min((LLm+1)/12,(MMm+1)/12)
       isp=(LLm+1)/12    !old version of the code



       do j=max(-1,jstr-1),jend
         do i=max(-1,istr-1),iend
          ibnd=isp

#   ifdef OBC_WEST
          ibnd=min(ibnd,i+coord(3))
#   endif
#   ifdef OBC_EAST
c          ibnd=min(ibnd,LLm+1-i-coord(3))
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
        enddo
      enddo


      do j=jstrR,jendR
        do i=istr,iendR
            dm_u(i,j) = 2./(pm(i,j)+pm(i-1,j))
            dn_u(i,j) = 2./(pn(i,j)+pn(i-1,j))
          enddo
        enddo

      do j=jstr,jendR
        do i=istrR,iendR
            dn_v(i,j) = 2./(pn(i,j)+pn(i,j-1))
            dm_v(i,j) = 2./(pm(i,j)+pm(i,j-1))
          enddo
        enddo


      do j=jstrR,jendR
        do i=istrR,iendR
              dm_r(i,j)=1./pm(i,j)
              dn_r(i,j)=1./pn(i,j)
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





!
! Compute horizontal harmonic viscosity along geopotential surfaces.
!--------------------------------------------------------------------
!
! Compute horizontal and vertical gradients.  Notice the recursive
! blocking sequence.  The vertical placement of the gradients is:
!
!    dZdx_r, dZde_r, dnUdx, dmVde(:,:,k1) k      rho-points
!    dZdx_r, dZde_r, dnUdx, dmVde(:,:,k2) k+1    rho-points
!    dZdx_p, dZde_p, dnVdx, dmUde(:,:,k1) k      psi-points
!    dZdx_p, dZde_p, dnVdx, dmUde(:,:,k2) k+1    psi-points
!                       UFs, dUdz(:,:,k1) k-1/2  WU-points
!                       UFs, dUdz(:,:,k2) k+1/2  WU-points
!                       VFs, dVdz(:,:,k1) k-1/2  WV-points
!                       VFs, dVdz(:,:,k2) k+1/2  WV-points
!
! Compute sigma-slopes (nondimensional) at RHO- and PSI-points.
! Compute momentum horizontal (1/m/s) and vertical (1/s) gradients.
!
      k2=1
      do k=0,N,+1  !--> irreversible
        k1=k2
        k2=3-k1
        if (k.lt.N) then
          do j=jstr-1,jend+1
            do i=istrU-1,iend+1
              UFx(i,j)=0.5*(z_r(i,j,k+1)-z_r(i-1,j,k+1))
     &                               *(pm(i-1,j)+pm(i,j))
# ifdef MASKING
     &                                        *umask(i,j)
# endif
            enddo
          enddo
          do j=jstrV-1,jend+1
            do i=istr-1,iend+1
              VFe(i,j)=0.5*(z_r(i,j,k+1)-z_r(i,j-1,k+1))
     &                               *(pn(i,j-1)+pn(i,j))
# ifdef MASKING
     &                                        *vmask(i,j)
# endif
            enddo
          enddo
          do j=jstrV-1,jend
            do i=istrU-1,iend
              dnUdx(i,j,k2)=0.5*pm(i,j)*(
     &               (pn(i  ,j)+pn(i+1,j))*u(i+1,j,k+1)
     &              -(pn(i-1,j)+pn(i  ,j))*u(i  ,j,k+1))
# ifdef MASKING
     &                                            *rmask(i,j)
# endif
              dmVde(i,j,k2)=0.5*pn(i,j)*(
     &               (pm(i,j  )+pm(i,j+1))*v(i,j+1,k+1)
     &              -(pm(i,j-1)+pm(i,j  ))*v(i,j  ,k+1))
# ifdef MASKING
     &                                            *rmask(i,j)
# endif
              dZdx_r(i,j,k2)=0.5*(UFx(i,j)+UFx(i+1,j))
              dZde_r(i,j,k2)=0.5*(VFe(i,j)+VFe(i,j+1))
            enddo
          enddo
          do j=jstr,jend+1
            do i=istr,iend+1
              dmUde(i,j,k2)=0.125*(pn(i,j)+pn(i-1,j)+pn(i,j-1)
     &                                            +pn(i-1,j-1))
     &            *( (pm(i-1,j  )+pm(i,j  ))*u(i,j  ,k+1)
     &              -(pm(i-1,j-1)+pm(i,j-1))*u(i,j-1,k+1))
# ifdef MASKING
     &                                             *pmask(i,j)
# endif
              dnVdx(i,j,k2)=0.125*(pm(i,j)+pm(i-1,j)+pm(i,j-1)
     &                                            +pm(i-1,j-1))
     &            *( (pn(i  ,j-1)+pn(i  ,j))*v(i  ,j,k+1)
     &              -(pn(i-1,j-1)+pn(i-1,j))*v(i-1,j,k+1))
# ifdef MASKING
     &                                             *pmask(i,j)
# endif
              dZde_p(i,j,k2)=0.5*(VFe(i-1,j)+VFe(i,j))
              dZdx_p(i,j,k2)=0.5*(UFx(i,j-1)+UFx(i,j))
            enddo
          enddo          !--> discard UFx,VFe, keep all others
        endif
 
 
c        do j=jstrV-1,jend+1       ! This
c          do i=istrU-1,iend+1
c            dZdx_r(i,j,k2)=0.
c            dZde_r(i,j,k2)=0.
c            dZdx_p(i,j,k2)=0.
c            dZde_p(i,j,k2)=0.
c          enddo
c        enddo
 
 
 
        if (k.eq.0 .or. k.eq.N) then
          do j=jstr-1,jend+1
            do i=istrU-1,iend+1
              dUdz(i,j,k2)=0.
              UFs(i,j,k2)=0.
            enddo
          enddo
          do j=jstrV-1,jend+1
            do i=istr-1,iend+1
              dVdz(i,j,k2)=0.
              VFs(i,j,k2)=0.
            enddo
          enddo
        else
          do j=jstr-1,jend+1
            do i=istrU-1,iend+1
              dUdz(i,j,k2)=2.*(u(i,j,k+1)-u(i,j,k))
     &                          /( z_r(i-1,j,k+1)-z_r(i-1,j,k)
     &                            +z_r(i  ,j,k+1)-z_r(i  ,j,k))
            enddo
          enddo
          do j=jstrV-1,jend+1
            do i=istr-1,iend+1
              dVdz(i,j,k2)=2.*(v(i,j,k+1)-v(i,j,k))
     &                          /( z_r(i,j-1,k+1)-z_r(i,j-1,k)
     &                            +z_r(i,j  ,k+1)-z_r(i,j  ,k))
            enddo
          enddo
        endif
!
! Compute components of the rotated viscous flux [m5/s2] along
! geopotential surfaces in the XI- and ETA-directions.
!
        if (k.gt.0) then
          do j=jstrV-1,jend
            do i=istrU-1,iend
              cff=visc2_r(i,j)*Hz(i,j,k)*(
     &         dn_r(i,j)*( dnUdx(i,j,k1) - 0.5*pn(i,j)*(
     &           min(dZdx_r(i,j,k1),0.)*(dUdz(i,j,k1)+dUdz(i+1,j,k2))
     &          +max(dZdx_r(i,j,k1),0.)*(dUdz(i,j,k2)+dUdz(i+1,j,k1))
     &                                                             ))
     &        -dm_r(i,j)*( dmVde(i,j,k1) - 0.5*pm(i,j)*(
     &           min(dZde_r(i,j,k1),0.)*(dVdz(i,j,k1)+dVdz(i,j+1,k2))
     &          +max(dZde_r(i,j,k1),0.)*(dVdz(i,j,k2)+dVdz(i,j+1,k1))
     &                                                            )))
# ifdef MASKING
     &                                                   *rmask(i,j)
# endif
              UFx(i,j)=dn_r(i,j)*dn_r(i,j)*cff
              VFe(i,j)=dm_r(i,j)*dm_r(i,j)*cff
            enddo
          enddo
          do j=jstr,jend+1
            do i=istr,iend+1
              cff=visc2_p(i,j)*0.25*( Hz(i,j,k) +Hz(i-1,j,k)
     &                               +Hz(i,j-1,k) +Hz(i-1,j-1,k))
     &       *( dn_p(i,j)*( dnVdx(i,j,k1)-0.125*( pn(i,j)+pn(i-1,j)
     &                                       +pn(i,j-1)+pn(i-1,j-1))
     &        *( min(dZdx_p(i,j,k1),0.)*(dVdz(i-1,j,k1)+dVdz(i,j,k2))
     &          +max(dZdx_p(i,j,k1),0.)*(dVdz(i-1,j,k2)+dVdz(i,j,k1))
 
     &     )) + dm_p(i,j)*( dmUde(i,j,k1)-0.125*( pm(i,j)+pm(i-1,j)
     &                                       +pm(i,j-1)+pm(i-1,j-1))
     &        *( min(dZde_p(i,j,k1),0.)*(dUdz(i,j-1,k1)+dUdz(i,j,k2))
     &          +max(dZde_p(i,j,k1),0.)*(dUdz(i,j-1,k2)+dUdz(i,j,k1))
     &                                                            )))
# ifdef MASKING
     &                                                   *pmask(i,j)
# endif
              UFe(i,j)=dm_p(i,j)*dm_p(i,j)*cff
              VFx(i,j)=dn_p(i,j)*dn_p(i,j)*cff
            enddo
          enddo
!



! Compute vertical flux [m^2/s^2] due to sloping terrain-following
! surfaces.
!
          if (k.lt.N) then
            do j=jstr,jend
              do i=istrU,iend
                cff1=0.5*(pn(i-1,j)+pn(i,j))
                cff2=0.5*(pm(i-1,j)+pm(i,j))
                cff=0.25*( dVdz(i,j,k2)+dVdz(i-1,j,k2)
     &                +dVdz(i,j+1,k2)+dVdz(i-1,j+1,k2))
                dnUdz=cff1*dUdz(i,j,k2)
                dmUdz=cff2*dUdz(i,j,k2)
                dnVdz=cff1*cff
                dmVdz=cff2*cff
 
                cff1=min(dZdx_r(i-1,j,k1),0.)
                cff2=min(dZdx_r(i  ,j,k2),0.)
                cff3=max(dZdx_r(i-1,j,k2),0.)
                cff4=max(dZdx_r(i  ,j,k1),0.)
                cff5=min(dZde_r(i-1,j,k1),0.)
                cff6=min(dZde_r(i  ,j,k2),0.)
                cff7=max(dZde_r(i-1,j,k2),0.)
                cff8=max(dZde_r(i  ,j,k1),0.)
 
                cff=dn_u(i,j)*( cff1*(cff1*dnUdz-dnUdx(i-1,j,k1))
     &                         +cff2*(cff2*dnUdz-dnUdx(i  ,j,k2))
     &                         +cff3*(cff3*dnUdz-dnUdx(i-1,j,k2))
     &                         +cff4*(cff4*dnUdz-dnUdx(i  ,j,k1))
     &                                                          )
     &             -dm_u(i,j)*( cff1*(cff5*dmVdz-dmVde(i-1,j,k1))
     &                         +cff2*(cff6*dmVdz-dmVde(i  ,j,k2))
     &                         +cff3*(cff7*dmVdz-dmVde(i-1,j,k2))
     &                         +cff4*(cff8*dmVdz-dmVde(i  ,j,k1))
     &                                                          )
                cff1=min(dZde_p(i,j  ,k1),0.)
                cff2=min(dZde_p(i,j+1,k2),0.)
                cff3=max(dZde_p(i,j  ,k2),0.)
                cff4=max(dZde_p(i,j+1,k1),0.)
                cff5=min(dZdx_p(i,j  ,k1),0.)
                cff6=min(dZdx_p(i,j+1,k2),0.)
                cff7=max(dZdx_p(i,j  ,k2),0.)
                cff8=max(dZdx_p(i,j+1,k1),0.)
 
                cff=cff + dm_u(i,j)*(
     &                          cff1*(cff1*dmUdz-dmUde(i,j  ,k1))
     &                         +cff2*(cff2*dmUdz-dmUde(i,j+1,k2))
     &                         +cff3*(cff3*dmUdz-dmUde(i,j  ,k2))
     &                         +cff4*(cff4*dmUdz-dmUde(i,j+1,k1))
     &                                                          )
     &             +dn_u(i,j)*( cff1*(cff5*dnVdz-dnVdx(i,j  ,k1))
     &                         +cff2*(cff6*dnVdz-dnVdx(i,j+1,k2))
     &                         +cff3*(cff7*dnVdz-dnVdx(i,j  ,k2))
     &                         +cff4*(cff8*dnVdz-dnVdx(i,j+1,k1))
     &                                                          )
 
                UFs(i,j,k2)=0.25*(visc2_r(i-1,j)+visc2_r(i,j))*cff
              enddo
            enddo
 
            do j=jstrV,jend
              do i=istr,iend
                cff1=0.5*(pn(i,j-1)+pn(i,j))
                cff2=0.5*(pm(i,j-1)+pm(i,j))
                cff=0.25*( dUdz(i,j,k2)+dUdz(i+1,j,k2)
     &                +dUdz(i,j-1,k2)+dUdz(i+1,j-1,k2))
                dnUdz=cff1*cff
                dmUdz=cff2*cff
                dnVdz=cff1*dVdz(i,j,k2)
                dmVdz=cff2*dVdz(i,j,k2)
 
                cff1=min(dZdx_p(i  ,j,k1),0.)
                cff2=min(dZdx_p(i+1,j,k2),0.)
                cff3=max(dZdx_p(i  ,j,k2),0.)
                cff4=max(dZdx_p(i+1,j,k1),0.)
                cff5=min(dZde_p(i  ,j,k1),0.)
                cff6=min(dZde_p(i+1,j,k2),0.)
                cff7=max(dZde_p(i  ,j,k2),0.)
                cff8=max(dZde_p(i+1,j,k1),0.)
 
                cff=dn_v(i,j)*( cff1*(cff1*dnVdz-dnVdx(i  ,j,k1))
     &                         +cff2*(cff2*dnVdz-dnVdx(i+1,j,k2))
     &                         +cff3*(cff3*dnVdz-dnVdx(i  ,j,k2))
     &                         +cff4*(cff4*dnVdz-dnVdx(i+1,j,k1))
     &                                                          )
     &             +dm_v(i,j)*( cff1*(cff5*dmUdz-dmUde(i  ,j,k1))
     &                         +cff2*(cff6*dmUdz-dmUde(i+1,j,k2))
     &                         +cff3*(cff7*dmUdz-dmUde(i  ,j,k2))
     &                         +cff4*(cff8*dmUdz-dmUde(i+1,j,k1))
     &                                                          )
                cff1=min(dZde_r(i,j-1,k1),0.)
                cff2=min(dZde_r(i,j  ,k2),0.)
                cff3=max(dZde_r(i,j-1,k2),0.)
                cff4=max(dZde_r(i,j  ,k1),0.)
                cff5=min(dZdx_r(i,j-1,k1),0.)
                cff6=min(dZdx_r(i,j  ,k2),0.)
                cff7=max(dZdx_r(i,j-1,k2),0.)
                cff8=max(dZdx_r(i,j  ,k1),0.)
 
                cff=cff+dm_v(i,j)*(
     &                          cff1*(cff1*dmVdz-dmVde(i,j-1,k1))
     &                         +cff2*(cff2*dmVdz-dmVde(i,j  ,k2))
     &                         +cff3*(cff3*dmVdz-dmVde(i,j-1,k2))
     &                         +cff4*(cff4*dmVdz-dmVde(i,j  ,k1))
     &                                                          )
     &             -dn_v(i,j)*( cff1*(cff5*dnUdz-dnUdx(i,j-1,k1))
     &                         +cff2*(cff6*dnUdz-dnUdx(i,j  ,k2))
     &                         +cff3*(cff7*dnUdz-dnUdx(i,j-1,k2))
     &                         +cff4*(cff8*dnUdz-dnUdx(i,j  ,k1))
     &                                                          )
 
                VFs(i,j,k2)=0.25*(visc2_r(i,j-1)+visc2_r(i,j))*cff
              enddo
            enddo
          endif
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
              MHmix(i,j,k,1) = (cff+UFs(i,j,k2)-UFs(i,j,k1))
     &                                  *dm_u(i,j)*dn_u(i,j)

# endif





            enddo
          enddo
          do j=jstrV,jend
            do i=istr,iend
              cff=0.125*(pm(i,j)+pm(i,j-1))*(pn(i,j) +pn(i,j-1))
     &               *( (pn(i,j-1)+pn(i,j))*(VFx(i+1,j)-VFx(i,j))
     &                 -(pm(i,j-1)+pm(i,j))*(VFe(i,j)-VFe(i,j-1))
     &                                                          )

# ifdef DIAGNOSTICS_UV
              MHmix(i,j,k,2) = (cff+VFs(i,j,k2)-VFs(i,j,k1))
     &                                  *dm_v(i,j)*dn_v(i,j)
# endif



            enddo
          enddo
        endif
      enddo
      return
      end

 
