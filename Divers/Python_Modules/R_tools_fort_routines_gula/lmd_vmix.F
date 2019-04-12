
# define SMOOTH_RIG  /* local Ri number smoothing */
 

      subroutine lmd_vmix (Lm,Mm,N,u,v,z_r
     &  ,bvf,Kv,Kt)
!
!  This subroutine computes vertical mixing coefficients for momentum
!  and tracers at the ocean interior using the Large, McWilliams and
!  Doney (1994) mixing scheme.
!
!  On Output:
!     Kv       vertical viscosity coefficient [m^2/s].
!     Kt       vertical diffusion coefficient for potential
!                 temperature [m^2/s].
!     Ks       vertical diffusion coefficient for salinity [m^2/s].
!
!  Reference:
!
!    Large, W.G., J.C. McWilliams, and S.C. Doney, 1994: A Review
!    and model with a nonlocal boundary layer parameterization,
!    Reviews of Geophysics, 32,363-403.
!
      implicit none

      integer Lm,Mm,N,imin,imax,jmin,jmax, i,j,k


      real*8 Rig(0:Lm+1,0:Mm+1,0:N),
     &      Kv(0:Lm+1,0:Mm+1,0:N),
     &      Kt(0:Lm+1,0:Mm+1,0:N)
!      &      ,Ks(0:Lm+1,0:Mm+1,0:N)

      real nu_sx, cff,dudz,dvdz

      real Ri0,  nuwm, nuws,  nu0m, nu0s, nu0c, lmd_nu, lmd_Rrho0,
     &     lmd_nuf,  lmd_fdd, lmd_tdd1, lmd_tdd2, lmd_tdd3, lmd_sdd1,
     &     lmd_sdd2, lmd_sdd3, eps

      parameter ( Ri0=0.7,
     &     nu0m=50.e-4,
     &     nu0s=50.e-4,
     &     nuwm=1.0e-4,
     &     nuws=0.1e-4,
     &     nu0c=0.1)

      parameter (eps=1.E-14)
 

      ! Variables IN
      real*8 z_r(0:Lm+1,0:Mm+1,N)
     &      ,u(1:Lm+1,0:Mm+1,N), v(0:Lm+1,1:Mm+1,N)
     &      ,bvf(0:Lm+1,0:Mm+1,0:N)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!# include "compute_tile_bounds.h"
# include "compute_auxiliary_bounds.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



Cf2py intent(in) Lm,Mm,N,u, v,z_r, bvf
Cf2py intent(inout) Kv,Kt

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# ifdef SMOOTH_RIG 
#  ifdef EW_PERIODIC
      imin=istr-1
      imax=iend+1
#  else
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
#  endif
#  ifdef NS_PERIODIC
      jmin=jstr-1
      jmax=jend+1
#  else
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
#  endif
#  define I_EXT_RANGE imin,imax
#  define J_EXT_RANGE jmin,jmax
# else
#  define I_EXT_RANGE istr,iend
#  define J_EXT_RANGE jstr,jend
# endif



!                                    ! Compute horizontal velocity
! Compute local Richardson number:   ! shear  (du/dz)^2+(dv/dz)^2  at
!-------- ----- ---------- -------   ! horizontal RHO- and vertical
!                                    ! W-points. Then compute gradient
      do k=1,N-1                     ! Richardson number (already
        do j=J_EXT_RANGE             ! divided by its critical value.
          do i=I_EXT_RANGE
            cff=0.5/(z_r(i,j,k+1)-z_r(i,j,k))
            dudz=cff*( u(i  ,j,k+1)-u(i  ,j,k)
     &                +u(i+1,j,k+1)-u(i+1,j,k))
            dvdz=cff*( v(i,j  ,k+1)-v(i,j  ,k)
     &                +v(i,j+1,k+1)-v(i,j+1,k))

            Rig(i,j,k)=bvf(i,j,k)/( Ri0*max(
     &                    dudz*dudz+dvdz*dvdz, 1.E-10 ))
          enddo
        enddo

# ifdef SMOOTH_RIG
#  ifndef EW_PERIODIC
        if (WESTERN_EDGE) then
          do j=J_EXT_RANGE
            Rig(istr-1,j,k)=Rig(istr,j,k)
          enddo
        endif
        if (EASTERN_EDGE) then
          do j=J_EXT_RANGE
            Rig(iend+1,j,k)=Rig(iend,j,k)
          enddo
        endif
#  endif
#  ifndef NS_PERIODIC
        if (SOUTHERN_EDGE) then
          do i=I_EXT_RANGE
            Rig(i,jstr-1,k)=Rig(i,jstr,k)
          enddo
        endif
        if (NORTHERN_EDGE) then
          do i=I_EXT_RANGE
            Rig(i,jend+1,k)=Rig(i,jend,k)
          enddo
        endif
#   ifndef EW_PERIODIC
        if (WESTERN_EDGE .and.
     &      SOUTHERN_EDGE) then
          Rig(istr-1,jstr-1,k)=Rig(istr,jstr,k)
        endif
        if (WESTERN_EDGE .and.
     &      NORTHERN_EDGE) then
          Rig(istr-1,jend+1,k)=Rig(istr,jend,k)
        endif
        if (EASTERN_EDGE .and.
     &      SOUTHERN_EDGE) then
          Rig(iend+1,jstr-1,k)=Rig(iend,jstr,k)
        endif
        if (EASTERN_EDGE .and.
     &      NORTHERN_EDGE) then
          Rig(iend+1,jend+1,k)=Rig(iend,jend,k)
        endif
#   endif
#  endif
                                   ! Smooth Rig horizontally: use
        do j=jstr-1,jend           ! array Rig(:,:,0) as scratch. 
          do i=istr-1,iend
            Rig(i,j,0)=0.25*(Rig(i,j  ,k)+Rig(i+1,j  ,k)
     &                      +Rig(i,j+1,k)+Rig(i+1,j+1,k))
          enddo
        enddo
        do j=jstr,jend
          do i=istr,iend
            Rig(i,j,k)=0.25*(Rig(i,j  ,0)+Rig(i-1,j  ,0)
     &                      +Rig(i,j-1,0)+Rig(i-1,j-1,0))
          enddo
        enddo    !--> discard  Rig(:,:,0)
# endif /* SMOOTH_RIG */

!
! Compute "interior" viscosities and diffusivities everywhere
! as the superposition of three processes: local Richardson number
! instability due to resolved vertical shear, internal wave breaking,
! and double diffusion.
!
        do j=jstr,jend
          do i=istr,iend
# ifdef LMD_RIMIX
            cff=min(1., max(0., Rig(i,j,k)))  ! Compute mixing die
            nu_sx=1. - cff*cff                ! to shear instability
            nu_sx=nu_sx*nu_sx*nu_sx           ! and internal wave 
                                              ! breaking.
            Kv(i,j,k)=nuwm + nu0m*nu_sx
            Kt(i,j,k)=nuws + nu0s*nu_sx

# else
            Kv(i,j,k)=nuwm                    ! internal waves
            Kt(i,j,k)=nuws                    ! only
# endif
          enddo
        enddo
      enddo     ! <-- k

!
! Pad out surface and bottom values for lmd_blmix calculations.
! The interior values used here may not be the best values to
! use for the padding.
!
!      do j=jstr,jend
!        do i=istr,iend
!          Kv(i,j,N)=Kv(i,j,N-1)
!           Ks(i,j,N)=Ks(i,j,N-1)
!          Kt(i,j,N)=Kt(i,j,N-1)
!          Kv(i,j,0)=Kv(i,j,  1)
!           Ks(i,j,0)=Ks(i,j,  1)
!          Kt(i,j,0)=Kt(i,j,  1)
!        enddo
!      enddo


 

      return
      end

 
 
