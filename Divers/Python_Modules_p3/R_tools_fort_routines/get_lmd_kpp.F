#include "cppdefs.h"
#if defined SOLVE3D && (defined LMD_KPP || defined LMD_BKPP)
#define INT_AT_W_POINTS

#define SMOOTH_HBL
#define LIMIT_UNSTABLE_ONLY

#define MERGE_OVERWRAP
c  MERGE_OVERWRAP is effective if both of SKPP and BKPP are defined.
c  By turning this ON, in case that surface and bottom mixed layer 
c  overwrap, hbl and hbbl are enforced to be a total depth. 


      subroutine get_lmd_kpp (istr,iend,jstr,jend, Kv,Kt,Ks,
     &                         ustar, Bo,Bosol,hbl, 
#ifdef LMD_BKPP
     &                         hbbl,    
#endif     
     &                          FX,FE,FE1,
     &                         Cr,FC, wrk1,wrk2,
     &                     Gm1,dGm1dS, Gt1,dGt1dS,  Gs1,dGs1dS,
     &                                                kbl
#ifdef LIMIT_MO_DEPTH     
     &    ,kmo
#endif      
     &                       )
      implicit none
# include "param.h"
# include "grid.h"
# include "ocean3d.h"
# include "forces.h"
# include "mixing.h"
# include "scalars.h"
# ifdef LMD_BKPP
# include "coupling.h"  
# endif
# define tind nstp

      integer istr,iend,jstr,jend, i,j,k
# ifdef SMOOTH_HBL
     &      , imin,imax,jmin,jmax
# endif 
      real Kv(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &     Kt(PRIVATE_2D_SCRATCH_ARRAY,0:N),
     &     Ks(PRIVATE_2D_SCRATCH_ARRAY,0:N),

     &  ustar(PRIVATE_2D_SCRATCH_ARRAY),  ustar3,
     &     Bo(PRIVATE_2D_SCRATCH_ARRAY),  Bfsfc,   zscale,
     &  Bosol(PRIVATE_2D_SCRATCH_ARRAY),  zetahat, ws,wm,
     &    hbl(PRIVATE_2D_SCRATCH_ARRAY),
     &    hbbl(PRIVATE_2D_SCRATCH_ARRAY),     
# ifdef SMOOTH_HBL
     &     FX(PRIVATE_2D_SCRATCH_ARRAY),
     &     FE(PRIVATE_2D_SCRATCH_ARRAY),
     &    FE1(PRIVATE_2D_SCRATCH_ARRAY),
# else
     &         FX, FE, FE1,
# endif

     &     Cr(PRIVATE_1D_SCRATCH_ARRAY,0:N),
     &     FC(PRIVATE_1D_SCRATCH_ARRAY,0:N),
     &   wrk1(PRIVATE_1D_SCRATCH_ARRAY,0:N),
     &   wrk2(PRIVATE_1D_SCRATCH_ARRAY,0:N),

     & Bfsfc_bl(PRIVATE_1D_SCRATCH_ARRAY),   z_bl,

     &    Gm1(PRIVATE_1D_SCRATCH_ARRAY),   Av_bl,
     & dGm1dS(PRIVATE_1D_SCRATCH_ARRAY),  dAv_bl, f1,
     &    Gt1(PRIVATE_1D_SCRATCH_ARRAY),   At_bl, a1, 
     & dGt1dS(PRIVATE_1D_SCRATCH_ARRAY),  dAt_bl, a2,
     &    Gs1(PRIVATE_1D_SCRATCH_ARRAY),   As_bl, a3,
     & dGs1dS(PRIVATE_1D_SCRATCH_ARRAY),  dAs_bl

      integer kmo(PRIVATE_1D_SCRATCH_ARRAY) 
      integer kbl(PRIVATE_1D_SCRATCH_ARRAY)

 
      real Kern, Vtc, Vtsq,  sigma, cff,cff1, cff_up,cff_dn
# ifdef LIMIT_MO_DEPTH
     &                                             , h_MO  
# ifdef INT_AT_RHO_POINTS
     &                                             , z_up
# endif
#endif

      real nubl, nu0c, Cv, Ricr, Ri_inv, betaT, epssfc, C_Ek, C_MO,
     &     Cstar, Cg, eps, zeta_m, a_m, c_m,   zeta_s, a_s, c_s,
     &                                                    r2,r3,r4
      parameter (
     &   nubl=0.01,   ! maximum allowed boundary layer

     &   nu0c=0.1,    ! convective adjustment for viscosity and diffusivity [m^2/s].                            
 
     &   Cv=1.8,      ! ratio of interior Brunt-Vaisala frequency
                      ! "N" to that at the entrainment depth "he".
     &   Ricr=0.45,   ! Critical bulk Richardson number (0.3).

     &   Ri_inv=1./Ricr,

     &   betaT=-0.2,  ! ratio of entrainment flux to surface
                      !         buoyancy flux.
 
     &   epssfc=0.1,  ! nondimensional extent of the surface layer.
 
     &   C_MO=1.,     ! constant for computaion Monin-Obukhov depth.

     &   C_Ek=258.,   ! constant for computating stabilization term 
                      ! due to Coriolis force (Ekman depth limit).
 
     &   Cstar=10.,   ! proportionality coefficient parameterizing
                      ! nonlocal transport.
     &   eps=1.E-20,

     &   zeta_m=-0.2, ! Maximum stability parameters "zeta"
     &   a_m=1.257,   ! value of the 1/3 power law regime of
     &   c_m=8.360,   ! flux profile for momentum and tracers
     &   zeta_s=-1.0, ! and coefficients of flux profile for
     &   a_s=-28.86,  ! momentum and tracers in their 1/3-power
     &   c_s=98.96,   ! law regime;

     &   r2=0.5, r3=1./3., r4=0.25)

#ifdef LMD_BKPP
      real  ustar2, Kv0, Kt0, Ks0, my_Akv_bak, 
     &                        my_Akt_bak, my_Aks_bak
#endif
 
      Cg=Cstar * vonKar * (c_s*vonKar*epssfc)**(1./3.)
      Vtc=Cv * sqrt(-betaT/(c_s*epssfc)) / (Ricr*vonKar**2)

# ifdef SMOOTH_HBL
#  define I_EXT_RANGE imin,imax
#  define J_EXT_RANGE jmin,jmax
# else
#  define I_EXT_RANGE istr,iend
#  define J_EXT_RANGE jstr,jend
# endif
!
! Compute thermal expansion coefficient "alpha" [kg/m^3/decC] and
! saline contraction coefficient "beta" [kg/m^3/PSU] at surface, then
! compute surface turbulent buoyancy forcing "Bo" [m^2/s^3] (in doing
! so remove incoming solar shortwave radiation component and save it
! separately as "Bosol"). Also get an approximation for ssurface
! layer depth using "epssfc" and boundary layer depth from previous
! time step (this is needed to estimate turbulent velocity scale
! in computation of "Vterm" in "Cr", before new hbl is found). Also
! compute turbulent  friction velocity "ustar" from wind stress at
! RHO-points. Finally, initialize boundary layer depth "hbl" and
! index "kbl" to the maximum (bottomed out) values.
!
# define alpha Bosol
# define beta  Bo
      call alfabeta_tile (istr,iend,jstr,jend, I_EXT_RANGE,
     &                             J_EXT_RANGE, alpha,beta)
      do j=J_EXT_RANGE
        do i=I_EXT_RANGE
          Bo(i,j)=g*( alpha(i,j)*(stflx(i,j,itemp)-srflx(i,j))
# ifdef SALINITY
     &                              -beta(i,j)*stflx(i,j,isalt)
# endif
     &                                                        )
          Bosol(i,j)=g*alpha(i,j)*srflx(i,j)
# undef beta
# undef alpha



# ifdef WND_AT_RHO_POINTS
          ustar(i,j)=sqrt(sqrt( sustr(i,j)**2+svstr(i,j)**2))
# else
          ustar(i,j)=sqrt(0.5*sqrt( (sustr(i,j)+sustr(i+1,j))**2
     &                             +(svstr(i,j)+svstr(i,j+1))**2))
# endif
          hbl(i,j)=hbls(i,j,tind) !<-- use old value as initial guess

          kbl(i)=0
# ifdef LIMIT_MO_DEPTH
#  define C_h_MO Gm1
          kmo(i)=0
          C_h_MO(i)=C_MO*ustar(i,j)**3/vonKar
# endif
          FC(i,N)=0.
          Cr(i,N)=0.
          Cr(i,0)=0.
        enddo  !--> discard alpha,beta; keep Bo,Bosol to the very end.

!
!======================================
!  Integral at W-points
!======================================
!

# if defined INT_AT_W_POINTS
!
!================================
!  Search for mixed layer depth
!================================
!
#  define du wrk1
#  define dv wrk2

      do k=1,N-1
        do i=I_EXT_RANGE
          cff=1./(Hz(i,j,k)+Hz(i,j,k+1))
          du(i,k)=cff*( u(i,j,k+1,tind)+u(i+1,j,k+1,tind)
     &                 -u(i,j,k  ,tind)-u(i+1,j,k  ,tind))
          dv(i,k)=cff*( v(i,j,k+1,tind)+v(i,j+1,k+1,tind)
     &                 -v(i,j,k  ,tind)-v(i,j+1,k  ,tind))
        enddo
      enddo
      do i=I_EXT_RANGE
        du(i,N)=du(i,N-1)
        dv(i,N)=dv(i,N-1)
        du(i,0)=du(i,  1)
        dv(i,0)=dv(i,  1)
      enddo

      do k=N,1,-1
        do i=I_EXT_RANGE
          zscale=z_w(i,j,N)-z_w(i,j,k-1)
          Kern=zscale/(zscale+epssfc*hbl(i,j))
          Bfsfc=Bo(i,j) +Bosol(i,j)*(1.-swr_frac(i,j,k-1))

#  include "lmd_wscale_ws_only.h"

          cff=bvf(i,j,k)*bvf(i,j,k-1)
          if (cff.gt.0.D0) then
            cff=cff/(bvf(i,j,k)+bvf(i,j,k-1))
          else
            cff=0.D0
          endif

          FC(i,k-1)=FC(i,k) + Kern*Hz(i,j,k)*( 0.375*( du(i,k)**2
     &                     +du(i,k-1)**2+dv(i,k)**2 +dv(i,k-1)**2 )
     &                 +0.25 *(du(i,k-1)*du(i,k)+dv(i,k-1)*dv(i,k))

     &              -Ri_inv*( cff + 0.25*(bvf(i,j,k)+bvf(i,j,k-1)))

     &                                        -C_Ek*f(i,j)*f(i,j) )

          Vtsq=Vtc*ws*sqrt(max(0., bvf(i,j,k-1)))

          Cr(i,k-1)=FC(i,k-1) +Vtsq
          if (kbl(i).eq.0 .and.  Cr(i,k-1).lt.0.) kbl(i)=k

#  ifdef LIMIT_MO_DEPTH
            if (kmo(i).eq.0 .and. Bfsfc*(z_w(i,j,N)-z_w(i,j,k-1))
     &                                  .gt. C_h_MO(i) ) kmo(i)=k
#  endif
        enddo
      enddo
#  undef dv
#  undef du

      do i=I_EXT_RANGE
c??     if (kbl(i).eq.N) then
c??       hbl(i,j)=z_w(i,j,N)-z_w(i,j,N-1)

        if (kbl(i).gt.0) then
          k=kbl(i)
          hbl(i,j)=z_w(i,j,N)-( z_w(i,j,k-1)*Cr(i,k)
     &                         -z_w(i,j,k)*Cr(i,k-1)
     &                         )/(Cr(i,k)-Cr(i,k-1))

c**          if (Cr(i,k)*Cr(i,k-1).gt.0.D0 ) write(*,*)
c**     &          '### ERROR', k, Cr(i,k), Cr(i,k-1), hbl(i,j)
          
        else
          hbl(i,j)=z_w(i,j,N)-z_w(i,j,0)+eps
        endif

#  ifdef LIMIT_MO_DEPTH
        if (kmo(i).gt.0) then
          k=kmo(i)
          cff_up=max(0., Bo(i,j)+Bosol(i,j)*(1.-swr_frac(i,j,k  )) )
          cff_dn=max(0., Bo(i,j)+Bosol(i,j)*(1.-swr_frac(i,j,k-1)) )

          h_MO=C_h_MO(i)*(cff_up*z_w(i,j,k) -cff_dn*z_w(i,j,k-1))
     &                    /( cff_dn*cff_up*(z_w(i,j,k)-z_w(i,j,k-1))
     &                                 + C_h_MO(i)*(cff_dn-cff_up) )

          hbl(i,j)=min(hbl(i,j), max(h_MO,0.))
        endif
#   undef C_h_MO
#  endif 
#  ifdef MASKING
        hbl(i,j)=hbl(i,j)*rmask(i,j)
#  endif
      enddo

# ifdef LMD_BKPP
!
!======================================
!  Search for bottom mixed layer depth
!======================================
!
        do i=I_EXT_RANGE
          kbl(i)  = 0       ! reset Cr at bottom and kbl for BKPP
          Cr(i,0) = 0.
        enddo
        do k=1,N,+1
          do i=I_EXT_RANGE
            Cr(i,k)=FC(i,k)-FC(i,0)
            if (kbl(i).eq.0 .and. Cr(i,k).gt.0.) kbl(i)=k
          enddo
        enddo
        do i=I_EXT_RANGE
          hbbl(i,j)=z_w(i,j,N)-z_w(i,j,0) !+eps ! total depth
          if (kbl(i).gt.0) then
            k=kbl(i)
            if (k.eq.1) then
              hbbl(i,j)=z_r(i,j,1)-z_w(i,j,0)   ! no BBL
            else
              hbbl(i,j)=( z_w(i,j,k-1)*Cr(i,k)-z_w(i,j,k)*Cr(i,k-1)
     &                            )/(Cr(i,k)-Cr(i,k-1) ) - z_w(i,j,0)
            endif
          endif
#  ifdef MASKING
          hbbl(i,j)=hbbl(i,j)*rmask(i,j)
#  endif
        enddo     !--> discard FC, Cr and kbl
# endif   /* LMD_BKPP */

!
!======================================
!  Integral at RHO-points
!======================================
!

# elif defined INT_AT_RHO_POINTS 
#  define swdk_r wrk1 
!
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
!/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\
!\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/
!
        do k=N-1,1,-1
          do i=I_EXT_RANGE
            zscale=z_w(i,j,N)-z_r(i,j,k)
            Kern=zscale/(zscale+epssfc*hbl(i,j))

            FC(i,k)=FC(i,k+1) + Kern*(
     &                 0.5*( ( u(i  ,j,k+1,tind)-u(i  ,j,k,tind)
     &                        +u(i+1,j,k+1,tind)-u(i+1,j,k,tind) )**2
     &                      +( v(i,j  ,k+1,tind)-v(i,j  ,k,tind)
     &                        +v(i,j+1,k+1,tind)-v(i,j+1,k,tind) )**2
     &                      )/(Hz(i,j,k)+Hz(i,j,k+1))

     &               -0.5*(Hz(i,j,k)+Hz(i,j,k+1))*( Ri_inv*bvf(i,j,k)

     &                                            +C_Ek*f(i,j)*f(i,j)
     &                                                             ))
          enddo
        enddo

        do k=N,1,-1
          do i=I_EXT_RANGE
            swdk_r(i,k)=sqrt(swr_frac(i,j,k)*swr_frac(i,j,k-1))

            zscale=z_w(i,j,N)-z_r(i,j,k)
            Bfsfc=Bo(i,j)+Bosol(i,j)*(1.-swdk_r(i,k))

#  include "lmd_wscale_ws_only.h"

c--         Vtsq=Vtc*ws*sqrt(max(0., 0.5*(bvf(i,j,k)+bvf(i,j,k-1)) ))

c fix by Gokhan

            Vtsq=Vtc*ws*sqrt(max(0., bvf(i,j,k-1) ))

            Cr(i,k)=FC(i,k)+Vtsq

            if (kbl(i).eq.0 .and. Cr(i,k).lt.0.) kbl(i)=k
#  ifdef LIMIT_MO_DEPTH
            if (kmo(i).eq.0 .and. Bfsfc*(z_w(i,j,N)-z_r(i,j,k))
     &                               .gt. C_h_MO(i) ) kmo(i)=k
#  endif
          enddo
        enddo   !<-- k, discard

        do i=I_EXT_RANGE
          hbl(i,j)=z_w(i,j,N)-z_w(i,j,0) +eps
          if (kbl(i).gt.0) then
            k=kbl(i)
            if (k.eq.N) then
              hbl(i,j)=z_w(i,j,N)-z_r(i,j,N)
            else
              hbl(i,j)=z_w(i,j,N)-( z_r(i,j,k)*Cr(i,k+1)
     &                              -z_r(i,j,k+1)*Cr(i,k)
     &                              )/(Cr(i,k+1)-Cr(i,k))
            endif
          endif
          
#  ifdef LIMIT_MO_DEPTH
          if (kmo(i).gt.0) then
            k=kmo(i)
            if (k.eq.N) then
              z_up=z_w(i,j,N)
              cff_up=max(0., Bo(i,j))
            else
              z_up=z_r(i,j,k+1)
              cff_up=max(0., Bo(i,j)+Bosol(i,j)*(1.-swdk_r(i,k+1)))
            endif
            cff_dn=max(0., Bo(i,j)+Bosol(i,j)*(1.-swdk_r(i,k)))

            h_MO=z_w(i,j,N)+C_h_MO(i)*(cff_up*z_up-cff_dn*z_r(i,j,k))
     &                             /( cff_up*cff_dn*(z_up-z_r(i,j,k))
     &                                  + C_h_MO(i)*(cff_dn-cff_up) )
            hbl(i,j)=min(hbl(i,j), max(h_MO,0.))
          endif
#   undef C_h_MO
#  endif /* LIMIT_MO_DEPTH */

#  ifdef MASKING
          hbl(i,j)=hbl(i,j)*rmask(i,j)
#  endif
        enddo     !--> discard kbl

# ifdef LMD_BKPP
!
!======================================
!  Search for bottom mixed layer depth
!======================================
!
        do i=I_EXT_RANGE
          kbl(i)=0       ! reset Cr at bottom and kbl for BKPP
          Cr(i,0)=0.D0
          FC(i,0)=1.5D0*FC(i,1)-0.5D0*FC(i,2)  ! linear extrapolation
        enddo
        do k=1,N,+1
          do i=I_EXT_RANGE
            Cr(i,k)=FC(i,k)-FC(i,0)

            if (kbl(i).eq.0 .and. Cr(i,k).gt.0.) kbl(i)=k
          enddo
        enddo
        do i=I_EXT_RANGE
          hbbl(i,j)=z_w(i,j,N)-z_w(i,j,0) !+eps ! total depth
          if (kbl(i).gt.0) then
            k=kbl(i)
            if (k.eq.1) then
              hbbl(i,j)=z_r(i,j,1)-z_w(i,j,0)   ! no BBL
            else
              hbbl(i,j)=( z_r(i,j,k-1)*Cr(i,k)-z_r(i,j,k)*Cr(i,k-1)
     &                            )/(Cr(i,k)-Cr(i,k-1)) -z_w(i,j,0)
            endif
          endif
#  ifdef MASKING
          hbbl(i,j)=hbbl(i,j)*rmask(i,j)
#  endif
        enddo     !--> discard FC, Cr and kbl
# endif   /* LMD_BKPP */


# else
      ERROR: Algorithm for Bulk Richardson number is not specified.
# endif 


!
!================================
!      Smoothing hbl/hbbl
!================================
!
# ifdef SMOOTH_HBL
      enddo  !<-- j terminate j-loop only if smothing takes place.
#  ifdef LMD_KPP
#   define wrk hbl
#   include "kpp_smooth.h"
#   undef wrk
#  endif
#  ifdef LMD_BKPP
#   define wrk hbbl
#   include "kpp_smooth.h"
#   undef wrk
#  endif
      do j=jstr,jend   !--> restart j-loop
# endif                     /* SMOOTH_HBL */
# undef I_EXT_RANGE
# undef J_EXT_RANGE
# define imin illegal
# define imax illegal
# define jmin illegal
# define jmax illegal


# if defined MERGE_OVERWRAP && (defined LMD_KPP && defined LMD_BKPP)
        do i=istr,iend
          cff=z_w(i,j,N)-z_w(i,j,0)
          if (hbl(i,j)+hbbl(i,j).ge.cff) then
            hbl(i,j)  = cff         ! if hbl and hbbl overwrap, replace
            hbbl(i,j) = cff         ! them with the total depth.
          endif
        enddo
# else
#  undef MERGE_OVERWRAP
# endif
!
!================================
!        Surface KPP
!================================
!
        do i=istr,iend
          kbl(i)=N             !<-- initialize search
        enddo
        do k=N-1,1,-1       ! find new boundary layer index "kbl".
          do i=istr,iend
            if (z_w(i,j,k) .gt. z_w(i,j,N)-hbl(i,j)) kbl(i)=k
          enddo
        enddo
!
! Find buoyancy forcing for final "hbl" values, and compute
! tubulent velocity scales (wm,ws) at "hbl".
! Then compute nondimensional shape function coefficients Gx( ) by
! matching values and vertical derivatives of  interior mixing
! coefficients at hbl (sigma=1).
!
        do i=istr,iend
          k=kbl(i)
          z_bl=z_w(i,j,N)-hbl(i,j)
          zscale=hbl(i,j)

          if (swr_frac(i,j,k-1).gt. 0.) then
            Bfsfc=Bo(i,j) +Bosol(i,j)*( 1. -swr_frac(i,j,k-1)
     &              *swr_frac(i,j,k)*(z_w(i,j,k)-z_w(i,j,k-1))
     &               /( swr_frac(i,j,k  )*(z_w(i,j,k)   -z_bl)
     &                 +swr_frac(i,j,k-1)*(z_bl -z_w(i,j,k-1))
     &                                                      ))
          else
            Bfsfc=Bo(i,j)+Bosol(i,j)
          endif

# include "lmd_wscale_wm_and_ws.h"

# ifdef LIMIT_UNSTABLE_ONLY
          f1=5.0 * max(0., Bfsfc) * vonKar/(ustar(i,j)**4+eps)
# else
          f1=0.
# endif
          cff=1./(z_w(i,j,k)-z_w(i,j,k-1))
          cff_up=cff*(z_bl -z_w(i,j,k-1))
          cff_dn=cff*(z_w(i,j,k)   -z_bl)

          Av_bl=cff_up*Kv(i,j,k)+cff_dn*Kv(i,j,k-1)
          dAv_bl=cff * (Kv(i,j,k)  -   Kv(i,j,k-1))
          Gm1(i)=Av_bl/(hbl(i,j)*wm+eps)
          dGm1dS(i)=min(0., Av_bl*f1-dAv_bl/(wm+eps))

          At_bl=cff_up*Kt(i,j,k)+cff_dn*Kt(i,j,k-1)
          dAt_bl=cff * (Kt(i,j,k)  -   Kt(i,j,k-1))
          Gt1(i)=At_bl/(hbl(i,j)*ws+eps)
          dGt1dS(i)=min(0., At_bl*f1-dAt_bl/(ws+eps))

# ifdef SALINITY
          As_bl=cff_up*Ks(i,j,k)+cff_dn*Ks(i,j,k-1)
          dAs_bl=cff * (Ks(i,j,k)  -   Ks(i,j,k-1))
          Gs1(i)=As_bl/(hbl(i,j)*ws+eps)
          dGs1dS(i)=min(0., As_bl*f1-dAs_bl/(ws+eps))
# endif
          Bfsfc_bl(i)=Bfsfc
        enddo

!
! Compute boundary layer mixing coefficients.
!--------- -------- ----- ------ -------------
! Compute turbulent velocity scales at vertical W-points.
!
        do i=istr,iend
          do k=N-1,kbl(i),-1
            Bfsfc=Bfsfc_bl(i)
            zscale=z_w(i,j,N)-z_w(i,j,k)

# include "lmd_wscale_wm_and_ws.h"
!
! Compute vertical mixing coefficients
!
            sigma=(z_w(i,j,N)-z_w(i,j,k))/max(hbl(i,j),eps)

            a1=sigma-2.
            a2=3.-2.*sigma
            a3=sigma-1.

            if (sigma.lt.0.07D0) then
              cff=0.5*(sigma-0.07D0)**2/0.07D0
            else
              cff=0.D0
            endif

            Kv(i,j,k)=wm*hbl(i,j)*( cff + sigma*( 1.+sigma*(
     &                           a1+a2*Gm1(i)+a3*dGm1dS(i) )))

            Kt(i,j,k)=ws*hbl(i,j)*( cff + sigma*( 1.+sigma*(
     &                           a1+a2*Gt1(i)+a3*dGt1dS(i) )))
# ifdef SALINITY
            Ks(i,j,k)=ws*hbl(i,j)*( cff + sigma*( 1.+sigma*(
     &                           a1+a2*Gs1(i)+a3*dGs1dS(i) ))) 
# endif
# ifdef LMD_NONLOCAL
            if (Bfsfc .lt. 0.) then
              ghat(i,j,k)=Cg * sigma*(1.-sigma)**2
            else
              ghat(i,j,k)=0.
            endif
# endif
          enddo
          do k=kbl(i)-1,1,-1
# ifdef LMD_NONLOCAL
            ghat(i,j,k)=0.
# endif 
# if  defined LMD_CONVEC && !defined LMD_BKPP
            if (bvf(i,j,k).lt.0.) then     ! Add convective
              Kv(i,j,k)=Kv(i,j,k) + nu0c   ! adjustment
              Kt(i,j,k)=Kt(i,j,k) + nu0c   ! if LMD_BKPP is defined, this will be done later !
#  ifdef SALINITY
              Ks(i,j,k)=Ks(i,j,k) + nu0c
#  endif
            endif
# endif
          enddo
        enddo
!
!================================
!        Bottom KPP
!================================
!
#if defined LMD_BKPP
        do i=istr,iend
          kbl(i)=N          !<-- initialize search
        enddo
        do k=N-1,1,-1       ! find new boundary layer index "kbl".
          do i=istr,iend
            if (z_r(i,j,k)-z_w(i,j,0).gt.hbbl(i,j)) kbl(i)=k
          enddo
        enddo
!
! Compute nondimensional shape function coefficients Gx( ) by
! matching values and vertical derivatives of interior mixing
! coefficients at hbbl (sigma=1).
!
        do i=istr,iend
          ustar2 = r_D(i,j)*sqrt(
     &             ( (0.5*(u(i,j,1,tind)+u(i+1,j,1,tind)))**2
     &              +(0.5*(v(i,j,1,tind)+v(i,j+1,1,tind)))**2 ) )
          wm=vonKar*sqrt(ustar2)  ! turbulent velocity scales (wm,ws)
          ws=wm                   ! with buoyancy effects neglected.
          k=kbl(i)
          z_bl=z_w(i,j,0)+hbbl(i,j)
          if (z_bl.lt.z_w(i,j,k-1)) k=k-1

          cff=1./(z_w(i,j,k)-z_w(i,j,k-1))
          cff_up=cff*(z_bl -z_w(i,j,k-1))
          cff_dn=cff*(z_w(i,j,k)   -z_bl)

          Av_bl=cff_up*Kv(i,j,k)+cff_dn*Kv(i,j,k-1)
          dAv_bl=cff * (Kv(i,j,k)  -   Kv(i,j,k-1))
          Gm1(i)=Av_bl/(hbbl(i,j)*wm+eps)
          dGm1dS(i)=min(0., -dAv_bl/(wm+eps))

          At_bl=cff_up*Kt(i,j,k)+cff_dn*Kt(i,j,k-1)
          dAt_bl=cff * (Kt(i,j,k)  -   Kt(i,j,k-1))
          Gt1(i)=At_bl/(hbbl(i,j)*ws+eps)
          dGt1dS(i)=min(0., -dAt_bl/(ws+eps))

#  ifdef SALINITY
          As_bl=cff_up*Ks(i,j,k)+cff_dn*Ks(i,j,k-1)
          dAs_bl=cff * (Ks(i,j,k)  -   Ks(i,j,k-1))
          Gs1(i)=As_bl/(hbbl(i,j)*ws+eps)
          dGs1dS(i)=min(0., -dAs_bl/(ws+eps))
#  endif
!
! Compute boundary layer mixing coefficients.
!--------- -------- ----- ------ -------------
! Compute turbulent velocity scales at vertical W-points.
!
          do k=1,N-1
            if (k.lt.kbl(i)) then
              sigma=min((z_w(i,j,k)-z_w(i,j,0)+Zob)/(hbbl(i,j)+Zob),1.)
              a1=sigma-2.
              a2=3.-2.*sigma
              a3=sigma-1.
              Kv0 =wm*hbbl(i,j)*( sigma*( 1.+sigma*(
     &                            a1+a2*Gm1(i)+a3*dGm1dS(i) )))
              Kt0 =ws*hbbl(i,j)*( sigma*( 1.+sigma*(
     &                            a1+a2*Gt1(i)+a3*dGt1dS(i) )))
#  ifdef SALINITY
              Ks0 =ws*hbbl(i,j)*( sigma*( 1.+sigma*(
     &                            a1+a2*Gs1(i)+a3*dGs1dS(i) )))
#  endif
!
! If BBL reaches into SBL, take the max of surface and bottom values.
!
              z_bl=z_w(i,j,N)-hbl(i,j)
              if (z_w(i,j,k).gt.z_bl) then
                Kv0=max(Kv(i,j,k),Kv0)
                Kt0=max(Kt(i,j,k),Kt0)
#  ifdef SALINITY
                Ks0=max(Ks(i,j,k),Ks0)
#  endif
              endif
              Kv(i,j,k)=Kv0
              Kt(i,j,k)=Kt0
#  ifdef SALINITY
              Ks(i,j,k)=Ks0
#  endif

#  ifdef LMD_CONVEC
            else  !<-- k > kbl(i)
              if (bvf(i,j,k).lt.0.) then
#   ifdef LMD_KPP
                z_bl=z_w(i,j,N)-hbl(i,j)
                if (z_w(i,j,k).lt.z_bl) then
#   endif
                  Kv(i,j,k)=Kv(i,j,k) + nu0c   ! Add convective
                  Kt(i,j,k)=Kt(i,j,k) + nu0c   ! adjustment outside
#   ifdef SALINITY
                  Ks(i,j,k)=Ks(i,j,k) + nu0c   ! of mixed layers.
#   endif
#   ifdef LMD_KPP
                endif
#   endif
              endif
#  endif
            endif  !<-- k < kbl(i)
          enddo    !<-- k
        enddo      !<-- i
# endif   /* LMD_BKPP */
!
!================================
!         Finalize
!================================
!
        do k=1,N-1                       ! Finalize: Copy into shared arrays
          do i=istr,iend
            Akv(i,j,k)=Kv(i,j,k)
            Akt(i,j,k,itemp)=Kt(i,j,k)
 
# ifdef SALINITY
            Akt(i,j,k,isalt)=Ks(i,j,k)
# endif
          enddo
        enddo
        
      enddo    !<-- j
      
      
      
      
      

      do j=jstr,jend                 ! Copy "hbl" into its shared
        do i=istr,iend               ! array and pad out ghost points
          hbls(i,j,3-nstp)=hbl(i,j)  ! at lateral-sideboundaries. 
#ifdef LMD_BKPP          
          hbbls(i,j,3-nstp)=hbbl(i,j)
#endif        
        enddo
      enddo
# ifndef EW_PERIODIC
      if (WESTERN_EDGE) then
        do j=jstr,jend
          hbls(istr-1,j,3-nstp)=hbls(istr,j,3-nstp)
#ifdef LMD_BKPP
          hbbls(istr-1,j,3-nstp)=hbbls(istr,j,3-nstp)
#endif          
        enddo
      endif
      if (EASTERN_EDGE) then
        do j=jstr,jend
          hbls(iend+1,j,3-nstp)=hbls(iend,j,3-nstp)
#ifdef LMD_BKPP
          hbbls(iend+1,j,3-nstp)=hbbls(iend,j,3-nstp)
#endif           
        enddo
      endif
# endif
# ifndef NS_PERIODIC
      if (SOUTHERN_EDGE) then
        do i=istr,iend
          hbls(i,jstr-1,3-nstp)=hbls(i,jstr,3-nstp)
#ifdef LMD_BKPP
          hbbls(i,jstr-1,3-nstp)=hbbls(i,jstr,3-nstp)
#endif           
        enddo
      endif
      if (NORTHERN_EDGE) then
        do i=istr,iend
          hbls(i,jend+1,3-nstp)=hbls(i,jend,3-nstp)
#ifdef LMD_BKPP
          hbbls(i,jend+1,3-nstp)=hbbls(i,jend,3-nstp)
#endif                     
        enddo
      endif
#  ifndef EW_PERIODIC
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        hbls(istr-1,jstr-1,3-nstp)=hbls(istr,jstr,3-nstp)
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        hbls(istr-1,jend+1,3-nstp)=hbls(istr,jend,3-nstp)
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        hbls(iend+1,jstr-1,3-nstp)=hbls(iend,jstr,3-nstp)
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        hbls(iend+1,jend+1,3-nstp)=hbls(iend,jend,3-nstp)
      endif
#ifdef LMD_BKPP
      if (WESTERN_EDGE .and. SOUTHERN_EDGE) then
        hbbls(istr-1,jstr-1,3-nstp)=hbbls(istr,jstr,3-nstp)
      endif
      if (WESTERN_EDGE .and. NORTHERN_EDGE) then
        hbbls(istr-1,jend+1,3-nstp)=hbbls(istr,jend,3-nstp)
      endif
      if (EASTERN_EDGE .and. SOUTHERN_EDGE) then
        hbbls(iend+1,jstr-1,3-nstp)=hbbls(iend,jstr,3-nstp)
      endif
      if (EASTERN_EDGE .and. NORTHERN_EDGE) then
        hbbls(iend+1,jend+1,3-nstp)=hbbls(iend,jend,3-nstp)
      endif
#endif      
#  endif
# endif

# if defined EW_PERIODIC || defined NS_PERIODIC || defined MPI
      call exchange_r2d_tile (istr,iend,jstr,jend,
     &                     hbls(START_2D_ARRAY,3-nstp))
# if defined LMD_BKPP
      call exchange_r2d_tile (istr,iend,jstr,jend,
     &                    hbbls(START_2D_ARRAY,3-nstp))
# endif     
      call exchange_w3d_tile (istr,iend,jstr,jend, Akv)
      call exchange_w3d_tile (istr,iend,jstr,jend,
     &                     Akt(START_2D_ARRAY,0,itemp))
      call exchange_w3d_tile (istr,iend,jstr,jend,
     &                     Akt(START_2D_ARRAY,0,isalt))

#endif
      return
      end


      subroutine check_kpp_switches (ierr)
!
! The purpose of this routine is to keep track of status of
! CPP-switches locally defined in this file. It must be placed
! here rather than stored in a separate file in order to be
! exposed to relevant CPP-settings. This routine does not affect
! any model results, other than signature in global attribute
! "CPPS" in output netCDF files.
!
      implicit none
      integer ierr, is,ie, lenstr
# include "param.h"
# include "strings.h"
      ie=lenstr(cpps)
      is=ie+2
      ie=is+10
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='<lmd_kpp.F>'
      MPI_master_only write(stdout,'(1x,A)') cpps(is:ie) 
      is=ie+2
# ifdef BULK_DIFFERENCE
      ie=is+14
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='BULK_DIFFERENCE'
      MPI_master_only write(stdout,'(10x,A)') cpps(is:ie)
      is=ie+2
# endif
# ifdef INT_AT_W_POINTS
      ie=is+14
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='INT_AT_W_POINTS'
      MPI_master_only write(stdout,'(10x,A)') cpps(is:ie)
      is=ie+2
# endif
# ifdef INT_AT_RHO_POINTS
      ie=is+16
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='INT_AT_RHO_POINTS'
      MPI_master_only write(stdout,'(10x,A)') cpps(is:ie)
      is=ie+2
# endif
# ifdef SMOOTH_HBL
      ie=is+9
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='SMOOTH_HBL'
      MPI_master_only write(stdout,'(10x,A)') cpps(is:ie)
      is=ie+2
# endif
# ifdef LIMIT_UNSTABLE_ONLY
      ie=is+18
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='LIMIT_UNSTABLE_ONLY'
      MPI_master_only write(stdout,'(10x,A)') cpps(is:ie)
      is=ie+2
# endif
# ifdef LIMIT_MO_DEPTH
      ie=is+13
      if (ie.gt.max_opt_size) goto 99
      cpps(is:ie)='LIMIT_MO_DEPTH'
      MPI_master_only write(stdout,'(10x,A)') cpps(is:ie)
      is=ie+2
# endif
      return
  99  MPI_master_only write(stdout,'(/1x,2A/12x,A/)') '### ERROR: ',
     &  'Unsufficient lenght of string "cpps" in file "strings.h".',
     &        'Increase parameter "max_opt_size" it and recompile.'
      ierr=ierr+1
      return
      end
#else
      subroutine lmd_kpp_empty
      end
#endif /* LMD_KPP */
 
