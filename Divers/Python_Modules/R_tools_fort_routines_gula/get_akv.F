

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subpart of the lmd_kpp.F routine (myroms version)
! used to compute the Kv
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


#define INT_AT_W_POINTS
#define SMOOTH_HBL
#define LIMIT_UNSTABLE_ONLY
c----#define WND_AT_RHO_POINTS

      subroutine get_akv (Lm,Mm,N,alpha,beta, z_r,z_w
     &  , stflx, srflx, swr_frac, sustr, svstr ,Ricr, hbls, f
     &  , u, v, bvf,  rmask, r_D
     &  , Kv) 

      implicit none

      integer Lm,Mm,N,NT, i,j,k
     &      ,itemp,isalt

      integer imin,imax,jmin,jmax

      real epsil

#include "scalars.h"

      parameter (NT=2)  
      parameter (itemp=1,isalt=2)  
      parameter (epsil=1.E-16)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      real ustar3,
     &     Bfsfc,   zscale,
     &     zetahat, ws,wm,
# ifdef LMD_BKPP
     &     zscaleb
# endif

      real*8 Bo(0:Lm+1,0:Mm+1), Bosol(0:Lm+1,0:Mm+1)
     &     ,Bfsfc_bl(0:Lm+1)
     &     ,z_bl
     &     ,ustar(0:Lm+1,0:Mm+1)
     &     ,Cr(0:Lm+1,0:N)
     &     ,FC(0:Lm+1,0:N)
     &     ,wrk1(0:Lm+1,0:N)
     &     ,wrk2(0:Lm+1,0:N)
     &     ,Hz(0:Lm+1,0:Mm+1,N)
# ifdef SMOOTH_HBL
     &     ,FX(0:Lm+1,0:Mm+1)
     &     ,FE(0:Lm+1,0:Mm+1)
     &     ,FE1(0:Lm+1,0:Mm+1)
# else
     &         ,FX, FE, FE1
# endif
     &    ,Gm1(0:Lm+1),   Av_bl,
     & dGm1dS(0:Lm+1),  dAv_bl, f1,
     &    Gt1(0:Lm+1),   At_bl, a1, 
     & dGt1dS(0:Lm+1),  dAt_bl, a2,
     &    Gs1(0:Lm+1),   As_bl, a3,
     & dGs1dS(0:Lm+1),  dAs_bl

      integer kbl(0:Lm+1)

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

      parameter (nubl=0.01,
     &   nu0c=0.1,Cv=1.8,
     &   betaT=-0.2,epssfc=0.1,C_MO=1.,C_Ek=258.,
     &   Cstar=10.,eps=1.E-20,zeta_m=-0.2,a_m=1.257,
     &   c_m=8.360,zeta_s=-1.0,a_s=-28.86,c_s=98.96,
     &   r2=0.5, r3=1./3., r4=0.25)



#ifdef LMD_BKPP
      real  ustar2, Kv0, Kt0, Ks0, my_Akv_bak, 
     &                        my_Akt_bak, my_Aks_bak
      real*8 hbbl(0:Lm+1,0:Mm+1)
#endif

      real*8 hbl(0:Lm+1,0:Mm+1)



      real*8 rdrg, Zob


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      real*8 ghat(0:Lm+1,0:Mm+1,N)
     &      ,r_D(0:Lm+1,0:Mm+1)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      ! Variables IN
      real*8 stflx(0:Lm+1,0:Mm+1,NT)
     &      ,srflx(0:Lm+1,0:Mm+1)
     &      ,swr_frac(0:Lm+1,0:Mm+1,0:N)
     &      ,alpha(0:Lm+1,0:Mm+1), beta(0:Lm+1,0:Mm+1)
     &      ,z_r(0:Lm+1,0:Mm+1,N), z_w(0:Lm+1,0:Mm+1,0:N)
     &      ,sustr(1:Lm+1,0:Mm+1), svstr(0:Lm+1,1:Mm+1)
     &      ,hbls(0:Lm+1,0:Mm+1), f(0:Lm+1,0:Mm+1)
     &      ,u(1:Lm+1,0:Mm+1,N), v(0:Lm+1,1:Mm+1,N)
     &      ,bvf(0:Lm+1,0:Mm+1,0:N),rmask(0:Lm+1,0:Mm+1)

      ! Variables OUT
      real*8 Kv(0:Lm+1,0:Mm+1,0:N)

      real*8 Kt(0:Lm+1,0:Mm+1,0:N)
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!# include "compute_tile_bounds.h"
# include "compute_auxiliary_bounds.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


Cf2py intent(in) Lm,Mm,N,alpha,beta ,z_r,z_w,stflx,srflx,  swr_frac, sustr, svstr ,Ricr,hbls, f, u, v, bvf, rmask, r_D
Cf2py intent(out) Kv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
 
      Ri_inv=1./Ricr


      Cg=Cstar * vonKar * (c_s*vonKar*epssfc)**(1./3.)
      Vtc=Cv * sqrt(-betaT/(c_s*epssfc)) / (Ricr*vonKar**2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

# ifdef SMOOTH_HBL
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
#  define I_EXT_RANGE_akv imin,imax
#  define J_EXT_RANGE_akv jmin,jmax
# else
#  define I_EXT_RANGE_akv istr,iend
#  define J_EXT_RANGE_akv jstr,jend
# endif




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

       CALL lmd_vmix (Lm,Mm,N,u,v,z_r
     &  ,bvf,Kv,Kt)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1




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

!       call alfabeta_tile (istr,iend,jstr,jend, I_EXT_RANGE_akv,
!      &                             J_EXT_RANGE_akv, alpha,beta)
      do j=J_EXT_RANGE_akv
        do i=I_EXT_RANGE_akv
            
          do k=1,N,+1
           Hz(i,j,k) = z_w(i,j,k) - z_w(i,j,k-1)
          enddo


          Bo(i,j)=g*( alpha(i,j)*(stflx(i,j,itemp)-srflx(i,j))
# ifdef SALINITY
     &                              -beta(i,j)*stflx(i,j,isalt)
# endif
     &                                                        )
          Bosol(i,j)=g*alpha(i,j)*srflx(i,j)


# ifdef WND_AT_RHO_POINTS
          ustar(i,j)=sqrt(sqrt( sustr(i,j)**2+svstr(i,j)**2))
# else
          ustar(i,j)=sqrt(0.5*sqrt( (sustr(i,j)+sustr(i+1,j))**2
     &                             +(svstr(i,j)+svstr(i,j+1))**2))
# endif
          hbl(i,j)=hbls(i,j) !<-- use old value as initial guess

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





!================================
!  Search for mixed layer depth
!================================
!
#  define du wrk1
#  define dv wrk2



      do k=1,N-1
        do i=I_EXT_RANGE_akv
          cff=1./(Hz(i,j,k)+Hz(i,j,k+1))
          du(i,k)=cff*( u(i,j,k+1)+u(i+1,j,k+1)
     &                 -u(i,j,k  )-u(i+1,j,k  ))
          dv(i,k)=cff*( v(i,j,k+1)+v(i,j+1,k+1)
     &                 -v(i,j,k  )-v(i,j+1,k  ))
        enddo
      enddo

      do i=I_EXT_RANGE_akv
        du(i,N)=du(i,N-1)
        dv(i,N)=dv(i,N-1)
        du(i,0)=du(i,  1)
        dv(i,0)=dv(i,  1)
      enddo


      do k=N,1,-1
        do i=I_EXT_RANGE_akv
          zscale=z_w(i,j,N)-z_w(i,j,k-1)
!          zscaleb=z_w(i,j,k)-z_w(i,j,0)
          Kern=zscale/(zscale+epssfc*hbl(i,j))
! # ifdef LMD_BKPP
!      &        *zscaleb/(zscaleb+epssfc*hbbls(i,j))
! #  endif
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



      do i=I_EXT_RANGE_akv
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
        do i=I_EXT_RANGE_akv
          kbl(i)  = 0       ! reset Cr at bottom and kbl for BKPP
          Cr(i,0) = 0.
        enddo
        do k=1,N,+1
          do i=I_EXT_RANGE_akv
            Cr(i,k)=FC(i,k)-FC(i,0)
            if (kbl(i).eq.0 .and. Cr(i,k).gt.0.) kbl(i)=k
          enddo
        enddo
        do i=I_EXT_RANGE_akv
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

!======================================



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

!           At_bl=cff_up*Kt(i,j,k)+cff_dn*Kt(i,j,k-1)
!           dAt_bl=cff * (Kt(i,j,k)  -   Kt(i,j,k-1))
!           Gt1(i)=At_bl/(hbl(i,j)*ws+eps)
!           dGt1dS(i)=min(0., At_bl*f1-dAt_bl/(ws+eps))
! 
! # ifdef SALINITY
!           As_bl=cff_up*Ks(i,j,k)+cff_dn*Ks(i,j,k-1)
!           dAs_bl=cff * (Ks(i,j,k)  -   Ks(i,j,k-1))
!           Gs1(i)=As_bl/(hbl(i,j)*ws+eps)
!           dGs1dS(i)=min(0., As_bl*f1-dAs_bl/(ws+eps))
! # endif
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

!             Kt(i,j,k)=ws*hbl(i,j)*( cff + sigma*( 1.+sigma*(
!      &                           a1+a2*Gt1(i)+a3*dGt1dS(i) )))
! # ifdef SALINITY
!             Ks(i,j,k)=ws*hbl(i,j)*( cff + sigma*( 1.+sigma*(
!      &                           a1+a2*Gs1(i)+a3*dGs1dS(i) ))) 
! # endif
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
!               Kt(i,j,k)=Kt(i,j,k) + nu0c   ! if LMD_BKPP is defined, this will be done later !
! #  ifdef SALINITY
!               Ks(i,j,k)=Ks(i,j,k) + nu0c
! #  endif
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



!================================

        Zob=0.01

! Compute nondimensional shape function coefficients Gx( ) by
! matching values and vertical derivatives of interior mixing
! coefficients at hbbl (sigma=1).
!
        do i=istr,iend
          ustar2 = r_D(i,j)*sqrt(
     &             ( (0.5*(u(i,j,1)+u(i+1,j,1)))**2
     &              +(0.5*(v(i,j,1)+v(i,j+1,1)))**2 ) )
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
! 
!           At_bl=cff_up*Kt(i,j,k)+cff_dn*Kt(i,j,k-1)
!           dAt_bl=cff * (Kt(i,j,k)  -   Kt(i,j,k-1))
!           Gt1(i)=At_bl/(hbbl(i,j)*ws+eps)
!           dGt1dS(i)=min(0., -dAt_bl/(ws+eps))
! 
! #  ifdef SALINITY
!           As_bl=cff_up*Ks(i,j,k)+cff_dn*Ks(i,j,k-1)
!           dAs_bl=cff * (Ks(i,j,k)  -   Ks(i,j,k-1))
!           Gs1(i)=As_bl/(hbbl(i,j)*ws+eps)
!           dGs1dS(i)=min(0., -dAs_bl/(ws+eps))
! #  endif
!
! Compute boundary layer mixing coefficients.
!--------- -------- ----- ------ -------------
! Compute turbulent velocity scales at vertical W-points.

          do k=1,N-1
            if (k.lt.kbl(i)) then
              sigma=min((z_w(i,j,k)-z_w(i,j,0)+Zob)/(hbbl(i,j)+Zob),1.)
              a1=sigma-2.
              a2=3.-2.*sigma
              a3=sigma-1.


              if (sigma.lt.0.07D0) then
                cff=0.5*(sigma-0.07D0)**2/0.07D0
              else
                cff=0.D0
              endif


              Kv0 =wm*hbbl(i,j)*(cff + sigma*( 1.+sigma*(
     &                            a1+a2*Gm1(i)+a3*dGm1dS(i) )))
!               Kt0 =ws*hbbl(i,j)*(cff + sigma*( 1.+sigma*(
!      &                            a1+a2*Gt1(i)+a3*dGt1dS(i) )))
! #  ifdef SALINITY
!               Ks0 =ws*hbbl(i,j)*(cff + sigma*( 1.+sigma*(
!      &                            a1+a2*Gs1(i)+a3*dGs1dS(i) )))
! #  endif
!
! If BBL reaches into SBL, take the max of surface and bottom values.
!
              z_bl=z_w(i,j,N)-hbl(i,j)
              if (z_w(i,j,k).gt.z_bl) then
                Kv0=max(Kv(i,j,k),Kv0)
!                 Kt0=max(Kt(i,j,k),Kt0)
! #  ifdef SALINITY
!                 Ks0=max(Ks(i,j,k),Ks0)
! #  endif
              endif
              Kv(i,j,k)=Kv0
!               Kt(i,j,k)=Kt0
! #  ifdef SALINITY
!               Ks(i,j,k)=Ks0
! #  endif

#  ifdef LMD_CONVEC
            else  !<-- k > kbl(i)
              if (bvf(i,j,k).lt.0.) then
#   ifdef LMD_KPP
                z_bl=z_w(i,j,N)-hbl(i,j)
                if (z_w(i,j,k).lt.z_bl) then
#   endif
                  Kv(i,j,k)=Kv(i,j,k) + nu0c   ! Add convective
!                   Kt(i,j,k)=Kt(i,j,k) + nu0c   ! adjustment outside
! #   ifdef SALINITY
!                   Ks(i,j,k)=Ks(i,j,k) + nu0c   ! of mixed layers.
! #   endif
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

        
      enddo    !<-- j
      
      

!======================================


      return
      end

