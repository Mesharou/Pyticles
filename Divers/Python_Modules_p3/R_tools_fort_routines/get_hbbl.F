

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Subpart of the lmd_kpp.F routine (myroms version)
! used to compute the new hbl
! (the part used to compute the new Kv, Kt has been removed)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

#define LMD_BKPP 

#define INT_AT_W_POINTS
#define SMOOTH_HBL
#define LIMIT_UNSTABLE_ONLY
c----#define WND_AT_RHO_POINTS

      subroutine get_hbbl (Lm,Mm,N,alpha,beta, z_r,z_w
     &  , stflx, srflx, swr_frac, sustr, svstr ,Ricr, hbls, f
     &  , u, v, bvf
     &  , hbbl, out1, out2, out3, out4) 

      implicit none

      integer Lm,Mm,N,NT, i,j,k
     &      ,istr,iend,jstr,jend
     &      ,itemp,isalt

      real epsil

      real nubl, nu0c, Cv, Ricr, Ri_inv, betaT, epssfc, C_Ek, C_MO,
     &     Cstar, Cg, eps, zeta_m, a_m, c_m,  zeta_s, a_s, c_s,
     &                                                    r2,r3,r4
      parameter (nubl=0.01,
     &   nu0c=0.1,Cv=1.8,
     &   betaT=-0.2,epssfc=0.1,C_MO=1.,C_Ek=258.,
     &   Cstar=10.,eps=1.E-20,zeta_m=-0.2,a_m=1.257,
     &   c_m=8.360,zeta_s=-1.0,a_s=-28.86,c_s=98.96,
     &   r2=0.5, r3=1./3., r4=0.25)

#include "scalars.h"

      parameter (NT=2)  
      parameter (itemp=1,isalt=2)  
      parameter (epsil=1.E-16)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      real*8 sigma
     &     ,Bo(0:Lm+1,0:Mm+1), Bosol(0:Lm+1,0:Mm+1)
     &     ,Bfsfc_bl(0:Lm+1)
     &     ,z_bl
     &     ,ustar(0:Lm+1,0:Mm+1)
     &     ,Cr(0:Lm+1,0:N)
     &     ,FC(0:Lm+1,0:N)
     &     ,wrk1(0:Lm+1,0:N)
     &     ,wrk2(0:Lm+1,0:N)
     &     ,cff, cff1
     &     ,Hz(0:Lm+1,0:Mm+1,N)
# ifdef SMOOTH_HBL
     &     ,FX(0:Lm+1,0:Mm+1)
     &     ,FE(0:Lm+1,0:Mm+1)
     &     ,FE1(0:Lm+1,0:Mm+1)
# else
     &         ,FX, FE, FE1
# endif

      integer kbl(0:Lm+1)

      real Kern, Vtc, Vtsq
     &     , Bfsfc,zscale
     &     , ustar3, zetahat, ws, wm


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      real*8 hbl(0:Lm+1,0:Mm+1)

      ! Variables IN
      real*8 stflx(0:Lm+1,0:Mm+1,NT)
     &      ,srflx(0:Lm+1,0:Mm+1)
     &      ,swr_frac(0:Lm+1,0:Mm+1,0:N)
     &      ,alpha(0:Lm+1,0:Mm+1), beta(0:Lm+1,0:Mm+1)
     &      ,z_r(0:Lm+1,0:Mm+1,N), z_w(0:Lm+1,0:Mm+1,0:N)
     &      ,sustr(1:Lm+1,0:Mm+1), svstr(0:Lm+1,1:Mm+1)
     &      ,hbls(0:Lm+1,0:Mm+1), f(0:Lm+1,0:Mm+1)
     &      ,u(1:Lm+1,0:Mm+1,N), v(0:Lm+1,1:Mm+1,N)
     &      ,bvf(0:Lm+1,0:Mm+1,0:N)

      ! Variables OUT
      real*8 hbbl(0:Lm+1,0:Mm+1)
     &      , out1(0:Lm+1,0:Mm+1,0:N), out2(0:Lm+1,0:Mm+1,0:N)
     &      , out3(0:Lm+1,0:Mm+1,0:N), out4(0:Lm+1,0:Mm+1,0:N)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


Cf2py intent(in) Lm,Mm,N,alpha,beta ,z_r,z_w,stflx,srflx,  swr_frac, sustr, svstr ,Ricr,hbls, f, u, v, bvf
Cf2py intent(out) hbbl, out1, out2, out3, out4

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
      Ri_inv=1./Ricr
      Cg=Cstar * vonKar * (c_s*vonKar*epssfc)**(1./3.)
      Vtc=Cv * sqrt(-betaT/(c_s*epssfc)) / (Ricr*vonKar**2)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


        istr=0+1
        iend=Lm
        jstr=0+1
        jend=Mm


#  define I_EXT_RANGE istr,iend
#  define J_EXT_RANGE jstr,jend



!================================

      do j=J_EXT_RANGE
        do i=I_EXT_RANGE
            
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

!!          hbl(i,j)=0. !<-- use old value as initial guess
          kbl(i)=0


# ifdef LIMIT_MO_DEPTH
#  define C_h_MO Gm1
          kmo(i)=0
          C_h_MO(i)=C_MO*ustar(i,j)**3/vonKar
# endif
          FC(i,N)=0.
          Cr(i,N)=0.
          Cr(i,0)=0.

          out1(i,j,N)=0.
          out2(i,j,N)=0.
          out3(i,j,N)=0.
          out4(i,j,N)=0.


        enddo 

!======================================
!  Integral at W-points
!======================================



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
          du(i,k)=cff*( u(i,j,k+1)+u(i+1,j,k+1)
     &                 -u(i,j,k  )-u(i+1,j,k  ))
          dv(i,k)=cff*( v(i,j,k+1)+v(i,j+1,k+1)
     &                 -v(i,j,k  )-v(i,j+1,k  ))
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          out1(i,j,k-1) = out1(i,j,k)+ Kern*Hz(i,j,k)*(
     &                     0.375*( du(i,k)**2 + du(i,k-1)**2
     &                            +dv(i,k)**2 + dv(i,k-1)**2 )
     &              +0.25 *(du(i,k-1)*du(i,k)+dv(i,k-1)*dv(i,k))
     &                                                )

          out2(i,j,k-1) = out2(i,j,k)+ Kern*Hz(i,j,k)*(
     &              -Ri_inv*( cff + 0.25*(bvf(i,j,k)+bvf(i,j,k-1)))
     &                 )

          out3(i,j,k-1) = out3(i,j,k)+ Kern*Hz(i,j,k)*(
     &                                        -C_Ek*f(i,j)*f(i,j)
     &                 )

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          FC(i,k-1)=FC(i,k) + Kern * Hz(i,j,k) * (
     &                       0.375*( du(i,k  )**2 +du(i,k-1)**2
     &                              +dv(i,k  )**2 +dv(i,k-1)**2  )
     &                      +0.25 *( du(i,k-1) * du(i,k)
     &                              +dv(i,k-1) * dv(i,k)         )
     &              -Ri_inv*( cff + 0.25*( bvf(i,j,k)+bvf(i,j,k-1)) )
     &              -C_Ek*f(i,j)*f(i,j)
     &                                           )


          Vtsq=Vtc*ws*sqrt(max(0., bvf(i,j,k-1)))

          out4(i,j,k-1) = Vtsq

          Cr(i,k)=FC(i,k)+Vtsq

        enddo
      enddo
#  undef dv
#  undef du





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
        enddo     !--> discard FC, Cr and kbl

# endif   /* LMD_BKPP */


!======================================

        
      enddo    !<-- j


# ifdef SMOOTH_HBL
#  ifdef LMD_KPP
#   define wrk hbbl
#   include "kpp_smooth.h"
#   undef wrk
#  endif
# endif 

!======================================


      return
      end

