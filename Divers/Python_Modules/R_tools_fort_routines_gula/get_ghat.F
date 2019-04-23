
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      subroutine get_ghat (Lm,Mm,N,alpha,beta, z_r,z_w
     &  ,stflx, srflx, hbl, swr_frac
     &  ,ghat) 

      implicit none

      integer Lm,Mm,N,NT, i,j,k
     &      ,imin,imax,jmin,jmax
     &      ,itemp,isalt

      real epsil

      real nubl, nu0c, Cv, Ricr, Ri_inv, betaT, epssfc, C_Ek, C_MO,
     &     Cstar, Cg, eps, zeta_m, a_m, c_m,  zeta_s, a_s, c_s,
     &                                                    r2,r3,r4
      parameter (nubl=0.01,
     &   nu0c=0.1,Cv=1.8,Ricr=0.45,Ri_inv=1./Ricr,
     &   betaT=-0.2,epssfc=0.1,C_MO=1.,C_Ek=258.,
     &   Cstar=10.,eps=1.E-20,zeta_m=-0.2,a_m=1.257,
     &   c_m=8.360,zeta_s=-1.0,a_s=-28.86,c_s=98.96,
     &   r2=0.5, r3=1./3., r4=0.25)

#include "scalars.h"

      parameter (NT=2)  
      parameter (itemp=1,isalt=2)  
      parameter (epsil=1.E-16)


      real*8 sigma, Bfsfc
     &     ,Bo(0:Lm+1,0:Mm+1), Bosol(0:Lm+1,0:Mm+1)
     &     , Bfsfc_bl(0:Lm+1)
     &     ,hbl(0:Lm+1,0:Mm+1), z_bl

      integer kbl(0:Lm+1)

      real*8 stflx(0:Lm+1,0:Mm+1,NT)
     &      ,srflx(0:Lm+1,0:Mm+1)
     &      ,ghat(0:Lm+1,0:Mm+1,N)
     &      ,swr_frac(0:Lm+1,0:Mm+1,0:N)
     &      ,alpha(0:Lm+1,0:Mm+1), beta(0:Lm+1,0:Mm+1)
     &      ,z_r(0:Lm+1,0:Mm+1,N), z_w(0:Lm+1,0:Mm+1,0:N)
!      &      ,Hz(0:Lm+1,0:Mm+1,N)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!# include "compute_tile_bounds.h"
# include "compute_auxiliary_bounds.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


Cf2py intent(in) Lm,Mm,N,alpha,beta ,z_r,z_w,stflx,srflx, hbl, swr_frac
Cf2py intent(out) ghat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 

      Cg=Cstar * vonKar * (c_s*vonKar*epssfc)**(1./3.)



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1



!       do j=jmin,jmax
!         do i=imin,imax
!           do k=1,N,+1
!            Hz(i,j,k) = z_w(i,j,k) - z_w(i,j,k-1)
!           enddo
!         enddo
!       enddo


       !!call get_swr_frac (Lm,Mm,N, Hz, swr_frac )


      do j=jmin,jmax
        do i=imin,imax
          Bo(i,j)=g*( alpha(i,j)*(stflx(i,j,itemp)-srflx(i,j))
# ifdef SALINITY
     &                              -beta(i,j)*stflx(i,j,isalt)
# endif
     &                                                        )
          Bosol(i,j)=g*alpha(i,j)*srflx(i,j)
        enddo 

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

        do i=istr,iend

          k=kbl(i)
          z_bl=z_w(i,j,N)-hbl(i,j)

          if (swr_frac(i,j,k-1).gt. 0.) then
            Bfsfc=Bo(i,j) +Bosol(i,j)*( 1. -swr_frac(i,j,k-1)
     &              *swr_frac(i,j,k)*(z_w(i,j,k)-z_w(i,j,k-1))
     &               /( swr_frac(i,j,k  )*(z_w(i,j,k)   -z_bl)
     &                 +swr_frac(i,j,k-1)*(z_bl -z_w(i,j,k-1))
     &                                                      ))
          else
            Bfsfc=Bo(i,j)+Bosol(i,j)
          endif


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

!
! Compute vertical mixing coefficients
!
            sigma=(z_w(i,j,N)-z_w(i,j,k))/max(hbl(i,j),eps)


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
          enddo


        enddo

       enddo



      return
      end

