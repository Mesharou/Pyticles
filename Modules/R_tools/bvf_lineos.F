!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Compute density anomaly (adapted from rho_eos.F)
!! without SALINITY and with linear EOS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine bvf_lineos(Lm,Mm,N, T,S,Tcoef, Scoef, z_r,R0,rho0,bvf)

      implicit none
      integer Lm,Mm,N, i,j,k,
     &      istrR,iendR,jstrR,jendR

      real*8 T(0:Lm+1,0:Mm+1,N),  S(0:Lm+1,0:Mm+1,N),
     &      z_r(0:Lm+1,0:Mm+1,N), 
     &      rho(0:Lm+1,0:Mm+1,N), bvf(0:Lm+1,0:Mm+1,0:N),
     &      cff, rho0, R0,Tcoef, Scoef

      real*8, parameter :: qp2=0.0000172, g=9.81

      integer numthreads, trd, chunk_size, margin, jstr,jend


Cf2py intent(in) Lm,Mm,N,T,S,Tcoef, Scoef, z_r,R0,rho0
Cf2py intent(out) bvf


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! R0         Coefficients for linear Equation of State (EOS)
! T0,Tcoef  
! S0,Scoef          rho = R0 - Tcoef*(T-T0) + Scoef*(S-S0)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      istrR=0
      iendR=Lm+1
      jstrR=0
      jendR=Mm+1

      do j=jstrR,jendR
        do k=1,N
          do i=istrR,iendR
            rho(i,j,k)=R0 -Tcoef*T(i,j,k)
#   ifdef SALINITY
     &                    +Scoef*S(i,j,k)
#   endif

          enddo
        enddo


        cff=g/rho0
        do k=1,N-1
          do i=istrR,iendR
            bvf(i,j,k)=cff*(rho(i,j,k)-rho(i,j,k+1))
     &                     /(z_r(i,j,k+1)-z_r(i,j,k))

          enddo
        enddo
        do i=istrR,iendR
          bvf(i,j,N)=bvf(i,j,N-1)
          bvf(i,j,0)=bvf(i,j,  1)
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      enddo    ! <-- j
      return
      end


