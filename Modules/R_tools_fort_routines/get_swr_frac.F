
 
      subroutine get_swr_frac (Lm,Mm,N, Hz, swr_frac )
!
! Compute fraction of solar shortwave flux penetrating to the
! specified depth due to exponential decay in Jerlov water type
! using Paulson and Simpson (1977) two-wavelength-band solar
! absorption model.
!
! Reference:  Paulson, C.A., and J.J. Simpson, 1977: Irradiance
! meassurements in the upper ocean, J. Phys. Oceanogr., 7, 952-956.
!
! This routine was adapted from Bill Large 1995 code.
!
! output: swr_frac (in "mixing.h")  shortwave radiation fraction
!
      implicit none

      integer Lm,Mm,N, i,j,k
      integer istr,iend,jstr,jend,  Jwt

      real*8  swdk1(0:Lm+1), swdk2(0:Lm+1)
      real*8 mu1(5),mu2(5), r1(5), attn1, attn2, xi1,xi2

      real*8 Hz(0:Lm+1,0:Mm+1,N)
     &      ,swr_frac(0:Lm+1,0:Mm+1,0:N)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


Cf2py intent(in) Lm,Mm,N, Hz
Cf2py intent(out) swr_frac

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      istr=0
      iend=Lm+1
      jstr=0
      jend=Mm+1


      mu1(1)=0.35    ! reciprocal of the absorption coefficient
      mu1(2)=0.6     ! for each of the two solar wavelength bands
      mu1(3)=1.0     ! as a function of Jerlov water type (Paulson
      mu1(4)=1.5     ! and Simpson, 1977) [dimensioned as length,
      mu1(5)=1.4     ! meters];
 
      mu2(1)=23.0
      mu2(2)=20.0
      mu2(3)=17.0
      mu2(4)=14.0
      mu2(5)=7.9
 
      r1(1)=0.58     ! fraction of the total radiance for
      r1(2)=0.62     ! wavelength band 1 as a function of Jerlov 
      r1(3)=0.67     ! water type (fraction for band 2 is always
      r1(4)=0.77     ! r2=1-r1);
      r1(5)=0.78
                     ! set Jerlov water type to assign everywhere
      Jwt=1          ! (an integer from 1 to 5).
 
      attn1=-1./mu1(Jwt)
      attn2=-1./mu2(Jwt) 


      do j=jstr,jend                     ! Algorithm: set fractions
        do i=istr,iend    
          swdk1(i)=r1(Jwt)               ! surface, then attenuate
          swdk2(i)=1.-swdk1(i)           ! them separately throughout
          swr_frac(i,j,N)=1.             ! the water column.
        enddo

        do k=N,1,-1
          do i=istr,iend
            xi1=attn1*Hz(i,j,k)
            if (xi1 .gt. -20.) then        ! this logic to avoid
              swdk1(i)=swdk1(i)*exp(xi1)   ! computing exponent for
            else                           ! a very large argument
              swdk1(i)=0.
            endif

            xi2=attn2*Hz(i,j,k)
            if (xi2 .gt. -20.) then
              swdk2(i)=swdk2(i)*exp(xi2)
            else
              swdk2(i)=0.
            endif
            swr_frac(i,j,k-1)=swdk1(i)+swdk2(i)
          enddo
        enddo
      enddo

      return
      end

