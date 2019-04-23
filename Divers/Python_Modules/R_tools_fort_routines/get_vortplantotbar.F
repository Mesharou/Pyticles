


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute Coriolis part of the barotropic vorticity balance equation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine get_vortplantotbar(Lm,Mm,u,v, Hz,pm,pn,f
     &  ,vrtCor)    

      implicit none

      integer Lm,Mm, imin,imax,jmin,jmax, i,j


      real*8 u(1:Lm+1,0:Mm+1),  v(0:Lm+1,1:Mm+1),
     &      FlxU(1:Lm+1,0:Mm+1),  FlxV(0:Lm+1,1:Mm+1),
     &      Hz(0:Lm+1,0:Mm+1), gamma,
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      dn_u(0:Lm+1,0:Mm+1), dm_v(0:Lm+1,0:Mm+1),
     &      f(0:Lm+1,0:Mm+1), fomn(0:Lm+1,0:Mm+1),
     &      dmde(0:Lm+1,0:Mm+1), dndx(0:Lm+1,0:Mm+1),
     &      var1, var2,var3, var4, cff, cff1, cff2


      real*8 wrkCor(0:Lm+1,0:Mm+1,2)

      real*8 wrk1(0:Lm+1,0:Mm+1), wrk2(0:Lm+1,0:Mm+1),
     &      UFx(0:Lm+1,0:Mm+1), UFe(0:Lm+1,0:Mm+1),
     &      VFx(0:Lm+1,0:Mm+1), VFe(0:Lm+1,0:Mm+1)

      real*8  vrtCor(1:Lm+1,1:Mm+1)


# include "scalars.h"

      parameter (gamma=0.25)  

Cf2py intent(in) Lm,Mm,N, u,v,z_r,z_w,pm,pn,f
Cf2py intent(out) vrtCor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      do j=jmin,jmax
        do i=imin,imax
           fomn(i,j)=f(i,j)/(pm(i,j)*pn(i,j))
        enddo
      enddo


!       do j=jmin+1,jmax-1
!         do i=imin+1,imax-1
!             dmde(i,j) = 0.5/pm(i,j+1)-0.5/pm(i,j-1)
!             dndx(i,j) = 0.5/pn(i+1,j)-0.5/pn(i-1,j)
!          enddo
!       enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





        do j=jmin+1,jmax-1
          do i=imin+1,imax-1
            cff=0.5*Hz(i,j)*(
     &              fomn(i,j)
!      &             +0.5*( (v(i,j)+v(i,j+1))*dndx(i,j)
!      &                   -(u(i,j)+u(i+1,j))*dmde(i,j))
     &                                                             )
            UFx(i,j)=cff*(v(i,j)+v(i,j+1))
            VFe(i,j)=cff*(u(i,j)+u(i+1,j))
          enddo
        enddo

        do j=jmin+1,jmax-1
          do i=imin+1,imax-1


                wrkCor(i,j,1) = 0.5*(UFx(i,j)+UFx(i-1,j))
 

          enddo
        enddo

        do j=jmin+1,jmax-1
          do i=imin+1,imax-1


                wrkCor(i,j,2) = -0.5*(VFe(i,j)+VFe(i,j-1))


          enddo
        enddo







! Divide all diagnostic terms by (pm*pn).
! There after the unit of these terms are :
!  s-2


      do j=jmin+1,jmax-1
        do i=imin+1,imax-1

            cff=0.25*(pm(i,j)+pm(i-1,j))
     &             *(pn(i,j)+pn(i-1,j))

            
            wrkCor(i,j,1)=wrkCor(i,j,1)*cff


            cff=0.25*(pm(i,j)+pm(i,j-1))
     &             *(pn(i,j)+pn(i,j-1))


            wrkCor(i,j,2)=wrkCor(i,j,2)*cff


          enddo
        enddo



      do j=jmin+1,jmax-1
        do i=imin+1,imax-1

           cff1 = 0.25*(pm(i,j) + pm(i-1,j) + pm(i,j-1) + pm(i-1,j-1))
           cff2 = 0.25*(pn(i,j) + pn(i-1,j) + pn(i,j-1) + pn(i-1,j-1))

            vrtCor(i,j) = (wrkCor(i,j,2) - wrkCor(i-1,j,2)) * cff1
     &                     - (wrkCor(i,j,1) - wrkCor(i,j-1,1)) * cff2



          enddo
         enddo





      return
      end











