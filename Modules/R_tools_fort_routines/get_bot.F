
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute PV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_bot(Lm,Mm, u,v, Hz
     &,rdrg,ubot,vbot)


      implicit none
      integer Lm,Mm, imin,imax,jmin,jmax, i,j,k
      real*8 ubot(1:Lm+1,0:Mm+1), vbot(0:Lm+1,1:Mm+1),
     &      u(1:Lm+1,0:Mm+1), v(0:Lm+1,1:Mm+1),
     &      Hz(0:Lm+1,0:Mm+1), rd(0:Lm+1,0:Mm+1),
     &      cff

      real*8 rdrg, Zob

# include "scalars.h"


Cf2py intent(in) Lm,Mm,u,v, Hz,rdrg
Cf2py intent(out) ubot,vbot

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bottom Drag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       Zob=0.01


       do j=jmin+1,jmax-1
         do i=imin+1,imax-1


# ifdef LINEAR_DRAG_ONLY
            rd(i,j)=rdrg
# else
            cff=sqrt( 0.333333333333*(
     &              u(i,j)**2 +u(i+1,j)**2
     &                    +u(i,j)*u(i+1,j)
     &              +v(i,j)**2+v(i,j+1)**2
     &                    +v(i,j)*v(i,j+1)
     &                                               ))
            rd(i,j)=rdrg + cff*(vonKar/log(Hz(i,j)/Zob))**2
# endif



          enddo
        enddo

       do j=jmin+1,jmax
            rd(imax ,j)=rd(imax-1 ,j)
            rd(imin ,j)=rd(imin+1 ,j)
        enddo

       do i=imin+1,imax
            rd(i ,jmax)=rd(i ,jmax-1)
            rd(i ,jmin)=rd(i ,jmin+1)
        enddo




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




        do j=jmin,jmax
          do i=imin+1,imax

            ubot(i,j) = 0.5 * (rd(i,j)+rd(i-1,j))*u(i,j)

        enddo !<- i
      enddo !<- j
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




        do j=jmin+1,jmax
          do i=imin,imax


            vbot(i,j) = 0.5 * (rd(i,j)+rd(i,j-1))*v(i,j)

        enddo !<- i
      enddo !<- j





      return
      end







 
 
