
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute PV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_windbot(Lm,Mm,N, ubot,vbot,uwind,vwind, Hz,
     &,rho0,pm,pn,rdrg,wind,bot)


      implicit none
      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k
      real*4 ubot(1:Lm+1,0:Mm+1), vbot(0:Lm+1,1:Mm+1),
     &      uwind(1:Lm+1,0:Mm+1), vwind(0:Lm+1,1:Mm+1),
     &      Hz(0:Lm+1,0:Mm+1), 
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1)


      real*4 rho0

Cf2py intent(in) Lm,Mm,N, ubot,vbot,uwind,vwind, Hz,rho0,pm,pn
Cf2py intent(out) wind,bot




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
! COMPUTE VORTICITY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        do j=jmin+1,jmax
          do i=imin+1,imax
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            ubot(i,j) = 0.5 * (rd(i,j)+rd(i-1,j))*u(i,j))
            vbot(i,j) = 0.5 * (rd(i,j)+rd(i,j-1))*v(i,j))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            dvbotdx = (vbot(i,j) - vbot(i-1,j)) * 0.25 * 
     &           (pm(i,j)+pm(i-1,j)+pm(i,j-1)+pm(i-1,j-1))



            dubotdy = (ubot(i,j) - ubot(i,j-1) )* 0.25 * 
     &           (pn(i,j)+pn(i-1,j)+pn(i,j-1)+pn(i-1,j-1))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            dvwinddx = (vwind(i,j) - vwind(i-1,j)) * 0.25 * 
     &           (pm(i,j)+pm(i-1,j)+pm(i,j-1)+pm(i-1,j-1))



            duwinddy = (uwind(i,j) - uwind(i,j-1) )* 0.25 * 
     &           (pn(i,j)+pn(i-1,j)+pn(i,j-1)+pn(i-1,j-1))


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





        enddo !<- i
      enddo !<- j



      return
      end







 
 
