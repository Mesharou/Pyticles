


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine get_intvortplanet (Lm,Mm,N,u,v, z_r,z_w,pm,pn,f
     &  ,depth1,depth2,vrtCor)    

      implicit none

      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k  


      real*8 u(1:Lm+1,0:Mm+1,N),  v(0:Lm+1,1:Mm+1,N),
     &      FlxU(1:Lm+1,0:Mm+1),  FlxV(0:Lm+1,1:Mm+1),
     &      z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N),
     &      zu(1:Lm+1,0:Mm+1,0:N),  zv(0:Lm+1,1:Mm+1,0:N),
     &      Hz(0:Lm+1,0:Mm+1,N), 
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      f(0:Lm+1,0:Mm+1),
     &      cff, cff1, cff2, depth1, depth2

      real*8  vrtCor(1:Lm+1,1:Mm+1)


# include "scalars.h"


Cf2py intent(in) Lm,Mm,N, u,v,z_r,z_w,pm,pn,f,depth1,depth2
Cf2py intent(out) vrtCor

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      depth1=abs(depth1)
      depth2=abs(depth2)

      if (depth2 .lt. depth1) then
            cff = depth1
            depth1 = depth2
            depth2 = cff
      end if

!       write(*,*) 'depths are', depth1, depth2

      do j=jmin,jmax
        do i=imin,imax
          do k=1,N,+1
           Hz(i,j,k) = z_w(i,j,k) - z_w(i,j,k-1)
          enddo

          do k=0,N,+1
           if (i.gt.0) zu(i,j,k) = -0.5*(z_w(i,j,k)+z_w(i-1,j,k))
           if (j.gt.0) zv(i,j,k) = -0.5*(z_w(i,j,k)+z_w(i,j-1,k))
          enddo
        enddo
      enddo


      do j=jmin,jmax
        do i=imin-1,imax

          FlxU(i,j) = 0
          do k=N,1,-1

            if (zu(i,j,k-1).le.depth2) then

                if ((zu(i,j,k-1).ge.depth1) .and.
     &                     (zu(i,j,k).le.depth1))  then
                    FlxU(i,j) =  (zu(i,j,k-1)-depth1) * u(i,j,k)

                else if (zu(i,j,k-1).ge.depth1)  then
                    FlxU(i,j) = FlxU(i,j) + 
     &             0.5*(Hz(i,j,k)+Hz(i-1,j,k)) * u(i,j,k)
                end if


            else if ((depth2-zu(i,j,k)).gt.0) then

                FlxU(i,j) = FlxU(i,j) + 
     &          (depth2-zu(i,j,k)) * u(i,j,k)

            end if

          enddo

        enddo !<- i
      enddo !<- j




      do j=jmin-1,jmax
        do i=imin,imax

          FlxV(i,j) = 0
          do k=N,1,-1

            if (zv(i,j,k-1).le.depth2) then

                if ((zv(i,j,k-1).ge.depth1) .and.
     &                      (zv(i,j,k).le.depth1))  then
                    FlxV(i,j) =  (zv(i,j,k-1)-depth1) * v(i,j,k)

                else if (zv(i,j,k-1).ge.depth1)  then
                    FlxV(i,j) = FlxV(i,j) + 
     &             0.5*(Hz(i,j,k)+Hz(i,j-1,k)) * v(i,j,k)
                end if


            else if ((depth2-zv(i,j,k)).gt.0) then

                FlxV(i,j) = FlxV(i,j) + 
     &          (depth2-zv(i,j,k)) * v(i,j,k)

            end if

          enddo
        enddo !<- i
      enddo !<- j




      do j=jmin+1,jmax
        do i=imin+1,imax

            cff1 = 0.5*FlxU(i,j) * 
     &          (f(i,j) - f(i-1,j)) * 0.5 * (pm(i,j)+ pm(i-1,j))
     &         + 0.5*FlxU(i,j-1) * 
     &          (f(i,j-1) - f(i-1,j-1)) * 0.5 * (pm(i,j-1)+ pm(i-1,j-1))


            cff2 = 0.5*FlxV(i,j) * 
     &          (f(i,j) - f(i,j-1)) * 0.5 * (pn(i,j)+ pn(i,j-1))
     &         + 0.5*FlxV(i-1,j) * 
     &          (f(i-1,j) - f(i-1,j-1)) * 0.5 * (pn(i-1,j)+ pn(i-1,j-1))


            vrtCor(i,j) = cff1 + cff2

        enddo !<- i
      enddo !<- j




















      return
      end











