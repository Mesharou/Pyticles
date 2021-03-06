
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute bottom drag [for croco version]
!!
!! - modified 18/01/03
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# define LIMIT_BSTRESS

      subroutine get_bot_croco(Lm,Mm, u,v, Hz, Hz2
     &,rdrg,rdrg2,Zob,Cdb_min,Cdb_max,dt,bustr,bvstr)

      implicit none
      integer Lm,Mm, imin,imax,jmin,jmax, i,j,k
      integer Istr,Iend,Jstr,Jend,IstrU,JstrV

      real*8 bustr(1:Lm+1,0:Mm+1), bvstr(0:Lm+1,1:Mm+1),
     &      u(1:Lm+1,0:Mm+1), v(0:Lm+1,1:Mm+1),
     &      Hz(0:Lm+1,0:Mm+1), wrk(0:Lm+1,0:Mm+1),
     &      Hz2(0:Lm+1,0:Mm+1),
     &      cff, cff1

      real*8 rdrg, rdrg2, Zob, Cdb_min, Cdb_max,dt

      !parameter (rdrg2 = 0.0, Cdb_max = 0.1, Cdb_min = 0.0001)


# include "scalars.h"


Cf2py intent(in) Lm,Mm,u,v, Hz,Hz2,rdrg,rdrg2,Zob,Cdb_min,Cdb_max,dt
Cf2py intent(out) bustr,bvstr

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1

      Istr = imin
      Iend = imax
      Jstr = jmin
      Jend = jmax

      IstrU = imin+1
      JstrV = jmin+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Bottom Drag
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!
!  Set bottom stress using logarithmic or linear 
!  and/or quadratic formulation. 
!
      if (Zob.ne.0.) then
        do j=JstrV-1,Jend
          do i=IstrU-1,Iend
            cff=vonKar/LOG((Hz2(i,j))/Zob)
            wrk(i,j)=MIN(Cdb_max,MAX(Cdb_min,cff*cff))
          enddo
        enddo
        do j=Jstr,Jend
          do i=IstrU,Iend
            cff=0.25*(v(i  ,j)+v(i  ,j+1)+
     &                v(i-1,j)+v(i-1,j+1))
            bustr(i,j)=0.5*(wrk(i-1,j)+wrk(i,j))*u(i,j)*
     &                 SQRT(u(i,j)*u(i,j)+cff*cff)
          enddo
        enddo
        do j=JstrV,Jend
          do i=Istr,Iend
            cff=0.25*(u(i,j  )+u(i+1,j)+
     &                u(i,j-1)+u(i+1,j-1))
            bvstr(i,j)=0.5*(wrk(i,j-1)+wrk(i,j))*v(i,j)*
     &                 SQRT(cff*cff+v(i,j)*v(i,j))
          enddo
        enddo
      elseif (rdrg2.gt.0.) then
        do j=JstrV,Jend
          do i=Istr,Iend

            cff=0.25*(v(i,j)+v(i,j+1)+v(i-1,j)+
     &                                            v(i-1,j+1))
            bustr(i,j)=u(i,j)*(rdrg2*sqrt(
     &                         u(i,j)*u(i,j)+cff*cff
     &                                                            ))
          enddo
        enddo
        do j=Jstr,Jend
          do i=IstrU,Iend
            cff=0.25*(u(i,j)+u(i+1,j)+u(i,j-1)+
     &                                            u(i+1,j-1))
            bvstr(i,j)=v(i,j)*(rdrg2*sqrt(
     &                         cff*cff+v(i,j)*v(i,j)
     &                                                            ))
          enddo
        enddo
      else
        do j=Jstr,Jend
          do i=Istr,Iend
            bustr(i,j)=rdrg*u(i,j)
          enddo
        enddo
        do j=Jstr,Jend
          do i=Istr,Iend
            bvstr(i,j)=rdrg*v(i,j)
          enddo
        enddo
      endif

#  ifdef LIMIT_BSTRESS
!
!  From J. Warner code:
!  Set limiting factor for bottom stress. The bottom stress is adjusted
!  to not change the direction of momentum.  It only should slow down
!  to zero.  The value of 0.75 is arbitrary limitation assigment.
!
      cff=0.75/dt
      do j=Jstr,Jend
        do i=IstrU,Iend
          cff1=cff*0.5*(Hz(i-1,j)+Hz(i,j))
          bustr(i,j)=SIGN(1.D0, bustr(i,j))* 
     &               MIN(ABS(bustr(i,j)),
     &                   ABS(u(i,j))*cff1)
        enddo
      enddo
      do j=JstrV,Jend
        do i=Istr,Iend
          cff1=cff*0.5*(Hz(i,j-1)+Hz(i,j))
          bvstr(i,j)=SIGN(1.D0, bvstr(i,j))*
     &               MIN(ABS(bvstr(i,j)),
     &                   ABS(v(i,j))*cff1)
        enddo
      enddo
#  endif


      return
      end







 
 
