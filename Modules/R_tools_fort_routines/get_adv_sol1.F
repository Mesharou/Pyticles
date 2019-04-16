!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute Advective part of the barotropic vorticity balance equation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






      subroutine get_adv_sol1 (Lm,Mm,N,u,v, w,z_r,z_w,pm,pn
     &  ,adv)    



      implicit none

      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k,
     &      istr,iend,jstr,jend,istrU,jstrV



      real*8 u(1:Lm+1,0:Mm+1,N),  v(0:Lm+1,1:Mm+1,N),
     &      Hz(0:Lm+1,0:Mm+1,N), 
     &      z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N),
     &      Hu,  Hv,
     &      pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     &      var1, var2,var3, var4

      real*8 adv(1:Lm+1,1:Mm+1)

      real*8 advux(1:Lm+1,0:Mm+1),advuy(1:Lm+1,0:Mm+1),
     &      advuz(1:Lm+1,0:Mm+1), 
     &      advvx(0:Lm+1,1:Mm+1), advvy(0:Lm+1,1:Mm+1),
     &      advvz(0:Lm+1,1:Mm+1)


      real*8 w(0:Lm+1,0:Mm+1,0:N)


# include "scalars.h"



Cf2py intent(in) Lm,Mm,N, u,v,w,z_r,z_w,pm,pn
Cf2py intent(out) adv

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! COMPUTE VORTICITY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!---------------------------------------------------------------------------------------


      do k=1,N
        do i=imin,imax
         do j=jmin,jmax

         Hz(i,j,k) = z_w(i,j,k) - z_w(i,j,k-1)

         enddo
        enddo
       enddo




        do i=imin+1,imax
         do j=jmin,jmax
         advux(i,j) = 0.
         advuy(i,j) = 0.
         advuz(i,j) = 0.
        enddo
       enddo

        do i=imin,imax
         do j=jmin+1,jmax
         advvx(i,j) = 0.
         advvy(i,j) = 0.
         advvz(i,j) = 0.
        enddo
       enddo



!---------------------------------------------------------------------------------------

      do k=1,N
       do i=imin+1,imax
        do j=jmin+1,jmax


         Hu = 0.5*(Hz(i,j,k) + Hz(i-1,j,k))

      if (i.eq.1) then
        var1 = u(i+1,j,k)-u(i,j,k)
      elseif (i.eq.imax) then
        var1 = u(i,j,k)-u(i-1,j,k)
      else
        var1 = 0.5* (u(i+1,j,k)-u(i-1,j,k))
      endif

         advux(i,j) = advux(i,j) + Hu *  u(i,j,k) * 
     &       var1 *
     &       0.5 *  (pm(i,j) + pm(i-1,j))


      if (j.eq.1) then
        var1 =  0.5 * ( v(i,j,k) + v(i-1,j,k))
      else
        var1 =  0.25 * ( v(i,j,k) + v(i,j-1,k)
     &   +  v(i-1,j,k) + v(i-1,j-1,k) )
      endif


      if (j.eq.jmax) then
         advuy(i,j) = advuy(i,j) + Hu *  
     &    var1*
     &    ( u(i,j,k)-u(i,j-1,k)
     &       )
     &    * 0.5 * (pn(i,j) + pn(i-1,j))
      else
         advuy(i,j) = advuy(i,j) + Hu *  
     &    var1 *
     &    0.5* ( u(i,j,k)-u(i,j-1,k)
     &       + u(i,j+1,k)-u(i,j,k) )
     &    * 0.5 * (pn(i,j) + pn(i-1,j))
      endif

      if (k.eq.1) then
         advuz(i,j) = advuz(i,j) +
     &    0.25 * ( w(i,j,k) + w(i-1,j,k)
     &    + w(i,j,k-1) + w(i-1,j,k-1))*
     &    (u(i,j,k+1) -  u(i,j,k)) 
      else
        if (k.eq.N) then
         advuz(i,j) = advuz(i,j) +
     &    0.25 * ( w(i,j,k) + w(i-1,j,k)
     &    + w(i,j,k-1) + w(i-1,j,k-1))*
     &    (u(i,j,k) -  u(i,j,k-1)) 
        else
         advuz(i,j) = advuz(i,j) +
     &    0.25 * ( w(i,j,k) + w(i-1,j,k)
     &    + w(i,j,k-1) + w(i-1,j,k-1))*
     &    0.5*(u(i,j,k+1) -  u(i,j,k-1)) 
        endif
      endif


         enddo
       enddo









       do i=imin+1,imax
        do j=jmin+1,jmax

         Hv = 0.5*(Hz(i,j,k) + Hz(i,j-1,k))





      if (i.eq.1) then
        var1 =  0.5 * ( u(i,j,k) + u(i,j-1,k))
      else
        var1 =  0.25 * ( u(i,j,k) + u(i,j-1,k)
     &   +  u(i-1,j,k) + u(i-1,j-1,k))
      endif





         if (i.eq.imax) then
         advvx(i,j) = advvx(i,j) + Hv *  
     &    var1*
     &    ( v(i,j,k)-v(i-1,j,k)
     &       )
     &    *0.5*(pm(i,j) + pm(i,j-1))
         else
         advvx(i,j) = advvx(i,j) + Hv *  
     &    var1*
     &    0.5* ( v(i,j,k)-v(i-1,j,k)
     &       + v(i+1,j,k)-v(i,j,k) )
     &    *0.5*(pm(i,j) + pm(i,j-1))
         endif



      if (j.eq.1) then
        var1 = v(i,j+1,k)-v(i,j,k)
      elseif (i.eq.imax) then
        var1 = v(i,j,k)-v(i,j-1,k)
      else
        var1 = 0.5* ( v(i,j+1,k)-v(i,j-1,k))
      endif


         advvy(i,j) = advvy(i,j) + Hv *  v(i,j,k) * 
     &       var1
     &   *0.5*(pn(i,j) + pn(i,j-1))


      if (k.eq.1) then
         advvz(i,j) = advvz(i,j) +
     &    0.25 * ( w(i,j,k) + w(i,j-1,k)
     &    + w(i,j,k-1) + w(i-1,j,k-1))*
     &    (v(i,j,k+1) -  v(i,j,k))
      else
        if (k.eq.N) then
         advvz(i,j) = advvz(i,j) +
     &    0.25 * ( w(i,j,k) + w(i,j-1,k)
     &    + w(i,j,k-1) + w(i-1,j,k-1))*
     &    (v(i,j,k) -  v(i,j,k-1))
        else
         advvz(i,j) = advvz(i,j) +
     &    0.25 * ( w(i,j,k) + w(i,j-1,k)
     &    + w(i,j,k-1) + w(i-1,j,k-1))*
     &    0.5*(v(i,j,k+1) -  v(i,j,k-1))
        endif
      endif

         enddo
       enddo
      enddo







        adv = 0.


        do j=jmin+2,jmax
          do i=imin+2,imax

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
          var1 = ( (advux(i,j) + advuy(i,j) + advuz(i,j)) -
     &        ( advux(i,j-1) + advuy(i,j-1) + advuz(i,j-1)))*
     &        0.25 * (pn(i,j)+pn(i-1,j)+pn(i,j-1)+pn(i-1,j-1))



          var2 = ((advvx(i,j) + advvy(i,j) + advvz(i,j)) -
     &        ( advvx(i-1,j) + advvy(i-1,j) + advvz(i-1,j) ))*
     &        0.25 * (pm(i,j)+pm(i-1,j)+pm(i,j-1)+pm(i-1,j-1))



         adv(i,j) = var2 - var1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


        enddo !<- i
       enddo !<- j



      return
      end

 
