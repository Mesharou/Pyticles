!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute z_r and z_w for VERT_COORD_TYPE_KAU
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      subroutine zlevs_kau(Lm,Mm,N, h, hc, Cs_r, Cs_w,z_r,z_w)


      implicit none

      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k

      real*8 Cs_w(0:N), Cs_r(N), cff_w, cff_r, Cs_r_k, Cs_w_k,
     &       hc, ds,
     &       z_r(0:Lm+1,0:Mm+1,N),z_w(0:Lm+1,0:Mm+1,0:N),
     &       h(0:Lm+1,0:Mm+1)

      real*8 hinv_r(0:Lm+1), hinv_w(0:Lm+1)
      real Gcoord_w,Gcoord_r,s_r, s_w
      real, parameter :: eps=1.E-20


Cf2py intent(in) Lm,Mm,N, h,hc, Cs_w, Cs_r
Cf2py intent(out) z_r,z_w


      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1


      ds=1.D0/dble(N)

      do j=jmin,jmax

        Cs_r_k=Cs_r(N)
        s_r = -0.5D0 * ds
        Gcoord_r= 0.16+abs(s_r)**0.8*(1+s_r)**0.2 + (0.3 *  exp(-1 /
     &                max(abs(s_r), eps)) / max(abs(s_r), eps)**1.7)
!       Gcoord_r= 1 !! SM09 CASE

        cff_r=hc*s_r*Gcoord_r

        do i=imin,imax

          hinv_r(i)=h(i,j)/(h(i,j)+hc*Gcoord_r)
          z_w(i,j,N)=0.
          z_r(i,j,N)=hinv_r(i)*( cff_r + Cs_r_k*h(i,j) )

          z_w(i,j,0)=-h(i,j)

        enddo
        
        do k=N-1,1,-1

          Cs_w_k = Cs_w(k)
          Cs_r_k = Cs_r(k)

          s_w = ds* dble(k-N)
          s_r = ds* (dble(k-N)-0.5D0)


! first define the Gcoord parameter
          Gcoord_w= 0.16+abs(s_w)**0.8*(1+s_w)**0.2 + (0.3 *  exp(-1 /
     &                max(abs(s_w), eps)) / max(abs(s_w), eps)**1.7)
!
          Gcoord_r= 0.16+abs(s_r)**0.8*(1+s_r)**0.2 + (0.3 *  exp(-1 /
     &                max(abs(s_r), eps)) / max(abs(s_r), eps)**1.7)

!              Gcoord_w= 1   !!! SM09 CASE
!              Gcoord_r= 1   !!! SM09 CASE

          cff_w=hc*s_w*Gcoord_w
          cff_r=hc*s_r*Gcoord_r

          do i=imin,imax

            hinv_w(i)=1./(h(i,j)+hc*Gcoord_w)
            z_w(i,j,k)= h(i,j) * (cff_w+Cs_w_k*h(i,j))*hinv_w(i)

            hinv_r(i)=1./(h(i,j)+hc*Gcoord_r)
            z_r(i,j,k)= h(i,j)*(cff_r+Cs_r_k*h(i,j))*hinv_r(i)


          enddo
        enddo
      enddo
      end


