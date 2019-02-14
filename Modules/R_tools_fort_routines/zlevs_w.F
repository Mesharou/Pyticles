!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute z_w for NEW_S_COORD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      subroutine zlevs_w(Lm,Mm,N, h,zeta, hc, Cs_w,z_w)


      implicit none

      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k

      real*8 Cs_w(0:N), cff_w, cff1_w,
!      &       Cs_r(N), cff_r, cff1_r,
     &       hc, ds,
     &       zeta(0:Lm+1,0:Mm+1),
     &       z_w(0:Lm+1,0:Mm+1,0:N),
     &       h(0:Lm+1,0:Mm+1),hinv(0:Lm+1,0:Mm+1)

     
Cf2py intent(in) Lm,Mm,N, h,zeta, hc, Cs_w
Cf2py intent(out) z_w


      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1


      ds=1.D0/dble(N)

      do j=jmin,jmax
        do i=imin,imax

          hinv(i,j)=1./(h(i,j)+hc)  
          z_w(i,j,0)=-h(i,j)

        enddo
        
        do k=1,N,+1  !--> irreversible because of recursion in Hz 

            
          cff_w=hc*ds* dble(k-N)
!           cff_r=hc*ds*(dble(k-N)-0.5)

          cff1_w=Cs_w(k)
!           cff1_r=Cs_r(k)


          do i=imin,imax

            z_w(i,j,k)=zeta(i,j) +(zeta(i,j)+h(i,j))
     &                        *(cff_w+cff1_w*h(i,j))*hinv(i,j) 

!             z_r(i,j,k)=zeta(i,j) +(zeta(i,j)+h(i,j))
!      &                        *(cff_r+cff1_r*h(i,j))*hinv(i,j)


          enddo
        enddo
      enddo
      end


