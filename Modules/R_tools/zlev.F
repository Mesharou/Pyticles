!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute z_r and z_w for NEW_S_COORD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      subroutine zlev(Lm,Mm,N, h,zeta, hc, Cs_r, Cs_w,z_r)


      implicit none

      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k

      real*8 Cs_w(0:N), Cs_r(N),  cff_r,cff1_r,
     &       hc, ds,
     &       zeta(0:Lm+1,0:Mm+1),
     &       z_r(0:Lm+1,0:Mm+1,N),
     &       h(0:Lm+1,0:Mm+1),hinv(0:Lm+1,0:Mm+1)

     
Cf2py intent(in) Lm,Mm,N, h,zeta, hc, Cs_w, Cs_r
Cf2py intent(out) z_r


      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1


      ds=1.D0/dble(N)

      do j=jmin,jmax
        do i=imin,imax


          hinv(i,j)=1./(h(i,j)+hc)  

        enddo
        
        do k=1,N,+1  !--> irreversible because of recursion in Hz 

            
          cff_r=hc*ds*(dble(k-N)-0.5)

          cff1_r=Cs_r(k)


          do i=imin,imax


            z_r(i,j,k)=zeta(i,j) +(zeta(i,j)+h(i,j))
     &                        *(cff_r+cff1_r*h(i,j))*hinv(i,j)


          enddo
        enddo
      enddo
      end


