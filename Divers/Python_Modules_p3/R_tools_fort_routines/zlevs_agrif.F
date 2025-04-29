!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute z_r and z_w for NEW_S_COORD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      subroutine zlevs_agrif(Lm,Mm,N, h,zeta, hc, Cs_r, Cs_w,
     &       sc_r,sc_w,z_r,z_w)


      implicit none

      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k

      real*8 Cs_w(0:N), Cs_r(N), cff_w, cff_r, cff1_w, cff1_r,
     &       sc_w(0:N), sc_r(N),
     &       hc, ds, cff0_w, cff0_r,
     &       zeta(0:Lm+1,0:Mm+1),
     &       z_r(0:Lm+1,0:Mm+1,N),z_w(0:Lm+1,0:Mm+1,0:N),
     &       h(0:Lm+1,0:Mm+1),hinv(0:Lm+1,0:Mm+1)

     
Cf2py intent(in) Lm,Mm,N, h,zeta, hc, Cs_w, Cs_r, sc_r,sc_w
Cf2py intent(out) z_r,z_w


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


          cff_w =hc*sc_w(k)
          cff_r =hc*sc_r(k)

          cff1_w=Cs_w(k)
          cff1_r=Cs_r(k)

          do i=imin,imax

            cff0_w=cff_w+cff1_w*h(i,j)
            cff0_r=cff_r+cff1_r*h(i,j)


            z_w(i,j,k)=cff0_w*h(i,j)*hinv(i,j)+zeta(i,j)
     &                     *(1.+cff0_w*hinv(i,j))
            z_r(i,j,k)=cff0_r*h(i,j)*hinv(i,j)+zeta(i,j)
     &                     *(1.+cff0_r*hinv(i,j))



          enddo
        enddo
      enddo
      end


