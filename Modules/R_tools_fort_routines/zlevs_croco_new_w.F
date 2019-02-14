!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute z_r and z_w for new or old S Coord
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# define NEW_S_COORD 


      subroutine zlevs_croco_new_w(Lm,Mm,N, h,zeta, hc, Cs_w,
     &                                      sc_w ,z_w)

      implicit none

      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k

      real*8 Cs_w(0:N), cff_w, cff1_w,
     &       sc_w(0:N),
     &       hc, ds, z_w0,
     &       zeta(0:Lm+1,0:Mm+1),
     &       z_w(0:Lm+1,0:Mm+1,0:N),
     &       h(0:Lm+1,0:Mm+1),hinv(0:Lm+1,0:Mm+1)

     
Cf2py intent(in) Lm,Mm,N, h,zeta, hc, Cs_w, sc_w
Cf2py intent(out) z_w


      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1


      ds=1.D0/dble(N)

      do j=jmin,jmax
        do i=imin,imax

#  ifdef NEW_S_COORD 
            hinv(i,j)=1./(h(i,j)+hc)      
#  else
            hinv(i,j)=1./h(i,j) 
#  endif

          z_w(i,j,0)=-h(i,j)

        enddo
        
        do k=1,N,+1  !--> irreversible because of recursion in Hz 

#  ifdef NEW_S_COORD
          cff_w =hc*sc_w(k)
          cff1_w=Cs_w(k)
#  else
          cff_w =hc*(sc_w(k)-Cs_w(k))
          cff1_w=Cs_w(k)

#  endif           


          do i=imin,imax

            z_w0=cff_w+cff1_w*h(i,j) 

#  ifdef NEW_S_COORD
            z_w(i,j,k)=z_w0*h(i,j)*hinv(i,j)+zeta(i,j)
     &                                       *(1.+z_w0*hinv(i,j))
#  else
            z_w(i,j,k)=z_w0+zeta(i,j)*(1.+z_w0*hinv(i,j))      
#  endif

          enddo
        enddo
      enddo
      end


