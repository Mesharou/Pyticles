
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!get_hbbls_from_AKt.F
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine get_hbbls_from_AKt(Lm,Mm,N1,N2,AKt,z_w,hbbls)


      implicit none
      integer Lm,Mm,N,N1,N2, imin,imax,jmin,jmax, i,j,k
      real*8 AKt(0:Lm+1,0:Mm+1,0:N1),  z_w(0:Lm+1,0:Mm+1,0:N2),
     &      hbbls(0:Lm+1,0:Mm+1)




Cf2py intent(in) Lm,Mm,N1,N2,AKt,z_w
Cf2py intent(out) hbbls




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1

      N = min(N1,N2)

      do j=jmin,jmax
        do i=imin,imax
            
            k=1

            do while((AKt(i,j,k)>1e-4).and.(k.lt.N))  
                k = k+1
            enddo

            hbbls(i,j) = z_w(i,j,max(1,k-1)) - z_w(i,j,0)


        enddo
      enddo






      return
      end







 
 
