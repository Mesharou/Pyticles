!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Z interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sigma_to_z_intr_bot_2d (Lm,N, nz, z_r, z_w, rmask, var,
     &                    z_lev, var_zlv, below,imin,kmin, FillValue)
!
! Interpolate field "var" defined in sigma-space to 3-D z_lev. 
!


      implicit none

      integer Lm,Mm,N, nz,  imin,imax, kmin, i,k,m

      integer km(0:Lm+1)

      real*8 var(imin:Lm+1,kmin:N), 
     &       z_r(0:Lm+1,N),   rmask(0:Lm+1),
     &       z_w(0:Lm+1,0:N), z_lev(0:Lm+1,nz),
     &       FillValue, var_zlv(imin:Lm+1,nz),
     &       zz(0:Lm+1,0:N+1), dpth, below
#ifdef CUBIC
     &     , dz(0:Lm+1,kmin-1:N), FC(0:Lm+1,kmin-1:N), p,q,cff
#endif
      integer numthreads, trd, chunk_size, margin, jstr,jend
C$    integer omp_get_num_threads, omp_get_thread_num


      imax=Lm+1


Cf2py intent(in) Lm,Mm,N, nz, z_r, z_w, rmask, var, z_lev, below, imin,kmin, FillValue
Cf2py intent(out) var_zlv



        if (kmin.eq.1) then
          if (imin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=z_r(i,k)
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=z_w(i,0)
              zz(i,N+1)=z_w(i,N)
            enddo
          elseif (imin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.5D0*(z_r(i,k)+z_r(i-1,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.5D0*(z_w(i-1,0)+z_w(i,0))
              zz(i,N+1)=0.5D0*(z_w(i-1,N)+z_w(i,N))
            enddo
          endif
        else
          if (imin.eq.0 ) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=z_w(i,k)
              enddo
            enddo
          elseif (imin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.5D0*(z_w(i,k)+z_w(i-1,k))
              enddo
            enddo
          endif
        endif
#ifdef CUBIC
        do k=kmin,N-1
          do i=imin,imax
            dz(i,k)=zz(i,k+1)-zz(i,k)
            FC(i,k)=var(i,k+1)-var(i,k)
          enddo
        enddo
        do i=imin,imax
          dz(i,kmin-1)=dz(i,kmin)
          FC(i,kmin-1)=FC(i,kmin)

          dz(i,N)=dz(i,N-1)
          FC(i,N)=FC(i,N-1)
        enddo
        do k=N,kmin,-1   !--> irreversible
          do i=imin,imax
            cff=FC(i,k)*FC(i,k-1)
            if (cff.gt.0.D0) then
              FC(i,k)=cff*(dz(i,k)+dz(i,k-1))/( (FC(i,k)+FC(i,k-1))
     &                                          *dz(i,k)*dz(i,k-1) )
            else
              FC(i,k)=0.D0
            endif
          enddo
        enddo
#endif
        do m=1,nz
       

          if (kmin.eq.0) then                      !
            do i=imin,imax                         !
              dpth=zz(i,N)-zz(i,0) 
              if (rmask(i).lt.0.5) then
                km(i)=-3     !--> masked out
              elseif (dpth*(z_lev(i,m)-zz(i,N)).gt.0.) then
                km(i)=N+2    !<-- above surface
              elseif (dpth*(zz(i,0)-z_lev(i,m)).gt.0.) then
                km(i)=-2     !<-- below bottom
              else
                km(i)=-1     !--> to search
              endif
            enddo
          else
            do i=imin,imax
              dpth=zz(i,N+1)-zz(i,0)     
              if (rmask(i).lt.0.5) then
                km(i)=-3     !--> masked out
              elseif (dpth*(z_lev(i,m)-zz(i,N+1)).gt.0.) then
                km(i)=N+2    !<-- above surface

              elseif (dpth*(z_lev(i,m)-zz(i,N)).gt.0.) then
                km(i)=N      !<-- below surface, but above z_r(N)
              elseif (dpth*(zz(i,0)-below-z_lev(i,m)).gt.0.) then
                km(i)=-3     !<-- below bottom 
              elseif (dpth*(zz(i,0)-z_lev(i,m)).gt.0.) then
                km(i)=-2     !<-- below bottom but close 
              elseif (dpth*(zz(i,1)-z_lev(i,m)).gt.0.) then
                km(i)=0      !<-- above bottom, but below z_r(1)
              else
                km(i)=-1     !--> to search
              endif
            enddo
          endif
          do k=N-1,kmin,-1
            do i=imin,imax
              if (km(i).eq.-1) then
                if((zz(i,k+1)-z_lev(i,m))*(z_lev(i,m)-zz(i,k))
     & .ge. 0.) km(i)=k
              endif
            enddo
          enddo

          do i=imin,imax
            if (km(i).eq.-3) then
              var_zlv(i,m)=FillValue             !<-- masked out
            elseif (km(i).eq.-2) then
#ifdef INTERP_BELOW
              var_zlv(i,m)=var(i,kmin)  !    
#ifdef CUBIC
     &           -FC(i,kmin)*(zz(i,kmin)-z_lev(i,m))
#else
     &          -(zz(i,kmin)-z_lev(i,m))*(var(i,kmin+1)-var(i,kmin))
     &                                  /(zz(i,kmin+1)-zz(i,kmin))
#endif
#else
              var_zlv(i,m)=FillValue      !<-- below bottom
#endif
            elseif (km(i).eq.N+2) then
#ifdef INTERP_ABOVE
              var_zlv(i,m)=var(i,N)     !-> R-point, above z_r(N)
#ifdef CUBIC
     &              +FC(i,N)*(z_lev(i,m)-zz(i,N))
#else
     &              +(z_lev(i,m)-zz(i,N))*(var(i,N)-var(i,N-1))
     &                                   /(zz(i,N)-zz(i,N-1))
#endif
#else
              var_zlv(i,m)=-FillValue     !<-- above surface
#endif
            elseif (km(i).eq.N) then
              var_zlv(i,m)=var(i,N)     !-> R-point, above z_r(N)
#ifdef CUBIC
     &              +FC(i,N)*(z_lev(i,m)-zz(i,N))
#else
     &              +(z_lev(i,m)-zz(i,N))*(var(i,N)-var(i,N-1))
     &                                   /(zz(i,N)-zz(i,N-1))
#endif
            elseif (km(i).eq.kmin-1) then   !-> R-point below z_r(1),
              var_zlv(i,m)=var(i,kmin)  !     but above bottom  
#ifdef CUBIC
     &           -FC(i,kmin)*(zz(i,kmin)-z_lev(i,m))
#else
     &          -(zz(i,kmin)-z_lev(i,m))*(var(i,kmin+1)-var(i,kmin))
     &                                  /(zz(i,kmin+1)-zz(i,kmin))
#endif
            else
              k=km(i)
              !write(*,*) k,km
#ifdef CUBIC
              cff=1.D0/(zz(i,k+1)-zz(i,k))
              p=z_lev(i,m)-zz(i,k)
              q=zz(i,k+1)-z_lev(i,m)

              var_zlv(i,m)=cff*( q*var(i,k) + p*var(i,k+1)
     &             -cff*p*q*( cff*(q-p)*(var(i,k+1)-var(i,k))
     &                                  +p*FC(i,k+1) -q*FC(i,k) )
     &                                                          )
#else
              var_zlv(i,m)=( var(i,k)*(zz(i,k+1)-z_lev(i,m))
     &                        +var(i,k+1)*(z_lev(i,m)-zz(i,k))
     &                               )/(zz(i,k+1)-zz(i,k))
#endif


            !write(*,*) 'bof',i,k,zz(i,k), zz(i,k+1), z_lev(i,m), m



            endif
          enddo
        enddo ! <-- m


      return
      end
 
