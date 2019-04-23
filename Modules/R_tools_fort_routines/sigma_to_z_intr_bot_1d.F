!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Z interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sigma_to_z_intr_bot_1d (N, nz, z_r, z_w, var,
     &                    z_lev,below, kmin, FillValue,var_zlv)
!
! Interpolate field "var" defined in sigma-space to z_lev. 
!


      implicit none

      integer N, nz,  kmin,k,m

      integer km

      real*8 var(kmin:N), 
     &       z_r(N),  
     &       z_w(0:N), z_lev(nz),
     &       FillValue, var_zlv(nz),
     &       zz(0:N+1), dpth, below
#ifdef CUBIC
     &     , dz(kmin-1:N), FC(kmin-1:N), p,q,cff
#endif


Cf2py intent(in) N, nz, z_r, z_w, var, z_lev, below,kmin, FillValue
Cf2py intent(out) var_zlv






        if (kmin.eq.1) then
            do k=1,N
                zz(k)=z_r(k)
            enddo
              zz(0)=z_w(0)
              zz(N+1)=z_w(N)
        else
            do k=0,N
                zz(k)=z_w(k)
            enddo
        endif
#ifdef CUBIC
        do k=kmin,N-1
            dz(k)=zz(k+1)-zz(k)
            FC(k)=var(k+1)-var(k)
        enddo
          dz(kmin-1)=dz(kmin)
          FC(kmin-1)=FC(kmin)
          dz(N)=dz(N-1)
          FC(N)=FC(N-1)

        do k=N,kmin,-1   !--> irreversible
            cff=FC(k)*FC(k-1)
            if (cff.gt.0.D0) then
              FC(k)=cff*(dz(k)+dz(k-1))/( (FC(k)+FC(k-1))
     &                                          *dz(k)*dz(k-1) )
            else
              FC(k)=0.D0
            endif
        enddo
#endif
        do m=1,nz
          if (kmin.eq.0) then                      !
              dpth=zz(N)-zz(0) 
              if (dpth*(z_lev(m)-zz(N)).gt.0.) then
                km=N+2    !<-- above surface
              elseif (dpth*(zz(0)-z_lev(m)).gt.0.) then
                km=-2     !<-- below bottom
              else
                km=-1     !--> to search
              endif
          else
              dpth=zz(N+1)-zz(0)     
              if (dpth*(z_lev(m)-zz(N+1)).gt.0.) then
                km=N+2    !<-- above surface
              elseif (dpth*(z_lev(m)-zz(N)).gt.0.) then
                km=N      !<-- below surface, but above z_r(N)
              elseif (dpth*(zz(0)-below-z_lev(m)).gt.0.) then
                km=-3     !<-- below bottom 
              elseif (dpth*(zz(0)-z_lev(m)).gt.0.) then
                km=-2     !<-- below bottom but close 
              elseif (dpth*(zz(1)-z_lev(m)).gt.0.) then
                km=0      !<-- above bottom, but below z_r(1)
              else
                km=-1     !--> to search
              endif
          endif
          do k=N-1,kmin,-1
              if (km.eq.-1) then
                if((zz(k+1)-z_lev(m))*(z_lev(m)-zz(k))
     & .ge. 0.) km=k
              endif
          enddo

            if (km.eq.-3) then
              var_zlv(m)=FillValue             !<-- masked out
            elseif (km.eq.-2) then
#ifdef INTERP_BELOW
              var_zlv(m)=var(kmin)  !    
#ifdef CUBIC
     &           -FC(kmin)*(zz(kmin)-z_lev(m))
#else
     &          -(zz(kmin)-z_lev(m))*(var(kmin+1)-var(kmin))
     &                                  /(zz(kmin+1)-zz(kmin))
#endif
#else
              var_zlv(m)=FillValue      !<-- below bottom
#endif
            elseif (km.eq.N+2) then
#ifdef INTERP_ABOVE
              var_zlv(m)=var(N)     !-> R-point, above z_r(N)
#ifdef CUBIC
     &              +FC(N)*(z_lev(m)-zz(N))
#else
     &              +(z_lev(m)-zz(N))*(var(N)-var(N-1))
     &                                   /(zz(N)-zz(N-1))
#endif
#else
              var_zlv(m)=-FillValue     !<-- above surface
#endif
            elseif (km.eq.N) then
              var_zlv(m)=var(N)     !-> R-point, above z_r(N)
#ifdef CUBIC
     &              +FC(N)*(z_lev(m)-zz(N))
#else
     &              +(z_lev(m)-zz(N))*(var(N)-var(N-1))
     &                                   /(zz(N)-zz(N-1))
#endif
            elseif (km.eq.kmin-1) then   !-> R-point below z_r(1),
              var_zlv(m)=var(kmin)  !     but above bottom  
#ifdef CUBIC
     &           -FC(kmin)*(zz(kmin)-z_lev(m))
#else
     &          -(zz(kmin)-z_lev(m))*(var(kmin+1)-var(kmin))
     &                                  /(zz(kmin+1)-zz(kmin))
#endif
            else
              k=km
              !write(*,*) k,km
#ifdef CUBIC
              cff=1.D0/(zz(k+1)-zz(k))
              p=z_lev(m)-zz(k)
              q=zz(k+1)-z_lev(m)

              var_zlv(m)=cff*( q*var(k) + p*var(k+1)
     &             -cff*p*q*( cff*(q-p)*(var(k+1)-var(k))
     &                                  +p*FC(k+1) -q*FC(k) )
     &                                                          )
#else
              var_zlv(m)=( var(k)*(zz(k+1)-z_lev(m))
     &                        +var(k+1)*(z_lev(m)-zz(k))
     &                               )/(zz(k+1)-zz(k))
#endif



            endif
        enddo ! <-- m


      return
      end
 
