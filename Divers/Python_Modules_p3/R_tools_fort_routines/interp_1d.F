

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Z interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine interp_1d(N, var,z_r,z_w,z_lev,var_zlv,kmin,verbose)

      implicit none

      integer k, N, km, verbose
      integer kmin

      real*8 var(kmin:N), z_r(N), z_w(0:N), z(0:N+1),
     &       z_lev, var_zlv, dpth,
     &       dz(kmin-1:N), FC(kmin-1:N),
     &       p,q,cff

      real*8, parameter :: FillValue = -9999.

Cf2py intent(in) N,var,z_r,z_w,z_lev,kmin,verbose
Cf2py intent(out) var_zlv



        !if (verbose.eq.1) write(*,*) 'depth was ', z_lev
        !if (z_lev.gt.0) write(*,*) 'depth is ', z_lev  !z_lev = -1*z_lev
        !if (verbose.eq.1) write(*,*) 'depth is ', z_lev


        if ((z_r(N).eq.0).and.(z_w(N).eq.0)) then
            !write(*,*) 'this point should be masked'
            var_zlv = FillValue

        else



!         if (verbose.eq.1) then
!             if (z_lev.gt.z_w(N)) then
!                 write(*,*) 'depths are ', z_lev, z_w(N)
!                 write(*,*) z_w
!                 write(*,*) z_r
!                 z_lev=z_w(N)
!             endif
!         endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            do k=1,N
                z(k)=z_r(k)
            enddo
              z(0)=z_w(0)
              z(N+1)=z_w(N)

        if (verbose.eq.1) then
           write(*,*)  ' '
           write(*,*) 'z',z
           write(*,*)  ' '
           write(*,*) 'z_r',z_r
           write(*,*)  ' '
           write(*,*) 'z_w',z_w
           write(*,*)  ' '

        endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



#ifdef CUBIC
        do k=kmin,N-1
            dz(k)=z(k+1)-z(k)
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if (kmin.eq.0) then     
                     
            dpth=z(N)-z(0) 

            if (dpth*(z_lev-z(N)).gt.0.) then
               km=N+2    !<-- above surface
            elseif (dpth*(z(0)-z_lev).gt.0.) then
               km=-2     !<-- below bottom
            else
               km=-1     !--> to search
            endif

        else

            dpth=z(N+1)-z(0)    

        if (verbose.eq.1) then
           write(*,*) 'dpth',dpth
        endif


            if (dpth*(z_lev-z(N+1)).gt.0.) then
                km=N+2    !<-- above surface
            elseif (dpth*(z_lev-z(N)).gt.0.) then
                km=N      !<-- below surface, but above z_r(N)
            elseif (dpth*(z(0)-z_lev).gt.0.) then
                km=-2     !<-- below bottom
            elseif (dpth*(z(1)-z_lev).gt.0.) then
                km=0      !<-- above bottom, but below z_r(1)
            else
                km=-1     !--> to search
            endif


        if (verbose.eq.1) then
           write(*,*) 'dpth*(z_lev-z(N+1)',dpth*(z_lev-z(N+1))
           write(*,*) 'dpth*(z_lev-z(N))',dpth*(z_lev-z(N))
           write(*,*) 'dpth*(z(0)-z_lev)',dpth*(z(0)-z_lev)
           write(*,*) 'dpth*(z(1)-z_lev)',dpth*(z(1)-z_lev)
           write(*,*) 'km',km
        endif

        endif

        do k=N-1,kmin,-1

!            write(*,*) 'test', k, z(k+1), z(k), z_lev
              if (km.eq.-1) then
                if ((z(k+1)-z_lev)*(z_lev-z(k)) .ge. 0.) then
                    km=k
                endif
              endif
        enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


          if (km.eq.-3) then
              var_zlv=0.             !<-- masked out
          elseif (km.eq.-2) then
#ifdef INTERP_BELOW
              var_zlv=var(kmin)  !     but above bottom  
#ifdef CUBIC
     &           -FC(kmin)*(z(kmin)-z_lev)
#else
     &          -(z(kmin)-z_lev)*(var(kmin+1)-var(kmin))
     &                                  /(z(kmin+1)-z(kmin))
#endif
#else
              var_zlv=FillValue       !<-- below bottom
#endif
            elseif (km.eq.N+2) then
#ifdef INTERP_ABOVE
              var_zlv=var(N)  !     but above bottom  
#ifdef CUBIC
     &              +FC(N)*(z_lev-z(N))
#else
     &              +(z_lev-z(N))*(var(N)-var(N-1))
     &                                   /(z(N)-z(N-1))
#endif
#else
              var_zlv=FillValue       !<-- above surface
#endif

            elseif (km.eq.N) then
              var_zlv=var(N)     !-> R-point, above z_r(N)
#ifdef CUBIC
     &              +FC(N)*(z_lev-z(N))
#else
     &              +(z_lev-z(N))*(var(N)-var(N-1))
     &                                   /(z(N)-z(N-1))
#endif
            elseif (km.eq.kmin-1) then   !-> R-point below z_r(1),
              var_zlv=var(kmin)  !     but above bottom  
#ifdef CUBIC
     &           -FC(kmin)*(z(kmin)-z_lev)
#else
     &          -(z(kmin)-z_lev)*(var(kmin+1)-var(kmin))
     &                                  /(z(kmin+1)-z(kmin))
#endif
            else
              k=km


#ifdef CUBIC
              cff=1.D0/(z(k+1)-z(k))
              p=z_lev-z(k)
              q=z(k+1)-z_lev

              !write(*,*) 'koko',cff,p,q,var(k),var(k+1)

              var_zlv=cff*( q*var(k) + p*var(k+1)
     &             -cff*p*q*( cff*(q-p)*(var(k+1)-var(k))
     &                                  +p*FC(k+1) -q*FC(k) )
     &                                                          )
#else
              var_zlv=( var(k)*(z(k+1)-z_lev)
     &                        +var(k+1)*(z_lev-z(k))
     &                               )/(z(k+1)-z(k))
#endif


        endif

        if (verbose.eq.1) then
           write(*,*) 'var_zlv', var_zlv
        endif


        endif


      return
      end


