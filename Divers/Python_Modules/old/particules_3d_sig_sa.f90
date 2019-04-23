# 1 "particules_3d_sig_sa.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "particules_3d_sig_sa.F"

!--------------------------------------------------------------------!!
! cpp particules_3d_sig_sa.F particules_3d_sig_sa.f90; f2py --fcompiler=intelem --compiler=intelem --f90flags="-extend_source" -DF2PY_REPORT_ON_ARRAY_COPY=1 -c -m particules_3d_sig_sa particules_3d_sig_sa.f90
!
!---------------------------------------------------------------------




!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Spatial interpolation routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# 1 "interp_3d_sig_sa.F" 1

!----------------------------------------------------------------------------------------------
! #define NEW_VERSION



!#define LINEAR_INTERPOLATION
! #define CUBIC_INTERPOLATION
! #define CRSPL_INTERPOLATION
! #define WENO_INTERPOLATION

!----------------------------------------------------------------------------------------------
# 405 "interp_3d_sig_sa.F"
!---------------------------------------------------------------------!
! Compute displacement given u,v,w at particule position
! with linear interpolation in space and time
!---------------------------------------------------------------------
!
!

       subroutine advance_3d(px,py,pz,u,v,w,itim,fct,pm,pn,
     & dz,dt,i0,j0,k0,nx,ny,nz,np,dpx,dpy,dpz)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz
       integer(kind=4) ,intent(in) :: np
       integer(kind=4) ,intent(in) :: i0,j0,k0
       integer(kind=4) ,dimension(0:1) ,intent(in) :: itim
       real(kind=8) ,dimension(nx-1,ny,nz,2),intent(in) :: u
       real(kind=8) ,dimension(nx,ny-1,nz,2),intent(in) :: v
       real(kind=8) ,dimension(nx,ny,nz+1,2),intent(in) :: w
       real(kind=8) ,dimension(nx,ny,nz,2),intent(in) :: dz
       real(kind=8) ,dimension(nx,ny),intent(in) :: pm,pn
       real(kind=8) ,intent(in) :: px,py,pz
       real(kind=8) ,intent(out) :: dpx,dpy,dpz
       real(kind=8) ,intent(in) :: dt,fct
! local
       integer(kind=8) :: i,j,k
       integer(kind=8) :: i_u,j_v,k_w
       real(kind=8) ,dimension(2,2,2,2) :: wt4
       real(kind=8) ,dimension(2,2,2) :: wt3
       real(kind=8) ,dimension(2,2) :: wt2
       real(kind=8) :: fcx,fcy,fcz,fctl
       real(kind=8) :: fcx_u,fcy_v,fcz_w
       real(kind=8) :: pu,pv,pw,pdz,ppm,ppn

!f2py intent(in) u,v,w
!f2py intent(in) px,py,pz
!f2py intent(out) dpx,dpy,dpz
!f2py intent(in) pm,pn,dt,fct
!f2py intent(in) dz
!f2py intent(in) nx,ny,nz
!f2py intent(in) i0,j0,k0
!f2py intent(in) itim
!f2py intent(in) np


           !---------------------------------------
           ! Old version where px,py base on rho-grid indices
           !---------------------------------------

! i = max(1,min(floor(px+1)-i0,nx-1))
! j = max(1,min(floor(py+1)-j0,ny-1))
! k = max(1,min(floor(pz+1-0.5)-k0,nz-1))
! i_u = max(1,min(floor(px+1-0.5)-i0,nx-2))
! j_v = max(1,min(floor(py+1-0.5)-j0,ny-2))
! k_w = max(1,min(floor(pz+1)-k0,nz))

           i = max(1,min(floor(px+1+0.5)-i0,nx-1))
           j = max(1,min(floor(py+1+0.5)-j0,ny-1))
           k = max(1,min(floor(pz+1-0.5)-k0,nz-1))

           i_u = max(1,min(floor(px+1)-i0,nx-2))
           j_v = max(1,min(floor(py+1)-j0,ny-2))
           k_w = max(1,min(floor(pz+1)-k0,nz))

           !---------------------------------------
           ! 1. Linear interpolation in space and time
           !---------------------------------------

           fcx = px+1+0.5 - i - i0;
           fcy = py+1+0.5 - j - j0;
           fcz = pz+1-0.5 - k - k0;

           fcx_u = px+1 - i - i0;
           fcy_v = py+1 - j - j0;
           fcz_w = pz+1 - k - k0;


           fctl = fct
           if (itim(0).eq.1) fctl = 1-fct

           !---------------------------------------
           ! Compute velocities and level depth at particle position
           !---------------------------------------

           CALL linear_4d(fcx_u,fcy,fcz,fctl,wt4)
           pu = sum(u(i_u:i_u+1,j:j+1,k:k+1,:)*wt4)

           CALL linear_4d(fcx,fcy_v,fcz,fctl,wt4)
           pv = sum(v(i:i+1,j_v:j_v+1,k:k+1,:)*wt4)

           CALL linear_4d(fcx,fcy,fcz_w,fctl,wt4)
           pw = sum(w(i:i+1,j:j+1,k_w:k_w+1,:)*wt4)

           CALL linear_4d(fcx,fcy,fcz,fctl,wt4)
           pdz = sum(dz(i:i+1,j:j+1,k:k+1,:)*wt4)

           if (pdz.gt.0) pdz=1./pdz

           CALL linear_2d(fcx,fcy,wt2)
           ppm = sum(pm(i:i+1,j:j+1)*wt2)
           ppn = sum(pn(i:i+1,j:j+1)*wt2)

           !---------------------------------------
           ! Update position
           !---------------------------------------

           dpx = dt*pu*ppm
           dpy = dt*pv*ppn
           dpz = dt*pw*pdz
! write(*,*) 'pxyz',px, py, pz
!
! write(*,*) 'pw, dpz', pw, dpz
!
       end






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get interpolation matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       subroutine linear_4d(fcx,fcy,fcz,fct,wt)
       implicit none
! import/export
       real(kind=8) ,intent(in) :: fcx,fcy,fcz,fct
       real(kind=8) ,dimension(2,2,2,2),intent(out) :: wt

!f2py intent(in) fcx,fcy,fcz,fct
!f2py intent(out) wt


           wt(1,1,1,1) = (1-fcz)*(1-fcy)*(1-fcx)*(1-fct);
           wt(1,1,2,1) = fcz *(1-fcy)*(1-fcx)*(1-fct);
           wt(1,2,1,1) = (1-fcz)* fcy *(1-fcx)*(1-fct);
           wt(1,2,2,1) = fcz * fcy *(1-fcx)*(1-fct);
           wt(2,1,1,1) = (1-fcz)*(1-fcy)* fcx *(1-fct);
           wt(2,1,2,1) = fcz *(1-fcy)* fcx *(1-fct);
           wt(2,2,1,1) = (1-fcz)* fcy * fcx *(1-fct);
           wt(2,2,2,1) = fcz * fcy * fcx *(1-fct);

           wt(1,1,1,2) = (1-fcz)*(1-fcy)*(1-fcx)* fct ;
           wt(1,1,2,2) = fcz *(1-fcy)*(1-fcx)* fct ;
           wt(1,2,1,2) = (1-fcz)* fcy *(1-fcx)* fct ;
           wt(1,2,2,2) = fcz * fcy *(1-fcx)* fct ;
           wt(2,1,1,2) = (1-fcz)*(1-fcy)* fcx * fct ;
           wt(2,1,2,2) = fcz *(1-fcy)* fcx * fct ;
           wt(2,2,1,2) = (1-fcz)* fcy * fcx * fct ;
           wt(2,2,2,2) = fcz * fcy * fcx * fct ;

       end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get interpolation matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       subroutine linear_3d(fcx,fcy,fcz,wt)
       implicit none
! import/export
       real(kind=8) ,intent(in) :: fcx,fcy,fcz
       real(kind=8) ,dimension(2,2,2),intent(out) :: wt

!f2py intent(in) fcx,fcy,fcz
!f2py intent(out) wt


           wt(1,1,1) = (1-fcz)*(1-fcy)*(1-fcx);
           wt(1,1,2) = fcz *(1-fcy)*(1-fcx);
           wt(1,2,1) = (1-fcz)* fcy *(1-fcx);
           wt(1,2,2) = fcz * fcy *(1-fcx);
           wt(2,1,1) = (1-fcz)*(1-fcy)* fcx ;
           wt(2,1,2) = fcz *(1-fcy)* fcx ;
           wt(2,2,1) = (1-fcz)* fcy * fcx ;
           wt(2,2,2) = fcz * fcy * fcx ;


       end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get interpolation matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       subroutine linear_2d(fcx,fcy,wt)
       implicit none
! import/export
       real(kind=8) ,intent(in) :: fcx,fcy
       real(kind=8) ,dimension(2,2),intent(out) :: wt

!f2py intent(in) fcx,fcy
!f2py intent(out) wt


           wt(1,1) = (1-fcy)*(1-fcx);
           wt(1,2) = (1-fcy)* fcx;
           wt(2,1) = fcy *(1-fcx);
           wt(2,2) = fcy * fcx;


       end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate T,S at each particle position (same than interp_3d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_3d_ts(pvar1,pvar2,px,py,pz,
     & var1,var2,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz
       integer(kind=4) ,intent(in) :: np
       integer(kind=4) ,intent(in) :: npmx
       integer(kind=4) ,intent(in) :: i0,j0,k0
       real(kind=8) ,dimension(nx,ny,nz),intent(in) :: var1,var2
       real(kind=8) ,dimension(np) ,intent(out):: pvar1,pvar2
       real(kind=8) ,dimension(np) ,intent(in):: px,py,pz
! local
       integer(kind=8) :: ip,jp,kp,i,j,k
       real(kind=8) ,dimension(2,2,2) :: wt3
       real(kind=8) :: fcx,fcy,fcz

!f2py intent(out) pvar1,pvar2
!f2py intent(in) var1,var2
!f2py intent(in) px,py,pz
!f2py intent(in) npmx
!f2py intent(in) nx,ny,nz
!f2py intent(in) i0,j0,k0
!f2py intent(in) np

       do ip = 1,np
         if (.not.(isnan(px(ip))) ) then

           i = max(1,min(floor(px(ip)+1+0.5)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+0.5)-j0,ny-1))
           k = max(1,min(floor(pz(ip)+1-0.5)-k0,nz-1))


           fcx = px(ip)+1+0.5 - i - i0;
           fcy = py(ip)+1+0.5 - j - j0;
           fcz = pz(ip)+1-0.5 - k - k0;

           CALL linear_3d(fcx,fcy,fcz,wt3)
           pvar1(ip) = sum(var1(i:i+1,j:j+1,k:k+1)*wt3)
           pvar2(ip) = sum(var2(i:i+1,j:j+1,k:k+1)*wt3)

         endif
       enddo
!
       end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate a 3D variable at each particle position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_3d(pvar1,px,py,pz,var1,npmx,i0,j0,k0,
     & nx,ny,nz,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz
       integer(kind=4) ,intent(in) :: np
       integer(kind=4) ,intent(in) :: npmx
       integer(kind=4) ,intent(in) :: i0,j0,k0
       real(kind=8) ,dimension(nx,ny,nz),intent(in) :: var1
       real(kind=8) ,dimension(np) ,intent(out):: pvar1
       real(kind=8) ,dimension(np) ,intent(in):: px,py,pz
! local
       integer(kind=8) :: ip,jp,kp,i,j,k
       real(kind=8) ,dimension(2,2,2) :: wt3
       real(kind=8) :: fcx,fcy,fcz

!f2py intent(out) pvar1
!f2py intent(in) var1
!f2py intent(in) px,py,pz
!f2py intent(in) npmx
!f2py intent(in) nx,ny,nz
!f2py intent(in) i0,j0,k0
!f2py intent(in) np

       do ip = 1,np
         if (.not.(isnan(px(ip))) ) then

           i = max(1,min(floor(px(ip)+1+0.5)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+0.5)-j0,ny-1))
           k = max(1,min(floor(pz(ip)+1-0.5)-k0,nz-1))

           fcx = px(ip)+1+0.5 - i - i0;
           fcy = py(ip)+1+0.5 - j - j0;
           fcz = pz(ip)+1-0.5 - k - k0;

           CALL linear_3d(fcx,fcy,fcz,wt3)
           pvar1(ip) = sum(var1(i:i+1,j:j+1,k:k+1)*wt3)

         endif
       enddo
!
       end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate a 3D variable at each particle position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_3d_psi(pvar1,px,py,pz,var1,npmx,i0,j0,k0,
     & nx,ny,nz,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz
       integer(kind=4) ,intent(in) :: np
       integer(kind=4) ,intent(in) :: npmx
       integer(kind=4) ,intent(in) :: i0,j0,k0
       real(kind=8) ,dimension(nx,ny,nz),intent(in) :: var1
       real(kind=8) ,dimension(np) ,intent(out):: pvar1
       real(kind=8) ,dimension(np) ,intent(in):: px,py,pz
! local
       integer(kind=8) :: ip,jp,kp,i,j,k
       real(kind=8) ,dimension(2,2,2) :: wt3
       real(kind=8) :: fcx,fcy,fcz

!f2py intent(out) pvar1
!f2py intent(in) var1
!f2py intent(in) px,py,pz
!f2py intent(in) npmx
!f2py intent(in) nx,ny,nz
!f2py intent(in) i0,j0,k0
!f2py intent(in) np

       do ip = 1,np
         if (.not.(isnan(px(ip))) ) then

           i = max(1,min(floor(px(ip)+1)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1)-j0,ny-1))
           k = max(1,min(floor(pz(ip)+1-0.5)-k0,nz-1))

           fcx = px(ip)+1 - i - i0;
           fcy = py(ip)+1 - j - j0;
           fcz = pz(ip)+1-0.5 - k - k0;

           CALL linear_3d(fcx,fcy,fcz,wt3)
           pvar1(ip) = sum(var1(i:i+1,j:j+1,k:k+1)*wt3)

         endif
       enddo
!
       end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate a 3D variable at each particle position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_3d_w(pvar1,px,py,pz,var1,npmx,i0,j0,k0,
     & nx,ny,nz,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz
       integer(kind=4) ,intent(in) :: np
       integer(kind=4) ,intent(in) :: npmx
       integer(kind=4) ,intent(in) :: i0,j0,k0
       real(kind=8) ,dimension(nx,ny,nz),intent(in) :: var1
       real(kind=8) ,dimension(np) ,intent(out):: pvar1
       real(kind=8) ,dimension(np) ,intent(in):: px,py,pz
! local
       integer(kind=8) :: ip,jp,kp,i,j,k
       real(kind=8) ,dimension(2,2,2) :: wt3
       real(kind=8) :: fcx,fcy,fcz

!f2py intent(out) pvar1
!f2py intent(in) var1
!f2py intent(in) px,py,pz
!f2py intent(in) npmx
!f2py intent(in) nx,ny,nz
!f2py intent(in) i0,j0,k0
!f2py intent(in) np

       do ip = 1,np
         if (.not.(isnan(px(ip))) ) then

           i = max(1,min(floor(px(ip)+1+0.5)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+0.5)-j0,ny-1))
           k = max(1,min(floor(pz(ip)+1)-k0,nz-1))

           fcx = px(ip)+1+0.5 - i - i0;
           fcy = py(ip)+1+0.5 - j - j0;
           fcz = pz(ip)+1 - k - k0;

           CALL linear_3d(fcx,fcy,fcz,wt3)
           pvar1(ip) = sum(var1(i:i+1,j:j+1,k:k+1)*wt3)

         endif
       enddo
!
       end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate a 3D variable at each particle position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_3d_psiw(pvar1,px,py,pz,var1,npmx,i0,j0,k0,
     & nx,ny,nz,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz
       integer(kind=4) ,intent(in) :: np
       integer(kind=4) ,intent(in) :: npmx
       integer(kind=4) ,intent(in) :: i0,j0,k0
       real(kind=8) ,dimension(nx,ny,nz),intent(in) :: var1
       real(kind=8) ,dimension(np) ,intent(out):: pvar1
       real(kind=8) ,dimension(np) ,intent(in):: px,py,pz
! local
       integer(kind=8) :: ip,jp,kp,i,j,k
       real(kind=8) ,dimension(2,2,2) :: wt3
       real(kind=8) :: fcx,fcy,fcz

!f2py intent(out) pvar1
!f2py intent(in) var1
!f2py intent(in) px,py,pz
!f2py intent(in) npmx
!f2py intent(in) nx,ny,nz
!f2py intent(in) i0,j0,k0
!f2py intent(in) np

       do ip = 1,np
         if (.not.(isnan(px(ip))) ) then

           i = max(1,min(floor(px(ip)+1)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1)-j0,ny-1))
           k = max(1,min(floor(pz(ip)+1)-k0,nz-1))

           fcx = px(ip)+1 - i - i0;
           fcy = py(ip)+1 - j - j0;
           fcz = pz(ip)+1 - k - k0;

           CALL linear_3d(fcx,fcy,fcz,wt3)
           pvar1(ip) = sum(var1(i:i+1,j:j+1,k:k+1)*wt3)

         endif
       enddo
!
       end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate a 2D variable at each horizontal particle position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_2d(pvar1,px,py,var1,npmx,i0,j0,nx,ny,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny
       integer(kind=4) ,intent(in) :: np
       integer(kind=4) ,intent(in) :: npmx
       integer(kind=4) ,intent(in) :: i0,j0
       real(kind=8) ,dimension(nx,ny),intent(in) :: var1
       real(kind=8) ,dimension(np) ,intent(out):: pvar1
       real(kind=8) ,dimension(np) ,intent(in):: px,py
! local
       integer(kind=8) :: ip,jp,i,j
       real(kind=8) ,dimension(2,2) :: wt
       real(kind=8) :: fcx,fcy

!f2py intent(out) pvar1
!f2py intent(in) var1
!f2py intent(in) px,py
!f2py intent(in) npmx
!f2py intent(in) nx,ny
!f2py intent(in) i0,j0
!f2py intent(in) np


       do ip = 1,np
         if (.not.(isnan(px(ip))) ) then

           i = max(1,min(floor(px(ip)+1+0.5)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+0.5)-j0,ny-1))

           fcx = px(ip)+1+0.5 - i - i0;
           fcy = py(ip)+1+0.5 - j - j0;

           CALL linear_2d(fcx,fcy,wt)

           pvar1(ip) = sum(var1(i:i+1,j:j+1)*wt)

         endif
       enddo
!
       end
# 17 "particules_3d_sig_sa.F" 2




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Time-stepping routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       !---------------------------------------------------------------------!
       ! Forward Euler time-stepping
       !---------------------------------------------------------------------

       subroutine timestep_FE(px,py,pz,u,v,w,itim,fct,dfct,
     & pm,pn,dz,dt,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz
       integer(kind=4) ,intent(in) :: np
       integer(kind=4) ,intent(in) :: npmx
       integer(kind=4) ,intent(in) :: i0,j0,k0
       integer(kind=4) ,dimension(0:1) ,intent(in) :: itim
       real(kind=8) ,dimension(nx-1,ny,nz,2),intent(in) :: u
       real(kind=8) ,dimension(nx,ny-1,nz,2),intent(in) :: v
       real(kind=8) ,dimension(nx,ny,nz+1,2),intent(in) :: w
       real(kind=8) ,dimension(nx,ny,nz,2),intent(in) :: dz
       real(kind=8) ,dimension(nx,ny),intent(in) :: pm,pn
       real(kind=8) ,dimension(np) ,intent(inout):: px,py,pz
       real(kind=8) ,intent(in) :: dt,fct,dfct
! local
       real(kind=8) :: dpx,dpy,dpz
       integer(kind=8) :: ip

!f2py intent(in) u,v,w
!f2py intent(inout) px,py,pz
!f2py intent(in) pm,pn,dt
!f2py intent(in) dz
!f2py intent(in) npmx
!f2py intent(in) nx,ny,nz
!f2py intent(in) i0,j0,k0
!f2py intent(in) itim
!f2py intent(in) np

       do ip = 1,min(npmx,np)

         if (.not.(isnan(px(ip))) ) then

           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,itim,fct,pm,pn
     & ,dz,dt,i0,j0,k0,nx,ny,nz,np,dpx,dpy,dpz)

           !---------------------------------------
           ! Update position
           !---------------------------------------

           px(ip) = px(ip) + dpx
           py(ip) = py(ip) + dpy
           pz(ip) = pz(ip) + dpz

         endif

       enddo


       end





       !---------------------------------------------------------------------!
       ! Runge-Kutta2 time-stepping
       !---------------------------------------------------------------------

       subroutine timestep_RK2(px,py,pz,u,v,w,itim,fct,dfct,
     & pm,pn,dz,dt,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz
       integer(kind=4) ,intent(in) :: np
       integer(kind=4) ,intent(in) :: npmx
       integer(kind=4) ,intent(in) :: i0,j0,k0
       integer(kind=4) ,dimension(0:1) ,intent(in) :: itim
       real(kind=8) ,dimension(nx-1,ny,nz,2),intent(in) :: u
       real(kind=8) ,dimension(nx,ny-1,nz,2),intent(in) :: v
       real(kind=8) ,dimension(nx,ny,nz+1,2),intent(in) :: w
       real(kind=8) ,dimension(nx,ny,nz,2),intent(in) :: dz
       real(kind=8) ,dimension(nx,ny),intent(in) :: pm,pn
       real(kind=8) ,dimension(np) ,intent(inout):: px,py,pz
       real(kind=8) ,intent(in) :: dt,fct,dfct
! local
       real(kind=8) :: dpx,dpy,dpz
       integer(kind=8) :: ip


!f2py intent(inout) px,py,pz
!f2py intent(in) u,v,w
!f2py intent(in) fct,dfct,dt
!f2py intent(in) pm,pn,dz
!f2py intent(in) npmx
!f2py intent(in) nx,ny,nz
!f2py intent(in) i0,j0,k0
!f2py intent(in) itim
!f2py intent(in) np

       do ip = 1,min(npmx,np)

         if (.not.(isnan(px(ip))) ) then

           ! Forward Euler / predictor
           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,itim,
     & fct,pm,pn,dz,dt,i0,j0,k0,nx,ny,nz,np,dpx,dpy,dpz)

           ! midpoint rule / corrector
           CALL advance_3d(px(ip)+0.5*dpx,py(ip)+0.5*dpy,pz(ip)+0.5*dpz,
     & u,v,w,itim,fct+0.5*dfct,pm,pn,dz,dt,i0,j0,k0,
     & nx,ny,nz,np,dpx,dpy,dpz)

           !---------------------------------------
           ! Update position
           !---------------------------------------

           px(ip) = px(ip) + dpx
           py(ip) = py(ip) + dpy
           pz(ip) = pz(ip) + dpz

         endif

       enddo


       end



       !---------------------------------------------------------------------!
       ! Runge-Kutta 4 time-stepping
       !---------------------------------------------------------------------

!

       subroutine timestep_RK4(px,py,pz,u,v,w,itim,fct,dfct,pm,pn,
     & dz,dt,npmx,i0,j0,k0,nx,ny,nz,np,dpxi,dpyi,dpzi)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz
       integer(kind=4) ,intent(in) :: np
       integer(kind=4) ,intent(in) :: npmx
       integer(kind=4) ,intent(in) :: i0,j0,k0
       integer(kind=4) ,dimension(0:1) ,intent(in) :: itim
       real(kind=8) ,dimension(nx-1,ny,nz,2),intent(in) :: u
       real(kind=8) ,dimension(nx,ny-1,nz,2),intent(in) :: v
       real(kind=8) ,dimension(nx,ny,nz+1,2),intent(in) :: w
       real(kind=8) ,dimension(nx,ny,nz,2),intent(in) :: dz
       real(kind=8) ,dimension(nx,ny),intent(in) :: pm,pn
       real(kind=8) ,dimension(np) ,intent(inout):: px,py,pz
       real(kind=8) ,dimension(np) ,intent(out) :: dpxi,dpyi,dpzi
       real(kind=8) ,intent(in) :: dt,fct,dfct
! local
       real(kind=8) ,dimension(0:3) :: dpx,dpy,dpz
       real(kind=8) :: coef
       integer(kind=8) :: ip

!f2py intent(inout) px,py,pz
!f2py intent(in) u,v,w
!f2py intent(in) fct,dfct,dt
!f2py intent(in) pm,pn,dz
!f2py intent(in) npmx
!f2py intent(in) nx,ny,nz
!f2py intent(in) i0,j0,k0
!f2py intent(in) itim
!f2py intent(in) np
!f2py intent(out) dpxi,dpyi,dpzi

       coef = 1./6.

       do ip = 1,min(npmx,np)

         if (.not.(isnan(px(ip))) ) then

           ! Forward Euler / predictor
           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,
     & itim,fct,pm,pn,dz,dt,i0,j0,k0,
     & nx,ny,nz,np,dpx(0),dpy(0),dpz(0))

           ! backward Euler / corrector
           CALL advance_3d(px(ip)+0.5*dpx(0),py(ip)+0.5*dpy(0),
     & pz(ip)+0.5*dpz(0),
     & u,v,w,itim,fct+0.5*dfct,pm,pn,dz,dt,i0,j0,k0,
     & nx,ny,nz,np,dpx(1),dpy(1),dpz(1))

           ! midpoint rule / predictor
           CALL advance_3d(px(ip)+0.5*dpx(1),py(ip)+0.5*dpy(1),
     & pz(ip)+0.5*dpz(1),
     & u,v,w,itim,fct+0.5*dfct,pm,pn,dz,dt,i0,j0,k0,
     & nx,ny,nz,np,dpx(2),dpy(2),dpz(2))

           ! Corrector
           CALL advance_3d(px(ip)+dpx(2),py(ip)+dpy(2),
     & pz(ip)+dpz(2),
     & u,v,w,itim,fct+1.*dfct,pm,pn,dz,dt,i0,j0,k0,
     & nx,ny,nz,np,dpx(3),dpy(3),dpz(3))

           !---------------------------------------
           ! Update position
           !---------------------------------------


           dpxi(ip) = coef * (dpx(0) + 2 * dpx(1) + 2 * dpx(2) + dpx(3))
           dpyi(ip) = coef * (dpy(0) + 2 * dpy(1) + 2 * dpy(2) + dpy(3))
           dpzi(ip) = coef * (dpz(0) + 2 * dpz(1) + 2 * dpz(2) + dpz(3))

           if (abs(dpyi(ip)).gt.200) then
                write(*,*) ' '
                write(*,*) 'before',px(ip),py(ip),pz(ip)
                write(*,*) dpx
                write(*,*) dpy
                write(*,*) dpz
                write(*,*) ' '
           endif

           px(ip) = px(ip) + dpxi(ip)
           py(ip) = py(ip) + dpyi(ip)
           pz(ip) = pz(ip) + dpzi(ip)

           if (abs(dpyi(ip)).gt.200) then
                write(*,*) 'after',dpxi(ip),dpyi(ip),dpzi(ip)
                write(*,*) 'after',px(ip),py(ip),pz(ip)
           endif


         endif

       enddo


       end


        !---------------------------------------------------------------------!
        ! Adams-Bashforth 2 time-stepping
        !---------------------------------------------------------------------
        ! For Adams-methods we need to keep the previous dpx,dpy,dpz values
        !

       subroutine timestep_AB2(px,py,pz,dpx,dpy,dpz,iab,
     & u,v,w,itim,fct,dfct,pm,pn,dz,dt,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none

! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz,np,npmx
       integer(kind=4) ,intent(in) :: i0,j0,k0
       integer(kind=4) ,dimension(0:1) ,intent(in) :: itim
       real(kind=8),dimension(nx-1,ny,nz,2),intent(in) :: u
       real(kind=8),dimension(nx,ny-1,nz,2),intent(in) :: v
       real(kind=8),dimension(nx,ny,nz+1,2),intent(in) :: w
       real(kind=8),dimension(nx,ny,nz,2) ,intent(in) :: dz
       real(kind=8),dimension(nx,ny) ,intent(in) :: pm,pn
       real(kind=8) ,dimension(np) ,intent(inout):: px,py,pz
       real(kind=8) ,intent(in) :: dt,fct,dfct

! import/export for AB scheme
       integer(kind=4),dimension(0:1) ,intent(in):: iab
       real(kind=8) ,dimension(np,0:1) ,intent(inout):: dpx,dpy,dpz

! local
       real(kind=8) :: coef
       integer(kind=8) :: ip
       real(kind=8) :: dpxi,dpyi,dpzi


!f2py intent(inout) px,py,pz
!f2py intent(in) u,v,w
!f2py intent(in) fct,dfct,dt
!f2py intent(in) pm,pn,dz
!f2py intent(in) npmx
!f2py intent(in) nx,ny,nz
!f2py intent(in) i0,j0,k0
!f2py intent(in) itim
!f2py intent(in) np
!f2py intent(in) iab
!f2py intent(inout) dpx,dpy,dpz


       coef = 0.5

       do ip = 1,min(npmx,np)
         if (.not.(isnan(px(ip))) ) then

           ! Forward Euler / predictor
           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,itim,
     & fct,pm,pn,dz,dt,i0,j0,k0,nx,ny,nz,np,
     & dpx(ip,iab(1)),dpy(ip,iab(1)),dpz(ip,iab(1)))

!
! if (.not.(isnan(dpx(ip,iab(0)))) ) then

            dpxi = coef * (3. * dpx(ip,iab(1)) - dpx(ip,iab(0)))
            dpyi = coef * (3. * dpy(ip,iab(1)) - dpy(ip,iab(0)))
            dpzi = coef * (3. * dpz(ip,iab(1)) - dpz(ip,iab(0)))
!
! else
!
! dpx(ip,iab(1)) = dpxi
! dpy(ip,iab(1)) = dpyi
! dpz(ip,iab(1)) = dpzi
!
! endif

            !---------------------------------------
            ! Update position
            !---------------------------------------

            px(ip) = px(ip) + dpxi
            py(ip) = py(ip) + dpyi
            pz(ip) = pz(ip) + dpzi

         endif

       enddo


       end


        !---------------------------------------------------------------------!
        ! Adams-Bashforth 3 time-stepping
        !---------------------------------------------------------------------
        ! For Adams-methods we need to keep the previous dpx,dpy,dpz values
        !

       subroutine timestep_AB3(px,py,pz,dpx,dpy,dpz,iab,
     & u,v,w,itim,fct,dfct,pm,pn,dz,dt,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none

! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz,np,npmx
       integer(kind=4) ,intent(in) :: i0,j0,k0
       integer(kind=4) ,dimension(0:1) ,intent(in) :: itim
       real(kind=8),dimension(nx-1,ny,nz,2),intent(in) :: u
       real(kind=8),dimension(nx,ny-1,nz,2),intent(in) :: v
       real(kind=8),dimension(nx,ny,nz+1,2),intent(in) :: w
       real(kind=8),dimension(nx,ny,nz,2) ,intent(in) :: dz
       real(kind=8),dimension(nx,ny) ,intent(in) :: pm,pn
       real(kind=8) ,dimension(np) ,intent(inout):: px,py,pz
       real(kind=8) ,intent(in) :: dt,fct,dfct

! import/export for AB scheme
       integer(kind=4),dimension(0:2) ,intent(inout):: iab
       real(kind=8) ,dimension(np,0:2) ,intent(inout):: dpx,dpy,dpz

! local
       real(kind=8) :: coef
       integer(kind=8) :: ip
       real(kind=8) :: dpxi,dpyi,dpzi


!f2py intent(inout) px,py,pz
!f2py intent(in) u,v,w
!f2py intent(in) fct,dfct,dt
!f2py intent(in) pm,pn,dz
!f2py intent(in) npmx
!f2py intent(in) nx,ny,nz
!f2py intent(in) i0,j0,k0
!f2py intent(in) itim
!f2py intent(in) np
!f2py intent(inout) iab
!f2py intent(inout) dpx,dpy,dpz


       coef = 1./12.

       do ip = 1,min(npmx,np)
         if (.not.(isnan(px(ip))) ) then

           ! Forward Euler / predictor
           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,itim,
     & fct,pm,pn,dz,dt,i0,j0,k0,nx,ny,nz,np,
     & dpx(ip,iab(2)),dpy(ip,iab(2)),dpz(ip,iab(2)))

            dpxi = coef * (23.*dpx(ip,iab(2)) - 16.*dpx(ip,iab(1))
     & + 5.*dpx(ip,iab(0)))
            dpyi = coef * (23.*dpy(ip,iab(2)) - 16.*dpy(ip,iab(1))
     & + 5.*dpy(ip,iab(0)))
            dpzi = coef * (23.*dpz(ip,iab(2)) - 16.*dpz(ip,iab(1))
     & + 5.*dpz(ip,iab(0)))


            !---------------------------------------
            ! Update position
            !---------------------------------------
            px(ip) = px(ip) + dpxi
            py(ip) = py(ip) + dpyi
            pz(ip) = pz(ip) + dpzi


         endif
       enddo



       end

        !---------------------------------------------------------------------!
        ! Adams-Bashforth 4 time-stepping
        !---------------------------------------------------------------------
        ! For Adams-methods we need to keep the previous dpx,dpy,dpz values
        !

       subroutine timestep_AB4(px,py,pz,dpx,dpy,dpz,iab,
     & u,v,w,itim,fct,dfct,pm,pn,dz,dt,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none

! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz,np,npmx
       integer(kind=4) ,intent(in) :: i0,j0,k0
       integer(kind=4) ,dimension(0:1) ,intent(in) :: itim
       real(kind=8),dimension(nx-1,ny,nz,2),intent(in) :: u
       real(kind=8),dimension(nx,ny-1,nz,2),intent(in) :: v
       real(kind=8),dimension(nx,ny,nz+1,2),intent(in) :: w
       real(kind=8),dimension(nx,ny,nz,2) ,intent(in) :: dz
       real(kind=8),dimension(nx,ny) ,intent(in) :: pm,pn
       real(kind=8) ,dimension(np) ,intent(inout):: px,py,pz
       real(kind=8) ,intent(in) :: dt,fct,dfct

! import/export for AB scheme
       integer(kind=4),dimension(0:3) ,intent(inout):: iab
       real(kind=8) ,dimension(np,0:3) ,intent(inout):: dpx,dpy,dpz

! local
       real(kind=8) :: coef
       integer(kind=8) :: ip
       real(kind=8) :: dpxi,dpyi,dpzi


!f2py intent(inout) px,py,pz
!f2py intent(in) u,v,w
!f2py intent(in) fct,dfct,dt
!f2py intent(in) pm,pn,dz
!f2py intent(in) npmx
!f2py intent(in) nx,ny,nz
!f2py intent(in) i0,j0,k0
!f2py intent(in) itim
!f2py intent(in) np
!f2py intent(inout) iab
!f2py intent(inout) dpx,dpy,dpz

       coef = 1./24.

       do ip = 1,min(npmx,np)
         if (.not.(isnan(px(ip))) ) then

           ! Forward Euler / predictor
           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,itim,
     & fct,pm,pn,dz,dt,i0,j0,k0,nx,ny,nz,np,
     & dpx(ip,iab(3)),dpy(ip,iab(3)),dpz(ip,iab(3)))
!
! if (.not.(isnan(dpx(ip,iab(0)))) ) then



           ! Adams-Bashforth 4
            dpxi = coef * (55.*dpx(ip,iab(3)) - 59.*dpx(ip,iab(2))
     & + 37.*dpx(ip,iab(1)) - 9.*dpx(ip,iab(0)))
            dpyi = coef * (55.*dpy(ip,iab(3)) - 59.*dpy(ip,iab(2))
     & + 37.*dpy(ip,iab(1)) - 9.*dpy(ip,iab(0)))
            dpzi = coef * (55.*dpz(ip,iab(3)) - 59.*dpz(ip,iab(2))
     & + 37.*dpz(ip,iab(1)) - 9.*dpz(ip,iab(0)))

!
! else if (.not.(isnan(dpx(ip,iab(1)))) ) then
!
! if (ip.eq.1) print *,'time step 2'
! dpxi = 1/12. * (23.*dpx(ip,iab(3)) - 16.*dpx(ip,iab(2))
! & + 5.*dpx(ip,iab(1)) )
! dpyi = 1/12. * (23.*dpy(ip,iab(3)) - 16.*dpy(ip,iab(2))
! & + 5.*dpy(ip,iab(1)) )
! dpzi = 1/12. * (23.*dpz(ip,iab(3)) - 16*dpz(ip,iab(2))
! & + 5.*dpz(ip,iab(1)) )
!
! else if (.not.(isnan(dpx(ip,iab(2)))) ) then
!
! if (ip.eq.1) print *,'time step 1'
! dpxi = 0.5 * (3. * dpx(ip,iab(3)) - dpx(ip,iab(2)))
! dpyi = 0.5 * (3. * dpy(ip,iab(3)) - dpy(ip,iab(2)))
! dpzi = 0.5 * (3. * dpz(ip,iab(3)) - dpz(ip,iab(2)))
!
! else
!
! if (ip.eq.1) print *,'time step 0'
! dpxi = dpx(ip,iab(3))
! dpyi = dpy(ip,iab(3))
! dpzi = dpz(ip,iab(3))
!
! endif


            !---------------------------------------
            ! Update position
            !---------------------------------------

            px(ip) = px(ip) + dpxi
            py(ip) = py(ip) + dpyi
            pz(ip) = pz(ip) + dpzi


         endif
       enddo



       end


        !---------------------------------------------------------------------!
        ! Adams-Bashforth 4 Predictor
        ! + Adams-Moulton 5 Corrector
        !---------------------------------------------------------------------
        ! For Adams-methods we need to keep the previous dpx,dpy,dpz values
        !

       subroutine timestep_ABM4(px,py,pz,dpx,dpy,dpz,iab,
     & u,v,w,itim,fct,dfct,pm,pn,dz,dt,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none

! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz,np,npmx
       integer(kind=4) ,intent(in) :: i0,j0,k0
       integer(kind=4) ,dimension(0:1) ,intent(in) :: itim
       real(kind=8),dimension(nx-1,ny,nz,2),intent(in) :: u
       real(kind=8),dimension(nx,ny-1,nz,2),intent(in) :: v
       real(kind=8),dimension(nx,ny,nz+1,2),intent(in) :: w
       real(kind=8),dimension(nx,ny,nz,2) ,intent(in) :: dz
       real(kind=8),dimension(nx,ny) ,intent(in) :: pm,pn
       real(kind=8) ,dimension(np) ,intent(inout):: px,py,pz
       real(kind=8) ,intent(in) :: dt,fct,dfct

! import/export for AB scheme
       integer(kind=4),dimension(0:3) ,intent(inout):: iab
       real(kind=8) ,dimension(np,0:3) ,intent(inout):: dpx,dpy,dpz

! local
       real(kind=8) :: coef
       integer(kind=8) :: ip
       real(kind=8) :: dpxc,dpyc,dpzc
       real(kind=8) :: dpxi,dpyi,dpzi

!f2py intent(inout) px,py,pz
!f2py intent(in) u,v,w
!f2py intent(in) fct,dfct,dt
!f2py intent(in) pm,pn,dz
!f2py intent(in) npmx
!f2py intent(in) nx,ny,nz
!f2py intent(in) i0,j0,k0
!f2py intent(in) itim
!f2py intent(in) np
!f2py intent(inout) iab
!f2py intent(inout) dpx,dpy,dpz



       coef = 1./24.

       do ip = 1,min(npmx,np)
         if (.not.(isnan(px(ip))) ) then

           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,itim,
     & fct,pm,pn,dz,dt,i0,j0,k0,nx,ny,nz,np,
     & dpx(ip,iab(3)),dpy(ip,iab(3)),dpz(ip,iab(3)))

           !---------------------------------------
           ! Predictor
           !---------------------------------------
!
! if (.not.(isnan(dpx(ip,iab(0)))) ) then

            dpxi = coef * (55.*dpx(ip,iab(3)) - 59.*dpx(ip,iab(2))
     & + 37.*dpx(ip,iab(1)) - 9.*dpx(ip,iab(0)))
            dpyi = coef * (55.*dpy(ip,iab(3)) - 59.*dpy(ip,iab(2))
     & + 37.*dpy(ip,iab(1)) - 9.*dpy(ip,iab(0)))
            dpzi = coef * (55.*dpz(ip,iab(3)) - 59.*dpz(ip,iab(2))
     & + 37.*dpz(ip,iab(1)) - 9.*dpz(ip,iab(0)))
!
! else if (.not.(isnan(dpx(ip,iab(1)))) ) then
!
! if (ip.eq.1) print *,'time step 2'
! dpxi = 1/12. * (23.*dpx(ip,iab(3)) - 16.*dpx(ip,iab(2))
! & + 5.*dpx(ip,iab(1)) )
! dpyi = 1/12. * (23.*dpy(ip,iab(3)) - 16.*dpy(ip,iab(2))
! & + 5.*dpy(ip,iab(1)) )
! dpzi = 1/12. * (23.*dpz(ip,iab(3)) - 16*dpz(ip,iab(2))
! & + 5.*dpz(ip,iab(1)) )
!
! else if (.not.(isnan(dpx(ip,iab(2)))) ) then
!
! if (ip.eq.1) print *,'time step 1'
! dpxi = 0.5 * (3. * dpx(ip,iab(3)) - dpx(ip,iab(2)))
! dpyi = 0.5 * (3. * dpy(ip,iab(3)) - dpy(ip,iab(2)))
! dpzi = 0.5 * (3. * dpz(ip,iab(3)) - dpz(ip,iab(2)))
!
! else
!
! if (ip.eq.1) print *,'time step 0'
! dpxi = dpxi
! dpyi = dpyi
! dpzi = dpzi
!
! endif



           !---------------------------------------
           ! Compute velocities at new time-step
           !---------------------------------------

           CALL advance_3d(px(ip)+dpxi,py(ip)+dpyi,pz(ip)+dpzi,
     & u,v,w,itim,
     & fct+dfct,pm,pn,dz,dt,i0,j0,k0,nx,ny,nz,np,dpxc,dpyc,dpzc)


           !---------------------------------------
           ! Corrector
           !---------------------------------------
!
! if (.not.(isnan(dpx(ip,iab(0)))) ) then

            dpxi = 1./270. * ( 19.* dpxi + 251.* ( 1./24. * (
     & 9. * dpxc + 19. * dpx(ip,iab(3))
     & - 5. * dpx(ip,iab(2)) + dpx(ip,iab(1)) )))

            dpyi = 1./270. * ( 19.* dpyi + 251.* ( 1./24. * (
     & 9. * dpyc + 19. * dpy(ip,iab(3))
     & - 5. * dpy(ip,iab(2)) + dpy(ip,iab(1)) )))

            dpzi = 1./270. * ( 19.* dpzi + 251.* ( 1./24. * (
     & 9. * dpzc + 19. * dpz(ip,iab(3))
     & - 5. * dpz(ip,iab(2)) + dpz(ip,iab(1)) )))
!
!
! else if (.not.(isnan(dpx(ip,iab(1)))) ) then
!
! print *,'time step 2'
! dpxi = 1./270. * ( 19.* dpxi + 251.* ( 1./12. * (
! & 5. * dpxc + 8. * dpx(ip,iab(3))
! & - 1. * dpx(ip,iab(2)) )))
! dpyi = 1./270. * ( 19.* dpyi + 251.* ( 1./12. * (
! & 5. * dpyc + 8. * dpy(ip,iab(3))
! & - 1. * dpy(ip,iab(2)) )))
! dpzi = 1./270. * ( 19.* dpzi + 251.* ( 1./12. * (
! & 5. * dpzc + 8. * dpz(ip,iab(3))
! & - 1. * dpz(ip,iab(2)) )))
!
! else if (.not.(isnan(dpx(ip,iab(2)))) ) then
!
! print *,'time step 1'
! dpxi = 1./270. * ( 19.* dpxi + 251.* ( 1./2. * (
! & 1. * dpxc + 1. * dpx(ip,iab(3)) )))
! dpyi = 1./270. * ( 19.* dpyi + 251.* ( 1./2. * (
! & 1. * dpyc + 1. * dpy(ip,iab(3)) )))
! dpzi = 1./270. * ( 19.* dpzi + 251.* ( 1./2. * (
! & 1. * dpzc + 1. * dpz(ip,iab(3)) )))
!
! else
!
! print *,'time step 0'
! pxi = 1./270. * ( 19.* pxi + 251.* ( px(ip) + dpxc ))
! pyi = 1./270. * ( 19.* pyi + 251.* ( py(ip) + dpyc ))
! pzi = 1./270. * ( 19.* pzi + 251.* ( pz(ip) + dpzc ))
!
! endif


            !---------------------------------------
            ! Update position
            !---------------------------------------

            px(ip) = px(ip) + dpxi
            py(ip) = py(ip) + dpyi
            pz(ip) = pz(ip) + dpzi


         endif
       enddo



       end





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Some soutines copied from the R_tools.F module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute z_r and z_w for NEW_S_COORD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# 1 "R_tools_fort_routines/zlevs.F" 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute z_r and z_w for NEW_S_COORD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      subroutine zlevs(Lm,Mm,N, h,zeta, hc, Cs_r, Cs_w,z_r,z_w)


      implicit none

      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k

      real*8 Cs_w(0:N), Cs_r(N), cff_w, cff_r, cff1_w, cff1_r,
     & hc, ds,
     & zeta(0:Lm+1,0:Mm+1),
     & z_r(0:Lm+1,0:Mm+1,N),z_w(0:Lm+1,0:Mm+1,0:N),
     & h(0:Lm+1,0:Mm+1),hinv(0:Lm+1,0:Mm+1)


Cf2py intent(in) Lm,Mm,N, h,zeta, hc, Cs_w, Cs_r
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

        do k=1,N,+1 !--> irreversible because of recursion in Hz


          cff_w=hc*ds* dble(k-N)
          cff_r=hc*ds*(dble(k-N)-0.5)

          cff1_w=Cs_w(k)
          cff1_r=Cs_r(k)


          do i=imin,imax

            z_w(i,j,k)=zeta(i,j) +(zeta(i,j)+h(i,j))
     & *(cff_w+cff1_w*h(i,j))*hinv(i,j)

            z_r(i,j,k)=zeta(i,j) +(zeta(i,j)+h(i,j))
     & *(cff_r+cff1_r*h(i,j))*hinv(i,j)


          enddo
        enddo
      enddo
      end
# 722 "particules_3d_sig_sa.F" 2

# 1 "R_tools_fort_routines/zlevs_w.F" 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute z_w for NEW_S_COORD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      subroutine zlevs_w(Lm,Mm,N, h,zeta, hc, Cs_w,z_w)


      implicit none

      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k

      real*8 Cs_w(0:N), cff_w, cff1_w,
! & Cs_r(N), cff_r, cff1_r,
     & hc, ds,
     & zeta(0:Lm+1,0:Mm+1),
     & z_w(0:Lm+1,0:Mm+1,0:N),
     & h(0:Lm+1,0:Mm+1),hinv(0:Lm+1,0:Mm+1)


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

        do k=1,N,+1 !--> irreversible because of recursion in Hz


          cff_w=hc*ds* dble(k-N)
! cff_r=hc*ds*(dble(k-N)-0.5)

          cff1_w=Cs_w(k)
! cff1_r=Cs_r(k)


          do i=imin,imax

            z_w(i,j,k)=zeta(i,j) +(zeta(i,j)+h(i,j))
     & *(cff_w+cff1_w*h(i,j))*hinv(i,j)

! z_r(i,j,k)=zeta(i,j) +(zeta(i,j)+h(i,j))
! & *(cff_r+cff1_r*h(i,j))*hinv(i,j)


          enddo
        enddo
      enddo
      end
# 724 "particules_3d_sig_sa.F" 2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Z interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# 1 "R_tools_fort_routines/sigma_to_z_intr_sfc.F" 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Z interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sigma_to_z_intr_sfc (Lm,Mm,N, nz, z_r, z_w, rmask, var,
     & z_lev, var_zlv, imin,jmin,kmin, FillValue)
!
! Interpolate field "var" defined in sigma-space to 3-D z_lev.
!


      implicit none

      integer Lm,Mm,N, nz, imin,imax,jmin,jmax, kmin, i,j,k,m

      integer km(0:Lm+1)

      real*8 var(imin:Lm+1,jmin:Mm+1,kmin:N),
     & z_r(0:Lm+1,0:Mm+1,N), rmask(0:Lm+1,0:Mm+1),
     & z_w(0:Lm+1,0:Mm+1,0:N), z_lev(0:Lm+1,0:Mm+1,nz),
     & FillValue, var_zlv(imin:Lm+1,jmin:Mm+1,nz),
     & zz(0:Lm+1,0:N+1), dpth

     & , dz(0:Lm+1,kmin-1:N), FC(0:Lm+1,kmin-1:N), p,q,cff

      integer numthreads, trd, chunk_size, margin, jstr,jend
C$ integer omp_get_num_threads, omp_get_thread_num


      imax=Lm+1
      jmax=Mm+1

      numthreads=1
C$ numthreads=omp_get_num_threads()
      trd=0
C$ trd=omp_get_thread_num()
      chunk_size=(jmax-jmin + numthreads)/numthreads
      margin=(chunk_size*numthreads -jmax+jmin-1)/2
      jstr=jmin !max( trd *chunk_size -margin, jmin )
      jend=jmax !min( (trd+1)*chunk_size-1-margin, jmax )


Cf2py intent(in) Lm,Mm,N, nz, z_r, z_w, rmask, var, z_lev, imin,jmin,kmin, FillValue
Cf2py intent(out) var_zlv
# 53 "R_tools_fort_routines/sigma_to_z_intr_sfc.F"
      do j=jstr,jend
        if (kmin.eq.1) then
          if (imin.eq.0 .and. jmin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=z_r(i,j,k)
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=z_w(i,j,0)
              zz(i,N+1)=z_w(i,j,N)
            enddo
          elseif (imin.eq.1 .and. jmin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.5D0*(z_r(i,j,k)+z_r(i-1,j,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.5D0*(z_w(i-1,j,0)+z_w(i,j,0))
              zz(i,N+1)=0.5D0*(z_w(i-1,j,N)+z_w(i,j,N))
            enddo
          elseif (imin.eq.0 .and. jmin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.5*(z_r(i,j,k)+z_r(i,j-1,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.5D0*(z_w(i,j,0)+z_w(i,j-1,0))
              zz(i,N+1)=0.5D0*(z_w(i,j,N)+z_w(i,j-1,N))
            enddo
          elseif (imin.eq.1 .and. jmin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.25D0*( z_r(i,j,k)+z_r(i-1,j,k)
     & +z_r(i,j-1,k)+z_r(i-1,j-1,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.25D0*( z_w(i,j,0)+z_w(i-1,j,0)
     & +z_w(i,j-1,0)+z_w(i-1,j-1,0))

              zz(i,N+1)=0.25D0*( z_w(i,j,N)+z_w(i-1,j,N)
     & +z_w(i,j-1,N)+z_w(i-1,j-1,N))
             enddo
          endif
        else
          if (imin.eq.0 .and. jmin.eq.0) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=z_w(i,j,k)
              enddo
            enddo
          elseif (imin.eq.1 .and. jmin.eq.0) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.5D0*(z_w(i,j,k)+z_w(i-1,j,k))
              enddo
            enddo
          elseif (imin.eq.0 .and. jmin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.5*(z_w(i,j,k)+z_w(i,j-1,k))
              enddo
            enddo
          elseif (imin.eq.1 .and. jmin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.25D0*( z_w(i,j,k)+z_w(i-1,j,k)
     & +z_w(i,j-1,k)+z_w(i-1,j-1,k))
              enddo
            enddo
          endif
        endif

        do k=kmin,N-1
          do i=imin,imax
            dz(i,k)=zz(i,k+1)-zz(i,k)
            FC(i,k)=var(i,j,k+1)-var(i,j,k)
          enddo
        enddo
        do i=imin,imax
          dz(i,kmin-1)=dz(i,kmin)
          FC(i,kmin-1)=FC(i,kmin)

          dz(i,N)=dz(i,N-1)
          FC(i,N)=FC(i,N-1)
        enddo
        do k=N,kmin,-1 !--> irreversible
          do i=imin,imax
            cff=FC(i,k)*FC(i,k-1)
            if (cff.gt.0.D0) then
              FC(i,k)=cff*(dz(i,k)+dz(i,k-1))/( (FC(i,k)+FC(i,k-1))
     & *dz(i,k)*dz(i,k-1) )
            else
              FC(i,k)=0.D0
            endif
          enddo
        enddo

        do m=1,nz


          if (kmin.eq.0) then !
            do i=imin,imax !
              dpth=zz(i,N)-zz(i,0)
              if (rmask(i,j).lt.0.5) then
                km(i)=-3 !--> masked out
              elseif (dpth*(z_lev(i,j,m)-zz(i,N)).gt.0.) then
                km(i)=N+2 !<-- above surface
              elseif (dpth*(zz(i,0)-z_lev(i,j,m)).gt.0.) then
                km(i)=-2 !<-- below bottom
              else
                km(i)=-1 !--> to search
              endif
            enddo
          else
            do i=imin,imax
              dpth=zz(i,N+1)-zz(i,0)
              if (rmask(i,j).lt.0.5) then
                km(i)=-3 !--> masked out
              elseif (dpth*(z_lev(i,j,m)-zz(i,N+1)).gt.0.) then
                km(i)=N+2 !<-- above surface

              elseif (dpth*(z_lev(i,j,m)-zz(i,N)).gt.0.) then
                km(i)=N !<-- below surface, but above z_r(N)
              elseif (dpth*(zz(i,0)-z_lev(i,j,m)).gt.0.) then
                km(i)=-2 !<-- below bottom
              elseif (dpth*(zz(i,1)-z_lev(i,j,m)).gt.0.) then
                km(i)=0 !<-- above bottom, but below z_r(1)
              else
                km(i)=-1 !--> to search
              endif
            enddo
          endif
          do k=N-1,kmin,-1
            do i=imin,imax
              if (km(i).eq.-1) then
                if((zz(i,k+1)-z_lev(i,j,m))*(z_lev(i,j,m)-zz(i,k))
     & .ge. 0.) km(i)=k
              endif
            enddo
          enddo

          do i=imin,imax
            if (km(i).eq.-3) then
              var_zlv(i,j,m)=0. !<-- masked out
            elseif (km(i).eq.-2) then

              var_zlv(i,j,m)=FillValue !<-- below bottom

            elseif (km(i).eq.N+2) then

              var_zlv(i,j,m)=var(i,j,N) !-> R-point, above z_r(N)

     & +FC(i,N)*(z_lev(i,j,m)-zz(i,N))







            elseif (km(i).eq.N) then
              var_zlv(i,j,m)=var(i,j,N) !-> R-point, above z_r(N)

     & +FC(i,N)*(z_lev(i,j,m)-zz(i,N))




            elseif (km(i).eq.kmin-1) then !-> R-point below z_r(1),
              var_zlv(i,j,m)=var(i,j,kmin) ! but above bottom

     & -FC(i,kmin)*(zz(i,kmin)-z_lev(i,j,m))




            else
              k=km(i)
              !write(*,*) k,km

              cff=1.D0/(zz(i,k+1)-zz(i,k))
              p=z_lev(i,j,m)-zz(i,k)
              q=zz(i,k+1)-z_lev(i,j,m)

              var_zlv(i,j,m)=cff*( q*var(i,j,k) + p*var(i,j,k+1)
     & -cff*p*q*( cff*(q-p)*(var(i,j,k+1)-var(i,j,k))
     & +p*FC(i,k+1) -q*FC(i,k) )
     & )







            !write(*,*) 'bof',i,j,k,zz(i,k), zz(i,k+1), z_lev(i,j,m), m
# 266 "R_tools_fort_routines/sigma_to_z_intr_sfc.F"
            endif
          enddo
        enddo ! <-- m
      enddo !<-- j

      return
      end
# 730 "particules_3d_sig_sa.F" 2

# 1 "R_tools_fort_routines/sigma_to_z_intr_bot.F" 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Z interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      subroutine sigma_to_z_intr_bot (Lm,Mm,N, nz, z_r, z_w, rmask, var,
     & z_lev, var_zlv, below, imin,jmin,kmin, FillValue)
!
! Interpolate field "var" defined in sigma-space to 3-D z_lev.
!


      implicit none

      integer Lm,Mm,N, nz, imin,imax,jmin,jmax, kmin, i,j,k,m

      integer km(0:Lm+1)

      real*8 var(imin:Lm+1,jmin:Mm+1,kmin:N),
     & z_r(0:Lm+1,0:Mm+1,N), rmask(0:Lm+1,0:Mm+1),
     & z_w(0:Lm+1,0:Mm+1,0:N), z_lev(0:Lm+1,0:Mm+1,nz),
     & FillValue, var_zlv(imin:Lm+1,jmin:Mm+1,nz),
     & zz(0:Lm+1,0:N+1), dpth, below

     & , dz(0:Lm+1,kmin-1:N), FC(0:Lm+1,kmin-1:N), p,q,cff

      integer numthreads, trd, chunk_size, margin, jstr,jend
C$ integer omp_get_num_threads, omp_get_thread_num


      imax=Lm+1
      jmax=Mm+1

      numthreads=1
C$ numthreads=omp_get_num_threads()
      trd=0
C$ trd=omp_get_thread_num()
      chunk_size=(jmax-jmin + numthreads)/numthreads
      margin=(chunk_size*numthreads -jmax+jmin-1)/2
      jstr=jmin !max( trd *chunk_size -margin, jmin )
      jend=jmax !min( (trd+1)*chunk_size-1-margin, jmax )


Cf2py intent(in) Lm,Mm,N, nz, z_r, z_w, rmask, var, z_lev, below, imin,jmin,kmin, FillValue
Cf2py intent(out) var_zlv
# 53 "R_tools_fort_routines/sigma_to_z_intr_bot.F"
      do j=jstr,jend
        if (kmin.eq.1) then
          if (imin.eq.0 .and. jmin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=z_r(i,j,k)
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=z_w(i,j,0)
              zz(i,N+1)=z_w(i,j,N)
            enddo
          elseif (imin.eq.1 .and. jmin.eq.0) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.5D0*(z_r(i,j,k)+z_r(i-1,j,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.5D0*(z_w(i-1,j,0)+z_w(i,j,0))
              zz(i,N+1)=0.5D0*(z_w(i-1,j,N)+z_w(i,j,N))
            enddo
          elseif (imin.eq.0 .and. jmin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.5*(z_r(i,j,k)+z_r(i,j-1,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.5D0*(z_w(i,j,0)+z_w(i,j-1,0))
              zz(i,N+1)=0.5D0*(z_w(i,j,N)+z_w(i,j-1,N))
            enddo
          elseif (imin.eq.1 .and. jmin.eq.1) then
            do k=1,N
              do i=imin,imax
                zz(i,k)=0.25D0*( z_r(i,j,k)+z_r(i-1,j,k)
     & +z_r(i,j-1,k)+z_r(i-1,j-1,k))
              enddo
            enddo
            do i=imin,imax
              zz(i,0)=0.25D0*( z_w(i,j,0)+z_w(i-1,j,0)
     & +z_w(i,j-1,0)+z_w(i-1,j-1,0))

              zz(i,N+1)=0.25D0*( z_w(i,j,N)+z_w(i-1,j,N)
     & +z_w(i,j-1,N)+z_w(i-1,j-1,N))
             enddo
          endif
        else
          if (imin.eq.0 .and. jmin.eq.0) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=z_w(i,j,k)
              enddo
            enddo
          elseif (imin.eq.1 .and. jmin.eq.0) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.5D0*(z_w(i,j,k)+z_w(i-1,j,k))
              enddo
            enddo
          elseif (imin.eq.0 .and. jmin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.5*(z_w(i,j,k)+z_w(i,j-1,k))
              enddo
            enddo
          elseif (imin.eq.1 .and. jmin.eq.1) then
            do k=0,N
              do i=imin,imax
                zz(i,k)=0.25D0*( z_w(i,j,k)+z_w(i-1,j,k)
     & +z_w(i,j-1,k)+z_w(i-1,j-1,k))
              enddo
            enddo
          endif
        endif

        do k=kmin,N-1
          do i=imin,imax
            dz(i,k)=zz(i,k+1)-zz(i,k)
            FC(i,k)=var(i,j,k+1)-var(i,j,k)
          enddo
        enddo
        do i=imin,imax
          dz(i,kmin-1)=dz(i,kmin)
          FC(i,kmin-1)=FC(i,kmin)

          dz(i,N)=dz(i,N-1)
          FC(i,N)=FC(i,N-1)
        enddo
        do k=N,kmin,-1 !--> irreversible
          do i=imin,imax
            cff=FC(i,k)*FC(i,k-1)
            if (cff.gt.0.D0) then
              FC(i,k)=cff*(dz(i,k)+dz(i,k-1))/( (FC(i,k)+FC(i,k-1))
     & *dz(i,k)*dz(i,k-1) )
            else
              FC(i,k)=0.D0
            endif
          enddo
        enddo

        do m=1,nz


          if (kmin.eq.0) then !
            do i=imin,imax !
              dpth=zz(i,N)-zz(i,0)
              if (rmask(i,j).lt.0.5) then
                km(i)=-3 !--> masked out
              elseif (dpth*(z_lev(i,j,m)-zz(i,N)).gt.0.) then
                km(i)=N+2 !<-- above surface
              elseif (dpth*(zz(i,0)-z_lev(i,j,m)).gt.0.) then
                km(i)=-2 !<-- below bottom
              else
                km(i)=-1 !--> to search
              endif
            enddo
          else
            do i=imin,imax
              dpth=zz(i,N+1)-zz(i,0)
              if (rmask(i,j).lt.0.5) then
                km(i)=-3 !--> masked out
              elseif (dpth*(z_lev(i,j,m)-zz(i,N+1)).gt.0.) then
                km(i)=N+2 !<-- above surface

              elseif (dpth*(z_lev(i,j,m)-zz(i,N)).gt.0.) then
                km(i)=N !<-- below surface, but above z_r(N)
              elseif (dpth*(zz(i,0)-below-z_lev(i,j,m)).gt.0.) then
                km(i)=-3 !<-- below bottom
              elseif (dpth*(zz(i,0)-z_lev(i,j,m)).gt.0.) then
                km(i)=-2 !<-- below bottom but close
              elseif (dpth*(zz(i,1)-z_lev(i,j,m)).gt.0.) then
                km(i)=0 !<-- above bottom, but below z_r(1)
              else
                km(i)=-1 !--> to search
              endif
            enddo
          endif
          do k=N-1,kmin,-1
            do i=imin,imax
              if (km(i).eq.-1) then
                if((zz(i,k+1)-z_lev(i,j,m))*(z_lev(i,j,m)-zz(i,k))
     & .ge. 0.) km(i)=k
              endif
            enddo
          enddo

          do i=imin,imax
            if (km(i).eq.-3) then
              var_zlv(i,j,m)=FillValue !<-- masked out
            elseif (km(i).eq.-2) then

              var_zlv(i,j,m)=var(i,j,kmin) !

     & -FC(i,kmin)*(zz(i,kmin)-z_lev(i,j,m))







            elseif (km(i).eq.N+2) then

              var_zlv(i,j,m)=var(i,j,N) !-> R-point, above z_r(N)

     & +FC(i,N)*(z_lev(i,j,m)-zz(i,N))







            elseif (km(i).eq.N) then
              var_zlv(i,j,m)=var(i,j,N) !-> R-point, above z_r(N)

     & +FC(i,N)*(z_lev(i,j,m)-zz(i,N))




            elseif (km(i).eq.kmin-1) then !-> R-point below z_r(1),
              var_zlv(i,j,m)=var(i,j,kmin) ! but above bottom

     & -FC(i,kmin)*(zz(i,kmin)-z_lev(i,j,m))




            else
              k=km(i)
              !write(*,*) k,km

              cff=1.D0/(zz(i,k+1)-zz(i,k))
              p=z_lev(i,j,m)-zz(i,k)
              q=zz(i,k+1)-z_lev(i,j,m)

              var_zlv(i,j,m)=cff*( q*var(i,j,k) + p*var(i,j,k+1)
     & -cff*p*q*( cff*(q-p)*(var(i,j,k+1)-var(i,j,k))
     & +p*FC(i,k+1) -q*FC(i,k) )
     & )







            !write(*,*) 'bof',i,j,k,zz(i,k), zz(i,k+1), z_lev(i,j,m), m
# 276 "R_tools_fort_routines/sigma_to_z_intr_bot.F"
            endif
          enddo
        enddo ! <-- m
      enddo !<-- j

      return
      end
# 732 "particules_3d_sig_sa.F" 2
# 740 "particules_3d_sig_sa.F"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Test if a (1-D) array if fortran formatted or not
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine test_arg_1d (Lm, T)

      implicit none
      integer Lm
      real*4 T(0:Lm+1)
Cf2py intent(in) Lm
Cf2py intent(inout) T

      write(*,*) T(2)

      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Test if a (1-D) array if fortran formatted or not
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine test_arg_1d_double (Lm, T)

      implicit none
      integer Lm
      real*8 T(0:Lm+1)
Cf2py intent(in) Lm
Cf2py intent(inout) T

      write(*,*) T(2)

      end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Test if a (2-D) array if fortran formatted or not
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine test_arg_2d (Lm,Mm, T)

      implicit none
      integer Lm,Mm
      real*4 T(0:Lm+1,0:Mm+1)
Cf2py intent(in) Lm,Mm
Cf2py intent(inout) T

      write(*,*) T(2,2)

      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Test if a (2-D) array if fortran formatted or not
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine test_arg_2d_double (Lm,Mm, T)

      implicit none
      integer Lm,Mm
      real*8 T(0:Lm+1,0:Mm+1)
Cf2py intent(in) Lm,Mm
Cf2py intent(inout) T

      write(*,*) T(2,2)

      end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Test if a 3-D array if fortran formatted or not
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine test_arg_3d (Lm,Mm,N, T)


      implicit none
      integer Lm,Mm,N
      real*4 T(0:Lm+1,0:Mm+1,N)
Cf2py intent(in) Lm,Mm,N
Cf2py intent(inout) T

      write(*,*) T(2,2,2)

      end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Test if a 3-D array if fortran formatted or not
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      subroutine test_arg_3d_double (Lm,Mm,N, T)


      implicit none
      integer Lm,Mm,N
      real*8 T(0:Lm+1,0:Mm+1,N)
Cf2py intent(in) Lm,Mm,N
Cf2py intent(inout) T

      write(*,*) T(2,2,2)

      end
