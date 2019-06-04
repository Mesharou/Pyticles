# 1 "pyticles_3d_sig_sa.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 31 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 32 "<command-line>" 2
# 1 "pyticles_3d_sig_sa.F"
!----------------------------------------------------------------------------------------------
! Fortran Routines for Pyticles
!
!----------------------------------------------------------------------------------------------
! 05/10/16:
! - add mask as a shared array and use it in the time_step routines
! - add the # define CHECK_MASK in pyticles_3d_sig_sa.F
! - add subroutines interp_3d_u / interp_3d_u
! 16/01/26:
! - Add parameter 'ng' corresponding to number of ghost points in the horizontal grid
!----------------------------------------------------------------------------------------------




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Spatial interpolation routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# 1 "interp_3d_for_pyticles.F" 1

!----------------------------------------------------------------------------------------------
! Spatial Interpolation Routines for Pyticles
!
!----------------------------------------------------------------------------------------------
! 05/10/16:
! - correction of advance_3d for ORIGINAL_VERSION (computation of fcx_u,fcx_v,fcx_w)
! - add subroutines interp_3d_u / interp_3d_u
! 16/05/03:
! - Correction of linear_2d [terms (2,1) and (1,2) inversion]
! 16/01/26:
! - Add parameter 'ng' corresponding to number of ghost points in the horizontal grid
! - interp_3d routines only for STAGGERED yet
!----------------------------------------------------------------------------------------------



!#undef NEW_VERSION
!#defin LINEAR_INTERPOLATION
!#define CUBIC_INTERPOLATION
!----------------------------------------------------------------------------------------------



! #define CRSPL_INTERPOLATION
! #define WENO_INTERPOLATION
!----------------------------------------------------------------------------------------------

!---------------------------------------------------------------------!
! if is defined, px,py,pz correpond to staggered positions (u,v,w grids)
! if is NOT defined, px,py,pz correspond to positions on the horizontal rho-grid and vertical w-grid
!---------------------------------------------------------------------!




!----------------------------------------------------------------------------------------------


       subroutine advance_3d(px,py,pz,u,v,w,itim,fct,pm,pn,
     & dz,dt,i0,j0,k0,nx,ny,nz,ng,np,dpx,dpy,dpz)

       !---------------------------------------------------------------------!
       ! Compute particle displacement with linear interpolation in time
       ! in: u,v,w; at their original staggered position
       ! px,py,pz; Particle position in index coordinates [0,nx] etc....
       ! pm,pn,dz; grid stretching at rho position (pm,pn 2d)
       !
       !---------------------------------------------------------------------

       implicit none
! import/export
       real(kind=8) ,intent(in) :: px,py,pz
       real(kind=8) ,dimension(0:nx-2,0:ny-1,nz,0:1) ,intent(in) :: u
       real(kind=8) ,dimension(0:nx-1,0:ny-2,nz,0:1) ,intent(in) :: v
       real(kind=8) ,dimension(0:nx-1,0:ny-1,0:nz,0:1),intent(in) :: w
       real(kind=8) ,dimension(0:nx-1,0:ny-1,nz,0:1) ,intent(in) :: dz
       real(kind=8) ,dimension(0:nx-1,0:ny-1) ,intent(in) :: pm,pn
       integer(kind=4) ,intent(in) :: nx,ny,nz,np
       integer(kind=4) ,intent(in) :: ng
       integer(kind=4) ,intent(in) :: i0,j0,k0
       real(kind=8) ,intent(out):: dpx,dpy,dpz
       real(kind=8) ,intent(in) :: dt,fct
       integer(kind=4) ,dimension(0:1) ,intent(in) :: itim
! local
       integer(kind=8) :: i,j,k
       integer(kind=8) :: iu,jv,kw
       real(kind=8) ,dimension(4,4,4) :: f
       real(kind=8) :: pxl,pyl,pzl
       real(kind=8) :: pxlu,pylv,pzlw
       real(kind=8) :: pu,pv,pw
       real(kind=8) :: ppm,ppn,pdz

!f2py intent(in) u,v,w
!f2py intent(in) px,py,pz
!f2py intent(out) dpx,dpy,dpz
!f2py intent(in) pm,pn,dt,fct
!f2py intent(in) dz
!f2py intent(in) nx,ny,nz,ng,np
!f2py intent(in) i0,j0,k0
!f2py intent(in) itim



           !! add ghost points to deal with boundaries
           !! for now, we ll add ghost points here at z = 0,nz if
           !! neccessary and rely on outside to keep a 2-point buffer in
           !! the horizontal




           iu= floor(px+ng-i0); pxlu= px+ng-i0-iu !! In x-dir u grid starts at 0
           jv= floor(py+ng-j0); pylv= py+ng-j0-jv !! In y-dir v grid starts at 0
           kw= floor(pz-k0); pzlw= pz-k0-kw !!
           i = floor(px+ng-i0+0.5); pxl = px+ng-i0-i+0.5 !! In x-dir rho grid starts at -0.5
           j = floor(py+ng-j0+0.5); pyl = py+ng-j0-j+0.5 !! In y-dir rho grid starts at -0.5
           k = floor(pz-k0+0.5); pzl = pz-k0-k+0.5 !!
# 109 "interp_3d_for_pyticles.F"
           !! for now, we do slip conditions near the bottom and
           !! extrapolation near the top.
           if (k .eq. 0) then
             f(:,:,3:4) = fct * u(iu-1:iu+2,j-1:j+2,k+1:k+2,itim(1))
     & + (1-fct) * u(iu-1:iu+2,j-1:j+2,k+1:k+2,itim(0))
             f(:,:, 2) = f(:,:,3)
             f(:,:, 1) = f(:,:,2)
           elseif (k.eq.1) then
             f(:,:,2:4) = fct * u(iu-1:iu+2,j-1:j+2, k:k+2,itim(1))
     & + (1-fct) * u(iu-1:iu+2,j-1:j+2, k:k+2,itim(0))
             f(:,:, 1) = f(:,:,2)
           elseif (k.eq.nz-1) then
             f(:,:,1:3) = fct * u(iu-1:iu+2,j-1:j+2,k-1:k+1,itim(1))
     & + (1-fct) * u(iu-1:iu+2,j-1:j+2,k-1:k+1,itim(0))
             f(:,:, 4) = 2*f(:,:,3)-f(:,:,2)
           elseif (k.eq.nz) then
             f(:,:,1:2) = fct * u(iu-1:iu+2,j-1:j+2,k-1:k,itim(1))
     & + (1-fct) * u(iu-1:iu+2,j-1:j+2,k-1:k,itim(0))
             f(:,:, 3) = 2*f(:,:,2)-f(:,:,1)
             f(:,:, 4) = 2*f(:,:,3)-f(:,:,2)
           else
             f = fct * u(iu-1:iu+2,j-1:j+2,k-1:k+2,itim(1))
     & + (1-fct) * u(iu-1:iu+2,j-1:j+2,k-1:k+2,itim(0))
           endif
           call interp3(f,pxlu,pyl,pzl,pu)

           if (k .eq. 0) then
             f(:,:,3:4) = fct * v(i-1:i+2,jv-1:jv+2,k+1:k+2,itim(1))
     & + (1-fct) * v(i-1:i+2,jv-1:jv+2,k+1:k+2,itim(0))
             f(:,:, 2) = f(:,:,3)
             f(:,:, 1) = f(:,:,2)
           elseif (k.eq.1) then
             f(:,:,2:4) = fct * v(i-1:i+2,jv-1:jv+2, k:k+2,itim(1))
     & + (1-fct) * v(i-1:i+2,jv-1:jv+2, k:k+2,itim(0))
             f(:,:, 1) = f(:,:,2)
           elseif (k.eq.nz-1) then
             f(:,:,1:3) = fct * v(i-1:i+2,jv-1:jv+2,k-1:k+1,itim(1))
     & + (1-fct) * v(i-1:i+2,jv-1:jv+2,k-1:k+1,itim(0))
             f(:,:, 4) = 2*f(:,:,3)-f(:,:,2)
           elseif (k.eq.nz) then
             f(:,:,1:2) = fct * v(i-1:i+2,jv-1:jv+2,k-1:k,itim(1))
     & + (1-fct) * v(i-1:i+2,jv-1:jv+2,k-1:k,itim(0))
             f(:,:, 3) = 2*f(:,:,2)-f(:,:,1)
             f(:,:, 4) = 2*f(:,:,3)-f(:,:,2)
           else
             f = fct * v(i-1:i+2,jv-1:jv+2,k-1:k+2,itim(1))
     & + (1-fct) * v(i-1:i+2,jv-1:jv+2,k-1:k+2,itim(0))
           endif
           call interp3(f,pxl,pylv,pzl,pv)

           if (kw .eq. 0) then
             f(:,:,2:4) = fct * w(i-1:i+2,j-1:j+2,kw:kw+2,itim(1))
     & + (1-fct) * w(i-1:i+2,j-1:j+2,kw:kw+2,itim(0))
             f(:,:, 1) = 0.
           elseif (kw.eq.nz-1) then
             f(:,:,1:3) = fct * w(i-1:i+2,j-1:j+2,kw-1:kw+1,itim(1))
     & + (1-fct) * w(i-1:i+2,j-1:j+2,kw-1:kw+1,itim(0))
             f(:,:, 4) = 0.
           else
             f = fct * w(i-1:i+2,j-1:j+2,kw-1:kw+2,itim(1))
     & + (1-fct) * w(i-1:i+2,j-1:j+2,kw-1:kw+2,itim(0))
           endif

           call interp3(f,pxl,pyl,pzlw,pw)

           if (k .eq. 0) then
             f(:,:,3:4) = fct * dz(i-1:i+2,j-1:j+2,k+1:k+2,itim(1))
     & + (1-fct) * dz(i-1:i+2,j-1:j+2,k+1:k+2,itim(0))
             f(:,:, 2) = f(:,:,3)
             f(:,:, 1) = f(:,:,2)
           elseif (k.eq.1) then
             f(:,:,2:4) = fct * dz(i-1:i+2,j-1:j+2, k:k+2,itim(1))
     & + (1-fct) * dz(i-1:i+2,j-1:j+2, k:k+2,itim(0))
             f(:,:, 1) = f(:,:,2)
           elseif (k.eq.nz-1) then
             f(:,:,1:3) = fct * dz(i-1:i+2,j-1:j+2,k-1:k+1,itim(1))
     & + (1-fct) * dz(i-1:i+2,j-1:j+2,k-1:k+1,itim(0))
             f(:,:, 4) = 2*f(:,:,3)-f(:,:,2)
           elseif (k.eq.nz) then
             f(:,:,1:2) = fct * dz(i-1:i+2,j-1:j+2,k-1:k,itim(1))
     & + (1-fct) * dz(i-1:i+2,j-1:j+2,k-1:k,itim(0))
             f(:,:, 3) = 2*f(:,:,2)-f(:,:,1)
             f(:,:, 4) = 2*f(:,:,3)-f(:,:,2)
           else
             f = fct * dz(i-1:i+2,j-1:j+2,k-1:k+2,itim(1))
     & + (1-fct) * dz(i-1:i+2,j-1:j+2,k-1:k+2,itim(0))
           endif
           call interp3(f,pxl,pyl,pzl,pdz)

           if (pdz.gt.0) pdz=1./pdz

         ppm = pm(i ,j)*(1-pxl)*(1-pyl) + pm(i ,j+1)*(1-pxl)*pyl +
     & pm(i+1,j)* pxl*(1-pyl) + pm(i+1,j+1)*pxl*pyl
         ppn = pn(i ,j)*(1-pxl)*(1-pyl) + pn(i ,j+1)*(1-pxl)*pyl +
     & pn(i+1,j)* pxl*(1-pyl) + pn(i+1,j+1)*pxl*pyl

         dpx = dt*pu*ppm
         dpy = dt*pv*ppn
         dpz = dt*pw*pdz



       end


!----------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------
# 246 "interp_3d_for_pyticles.F"
!----------------------------------------------------------------------------------------------

       subroutine interp3(f,x,y,z,fi)
       !---------------------------------------------------------------------!
       ! Compute cubic 3d interpolant
       ! in: f(4,4,4)
       ! x,y,z; local coordinates (should be between 0 and 1)
       ! out: fi
       !
       !---------------------------------------------------------------------
       implicit none
! import/export
       real(kind=8),dimension(4,4,4),intent(in) :: f
       real(kind=8), intent(in) :: x,y,z
       real(kind=8), intent(out):: fi
! local
       real(kind=8) :: x2,x3,y2,y3,z2,z3
       real(kind=8),dimension(4,4) :: fx,a1,a2,a3,a4
       real(kind=8),dimension(4) :: fy,b1,b2,b3,b4
       real(kind=8) :: c1,c2,c3,c4


       x2 = x*x; x3 = x*x2
       y2 = y*y; y3 = y*y2
       z2 = z*z; z3 = z*z2

       a4 = -1./6*f(1,:,:) + 0.5*f(2,:,:) - 0.5*f(3,:,:) + 1./6*f(4,:,:);
       a3 = 0.5* f(1,:,:) - f(2,:,:) + 0.5*f(3,:,:);
       a2 = -1./3*f(1,:,:) - 0.5*f(2,:,:) + f(3,:,:) - 1./6*f(4,:,:);
       a1 = f(2,:,:);

       fx= a4*x3+ a3*x2 + a2*x + a1;

       b4 = -1./6*fx(1,:) + 0.5*fx(2,:) - 0.5*fx(3,:) + 1./6*fx(4,:);
       b3 = 0.5* fx(1,:) - fx(2,:) + 0.5*fx(3,:);
       b2 = -1./3*fx(1,:) - 0.5*fx(2,:) + fx(3,:) - 1./6*fx(4,:);
       b1 = fx(2,:);

       fy= b4*y3+ b3*y2 + b2*y + b1;

       c4 = -1./6*fy(1) + 0.5*fy(2) - 0.5*fy(3) + 1./6*fy(4);
       c3 = 0.5* fy(1) - fy(2) + 0.5*fy(3);
       c2 = -1./3*fy(1) - 0.5*fy(2) + fy(3) - 1./6*fy(4);
       c1 = fy(2);

       fi= c4*z3+ c3*z2 + c2*z + c1;


       end
!----------------------------------------------------------------------------------------------
# 576 "interp_3d_for_pyticles.F"
!---------------------------------------------------------------------!
! Compute 2d displacement given u,v at particule position
! with linear interpolation in space and time
!---------------------------------------------------------------------



       subroutine advance_2d(px,py,u,v,itim,fct,pm,pn,
     & dt,i0,j0,nx,ny,ng,np,dpx,dpy)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny
       integer(kind=4) ,intent(in) :: ng
       integer(kind=4) ,intent(in) :: np
       integer(kind=4) ,intent(in) :: i0,j0
       integer(kind=4) ,dimension(0:1) ,intent(in) :: itim
       real(kind=8) ,dimension(nx-1,ny,2),intent(in) :: u
       real(kind=8) ,dimension(nx,ny-1,2),intent(in) :: v
       real(kind=8) ,dimension(nx,ny),intent(in) :: pm,pn
       real(kind=8) ,intent(in) :: px,py
       real(kind=8) ,intent(out) :: dpx,dpy
       real(kind=8) ,intent(in) :: dt,fct
! local
       integer(kind=8) :: i,j
       integer(kind=8) :: i_u,j_v
       real(kind=8) ,dimension(2,2,2) :: wt3
       real(kind=8) ,dimension(2,2) :: wt2
       real(kind=8) :: fcx,fcy,fctl
       real(kind=8) :: fcx_u,fcy_v
       real(kind=8) :: pu,pv,ppm,ppn

!f2py intent(in) u,v
!f2py intent(in) px,py
!f2py intent(out) dpx,dpy
!f2py intent(in) pm,pn,dt,fct
!f2py intent(in) nx,ny
!f2py intent(in) i0,j0
!f2py intent(in) itim
!f2py intent(in) ng,np


           !---------------------------------------
           ! 1. Linear interpolation in space and time
           !---------------------------------------


           i = floor(px+1+0.5+ng)-i0
           j = floor(py+1+0.5+ng)-j0

           i_u = floor(px+1+ng)-i0
           j_v = floor(py+1+ng)-j0

           fcx = px+1+0.5+ng - i - i0;
           fcy = py+1+0.5+ng - j - j0;

           fcx_u = px+1+ng - i_u - i0;
           fcy_v = py+1+ng - j_v - j0;
# 649 "interp_3d_for_pyticles.F"
           fctl = fct
           if (itim(0).eq.1) fctl = 1-fct

           !---------------------------------------
           ! Compute velocities and level depth at particle position
           !---------------------------------------


           CALL linear_3d(fcx_u,fcy,fctl,wt3)
           pu = sum(u(i_u:i_u+1,j:j+1,:)*wt3)

           CALL linear_3d(fcx,fcy_v,fctl,wt3)
           pv = sum(v(i:i+1,j_v:j_v+1,:)*wt3)

           CALL linear_2d(fcx,fcy,wt2)
           ppm = sum(pm(i:i+1,j:j+1)*wt2)
           ppn = sum(pn(i:i+1,j:j+1)*wt2)

           !---------------------------------------
           ! Update position
           !---------------------------------------

           dpx = dt*pu*ppm
           dpy = dt*pv*ppn

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


           !wt(1,1) = (1-fcy)*(1-fcx);
           !wt(1,2) = (1-fcy)* fcx;
           !wt(2,1) = fcy *(1-fcx);
           !wt(2,2) = fcy * fcx;

           wt(1,1) = (1-fcy)*(1-fcx);
           wt(1,2) = fcy *(1-fcx);
           wt(2,1) = (1-fcy)* fcx;
           wt(2,2) = fcy * fcx;

       end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate T,S at each particle position (same than interp_3d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_3d_ts(pvar1,pvar2,px,py,pz,
     & var1,var2,ng,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz
       integer(kind=4) ,intent(in) :: np,ng
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
!f2py intent(in) ng,np

       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           i = max(1,min(floor(px(ip)+1+0.5+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+0.5+ng)-j0,ny-1))
           k = max(1,min(floor(pz(ip)+1-0.5)-k0,nz-1))

           fcx = px(ip)+1+0.5+ng - i - i0;
           fcy = py(ip)+1+0.5+ng - j - j0;
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


       subroutine interp_3d(pvar1,px,py,pz,var1,ng,npmx,i0,j0,k0,
     & nx,ny,nz,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz
       integer(kind=4) ,intent(in) :: np,ng
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
!f2py intent(in) ng,np

       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then


           i = max(1,min(floor(px(ip)+1+0.5+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+0.5+ng)-j0,ny-1))
           k = max(1,min(floor(pz(ip)+1-0.5)-k0,nz-1))

           fcx = px(ip)+1+0.5+ng - i - i0;
           fcy = py(ip)+1+0.5+ng - j - j0;
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


       subroutine interp_3d_u(pvar1,px,py,pz,var1,ng,npmx,i0,j0,k0,
     & nx,ny,nz,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz
       integer(kind=4) ,intent(in) :: np,ng
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
!f2py intent(in) ng,np

       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           i = max(1,min(floor(px(ip)+1+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+0.5+ng)-j0,ny-1))
           k = max(1,min(floor(pz(ip)+1-0.5)-k0,nz-1))

           fcx = px(ip)+1+ng - i - i0;
           fcy = py(ip)+1+0.5+ng - j - j0;
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


       subroutine interp_3d_v(pvar1,px,py,pz,var1,ng,npmx,i0,j0,k0,
     & nx,ny,nz,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz
       integer(kind=4) ,intent(in) :: np,ng
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
!f2py intent(in) ng,np

       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           i = max(1,min(floor(px(ip)+1+0.5+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+ng)-j0,ny-1))
           k = max(1,min(floor(pz(ip)+1-0.5)-k0,nz-1))

           fcx = px(ip)+1+0.5+ng - i - i0;
           fcy = py(ip)+1+ng - j - j0;
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


       subroutine interp_3d_psi(pvar1,px,py,pz,var1,ng,npmx,i0,j0,k0,
     & nx,ny,nz,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz
       integer(kind=4) ,intent(in) :: np,ng
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
!f2py intent(in) ng,np

       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           i = max(1,min(floor(px(ip)+1+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+ng)-j0,ny-1))
           k = max(1,min(floor(pz(ip)+1-0.5)-k0,nz-1))

           fcx = px(ip)+1+ng - i - i0;
           fcy = py(ip)+1+ng - j - j0;
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


       subroutine interp_3d_w(pvar1,px,py,pz,var1,ng,npmx,i0,j0,k0,
     & nx,ny,nz,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz
       integer(kind=4) ,intent(in) :: np,ng
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
!f2py intent(in) ng,np

       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           i = max(1,min(floor(px(ip)+1+0.5+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+0.5+ng)-j0,ny-1))
           k = max(1,min(floor(pz(ip)+1)-k0,nz-1))

           fcx = px(ip)+1+0.5+ng - i - i0;
           fcy = py(ip)+1+0.5+ng - j - j0;
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


       subroutine interp_3d_psiw(pvar1,px,py,pz,var1,ng,npmx,i0,j0,k0,
     & nx,ny,nz,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz
       integer(kind=4) ,intent(in) :: np,ng
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
!f2py intent(in) ng,np

       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           i = max(1,min(floor(px(ip)+1+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+ng)-j0,ny-1))
           k = max(1,min(floor(pz(ip)+1)-k0,nz-1))

           fcx = px(ip)+1+ng - i - i0;
           fcy = py(ip)+1+ng - j - j0;
           fcz = pz(ip)+1 - k - k0;

           CALL linear_3d(fcx,fcy,fcz,wt3)
           pvar1(ip) = sum(var1(i:i+1,j:j+1,k:k+1)*wt3)

         endif
       enddo
!
       end




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate T,S at each particle position (same than interp_2d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_2d_ts(pvar1,pvar2,px,py,
     & var1,var2,ng,npmx,i0,j0,nx,ny,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny
       integer(kind=4) ,intent(in) :: np,ng
       integer(kind=4) ,intent(in) :: npmx
       integer(kind=4) ,intent(in) :: i0,j0
       real(kind=8) ,dimension(nx,ny),intent(in) :: var1,var2
       real(kind=8) ,dimension(np) ,intent(out):: pvar1,pvar2
       real(kind=8) ,dimension(np) ,intent(in):: px,py
! local
       integer(kind=8) :: ip,jp,i,j
       real(kind=8) ,dimension(2,2) :: wt2
       real(kind=8) :: fcx,fcy

!f2py intent(out) pvar1,pvar2
!f2py intent(in) var1,var2
!f2py intent(in) px,py
!f2py intent(in) npmx
!f2py intent(in) nx,ny
!f2py intent(in) i0,j0
!f2py intent(in) ng,np

       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip))) ) then

           i = max(1,min(floor(px(ip)+1+0.5+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+0.5+ng)-j0,ny-1))

           fcx = px(ip)+1+0.5+ng - i - i0;
           fcy = py(ip)+1+0.5+ng - j - j0;

           CALL linear_2d(fcx,fcy,wt2)
           pvar1(ip) = sum(var1(i:i+1,j:j+1)*wt2)
           pvar2(ip) = sum(var2(i:i+1,j:j+1)*wt2)

         endif
       enddo
!
       end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate u at each horizontal particle position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_2d_u(pvar1,px,py,var1,ng,npmx,i0,j0,nx,ny,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny
       integer(kind=4) ,intent(in) :: ng,np
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
!f2py intent(in) ng,np


       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip))) ) then

           i = max(1,min(floor(px(ip)+1+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+0.5+ng)-j0,ny-1))

           fcx = px(ip)+1+ng - i - i0;
           fcy = py(ip)+1+0.5+ng - j - j0;

           CALL linear_2d(fcx,fcy,wt)
           pvar1(ip) = sum(var1(i:i+1,j:j+1)*wt)

         endif
       enddo
!
       end





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate v at each horizontal particle position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_2d_v(pvar1,px,py,var1,ng,npmx,i0,j0,nx,ny,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny
       integer(kind=4) ,intent(in) :: ng,np
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
!f2py intent(in) ng,np


       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip))) ) then

           i = max(1,min(floor(px(ip)+1+0.5+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+ng)-j0,ny-1))

           fcx = px(ip)+1+0.5+ng - i - i0;
           fcy = py(ip)+1+ng - j - j0;

           CALL linear_2d(fcx,fcy,wt)
           pvar1(ip) = sum(var1(i:i+1,j:j+1)*wt)

         endif
       enddo
!
       end





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate a 2D variable at each horizontal particle position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_2d(pvar1,px,py,var1,ng,npmx,i0,j0,nx,ny,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny
       integer(kind=4) ,intent(in) :: ng,np
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
!f2py intent(in) ng,np


       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip))) ) then

           i = max(1,min(floor(px(ip)+1+0.5+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+0.5+ng)-j0,ny-1))

           fcx = px(ip)+1+0.5+ng - i - i0;
           fcy = py(ip)+1+0.5+ng - j - j0;

           CALL linear_2d(fcx,fcy,wt)
           pvar1(ip) = sum(var1(i:i+1,j:j+1)*wt)

         endif
       enddo
!
       end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate a 2D variable at each horizontal particle position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_2d_psi(pvar1,px,py,var1,ng,npmx,i0,j0,nx,ny,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny
       integer(kind=4) ,intent(in) :: ng,np
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
!f2py intent(in) ng,np


       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip))) ) then

           i = max(1,min(floor(px(ip)+1+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+ng)-j0,ny-1))

           fcx = px(ip)+1+ng - i - i0;
           fcy = py(ip)+1+ng - j - j0;

           CALL linear_2d(fcx,fcy,wt)
           pvar1(ip) = sum(var1(i:i+1,j:j+1)*wt)

         endif
       enddo
!
       end
# 21 "pyticles_3d_sig_sa.F" 2



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Time-stepping routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       !---------------------------------------------------------------------!
       ! Forward Euler time-stepping
       !---------------------------------------------------------------------

       subroutine timestep_FE(px,py,pz,u,v,w,itim,fct,dfct,

     & mask,

     & pm,pn,dz,dt,ng,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz
       integer(kind=4) ,intent(in) :: ng,np
       integer(kind=4) ,intent(in) :: npmx
       integer(kind=4) ,intent(in) :: i0,j0,k0
       integer(kind=4) ,dimension(0:1) ,intent(in) :: itim
       real(kind=8) ,dimension(nx-1,ny,nz,2),intent(in) :: u
       real(kind=8) ,dimension(nx,ny-1,nz,2),intent(in) :: v
       real(kind=8) ,dimension(nx,ny,nz+1,2),intent(in) :: w
       real(kind=8) ,dimension(nx,ny,nz,2),intent(in) :: dz
       real(kind=8) ,dimension(nx,ny),intent(in) :: pm,pn

       real(kind=8) ,dimension(nx,ny),intent(in) :: mask

       real(kind=8) ,dimension(np) ,intent(inout):: px,py,pz
       real(kind=8) ,intent(in) :: dt,fct,dfct
! local
       real(kind=8) :: dpx,dpy,dpz
       integer(kind=8) :: ip

!f2py intent(inout) px,py,pz
!f2py intent(in) u,v,w
!f2py intent(in) itim
!f2py intent(in) fct,dfct,dt
!f2py intent(in) pm,pn
!f2py intent(in) dz
!f2py intent(in) ng,npmx
!f2py intent(in) i0,j0,k0
!f2py intent(in) nx,ny,nz
!f2py intent(in) np

!f2py intent(in) mask


       do ip = 1,min(npmx,np)

         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,itim,fct,pm,pn
     & ,dz,dt,i0,j0,k0,nx,ny,nz,ng,np,dpx,dpy,dpz)

           !---------------------------------------
           ! Update position
           !---------------------------------------

           call check_mask(mask,px(ip),py(ip),dpx,dpy
     & ,ng,npmx,i0,j0,nx,ny)

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

     & mask,

     & pm,pn,dz,dt,ng,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz
       integer(kind=4) ,intent(in) :: ng,np
       integer(kind=4) ,intent(in) :: npmx
       integer(kind=4) ,intent(in) :: i0,j0,k0
       integer(kind=4) ,dimension(0:1) ,intent(in) :: itim
       real(kind=8) ,dimension(nx-1,ny,nz,2),intent(in) :: u
       real(kind=8) ,dimension(nx,ny-1,nz,2),intent(in) :: v
       real(kind=8) ,dimension(nx,ny,nz+1,2),intent(in) :: w
       real(kind=8) ,dimension(nx,ny,nz,2),intent(in) :: dz
       real(kind=8) ,dimension(nx,ny),intent(in) :: pm,pn

       real(kind=8) ,dimension(nx,ny),intent(in) :: mask

       real(kind=8) ,dimension(np) ,intent(inout):: px,py,pz
       real(kind=8) ,intent(in) :: dt,fct,dfct
! local
       real(kind=8) :: dpx,dpy,dpz
       integer(kind=8) :: ip



!f2py intent(inout) px,py,pz
!f2py intent(in) u,v,w
!f2py intent(in) itim
!f2py intent(in) fct,dfct,dt
!f2py intent(in) pm,pn
!f2py intent(in) dz
!f2py intent(in) ng,npmx
!f2py intent(in) i0,j0,k0
!f2py intent(in) nx,ny,nz
!f2py intent(in) np

!f2py intent(in) mask


       do ip = 1,min(npmx,np)

         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           ! Forward Euler / predictor
           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,itim,
     & fct,pm,pn,dz,dt,i0,j0,k0,nx,ny,nz,ng,np,dpx,dpy,dpz)

           if (.not.(isnan(dpx*dpy*dpz)) ) then

             ! midpoint rule / corrector
           CALL advance_3d(px(ip)+0.5*dpx,py(ip)+0.5*dpy,pz(ip)+0.5*dpz,
     & u,v,w,itim,fct+0.5*dfct,pm,pn,dz,dt,i0,j0,k0,
     & nx,ny,nz,ng,np,dpx,dpy,dpz)

           endif

             !---------------------------------------
             ! Update position
             !---------------------------------------

           call check_mask(mask,px(ip),py(ip),dpx,dpy
     & ,ng,npmx,i0,j0,nx,ny)

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

     & mask,

     & dz,dt,ng,npmx,i0,j0,k0,nx,ny,nz,np,dpxi,dpyi,dpzi)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz
       integer(kind=4) ,intent(in) :: ng,np
       integer(kind=4) ,intent(in) :: npmx
       integer(kind=4) ,intent(in) :: i0,j0,k0
       integer(kind=4) ,dimension(0:1) ,intent(in) :: itim
       real(kind=8) ,dimension(nx-1,ny,nz,2),intent(in) :: u
       real(kind=8) ,dimension(nx,ny-1,nz,2),intent(in) :: v
       real(kind=8) ,dimension(nx,ny,nz+1,2),intent(in) :: w
       real(kind=8) ,dimension(nx,ny,nz,2),intent(in) :: dz
       real(kind=8) ,dimension(nx,ny),intent(in) :: pm,pn

       real(kind=8) ,dimension(nx,ny),intent(in) :: mask

       real(kind=8) ,dimension(np) ,intent(inout):: px,py,pz
       real(kind=8) ,dimension(np) ,intent(out) :: dpxi,dpyi,dpzi
       real(kind=8) ,intent(in) :: dt,fct,dfct
! local
       real(kind=8) ,dimension(0:3) :: dpx,dpy,dpz
       real(kind=8) :: coef
       integer(kind=8) :: ip

!f2py intent(inout) px,py,pz
!f2py intent(in) u,v,w
!f2py intent(in) itim
!f2py intent(in) fct,dfct,dt
!f2py intent(in) pm,pn
!f2py intent(in) dz
!f2py intent(in) ng,npmx
!f2py intent(in) i0,j0,k0
!f2py intent(in) nx,ny,nz
!f2py intent(in) np

!f2py intent(in) mask

!f2py intent(out) dpxi,dpyi,dpzi

       coef = 1./6.

       do ip = 1,min(npmx,np)

         if (.not.(isnan(px(ip)*py(ip)*pz(ip)))) then

           ! Forward Euler / predictor
           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,
     & itim,fct,pm,pn,dz,dt,i0,j0,k0,
     & nx,ny,nz,ng,np,dpx(0),dpy(0),dpz(0))

           if (.not.(isnan(dpx(0)*dpy(0)*dpz(0)))) then

           ! backward Euler / corrector
           CALL advance_3d(px(ip)+0.5*dpx(0),py(ip)+0.5*dpy(0),
     & pz(ip)+0.5*dpz(0),
     & u,v,w,itim,fct+0.5*dfct,pm,pn,dz,dt,i0,j0,k0,
     & nx,ny,nz,ng,np,dpx(1),dpy(1),dpz(1))

            if (.not.(isnan(dpx(1)*dpy(1)*dpz(1)))) then

           ! midpoint rule / predictor
           CALL advance_3d(px(ip)+0.5*dpx(1),py(ip)+0.5*dpy(1),
     & pz(ip)+0.5*dpz(1),
     & u,v,w,itim,fct+0.5*dfct,pm,pn,dz,dt,i0,j0,k0,
     & nx,ny,nz,ng,np,dpx(2),dpy(2),dpz(2))

             if (.not.(isnan(dpx(2)*dpy(2)*dpz(2)))) then

           ! Corrector
           CALL advance_3d(px(ip)+dpx(2),py(ip)+dpy(2),
     & pz(ip)+dpz(2),
     & u,v,w,itim,fct+1.*dfct,pm,pn,dz,dt,i0,j0,k0,
     & nx,ny,nz,ng,np,dpx(3),dpy(3),dpz(3))

             !else
              !write(*,*) 'fail MP 3',ip,px(ip),py(ip),pz(ip)
             endif

            !else
             !write(*,*) 'fail BE 2',ip,px(ip),py(ip),pz(ip)
            endif

           !else
             !write(*,*) '_________________________________'
             !write(*,*) 'fail FE 1',ip,px(ip),py(ip),pz(ip)
             !write(*,*) 'fail FE 1',dpx(0),dpy(0),dpz(0)
             !write(*,*) 'fail FE 1',fct,dfct
             !write(*,*) '_________________________________'
           endif




           !---------------------------------------
           ! Update position
           !---------------------------------------


           dpxi(ip) = coef * (dpx(0) + 2 * dpx(1) + 2 * dpx(2) + dpx(3))
           dpyi(ip) = coef * (dpy(0) + 2 * dpy(1) + 2 * dpy(2) + dpy(3))
           dpzi(ip) = coef * (dpz(0) + 2 * dpz(1) + 2 * dpz(2) + dpz(3))


           !write(*,*) 'in RK4',px(ip),py(ip),dpxi(ip),dpyi(ip)
           call check_mask(mask,px(ip),py(ip),dpxi(ip),dpyi(ip)
     & ,ng,npmx,i0,j0,nx,ny)


           px(ip) = px(ip) + dpxi(ip)
           py(ip) = py(ip) + dpyi(ip)
           pz(ip) = pz(ip) + dpzi(ip)


          endif

       enddo


       end


        !---------------------------------------------------------------------!
        ! Adams-Bashforth 2 time-stepping
        !---------------------------------------------------------------------
        ! For Adams-methods we need to keep the previous dpx,dpy,dpz values
        !

       subroutine timestep_AB2(px,py,pz,dpx,dpy,dpz,iab,

     & mask,

     & u,v,w,itim,fct,dfct,pm,pn,dz,dt,ng,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none

! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz,ng,np,npmx
       integer(kind=4) ,intent(in) :: i0,j0,k0
       integer(kind=4) ,dimension(0:1) ,intent(in) :: itim
       real(kind=8),dimension(nx-1,ny,nz,2),intent(in) :: u
       real(kind=8),dimension(nx,ny-1,nz,2),intent(in) :: v
       real(kind=8),dimension(nx,ny,nz+1,2),intent(in) :: w
       real(kind=8),dimension(nx,ny,nz,2) ,intent(in) :: dz
       real(kind=8),dimension(nx,ny) ,intent(in) :: pm,pn

       real(kind=8) ,dimension(nx,ny),intent(in) :: mask

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
!f2py intent(in) ng,npmx
!f2py intent(in) nx,ny,nz
!f2py intent(in) i0,j0,k0
!f2py intent(in) itim
!f2py intent(in) np
!f2py intent(in) iab

!f2py intent(in) mask

!f2py intent(inout) dpx,dpy,dpz


       coef = 0.5

       do ip = 1,min(npmx,np)
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           ! Forward Euler / predictor
           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,itim,
     & fct,pm,pn,dz,dt,i0,j0,k0,nx,ny,nz,ng,np,
     & dpx(ip,iab(1)),dpy(ip,iab(1)),dpz(ip,iab(1)))



            dpxi = coef * (3. * dpx(ip,iab(1)) - dpx(ip,iab(0)))
            dpyi = coef * (3. * dpy(ip,iab(1)) - dpy(ip,iab(0)))
            dpzi = coef * (3. * dpz(ip,iab(1)) - dpz(ip,iab(0)))

            !---------------------------------------
            ! Update position
            !---------------------------------------


           call check_mask(mask,px(ip),py(ip),dpxi,dpyi
     & ,ng,npmx,i0,j0,nx,ny)

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

     & mask,

     & u,v,w,itim,fct,dfct,pm,pn,dz,dt,ng,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none

! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz,ng,np,npmx
       integer(kind=4) ,intent(in) :: i0,j0,k0
       integer(kind=4) ,dimension(0:1) ,intent(in) :: itim
       real(kind=8),dimension(nx-1,ny,nz,2),intent(in) :: u
       real(kind=8),dimension(nx,ny-1,nz,2),intent(in) :: v
       real(kind=8),dimension(nx,ny,nz+1,2),intent(in) :: w
       real(kind=8),dimension(nx,ny,nz,2) ,intent(in) :: dz
       real(kind=8),dimension(nx,ny) ,intent(in) :: pm,pn

       real(kind=8) ,dimension(nx,ny),intent(in) :: mask

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
!f2py intent(in) ng,npmx
!f2py intent(in) nx,ny,nz
!f2py intent(in) i0,j0,k0
!f2py intent(in) itim
!f2py intent(in) np

!f2py intent(in) mask

!f2py intent(inout) iab
!f2py intent(inout) dpx,dpy,dpz


       coef = 1./12.

       do ip = 1,min(npmx,np)
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           ! Forward Euler / predictor
           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,itim,
     & fct,pm,pn,dz,dt,i0,j0,k0,nx,ny,nz,ng,np,
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


           call check_mask(mask,px(ip),py(ip),dpxi,dpyi
     & ,ng,npmx,i0,j0,nx,ny)


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

     & mask,

     & u,v,w,itim,fct,dfct,pm,pn,dz,dt,ng,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none

! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz,ng,np,npmx
       integer(kind=4) ,intent(in) :: i0,j0,k0
       integer(kind=4) ,dimension(0:1) ,intent(in) :: itim
       real(kind=8),dimension(nx-1,ny,nz,2),intent(in) :: u
       real(kind=8),dimension(nx,ny-1,nz,2),intent(in) :: v
       real(kind=8),dimension(nx,ny,nz+1,2),intent(in) :: w
       real(kind=8),dimension(nx,ny,nz,2) ,intent(in) :: dz
       real(kind=8),dimension(nx,ny) ,intent(in) :: pm,pn

       real(kind=8) ,dimension(nx,ny),intent(in) :: mask

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
!f2py intent(in) ng,npmx
!f2py intent(in) nx,ny,nz
!f2py intent(in) i0,j0,k0
!f2py intent(in) itim
!f2py intent(in) np

!f2py intent(in) mask

!f2py intent(inout) iab
!f2py intent(inout) dpx,dpy,dpz

       coef = 1./24.

       do ip = 1,min(npmx,np)
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           ! Forward Euler / predictor
           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,itim,
     & fct,pm,pn,dz,dt,i0,j0,k0,nx,ny,nz,ng,np,
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


            !---------------------------------------
            ! Update position
            !---------------------------------------


           call check_mask(mask,px(ip),py(ip),dpxi,dpyi
     & ,ng,npmx,i0,j0,nx,ny)


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

     & mask,

     & u,v,w,itim,fct,dfct,pm,pn,dz,dt,ng,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none

! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz,ng,np,npmx
       integer(kind=4) ,intent(in) :: i0,j0,k0
       integer(kind=4) ,dimension(0:1) ,intent(in) :: itim
       real(kind=8),dimension(nx-1,ny,nz,2),intent(in) :: u
       real(kind=8),dimension(nx,ny-1,nz,2),intent(in) :: v
       real(kind=8),dimension(nx,ny,nz+1,2),intent(in) :: w
       real(kind=8),dimension(nx,ny,nz,2) ,intent(in) :: dz
       real(kind=8),dimension(nx,ny) ,intent(in) :: pm,pn

       real(kind=8) ,dimension(nx,ny),intent(in) :: mask

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
!f2py intent(in) ng,np

!f2py intent(in) mask

!f2py intent(inout) iab
!f2py intent(inout) dpx,dpy,dpz



       coef = 1./24.

       do ip = 1,min(npmx,np)
         if (.not.(isnan(px(ip))) ) then

           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,itim,
     & fct,pm,pn,dz,dt,i0,j0,k0,nx,ny,nz,ng,np,
     & dpx(ip,iab(3)),dpy(ip,iab(3)),dpz(ip,iab(3)))

           !---------------------------------------
           ! Predictor
           !---------------------------------------


            dpxi = coef * (55.*dpx(ip,iab(3)) - 59.*dpx(ip,iab(2))
     & + 37.*dpx(ip,iab(1)) - 9.*dpx(ip,iab(0)))
            dpyi = coef * (55.*dpy(ip,iab(3)) - 59.*dpy(ip,iab(2))
     & + 37.*dpy(ip,iab(1)) - 9.*dpy(ip,iab(0)))
            dpzi = coef * (55.*dpz(ip,iab(3)) - 59.*dpz(ip,iab(2))
     & + 37.*dpz(ip,iab(1)) - 9.*dpz(ip,iab(0)))


           !---------------------------------------
           ! Compute velocities at new time-step
           !---------------------------------------

           CALL advance_3d(px(ip)+dpxi,py(ip)+dpyi,pz(ip)+dpzi,
     & u,v,w,itim,
     & fct+dfct,pm,pn,dz,dt,i0,j0,k0,nx,ny,nz,ng,np,dpxc,dpyc,dpzc)


           !---------------------------------------
           ! Corrector
           !---------------------------------------

            dpxi = 1./270. * ( 19.* dpxi + 251.* ( 1./24. * (
     & 9. * dpxc + 19. * dpx(ip,iab(3))
     & - 5. * dpx(ip,iab(2)) + dpx(ip,iab(1)) )))

            dpyi = 1./270. * ( 19.* dpyi + 251.* ( 1./24. * (
     & 9. * dpyc + 19. * dpy(ip,iab(3))
     & - 5. * dpy(ip,iab(2)) + dpy(ip,iab(1)) )))

            dpzi = 1./270. * ( 19.* dpzi + 251.* ( 1./24. * (
     & 9. * dpzc + 19. * dpz(ip,iab(3))
     & - 5. * dpz(ip,iab(2)) + dpz(ip,iab(1)) )))


            !---------------------------------------
            ! Update position
            !---------------------------------------


           call check_mask(mask,px(ip),py(ip),dpxi,dpyi
     & ,ng,npmx,i0,j0,nx,ny)


            px(ip) = px(ip) + dpxi
            py(ip) = py(ip) + dpyi
            pz(ip) = pz(ip) + dpzi


         endif
       enddo



       end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Else
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       subroutine check_mask(mask,px,py,dpx,dpy,ng,npmx,i0,j0,nx,ny)

       !---------------------------------------------------------------------!
       ! Check mask routine (check if particle is overshooting into mask)
       !---------------------------------------------------------------------


       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny
       integer(kind=4) ,intent(in) :: ng,npmx
       real(kind=8),dimension(nx,ny) ,intent(in) :: mask
       real(kind=8) ,intent(in) :: px,py
       real(kind=8) ,intent(inout):: dpx,dpy
       integer(kind=4) ,intent(in) :: i0,j0
! local
       integer(kind=8) :: i,j
       real(kind=8),dimension(2,2) :: wt
       real(kind=8) :: fcx,fcy,pmask
       real(kind=8),dimension(4) :: dmask


!f2py intent(inout) dpx,dpy
!f2py intent(in) mask,px,py
!f2py intent(in) npmx
!f2py intent(in) i0,j0
!f2py intent(in) ng
!f2py intent(in) nx,ny

        !---------------------------------------------------------------------
        ! compute pmask
        !---------------------------------------------------------------------

        i = max(1,min(floor(px+dpx+1+0.5+ng)-i0,nx-1))
        j = max(1,min(floor(py+dpy+1+0.5+ng)-j0,ny-1))

        fcx = px+dpx+1+0.5+ng - i - i0;
        fcy = py+dpy+1+0.5+ng - j - j0;

        CALL linear_2d(fcx,fcy,wt)
        pmask = sum(mask(i:i+1,j:j+1)*wt)


        !write(*,*) 'pmask ', pmask

        if (pmask.lt.1.) then

          !write(*,*) ' '
          !write(*,*) 'correcting mask for ', px,py
          !write(*,*) 'correcting ', dpx,dpy

          dmask(1) = (1-fcy)*mask(i,j) + fcy*mask(i,j+1)
          dmask(2) = (1-fcy)*mask(i+1,j) + fcy*mask(i+1,j+1)
          dmask(3) = (1-fcx)*mask(i,j) + fcx*mask(i+1,j)
          dmask(4) = (1-fcx)*mask(i,j+1) + fcx*mask(i+1,j+1)

          !write(*,*) 'dmask(1),dmask(2) ', dmask(1),dmask(2)
          !write(*,*) 'dmask(3),dmask(4) ', dmask(3),dmask(4)

          if ( ((dpx.gt.0).and.((dmask(1)-dmask(2)).ge.0.5)) .or.
     & ((dpx.lt.0).and.((dmask(2)-dmask(1)).ge.0.5)) ) dpx = 0.

          if ( ((dpy.gt.0).and.((dmask(3)-dmask(4)).ge.0.5)) .or.
     & ((dpy.lt.0).and.((dmask(4)-dmask(3)).ge.0.5)) ) dpy = 0.

          !write(*,*) 'correcting mask for ', px,py
          !write(*,*) 'corrected ', dpx,dpy

        endif

       end







       !---------------------------------------------------------------------!
       ! Runge-Kutta 4 time-stepping for 2D advection
       !---------------------------------------------------------------------

!


       subroutine timestep2d_RK4(px,py,u,v,itim,fct,dfct,pm,pn,

     & mask,

     & dt,ng,npmx,i0,j0,nx,ny,np,dpxi,dpyi)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny
       integer(kind=4) ,intent(in) :: ng,np
       integer(kind=4) ,intent(in) :: npmx
       integer(kind=4) ,intent(in) :: i0,j0
       integer(kind=4) ,dimension(0:1) ,intent(in) :: itim
       real(kind=8) ,dimension(nx-1,ny,2),intent(in) :: u
       real(kind=8) ,dimension(nx,ny-1,2),intent(in) :: v
       real(kind=8) ,dimension(nx,ny),intent(in) :: pm,pn

       real(kind=8) ,dimension(nx,ny),intent(in) :: mask

       real(kind=8) ,dimension(np) ,intent(inout):: px,py
       real(kind=8) ,dimension(np) ,intent(out) :: dpxi,dpyi
       real(kind=8) ,intent(in) :: dt,fct,dfct
! local
       real(kind=8) ,dimension(0:3) :: dpx,dpy
       real(kind=8) :: coef
       integer(kind=8) :: ip

!f2py intent(inout) px,py
!f2py intent(in) u,v
!f2py intent(in) itim
!f2py intent(in) fct,dfct,dt
!f2py intent(in) pm,pn
!f2py intent(in) ng,npmx
!f2py intent(in) i0,j0
!f2py intent(in) nx,ny
!f2py intent(in) np

!f2py intent(in) mask

!f2py intent(out) dpxi,dpyi

       coef = 1./6.

       do ip = 1,min(npmx,np)

         if (.not.(isnan(px(ip)*py(ip)))) then

           ! Forward Euler / predictor
           CALL advance_2d(px(ip),py(ip),u,v,
     & itim,fct,pm,pn,dt,i0,j0,
     & nx,ny,ng,np,dpx(0),dpy(0))

           if (.not.(isnan(dpx(0)*dpy(0)))) then

           ! backward Euler / corrector
           CALL advance_2d(px(ip)+0.5*dpx(0),py(ip)+0.5*dpy(0),
     & u,v,itim,fct+0.5*dfct,pm,pn,dt,i0,j0,
     & nx,ny,ng,np,dpx(1),dpy(1))

            if (.not.(isnan(dpx(1)*dpy(1)))) then

           ! midpoint rule / predictor
           CALL advance_2d(px(ip)+0.5*dpx(1),py(ip)+0.5*dpy(1),
     & u,v,itim,fct+0.5*dfct,pm,pn,dt,i0,j0,
     & nx,ny,ng,np,dpx(2),dpy(2))

             if (.not.(isnan(dpx(2)*dpy(2)))) then

           ! Corrector
           CALL advance_2d(px(ip)+dpx(2),py(ip)+dpy(2),
     & u,v,itim,fct+1.*dfct,pm,pn,dt,i0,j0,
     & nx,ny,ng,np,dpx(3),dpy(3))

             endif

            endif

           endif




           !---------------------------------------
           ! Update position
           !---------------------------------------

           dpxi(ip) = coef * (dpx(0) + 2 * dpx(1) + 2 * dpx(2) + dpx(3))
           dpyi(ip) = coef * (dpy(0) + 2 * dpy(1) + 2 * dpy(2) + dpy(3))

           px(ip) = px(ip) + dpxi(ip)
           py(ip) = py(ip) + dpyi(ip)


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
# 928 "pyticles_3d_sig_sa.F" 2

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
# 930 "pyticles_3d_sig_sa.F" 2

# 1 "R_tools_fort_routines/zlevs_croco_new.F" 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute z_r and z_w for new or old S Coord
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





      subroutine zlevs_croco_new(Lm,Mm,N, h,zeta, hc, Cs_r, Cs_w,
     & sc_r,sc_w,z_r,z_w)

      implicit none

      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k

      real*8 Cs_w(0:N), Cs_r(N), cff_w, cff_r, cff1_w, cff1_r,
     & sc_w(0:N), sc_r(N),
     & hc, ds, z_w0, z_r0,
     & zeta(0:Lm+1,0:Mm+1),
     & z_r(0:Lm+1,0:Mm+1,N),z_w(0:Lm+1,0:Mm+1,0:N),
     & h(0:Lm+1,0:Mm+1),hinv(0:Lm+1,0:Mm+1)


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

        do k=1,N,+1 !--> irreversible because of recursion in Hz


          cff_w =hc*sc_w(k)
          cff_r =hc*sc_r(k)
          cff1_w=Cs_w(k)
          cff1_r=Cs_r(k)
# 64 "R_tools_fort_routines/zlevs_croco_new.F"
          do i=imin,imax

            z_w0=cff_w+cff1_w*h(i,j)
            z_r0=cff_r+cff1_r*h(i,j)


            z_w(i,j,k)=z_w0*h(i,j)*hinv(i,j)+zeta(i,j)
     & *(1.+z_w0*hinv(i,j))
            z_r(i,j,k)=z_r0*h(i,j)*hinv(i,j)+zeta(i,j)
     & *(1.+z_r0*hinv(i,j))





          enddo
        enddo
      enddo
      end
# 932 "pyticles_3d_sig_sa.F" 2

# 1 "R_tools_fort_routines/zlevs_croco_new_w.F" 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute z_r and z_w for new or old S Coord
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





      subroutine zlevs_croco_new_w(Lm,Mm,N, h,zeta, hc, Cs_w,
     & sc_w ,z_w)

      implicit none

      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k

      real*8 Cs_w(0:N), cff_w, cff1_w,
     & sc_w(0:N),
     & hc, ds, z_w0,
     & zeta(0:Lm+1,0:Mm+1),
     & z_w(0:Lm+1,0:Mm+1,0:N),
     & h(0:Lm+1,0:Mm+1),hinv(0:Lm+1,0:Mm+1)


Cf2py intent(in) Lm,Mm,N, h,zeta, hc, Cs_w, sc_w
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


          cff_w =hc*sc_w(k)
          cff1_w=Cs_w(k)







          do i=imin,imax

            z_w0=cff_w+cff1_w*h(i,j)


            z_w(i,j,k)=z_w0*h(i,j)*hinv(i,j)+zeta(i,j)
     & *(1.+z_w0*hinv(i,j))




          enddo
        enddo
      enddo
      end
# 934 "pyticles_3d_sig_sa.F" 2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute z_r and z_w for OLD_S_COORD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# 1 "R_tools_fort_routines/zlevs_croco_old.F" 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute z_r and z_w for new or old S Coord
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





      subroutine zlevs_croco_old(Lm,Mm,N, h,zeta, hc, Cs_r, Cs_w,
     & sc_r,sc_w,z_r,z_w)

      implicit none

      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k

      real*8 Cs_w(0:N), Cs_r(N), cff_w, cff_r, cff1_w, cff1_r,
     & sc_w(0:N), sc_r(N),
     & hc, ds, z_w0, z_r0,
     & zeta(0:Lm+1,0:Mm+1),
     & z_r(0:Lm+1,0:Mm+1,N),z_w(0:Lm+1,0:Mm+1,0:N),
     & h(0:Lm+1,0:Mm+1),hinv(0:Lm+1,0:Mm+1)


Cf2py intent(in) Lm,Mm,N, h,zeta, hc, Cs_w, Cs_r, sc_r,sc_w
Cf2py intent(out) z_r,z_w


      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1


      ds=1.D0/dble(N)

      do j=jmin,jmax
        do i=imin,imax




            hinv(i,j)=1./h(i,j)


          z_w(i,j,0)=-h(i,j)

        enddo

        do k=1,N,+1 !--> irreversible because of recursion in Hz







          cff_w =hc*(sc_w(k)-Cs_w(k))
          cff_r =hc*(sc_r(k)-Cs_r(k))
          cff1_w=Cs_w(k)
          cff1_r=Cs_r(k)



          do i=imin,imax

            z_w0=cff_w+cff1_w*h(i,j)
            z_r0=cff_r+cff1_r*h(i,j)







            z_w(i,j,k)=z_w0+zeta(i,j)*(1.+z_w0*hinv(i,j))
            z_r(i,j,k)=z_r0+zeta(i,j)*(1.+z_r0*hinv(i,j))


          enddo
        enddo
      enddo
      end
# 940 "pyticles_3d_sig_sa.F" 2

# 1 "R_tools_fort_routines/zlevs_croco_old_w.F" 1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute z_r and z_w for new or old S Coord
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!





      subroutine zlevs_croco_old_w(Lm,Mm,N, h,zeta, hc, Cs_w,
     & sc_w,z_w)

      implicit none

      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k

      real*8 Cs_w(0:N), cff_w, cff1_w,
     & sc_w(0:N),
     & hc, ds, z_w0,
     & zeta(0:Lm+1,0:Mm+1),
     & z_w(0:Lm+1,0:Mm+1,0:N),
     & h(0:Lm+1,0:Mm+1),hinv(0:Lm+1,0:Mm+1)


Cf2py intent(in) Lm,Mm,N, h,zeta, hc, Cs_w,sc_w
Cf2py intent(out) z_w


      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1


      ds=1.D0/dble(N)

      do j=jmin,jmax
        do i=imin,imax




            hinv(i,j)=1./h(i,j)


          z_w(i,j,0)=-h(i,j)

        enddo

        do k=1,N,+1 !--> irreversible because of recursion in Hz





          cff_w =hc*(sc_w(k)-Cs_w(k))
          cff1_w=Cs_w(k)



          do i=imin,imax

            z_w0=cff_w+cff1_w*h(i,j)





            z_w(i,j,k)=z_w0+zeta(i,j)*(1.+z_w0*hinv(i,j))


          enddo
        enddo
      enddo
      end
# 942 "pyticles_3d_sig_sa.F" 2

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
# 948 "pyticles_3d_sig_sa.F" 2

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
# 950 "pyticles_3d_sig_sa.F" 2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Omega
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# 1 "R_tools_fort_routines/get_omega.F" 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute vertical velocity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






      subroutine get_omega(Lm,Mm,N,u,v, z_r,z_w,pm,pn
     & ,W)


      implicit none

      integer Lm,Mm,N, imin,imax,jmin,jmax, i,j,k,
     & istr,iend,jstr,jend,istrU,jstrV


      real*8 u(1:Lm+1,0:Mm+1,N), v(0:Lm+1,1:Mm+1,N),
     & FlxU(1:Lm+1,0:Mm+1,N), FlxV(0:Lm+1,1:Mm+1,N),
     & z_r(0:Lm+1,0:Mm+1,N), z_w(0:Lm+1,0:Mm+1,0:N),
     & Hz(0:Lm+1,0:Mm+1,N),
     & pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1),
     & dn_u(0:Lm+1,0:Mm+1), dm_v(0:Lm+1,0:Mm+1),
     & var1, var2,var3, var4


      real*8 wrk(0:Lm+1),W(0:Lm+1,0:Mm+1,0:N)

# 1 "R_tools_fort_routines/scalars.h" 1
! This is include file "scalars.h"
!----- -- ------- ---- -----------
! The following common block contains time variables and indices
! for 2D (k-indices) and 3D (n-indices) computational engines. Since
! they are changed together, they are placed into the same cache line
! despite their mixed type, so that only one cachene is being
! invalidated and has to be propagated accross the cluster.
! Additionally, variables proc and CPU_time are to hold process ID
! numbers of individual threads and to measure CPU time consumed by
! each of them during the whole model run (these are for purely
! diagnostic/performance measurements and do not affect the model
! results.)
!
! Note that real values are placed first into the common block before
! integers. This is done to prevent misallignment of the 8-byte
! objects in the case when an uneven number of 4-byte integers is
! placed before a 8-byte real (in the case when default real size is
! set to 8 Bytes). Although misallignment is not formally a violation
! of fortran standard, it may cause performance degradation and/or
! make compiler issue a warning message (Sun, DEC Alpha) or even
! crash (Alpha).
!

!
! Physical constants: Earth radius [m]; Aceleration of gravity
!--------- ---------- duration of the day in seconds; Specific
! heat [Joules/kg/degC] for seawater (it is approximately 4000,
! and varies only slightly, see Gill, 1982, Appendix 3); von
! Karman constant.
!
      real pi, Eradius,g, Cp,vonKar, deg2rad,rad2deg,day2sec,sec2day
      parameter (pi=3.14159265358979323, Eradius=6371315.,
     & deg2rad=pi/180., rad2deg=180./pi, day2sec=86400.,
     & sec2day=1./86400., Cp=3985., vonKar=0.41)
      parameter (g=9.81)
# 35 "R_tools_fort_routines/get_omega.F" 2




Cf2py intent(in) Lm,Mm,N, u,v,z_r,z_w,pm,pn
Cf2py intent(out) W

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1


      istr=1
      iend=Lm
      jstr=1
      jend=Mm

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!
! Compute "omega" vertical velocity by means of integration of mass
! divergence of mass fluxes from bottom up. In this computation,
! unlike that in omega.F, there is (1) immediate multiplication by
! pm*pn so that the result has meaning of velocity, rather than
! finite volume mass flux through vertical facet of tracer grid box;
! and (2, also unlike omega.F) no subtraction of vertical velocity
! of moving grid-box interface (the effect of "breething" of vertical
! grid system due to evolving free surface) is made now.
! Consequently, Wrk(:,N).ne.0, unlike its counterpart W(:,:,N).eqv.0
! in omega.F. Once omega vertical velocity is computed, interpolate
! it to vertical RHO-points.
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



      do j=jmin,jmax
        do i=imin,imax
          do k=1,N,+1
           Hz(i,j,k) = z_w(i,j,k) - z_w(i,j,k-1)
          enddo
        enddo
      enddo




      do j=jmin,jmax
        do i=imin+1,imax
            dn_u(i,j) = 2./(pn(i,j)+pn(i-1,j))
            do k=1,N,+1
              FlxU(i,j,k) = 0.5*(Hz(i,j,k)+Hz(i-1,j,k))*dn_u(i,j)
     & * u(i,j,k)
            enddo
          enddo
      enddo



      do j=jmin+1,jmax
        do i=imin,imax
            dm_v(i,j) = 2./(pm(i,j)+pm(i,j-1))
            do k=1,N,+1
              FlxV(i,j,k) = 0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)
     & * v(i,j,k)
            enddo
          enddo
      enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=jstr,jend
        do i=istr,iend
          W(i,j,0)=0.
        enddo

        do k=1,N,+1 !--> recursive
          do i=istr,iend
            W(i,j,k)=W(i,j,k-1) -FlxU(i+1,j,k)+FlxU(i,j,k)
     & -FlxV(i,j+1,k)+FlxV(i,j,k)
          enddo
        enddo

        do i=istr,iend
          wrk(i)=W(i,j,N)/(z_w(i,j,N)-z_w(i,j,0))
          W(i,j,N)=0.
        enddo

        do k=N-1,1,-1
          do i=istr,iend
            W(i,j,k)=W(i,j,k)-wrk(i)*(z_w(i,j,k)-z_w(i,j,0))
          enddo
        enddo
      enddo



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      do j=jstr,jend
        do i=istr,iend
           do k=1,N,+1

             W(i,j,k)=W(i,j,k)*pm(i,j)*pn(i,j)

          enddo
        enddo
      enddo


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

                       ! Set lateral
        do k=0,N ! boundary
          do j=jstr,jend ! conditions
            W(istr-1,j,k)=W(istr,j,k)
            W(iend+1,j,k)=W(iend,j,k)
          enddo
        enddo

        do k=0,N
          do i=istr,iend
            W(i,jstr-1,k)=W(i,jstr,k)
            W(i,jend+1,k)=W(i,jend,k)
          enddo
        enddo

        do k=0,N
          W(istr-1,jstr-1,k)=W(istr,jstr,k)
        enddo

        do k=0,N
          W(istr-1,jend+1,k)=W(istr,jend,k)
        enddo

        do k=0,N
          W(iend+1, jstr-1,k)=W(iend,jstr,k)
        enddo

        do k=0,N
          W(iend+1,jend+1,k)=W(iend,jend,k)
        enddo

      return
      end
# 957 "pyticles_3d_sig_sa.F" 2
# 966 "pyticles_3d_sig_sa.F"
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
