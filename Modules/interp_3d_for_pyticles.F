
!----------------------------------------------------------------------------------------------
! Spatial Interpolation Routines for Pyticles
!
!----------------------------------------------------------------------------------------------
! 05/10/16:
!     - correction of advance_3d for ORIGINAL_VERSION (computation of fcx_u,fcx_v,fcx_w)
!     - add subroutines interp_3d_u / interp_3d_u
! 16/05/03:
!     - Correction of linear_2d [terms (2,1) and (1,2) inversion]
! 16/01/26:
!     - Add parameter 'ng' corresponding to number of ghost points in the horizontal grid
!     - interp_3d routines only for STAGGERED yet
!----------------------------------------------------------------------------------------------



!----------------------------------------------------------------------------------------------
#undef NEW_VERSION
#define LINEAR_INTERPOLATION
#undef CUBIC_INTERPOLATION
! #define CRSPL_INTERPOLATION
! #define WENO_INTERPOLATION
!----------------------------------------------------------------------------------------------
#define STAGGERED
!---------------------------------------------------------------------!
! if STAGGERED is defined, px,py,pz correpond to staggered positions (u,v,w grids)
! if STAGGERED is NOT defined, px,py,pz correspond to positions on the horizontal rho-grid and vertical w-grid
!---------------------------------------------------------------------!



#ifdef NEW_VERSION 
!----------------------------------------------------------------------------------------------


       subroutine advance_3d(px,py,pz,u,v,w,itim,fct,pm,pn,
     &          dz,dt,i0,j0,k0,nx,ny,nz,ng,np,dpx,dpy,dpz)

       !---------------------------------------------------------------------!
       !  Compute particle displacement with linear interpolation in time
       !  in: u,v,w;  at their original staggered position
       !      px,py,pz;  Particle position in index coordinates [0,nx] etc.... 
       !      pm,pn,dz;  grid stretching at rho position (pm,pn 2d)
       ! 
       !---------------------------------------------------------------------

       implicit none 
!      import/export
       real(kind=8)                                 ,intent(in) :: px,py,pz
       real(kind=8) ,dimension(0:nx-2,0:ny-1,nz,0:1)    ,intent(in) :: u
       real(kind=8) ,dimension(0:nx-1,0:ny-2,nz,0:1)    ,intent(in) :: v
       real(kind=8) ,dimension(0:nx-1,0:ny-1,0:nz,0:1),intent(in) :: w
       real(kind=8) ,dimension(0:nx-1,0:ny-1,nz,0:1)  ,intent(in) :: dz
       real(kind=8) ,dimension(0:nx-1,0:ny-1)       ,intent(in) :: pm,pn
       integer(kind=4)                              ,intent(in) :: nx,ny,nz,np
       integer(kind=4)                              ,intent(in) :: ng
       integer(kind=4)                              ,intent(in) :: i0,j0,k0
       real(kind=8)                                 ,intent(out):: dpx,dpy,dpz
       real(kind=8)                                 ,intent(in) :: dt,fct
       integer(kind=4) ,dimension(0:1)              ,intent(in) :: itim
!      local
       integer(kind=8)                :: i,j,k
       integer(kind=8)                :: iu,jv,kw
       real(kind=8) ,dimension(4,4,4) :: f
       real(kind=8)                   :: pxl,pyl,pzl
       real(kind=8)                   :: pxlu,pylv,pzlw
       real(kind=8)                   :: pu,pv,pw
       real(kind=8)                   :: ppm,ppn,pdz

!f2py intent(in)   u,v,w
!f2py intent(in)   px,py,pz
!f2py intent(out)  dpx,dpy,dpz
!f2py intent(in)   pm,pn,dt,fct
!f2py intent(in)   dz
!f2py intent(in)   nx,ny,nz,ng,np
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   itim



           !! add ghost points to deal with boundaries
           !! for now, we ll add ghost points here at z = 0,nz if
           !! neccessary and rely on outside to keep a 2-point buffer in
           !! the horizontal



# ifdef STAGGERED 
           iu= floor(px+ng-i0);     pxlu= px+ng-i0-iu     !! In x-dir u grid starts at 0
           jv= floor(py+ng-j0);     pylv= py+ng-j0-jv     !! In y-dir v grid starts at 0
           kw= floor(pz-k0);        pzlw= pz-k0-kw     !! 
           i = floor(px+ng-i0+0.5); pxl = px+ng-i0-i+0.5  !! In x-dir rho grid starts at -0.5
           j = floor(py+ng-j0+0.5); pyl = py+ng-j0-j+0.5  !! In y-dir rho grid starts at -0.5
           k = floor(pz-k0+0.5);    pzl = pz-k0-k+0.5  !! 
# else
           iu= floor(px+ng-i0-0.5); pxlu= px+ng-i0-iu-0.5 !! In x-dir u grid starts at 0.5
           jv= floor(py+ng-j0-0.5); pylv= py+ng-j0-jv-0.5 !! In y-dir v grid starts at 0.5
           kw= floor(pz-k0);        pzlw= pz-k0-kw     !! 
           i = floor(px+ng-i0);     pxl = px+ng-i0-i      !! In x-dir rho grid starts at 0
           j = floor(py+ng-j0);     pyl = py+ng-j0-j      !! In y-dir rho grid starts at 0
           k = floor(pz-k0+0.5);    pzl = pz-k0-k+0.5  !! 
# endif


           !! for now, we do slip conditions near the bottom and 
           !! extrapolation near the top.
           if (k .eq. 0) then  
             f(:,:,3:4) = fct * u(iu-1:iu+2,j-1:j+2,k+1:k+2,itim(1)) 
     &              + (1-fct) * u(iu-1:iu+2,j-1:j+2,k+1:k+2,itim(0)) 
             f(:,:,  2) = f(:,:,3)
             f(:,:,  1) = f(:,:,2)
           elseif (k.eq.1) then
             f(:,:,2:4) = fct * u(iu-1:iu+2,j-1:j+2, k:k+2,itim(1)) 
     &              + (1-fct) * u(iu-1:iu+2,j-1:j+2, k:k+2,itim(0)) 
             f(:,:,  1) = f(:,:,2)
           elseif (k.eq.nz-1) then
             f(:,:,1:3) = fct * u(iu-1:iu+2,j-1:j+2,k-1:k+1,itim(1)) 
     &              + (1-fct) * u(iu-1:iu+2,j-1:j+2,k-1:k+1,itim(0)) 
             f(:,:,  4) = 2*f(:,:,3)-f(:,:,2)
           elseif (k.eq.nz) then
             f(:,:,1:2) = fct * u(iu-1:iu+2,j-1:j+2,k-1:k,itim(1)) 
     &              + (1-fct) * u(iu-1:iu+2,j-1:j+2,k-1:k,itim(0)) 
             f(:,:,  3) = 2*f(:,:,2)-f(:,:,1)
             f(:,:,  4) = 2*f(:,:,3)-f(:,:,2)
           else
             f = fct * u(iu-1:iu+2,j-1:j+2,k-1:k+2,itim(1)) 
     &     + (1-fct) * u(iu-1:iu+2,j-1:j+2,k-1:k+2,itim(0)) 
           endif
           call interp3(f,pxlu,pyl,pzl,pu)

           if (k .eq. 0) then  
             f(:,:,3:4) = fct * v(i-1:i+2,jv-1:jv+2,k+1:k+2,itim(1)) 
     &              + (1-fct) * v(i-1:i+2,jv-1:jv+2,k+1:k+2,itim(0)) 
             f(:,:,  2) = f(:,:,3)
             f(:,:,  1) = f(:,:,2)
           elseif (k.eq.1) then
             f(:,:,2:4) = fct * v(i-1:i+2,jv-1:jv+2, k:k+2,itim(1)) 
     &              + (1-fct) * v(i-1:i+2,jv-1:jv+2, k:k+2,itim(0)) 
             f(:,:,  1) = f(:,:,2)
           elseif (k.eq.nz-1) then
             f(:,:,1:3) = fct * v(i-1:i+2,jv-1:jv+2,k-1:k+1,itim(1)) 
     &              + (1-fct) * v(i-1:i+2,jv-1:jv+2,k-1:k+1,itim(0)) 
             f(:,:,  4) = 2*f(:,:,3)-f(:,:,2)
           elseif (k.eq.nz) then
             f(:,:,1:2) = fct * v(i-1:i+2,jv-1:jv+2,k-1:k,itim(1)) 
     &              + (1-fct) * v(i-1:i+2,jv-1:jv+2,k-1:k,itim(0)) 
             f(:,:,  3) = 2*f(:,:,2)-f(:,:,1)
             f(:,:,  4) = 2*f(:,:,3)-f(:,:,2)
           else
             f = fct * v(i-1:i+2,jv-1:jv+2,k-1:k+2,itim(1)) 
     &     + (1-fct) * v(i-1:i+2,jv-1:jv+2,k-1:k+2,itim(0))
           endif
           call interp3(f,pxl,pylv,pzl,pv)

           if (kw .eq. 0) then
             f(:,:,2:4) = fct * w(i-1:i+2,j-1:j+2,kw:kw+2,itim(1)) 
     &              + (1-fct) * w(i-1:i+2,j-1:j+2,kw:kw+2,itim(0)) 
             f(:,:,  1) = 0.
           elseif (kw.eq.nz-1) then
             f(:,:,1:3) = fct * w(i-1:i+2,j-1:j+2,kw-1:kw+1,itim(1)) 
     &             + (1-fct) *  w(i-1:i+2,j-1:j+2,kw-1:kw+1,itim(0)) 
             f(:,:,  4) = 0.
           else
             f = fct * w(i-1:i+2,j-1:j+2,kw-1:kw+2,itim(1)) 
     &     + (1-fct) * w(i-1:i+2,j-1:j+2,kw-1:kw+2,itim(0)) 
           endif

           call interp3(f,pxl,pyl,pzlw,pw)

           if (k .eq. 0) then  
             f(:,:,3:4) = fct * dz(i-1:i+2,j-1:j+2,k+1:k+2,itim(1)) 
     &              + (1-fct) * dz(i-1:i+2,j-1:j+2,k+1:k+2,itim(0)) 
             f(:,:,  2) = f(:,:,3)
             f(:,:,  1) = f(:,:,2)
           elseif (k.eq.1) then
             f(:,:,2:4) = fct * dz(i-1:i+2,j-1:j+2, k:k+2,itim(1)) 
     &              + (1-fct) * dz(i-1:i+2,j-1:j+2, k:k+2,itim(0)) 
             f(:,:,  1) = f(:,:,2)
           elseif (k.eq.nz-1) then
             f(:,:,1:3) = fct * dz(i-1:i+2,j-1:j+2,k-1:k+1,itim(1)) 
     &              + (1-fct) * dz(i-1:i+2,j-1:j+2,k-1:k+1,itim(0)) 
             f(:,:,  4) = 2*f(:,:,3)-f(:,:,2)
           elseif (k.eq.nz) then
             f(:,:,1:2) = fct * dz(i-1:i+2,j-1:j+2,k-1:k,itim(1)) 
     &              + (1-fct) * dz(i-1:i+2,j-1:j+2,k-1:k,itim(0)) 
             f(:,:,  3) = 2*f(:,:,2)-f(:,:,1)
             f(:,:,  4) = 2*f(:,:,3)-f(:,:,2)
           else
             f = fct * dz(i-1:i+2,j-1:j+2,k-1:k+2,itim(1)) 
     &         + (1-fct) * dz(i-1:i+2,j-1:j+2,k-1:k+2,itim(0)) 
           endif
           call interp3(f,pxl,pyl,pzl,pdz)

           if (pdz.gt.0) pdz=1./pdz

         ppm = pm(i  ,j)*(1-pxl)*(1-pyl) + pm(i  ,j+1)*(1-pxl)*pyl + 
     &         pm(i+1,j)*    pxl*(1-pyl) + pm(i+1,j+1)*pxl*pyl
         ppn = pn(i  ,j)*(1-pxl)*(1-pyl) + pn(i  ,j+1)*(1-pxl)*pyl + 
     &         pn(i+1,j)*    pxl*(1-pyl) + pn(i+1,j+1)*pxl*pyl

         dpx = dt*pu*ppm
         dpy = dt*pv*ppn
         dpz = dt*pw*pdz



       end


!----------------------------------------------------------------------------------------------

!----------------------------------------------------------------------------------------------

# ifdef LINEAR_INTERPOLATION

       subroutine interp3(f,x,y,z,fi)

       !---------------------------------------------------------------------!
       !  Compute linear 3d interpolant
       !  in: f(4,4,4) 
       !      x,y,z;  local coordinates (should be between 0 and 1) 
       !  out: fi 
       ! 
       !---------------------------------------------------------------------

       implicit none 
!      import/export
       real(kind=8),dimension(4,4,4),intent(in) :: f
       real(kind=8),                 intent(in) :: x,y,z
       real(kind=8),                 intent(out):: fi
!      local
       real(kind=8),dimension(2,2,2)            :: wt

       CALL linear_3d(x,y,z,wt)
       fi= sum(f(2:3,2:3,2:3)*wt)

       end



# elif defined CUBIC_INTERPOLATION
!----------------------------------------------------------------------------------------------

       subroutine interp3(f,x,y,z,fi)
       !---------------------------------------------------------------------!
       !  Compute cubic 3d interpolant
       !  in: f(4,4,4) 
       !      x,y,z;  local coordinates (should be between 0 and 1) 
       !  out: fi 
       ! 
       !---------------------------------------------------------------------
       implicit none 
!      import/export
       real(kind=8),dimension(4,4,4),intent(in) :: f
       real(kind=8),                 intent(in) :: x,y,z
       real(kind=8),                 intent(out):: fi
!      local
       real(kind=8)  :: x2,x3,y2,y3,z2,z3
       real(kind=8),dimension(4,4) :: fx,a1,a2,a3,a4
       real(kind=8),dimension(4)   :: fy,b1,b2,b3,b4
       real(kind=8)                :: c1,c2,c3,c4


       x2 = x*x; x3 = x*x2 
       y2 = y*y; y3 = y*y2 
       z2 = z*z; z3 = z*z2 

       a4 = -1./6*f(1,:,:) + 0.5*f(2,:,:) - 0.5*f(3,:,:) + 1./6*f(4,:,:);
       a3 =  0.5* f(1,:,:) -     f(2,:,:) + 0.5*f(3,:,:);
       a2 = -1./3*f(1,:,:) - 0.5*f(2,:,:) +     f(3,:,:) - 1./6*f(4,:,:);
       a1 =                      f(2,:,:);

       fx= a4*x3+ a3*x2 + a2*x + a1;

       b4 = -1./6*fx(1,:) + 0.5*fx(2,:) - 0.5*fx(3,:) + 1./6*fx(4,:);
       b3 =  0.5* fx(1,:) -     fx(2,:) + 0.5*fx(3,:);
       b2 = -1./3*fx(1,:) - 0.5*fx(2,:) +     fx(3,:) - 1./6*fx(4,:);
       b1 =                     fx(2,:);

       fy= b4*y3+ b3*y2 + b2*y + b1;

       c4 = -1./6*fy(1) + 0.5*fy(2) - 0.5*fy(3) + 1./6*fy(4);
       c3 =  0.5* fy(1) -     fy(2) + 0.5*fy(3);
       c2 = -1./3*fy(1) - 0.5*fy(2) +     fy(3) - 1./6*fy(4);
       c1 =                   fy(2);

       fi= c4*z3+ c3*z2 + c2*z + c1;


       end
!----------------------------------------------------------------------------------------------


# elif defined CRSPL_INTERPOLATION

       subroutine interp3(f,x,y,z,fi)
       !---------------------------------------------------------------------!
       !  Compute cartmull-rom spline 3d interpolant
       !  in: f(4,4,4) 
       !      x,y,z;  local coordinates (should be between 0 and 1) 
       !  out: fi 
       ! 
       !---------------------------------------------------------------------
       implicit none 
!      import/export
       real(kind=8),dimension(4,4,4),intent(in) :: f
       real(kind=8),                 intent(in) :: x,y,z
       real(kind=8),                 intent(out):: fi
!      local
       real(kind=8)  :: x2,x3,y2,y3,z2,z3
       real(kind=8),dimension(4,4) :: fx,a1,a2,a3,a4
       real(kind=8),dimension(4)   :: fy,b1,b2,b3,b4
       real(kind=8)                :: c1,c2,c3,c4


       x2 = x*x; x3 = x*x2 
       y2 = y*y; y3 = y*y2 
       z2 = z*z; z3 = z*z2 

       a4 = -0.5*f(1,:,:) + 1.5*f(2,:,:) - 1.5*f(3,:,:) + 0.5*f(4,:,:) 
       a3 =      f(1,:,:) - 2.5*f(2,:,:) + 2.0*f(3,:,:) - 0.5*f(4,:,:) 
       a2 = -0.5*f(1,:,:)                + 0.5*f(3,:,:)
       a1 =                     f(2,:,:);
       fx= a4*x3+ a3*x2 + a2*x + a1;

       b4 = -0.5*fx(1,:) + 1.5*fx(2,:) - 1.5*fx(3,:) + 0.5*fx(4,:);
       b3 =      fx(1,:) - 2.5*fx(2,:) + 2.0*fx(3,:) - 0.5*fx(4,:);
       b2 = -0.5*fx(1,:)               + 0.5*fx(3,:);
       b1 =                    fx(2,:);
       fy= b4*y3+ b3*y2 + b2*y + b1;

       c4 = -0.5*fy(1) + 1.5*fy(2) - 1.5*fy(3) + 0.5*fy(4)
       c3 =      fy(1) - 2.5*fy(2) + 2.0*fy(3) - 0.5*fy(4)
       c2 = -0.5*fy(1)             + 0.5*fy(3) 
       c1 =                  fy(2);
       fi= c4*z3+ c3*z2 + c2*z + c1;

       end
!----------------------------------------------------------------------------------------------

# elif defined WENO_INTERPOLATION

       subroutine interp3(f,x,y,z,fi)
       !---------------------------------------------------------------------!
       !  Compute 4th order weno 3d interpolant
       !  in: f(4,4,4) 
       !      x,y,z;  local coordinates (should be between 0 and 1) 
       !  out: fi 
       ! 
       !---------------------------------------------------------------------
       implicit none 
!      import/export
       real(kind=8),dimension(4,4,4),intent(in) :: f
       real(kind=8),                 intent(in) :: x,y,z
       real(kind=8),                 intent(out):: fi
!      local
       real(kind=8),dimension(4,4) :: fx
       real(kind=8),dimension(4)   :: fy
       integer(kind=4)             :: j,k


       do j = 1,4
         do k = 1,4
           call weno(f(1,j,k),f(2,j,k),f(3,j,k),f(4,j,k),x,fx(j,k));
         enddo
       enddo

       do k = 1,4
         call weno(fx(1,k),fx(2,k),fx(3,k),fx(4,k),y,fy(k));
       enddo

       call  weno(fy(1),fy(2),fy(3),fy(4),z,fi);

       end
!----------------------------------------------------------------------------------------------

       subroutine weno(f1,f2,f3,f4,x,fi)
       !---------------------------------------------------------------------!
       !  Compute 1-dimensional 4th order weno interpolant
       !  in: f1,f2,f3,f4
       !      x  local coordinate
       !  out: fi 
       ! 
       !---------------------------------------------------------------------
       implicit none 
!      import/export
       real(kind=8),intent(in) :: f1,f2,f3,f4
       real(kind=8),intent(in) :: x
       real(kind=8),intent(out):: fi
!      local
       integer(kind=4)     :: j,k
       real(kind=8)        :: w1,w2,p1,p2

        call wenoPara(f1,f2,f3,  x,w1,p1);
        call wenoPara(f4,f3,f2,1-x,w2,p2);

       fi = (w1*p1 + w2*p2)/(w1+w2);

       end

!----------------------------------------------------------------------------------------------
       subroutine wenoPara(f1,f2,f3,x,w,p)
       !---------------------------------------------------------------------!
       !  Compute 1-dimensional 4th order weno interpolant
       !  in: f1,f2,f3,f4
       !      x  local coordinate
       !  out: fi 
       ! 
       !---------------------------------------------------------------------
       implicit none 
!      import/export
       real(kind=8),intent(in) :: f1,f2,f3
       real(kind=8),intent(in) :: x
       real(kind=8),intent(out):: w,p
!      local
       real(kind=8)    :: d,dd,S

       d = (f3-f1)*0.5;                  !! 1st deriv at 0
       dd= (f1-2*f2+f3);                 !! 2nd derivative
       S = d*(d+dd) + 4/3*dd*dd;         !! Smoothness over x=[0,1]
       w = (2-x)/(1e-6+S*S);             !! relative weigth
       p = f2 + x*(d + 0.5*x*dd);        !! value of parabola at x


       end

#endif /* Choice of CUBIC, CRSPL, WENO */


#else /* ORIGINAL_VERSION */

!---------------------------------------------------------------------!
!  Compute displacement given u,v,w at particule position 
!  with linear interpolation in space and time
!---------------------------------------------------------------------

! 16/01/26:
! Modified sizes of u,v,w,dz,pm,pn to add ghost points (ng*2)



       subroutine advance_3d(px,py,pz,u,v,w,itim,fct,pm,pn,
     &          dz,dt,i0,j0,k0,nx,ny,nz,ng,np,dpx,dpy,dpz)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny,nz
       integer(kind=4)                    ,intent(in)   :: ng
       integer(kind=4)                    ,intent(in)   :: np
       integer(kind=4)                    ,intent(in)   :: i0,j0,k0
       integer(kind=4) ,dimension(0:1)    ,intent(in)   :: itim
       real(kind=8) ,dimension(nx-1,ny,nz,2),intent(in) :: u
       real(kind=8) ,dimension(nx,ny-1,nz,2),intent(in) :: v
       real(kind=8) ,dimension(nx,ny,nz+1,2),intent(in) :: w
       real(kind=8) ,dimension(nx,ny,nz,2),intent(in)   :: dz
       real(kind=8) ,dimension(nx,ny),intent(in)        :: pm,pn
       real(kind=8)                       ,intent(in)   :: px,py,pz
       real(kind=8)                       ,intent(out)  :: dpx,dpy,dpz
       real(kind=8)                       ,intent(in)   :: dt,fct
!      local
       integer(kind=8)                  :: i,j,k
       integer(kind=8)                  :: i_u,j_v,k_w
       real(kind=8) ,dimension(2,2,2,2) :: wt4
       real(kind=8) ,dimension(2,2,2)   :: wt3
       real(kind=8) ,dimension(2,2)     :: wt2
       real(kind=8)                     :: fcx,fcy,fcz,fctl
       real(kind=8)                     :: fcx_u,fcy_v,fcz_w
       real(kind=8)                     :: pu,pv,pw,pdz,ppm,ppn

!f2py intent(in)   u,v,w
!f2py intent(in)   px,py,pz
!f2py intent(out)  dpx,dpy,dpz
!f2py intent(in)   pm,pn,dt,fct
!f2py intent(in)   dz
!f2py intent(in)   nx,ny,nz
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   itim
!f2py intent(in)   ng,np


           !---------------------------------------
           ! 1. Linear interpolation in space and time
           !---------------------------------------        

# ifdef STAGGERED 
           i = floor(px+1+0.5+ng)-i0
           j = floor(py+1+0.5+ng)-j0
           k = max(1,min(floor(pz+1-0.5)-k0,nz-1))

           i_u = floor(px+1+ng)-i0
           j_v = floor(py+1+ng)-j0
           k_w = max(1,min(floor(pz+1)-k0,nz))

           fcx = px+1+0.5+ng - i - i0;
           fcy = py+1+0.5+ng - j - j0;
           fcz = pz+1-0.5 - k - k0;

           fcx_u = px+1+ng - i_u - i0;
           fcy_v = py+1+ng - j_v - j0;
           fcz_w = pz+1 - k_w - k0;

# else
           i = floor(px+1+ng)-i0
           j = floor(py+1+ng)-j0
           k = max(1,min(floor(pz+1-0.5)-k0,nz-1))

           i_u = floor(px+1-0.5+ng)-i0
           j_v = floor(py+1-0.5+ng)-j0
           k_w = max(1,min(floor(pz+1)-k0,nz))

           fcx = px+1+ng - i - i0;
           fcy = py+1+ng - j - j0;
           fcz = pz+1-0.5 - k - k0;

           fcx_u = px+1-0.5+ng - i_u - i0;
           fcy_v = py+1-0.5+ng - j_v - j0;
           fcz_w = pz+1 - k_w - k0;
# endif

           fctl = fct
           if (itim(0).eq.1)  fctl = 1-fct

           !---------------------------------------
           ! Compute velocities and level depth at particle position
           !---------------------------------------       

           if (i_u+1 .gt. nx-1) then
               write(*,*) 'px,nx-1', px,nx-1
               write(*,*) 'i_u,i_u+1', i_u,i_u+1
               write(*,*) 'fcx_u', fcx_u
           endif

           CALL linear_4d(fcx_u,fcy,fcz,fctl,wt4)
           pu   = sum(u(i_u:i_u+1,j:j+1,k:k+1,:)*wt4)

           !write(*,*) 'pu', pu
           !write(*,*) 'px', px
           !write(*,*) 'i_u,j', i_u,j
           !write(*,*) 'fcx_u,fcy', fcx_u,fcy

           CALL linear_4d(fcx,fcy_v,fcz,fctl,wt4)
           pv   = sum(v(i:i+1,j_v:j_v+1,k:k+1,:)*wt4)

           CALL linear_4d(fcx,fcy,fcz_w,fctl,wt4)
           pw   = sum(w(i:i+1,j:j+1,k_w:k_w+1,:)*wt4)


           CALL linear_4d(fcx,fcy,fcz,fctl,wt4)
           pdz   = sum(dz(i:i+1,j:j+1,k:k+1,:)*wt4)


           if (pdz.gt.0) pdz=1./pdz

           CALL linear_2d(fcx,fcy,wt2)
           ppm    = sum(pm(i:i+1,j:j+1)*wt2)
           ppn    = sum(pn(i:i+1,j:j+1)*wt2)

           !---------------------------------------
           ! Update position
           !---------------------------------------        

           dpx = dt*pu*ppm
           dpy = dt*pv*ppn
           dpz = dt*pw*pdz


       end

#endif



!---------------------------------------------------------------------!
!  Compute 2d displacement given u,v at particule position 
!  with linear interpolation in space and time
!---------------------------------------------------------------------



       subroutine advance_2d(px,py,u,v,itim,fct,pm,pn,
     &          dt,i0,j0,nx,ny,ng,np,dpx,dpy)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny
       integer(kind=4)                    ,intent(in)   :: ng
       integer(kind=4)                    ,intent(in)   :: np
       integer(kind=4)                    ,intent(in)   :: i0,j0
       integer(kind=4) ,dimension(0:1)    ,intent(in)   :: itim
       real(kind=8) ,dimension(nx-1,ny,2),intent(in) :: u
       real(kind=8) ,dimension(nx,ny-1,2),intent(in) :: v
       real(kind=8) ,dimension(nx,ny),intent(in)        :: pm,pn
       real(kind=8)                       ,intent(in)   :: px,py
       real(kind=8)                       ,intent(out)  :: dpx,dpy
       real(kind=8)                       ,intent(in)   :: dt,fct
!      local
       integer(kind=8)                  :: i,j
       integer(kind=8)                  :: i_u,j_v
       real(kind=8) ,dimension(2,2,2)   :: wt3
       real(kind=8) ,dimension(2,2)     :: wt2
       real(kind=8)                     :: fcx,fcy,fctl
       real(kind=8)                     :: fcx_u,fcy_v
       real(kind=8)                     :: pu,pv,ppm,ppn

!f2py intent(in)   u,v
!f2py intent(in)   px,py
!f2py intent(out)  dpx,dpy
!f2py intent(in)   pm,pn,dt,fct
!f2py intent(in)   nx,ny
!f2py intent(in)   i0,j0
!f2py intent(in)   itim
!f2py intent(in)   ng,np


           !---------------------------------------
           ! 1. Linear interpolation in space and time
           !---------------------------------------        

# ifdef STAGGERED 
           i = floor(px+1+0.5+ng)-i0
           j = floor(py+1+0.5+ng)-j0

           i_u = floor(px+1+ng)-i0
           j_v = floor(py+1+ng)-j0

           fcx = px+1+0.5+ng - i - i0;
           fcy = py+1+0.5+ng - j - j0;

           fcx_u = px+1+ng - i_u - i0;
           fcy_v = py+1+ng - j_v - j0;

# else
           i = floor(px+1+ng)-i0
           j = floor(py+1+ng)-j0

           i_u = floor(px+1-0.5+ng)-i0
           j_v = floor(py+1-0.5+ng)-j0

           fcx = px+1+ng - i - i0;
           fcy = py+1+ng - j - j0;

           fcx_u = px+1-0.5+ng - i_u - i0;
           fcy_v = py+1-0.5+ng - j_v - j0;

# endif

           fctl = fct
           if (itim(0).eq.1)  fctl = 1-fct

           !---------------------------------------
           ! Compute velocities and level depth at particle position
           !---------------------------------------       
 

           CALL linear_3d(fcx_u,fcy,fctl,wt3)
           pu   = sum(u(i_u:i_u+1,j:j+1,:)*wt3)

           CALL linear_3d(fcx,fcy_v,fctl,wt3)
           pv   = sum(v(i:i+1,j_v:j_v+1,:)*wt3)

           CALL linear_2d(fcx,fcy,wt2)
           ppm    = sum(pm(i:i+1,j:j+1)*wt2)
           ppn    = sum(pn(i:i+1,j:j+1)*wt2)

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
!      import/export
       real(kind=8) ,intent(in)                    :: fcx,fcy,fcz,fct
       real(kind=8) ,dimension(2,2,2,2),intent(out) :: wt

!f2py intent(in)   fcx,fcy,fcz,fct
!f2py intent(out) wt


           wt(1,1,1,1) = (1-fcz)*(1-fcy)*(1-fcx)*(1-fct);
           wt(1,1,2,1) =    fcz *(1-fcy)*(1-fcx)*(1-fct);
           wt(1,2,1,1) = (1-fcz)*   fcy *(1-fcx)*(1-fct);
           wt(1,2,2,1) =    fcz *   fcy *(1-fcx)*(1-fct);
           wt(2,1,1,1) = (1-fcz)*(1-fcy)*   fcx *(1-fct);
           wt(2,1,2,1) =    fcz *(1-fcy)*   fcx *(1-fct);
           wt(2,2,1,1) = (1-fcz)*   fcy *   fcx *(1-fct);
           wt(2,2,2,1) =    fcz *   fcy *   fcx *(1-fct);

           wt(1,1,1,2) = (1-fcz)*(1-fcy)*(1-fcx)*   fct ;
           wt(1,1,2,2) =    fcz *(1-fcy)*(1-fcx)*   fct ;
           wt(1,2,1,2) = (1-fcz)*   fcy *(1-fcx)*   fct ;
           wt(1,2,2,2) =    fcz *   fcy *(1-fcx)*   fct ;
           wt(2,1,1,2) = (1-fcz)*(1-fcy)*   fcx *   fct ;
           wt(2,1,2,2) =    fcz *(1-fcy)*   fcx *   fct ;
           wt(2,2,1,2) = (1-fcz)*   fcy *   fcx *   fct ;
           wt(2,2,2,2) =    fcz *   fcy *   fcx *   fct ;

       end



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get interpolation matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       subroutine linear_3d(fcx,fcy,fcz,wt)
       implicit none 
!      import/export
       real(kind=8) ,intent(in)                    :: fcx,fcy,fcz
       real(kind=8) ,dimension(2,2,2),intent(out) :: wt

!f2py intent(in)  fcx,fcy,fcz
!f2py intent(out) wt


           wt(1,1,1) = (1-fcz)*(1-fcy)*(1-fcx);
           wt(1,1,2) =    fcz *(1-fcy)*(1-fcx);
           wt(1,2,1) = (1-fcz)*   fcy *(1-fcx);
           wt(1,2,2) =    fcz *   fcy *(1-fcx);
           wt(2,1,1) = (1-fcz)*(1-fcy)*   fcx ;
           wt(2,1,2) =    fcz *(1-fcy)*   fcx ;
           wt(2,2,1) = (1-fcz)*   fcy *   fcx ;
           wt(2,2,2) =    fcz *   fcy *   fcx ;


       end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Get interpolation matrix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       subroutine linear_2d(fcx,fcy,wt)
       implicit none 
!      import/export
       real(kind=8) ,intent(in)                    :: fcx,fcy
       real(kind=8) ,dimension(2,2),intent(out) :: wt

!f2py intent(in)  fcx,fcy
!f2py intent(out) wt


           !wt(1,1) = (1-fcy)*(1-fcx);
           !wt(1,2) = (1-fcy)*   fcx;
           !wt(2,1) =    fcy *(1-fcx);
           !wt(2,2) =    fcy *   fcx;

           wt(1,1) = (1-fcy)*(1-fcx);
           wt(1,2) =    fcy *(1-fcx);
           wt(2,1) = (1-fcy)*   fcx;
           wt(2,2) =    fcy *   fcx;

       end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Interpolate T,S at each particle position (same than interp_3d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_3d_ts(pvar1,pvar2,px,py,pz,
     &          var1,var2,ng,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny,nz
       integer(kind=4)                    ,intent(in)   :: np,ng
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0,k0
       real(kind=8)   ,dimension(nx,ny,nz),intent(in)   :: var1,var2
       real(kind=8)   ,dimension(np)      ,intent(out):: pvar1,pvar2
       real(kind=8)   ,dimension(np)      ,intent(in):: px,py,pz
!      local
       integer(kind=8)                  :: ip,jp,kp,i,j,k
       real(kind=8)   ,dimension(2,2,2) :: wt3
       real(kind=8)                     :: fcx,fcy,fcz

!f2py intent(out)   pvar1,pvar2
!f2py intent(in)   var1,var2
!f2py intent(in)   px,py,pz
!f2py intent(in)   npmx
!f2py intent(in)   nx,ny,nz
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   ng,np

       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           i = max(1,min(floor(px(ip)+1+0.5+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+0.5+ng)-j0,ny-1))
           k = max(1,min(floor(pz(ip)+1-0.5)-k0,nz-1))

           fcx = px(ip)+1+0.5+ng - i - i0;
           fcy = py(ip)+1+0.5+ng - j - j0;
           fcz = pz(ip)+1-0.5 - k - k0;

           CALL linear_3d(fcx,fcy,fcz,wt3)
           pvar1(ip)   = sum(var1(i:i+1,j:j+1,k:k+1)*wt3)
           pvar2(ip)   = sum(var2(i:i+1,j:j+1,k:k+1)*wt3)

         endif
       enddo
!
       end    


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate a 3D variable at each particle position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_3d(pvar1,px,py,pz,var1,ng,npmx,i0,j0,k0,
     &          nx,ny,nz,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny,nz
       integer(kind=4)                    ,intent(in)   :: np,ng
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0,k0
       real(kind=8)   ,dimension(nx,ny,nz),intent(in)   :: var1
       real(kind=8)   ,dimension(np)      ,intent(out):: pvar1
       real(kind=8)   ,dimension(np)      ,intent(in):: px,py,pz
!      local
       integer(kind=8)                  :: ip,jp,kp,i,j,k
       real(kind=8)   ,dimension(2,2,2) :: wt3
       real(kind=8)                     :: fcx,fcy,fcz

!f2py intent(out)   pvar1
!f2py intent(in)   var1
!f2py intent(in)   px,py,pz
!f2py intent(in)   npmx
!f2py intent(in)   nx,ny,nz
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   ng,np

       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then


           i = max(1,min(floor(px(ip)+1+0.5+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+0.5+ng)-j0,ny-1))
           k = max(1,min(floor(pz(ip)+1-0.5)-k0,nz-1))

           fcx = px(ip)+1+0.5+ng - i - i0;
           fcy = py(ip)+1+0.5+ng - j - j0;
           fcz = pz(ip)+1-0.5 - k - k0;

           CALL linear_3d(fcx,fcy,fcz,wt3)
           pvar1(ip)   = sum(var1(i:i+1,j:j+1,k:k+1)*wt3)

         endif
       enddo
!
       end    


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate a 3D variable at each particle position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_3d_u(pvar1,px,py,pz,var1,ng,npmx,i0,j0,k0,
     &          nx,ny,nz,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny,nz
       integer(kind=4)                    ,intent(in)   :: np,ng
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0,k0
       real(kind=8)   ,dimension(nx,ny,nz),intent(in)   :: var1
       real(kind=8)   ,dimension(np)      ,intent(out):: pvar1
       real(kind=8)   ,dimension(np)      ,intent(in):: px,py,pz
!      local
       integer(kind=8)                  :: ip,jp,kp,i,j,k
       real(kind=8)   ,dimension(2,2,2) :: wt3
       real(kind=8)                     :: fcx,fcy,fcz

!f2py intent(out)   pvar1
!f2py intent(in)   var1
!f2py intent(in)   px,py,pz
!f2py intent(in)   npmx
!f2py intent(in)   nx,ny,nz
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   ng,np

       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           i = max(1,min(floor(px(ip)+1+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+0.5+ng)-j0,ny-1))
           k = max(1,min(floor(pz(ip)+1-0.5)-k0,nz-1))

           fcx = px(ip)+1+ng - i - i0;
           fcy = py(ip)+1+0.5+ng - j - j0;
           fcz = pz(ip)+1-0.5 - k - k0;

           CALL linear_3d(fcx,fcy,fcz,wt3)
           pvar1(ip)   = sum(var1(i:i+1,j:j+1,k:k+1)*wt3)


         endif
       enddo
!
       end    


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate a 3D variable at each particle position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_3d_v(pvar1,px,py,pz,var1,ng,npmx,i0,j0,k0,
     &          nx,ny,nz,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny,nz
       integer(kind=4)                    ,intent(in)   :: np,ng
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0,k0
       real(kind=8)   ,dimension(nx,ny,nz),intent(in)   :: var1
       real(kind=8)   ,dimension(np)      ,intent(out):: pvar1
       real(kind=8)   ,dimension(np)      ,intent(in):: px,py,pz
!      local
       integer(kind=8)                  :: ip,jp,kp,i,j,k
       real(kind=8)   ,dimension(2,2,2) :: wt3
       real(kind=8)                     :: fcx,fcy,fcz

!f2py intent(out)   pvar1
!f2py intent(in)   var1
!f2py intent(in)   px,py,pz
!f2py intent(in)   npmx
!f2py intent(in)   nx,ny,nz
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   ng,np

       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           i = max(1,min(floor(px(ip)+1+0.5+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+ng)-j0,ny-1))
           k = max(1,min(floor(pz(ip)+1-0.5)-k0,nz-1))

           fcx = px(ip)+1+0.5+ng - i - i0;
           fcy = py(ip)+1+ng - j - j0;
           fcz = pz(ip)+1-0.5 - k - k0;

           CALL linear_3d(fcx,fcy,fcz,wt3)
           pvar1(ip)   = sum(var1(i:i+1,j:j+1,k:k+1)*wt3)

         endif
       enddo
!
       end    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate a 3D variable at each particle position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_3d_psi(pvar1,px,py,pz,var1,ng,npmx,i0,j0,k0,
     &          nx,ny,nz,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny,nz
       integer(kind=4)                    ,intent(in)   :: np,ng
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0,k0
       real(kind=8)   ,dimension(nx,ny,nz),intent(in)   :: var1
       real(kind=8)   ,dimension(np)      ,intent(out):: pvar1
       real(kind=8)   ,dimension(np)      ,intent(in):: px,py,pz
!      local
       integer(kind=8)                  :: ip,jp,kp,i,j,k
       real(kind=8)   ,dimension(2,2,2) :: wt3
       real(kind=8)                     :: fcx,fcy,fcz

!f2py intent(out)   pvar1
!f2py intent(in)   var1
!f2py intent(in)   px,py,pz
!f2py intent(in)   npmx
!f2py intent(in)   nx,ny,nz
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   ng,np

       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           i = max(1,min(floor(px(ip)+1+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+ng)-j0,ny-1))
           k = max(1,min(floor(pz(ip)+1-0.5)-k0,nz-1))

           fcx = px(ip)+1+ng - i - i0;
           fcy = py(ip)+1+ng - j - j0;
           fcz = pz(ip)+1-0.5 - k - k0;

           CALL linear_3d(fcx,fcy,fcz,wt3)
           pvar1(ip)   = sum(var1(i:i+1,j:j+1,k:k+1)*wt3)

         endif
       enddo
!
       end    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate a 3D variable at each particle position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_3d_w(pvar1,px,py,pz,var1,ng,npmx,i0,j0,k0,
     &          nx,ny,nz,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny,nz
       integer(kind=4)                    ,intent(in)   :: np,ng
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0,k0
       real(kind=8)   ,dimension(nx,ny,nz),intent(in)   :: var1
       real(kind=8)   ,dimension(np)      ,intent(out):: pvar1
       real(kind=8)   ,dimension(np)      ,intent(in):: px,py,pz
!      local
       integer(kind=8)                  :: ip,jp,kp,i,j,k
       real(kind=8)   ,dimension(2,2,2) :: wt3
       real(kind=8)                     :: fcx,fcy,fcz

!f2py intent(out)   pvar1
!f2py intent(in)   var1
!f2py intent(in)   px,py,pz
!f2py intent(in)   npmx
!f2py intent(in)   nx,ny,nz
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   ng,np

       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           i = max(1,min(floor(px(ip)+1+0.5+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+0.5+ng)-j0,ny-1))
           k = max(1,min(floor(pz(ip)+1)-k0,nz-1))

           fcx = px(ip)+1+0.5+ng - i - i0;
           fcy = py(ip)+1+0.5+ng - j - j0;
           fcz = pz(ip)+1 - k - k0;

           CALL linear_3d(fcx,fcy,fcz,wt3)
           pvar1(ip)   = sum(var1(i:i+1,j:j+1,k:k+1)*wt3)

         endif
       enddo
!
       end    


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate a 3D variable at each particle position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_3d_psiw(pvar1,px,py,pz,var1,ng,npmx,i0,j0,k0,
     &          nx,ny,nz,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny,nz
       integer(kind=4)                    ,intent(in)   :: np,ng
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0,k0
       real(kind=8)   ,dimension(nx,ny,nz),intent(in)   :: var1
       real(kind=8)   ,dimension(np)      ,intent(out):: pvar1
       real(kind=8)   ,dimension(np)      ,intent(in):: px,py,pz
!      local
       integer(kind=8)                  :: ip,jp,kp,i,j,k
       real(kind=8)   ,dimension(2,2,2) :: wt3
       real(kind=8)                     :: fcx,fcy,fcz

!f2py intent(out)   pvar1
!f2py intent(in)   var1
!f2py intent(in)   px,py,pz
!f2py intent(in)   npmx
!f2py intent(in)   nx,ny,nz
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   ng,np

       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           i = max(1,min(floor(px(ip)+1+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+ng)-j0,ny-1))
           k = max(1,min(floor(pz(ip)+1)-k0,nz-1))

           fcx = px(ip)+1+ng - i - i0;
           fcy = py(ip)+1+ng - j - j0;
           fcz = pz(ip)+1 - k - k0;

           CALL linear_3d(fcx,fcy,fcz,wt3)
           pvar1(ip)   = sum(var1(i:i+1,j:j+1,k:k+1)*wt3)

         endif
       enddo
!
       end    




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Interpolate T,S at each particle position (same than interp_2d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_2d_ts(pvar1,pvar2,px,py,
     &          var1,var2,ng,npmx,i0,j0,nx,ny,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny
       integer(kind=4)                    ,intent(in)   :: np,ng
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0
       real(kind=8)   ,dimension(nx,ny),intent(in)   :: var1,var2
       real(kind=8)   ,dimension(np)      ,intent(out):: pvar1,pvar2
       real(kind=8)   ,dimension(np)      ,intent(in):: px,py
!      local
       integer(kind=8)                  :: ip,jp,i,j
       real(kind=8)   ,dimension(2,2) :: wt2
       real(kind=8)                     :: fcx,fcy

!f2py intent(out)   pvar1,pvar2
!f2py intent(in)   var1,var2
!f2py intent(in)   px,py
!f2py intent(in)   npmx
!f2py intent(in)   nx,ny
!f2py intent(in)   i0,j0
!f2py intent(in)   ng,np

       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip))) ) then

           i = max(1,min(floor(px(ip)+1+0.5+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+0.5+ng)-j0,ny-1))

           fcx = px(ip)+1+0.5+ng - i - i0;
           fcy = py(ip)+1+0.5+ng - j - j0;

           CALL linear_2d(fcx,fcy,wt2)
           pvar1(ip)   = sum(var1(i:i+1,j:j+1)*wt2)
           pvar2(ip)   = sum(var2(i:i+1,j:j+1)*wt2)

         endif
       enddo
!
       end    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate u at each horizontal particle position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_2d_u(pvar1,px,py,var1,ng,npmx,i0,j0,nx,ny,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny
       integer(kind=4)                    ,intent(in)   :: ng,np
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0
       real(kind=8)   ,dimension(nx,ny),intent(in)   :: var1
       real(kind=8)   ,dimension(np)      ,intent(out):: pvar1
       real(kind=8)   ,dimension(np)      ,intent(in):: px,py
!      local
       integer(kind=8)                  :: ip,jp,i,j
       real(kind=8)   ,dimension(2,2) :: wt
       real(kind=8)                     :: fcx,fcy

!f2py intent(out)   pvar1
!f2py intent(in)   var1
!f2py intent(in)   px,py
!f2py intent(in)   npmx
!f2py intent(in)   nx,ny
!f2py intent(in)   i0,j0
!f2py intent(in)   ng,np


       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip))) ) then

           i = max(1,min(floor(px(ip)+1+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+0.5+ng)-j0,ny-1))

           fcx = px(ip)+1+ng - i - i0;
           fcy = py(ip)+1+0.5+ng - j - j0;

           CALL linear_2d(fcx,fcy,wt)
           pvar1(ip)   = sum(var1(i:i+1,j:j+1)*wt)

         endif
       enddo
!
       end   





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate v at each horizontal particle position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_2d_v(pvar1,px,py,var1,ng,npmx,i0,j0,nx,ny,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny
       integer(kind=4)                    ,intent(in)   :: ng,np
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0
       real(kind=8)   ,dimension(nx,ny),intent(in)   :: var1
       real(kind=8)   ,dimension(np)      ,intent(out):: pvar1
       real(kind=8)   ,dimension(np)      ,intent(in):: px,py
!      local
       integer(kind=8)                  :: ip,jp,i,j
       real(kind=8)   ,dimension(2,2) :: wt
       real(kind=8)                     :: fcx,fcy

!f2py intent(out)   pvar1
!f2py intent(in)   var1
!f2py intent(in)   px,py
!f2py intent(in)   npmx
!f2py intent(in)   nx,ny
!f2py intent(in)   i0,j0
!f2py intent(in)   ng,np


       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip))) ) then

           i = max(1,min(floor(px(ip)+1+0.5+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+ng)-j0,ny-1))

           fcx = px(ip)+1+0.5+ng - i - i0;
           fcy = py(ip)+1+ng - j - j0;

           CALL linear_2d(fcx,fcy,wt)
           pvar1(ip)   = sum(var1(i:i+1,j:j+1)*wt)

         endif
       enddo
!
       end   





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate a 2D variable at each horizontal particle position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_2d(pvar1,px,py,var1,ng,npmx,i0,j0,nx,ny,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny
       integer(kind=4)                    ,intent(in)   :: ng,np
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0
       real(kind=8)   ,dimension(nx,ny),intent(in)   :: var1
       real(kind=8)   ,dimension(np)      ,intent(out):: pvar1
       real(kind=8)   ,dimension(np)      ,intent(in):: px,py
!      local
       integer(kind=8)                  :: ip,jp,i,j
       real(kind=8)   ,dimension(2,2) :: wt
       real(kind=8)                     :: fcx,fcy

!f2py intent(out)   pvar1
!f2py intent(in)   var1
!f2py intent(in)   px,py
!f2py intent(in)   npmx
!f2py intent(in)   nx,ny
!f2py intent(in)   i0,j0
!f2py intent(in)   ng,np


       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip))) ) then

           i = max(1,min(floor(px(ip)+1+0.5+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+0.5+ng)-j0,ny-1))

           fcx = px(ip)+1+0.5+ng - i - i0;
           fcy = py(ip)+1+0.5+ng - j - j0;

           CALL linear_2d(fcx,fcy,wt)
           pvar1(ip)   = sum(var1(i:i+1,j:j+1)*wt)

         endif
       enddo
!
       end   


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate a 2D variable at each horizontal particle position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_2d_psi(pvar1,px,py,var1,ng,npmx,i0,j0,nx,ny,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny
       integer(kind=4)                    ,intent(in)   :: ng,np
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0
       real(kind=8)   ,dimension(nx,ny),intent(in)   :: var1
       real(kind=8)   ,dimension(np)      ,intent(out):: pvar1
       real(kind=8)   ,dimension(np)      ,intent(in):: px,py
!      local
       integer(kind=8)                  :: ip,jp,i,j
       real(kind=8)   ,dimension(2,2) :: wt
       real(kind=8)                     :: fcx,fcy

!f2py intent(out)   pvar1
!f2py intent(in)   var1
!f2py intent(in)   px,py
!f2py intent(in)   npmx
!f2py intent(in)   nx,ny
!f2py intent(in)   i0,j0
!f2py intent(in)   ng,np


       do ip = 1,np
         if (.not.(isnan(px(ip)*py(ip))) ) then

           i = max(1,min(floor(px(ip)+1+ng)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1+ng)-j0,ny-1))

           fcx = px(ip)+1+ng - i - i0;
           fcy = py(ip)+1+ng - j - j0;

           CALL linear_2d(fcx,fcy,wt)
           pvar1(ip)   = sum(var1(i:i+1,j:j+1)*wt)

         endif
       enddo
!
       end   


