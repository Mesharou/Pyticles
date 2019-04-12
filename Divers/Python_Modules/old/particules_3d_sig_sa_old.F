
!---------------------------------------------------------------------!! 
!! 
!! cpp particules_3d_sig_sa.F particules_3d_sig_sa.f90; f2py --f90flags="-extend_source -O1 " -DF2PY_REPORT_ON_ARRAY_COPY=1 -c -m particules_3d_sig_sa particules_3d_sig_sa.f90
!! 
!---------------------------------------------------------------------

!---------------------------------------------------------------------!
!  Forward Euler time-stepping
!---------------------------------------------------------------------

!

       subroutine timestep_FE(px,py,pz,u,v,w,fct,dfct,pm,pn,dz,dt,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny,nz
       integer(kind=4)                    ,intent(in)   :: np
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0,k0
       real(kind=8) ,dimension(nx-1,ny,nz,2),intent(in)   :: u
       real(kind=8) ,dimension(nx,ny-1,nz,2),intent(in)   :: v
       real(kind=8) ,dimension(nx,ny,nz+1,2),intent(in)   :: w
       real(kind=8) ,dimension(nx,ny,nz,2),intent(in)   :: dz
       real(kind=8) ,dimension(nx,ny),intent(in)   :: pm,pn
       real(kind=8)   ,dimension(np)      ,intent(inout):: px,py,pz
       real(kind=8)                       ,intent(in)   :: dt,fct,dfct
!      local
       real(kind=8)           :: dpx,dpy,dpz
       integer(kind=8)                  :: ip

!f2py intent(in)   u,v,w
!f2py intent(inout) px,py,pz
!f2py intent(in)   pm,pn,dt
!f2py intent(in)   dz
!f2py intent(in)   npmx
!f2py intent(in)   nx,ny,nz
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   np


       do ip = 1,min(npmx,np)
        
         if (.not.(isnan(px(ip))) ) then       

           CALL advance_3d_lin(px(ip),py(ip),pz(ip),u,v,w,fct,pm,pn,dz,dt,npmx,i0,j0,k0,nx,ny,nz,np,dpx,dpy,dpz)

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
!  Runge-Kutta2 time-stepping
!---------------------------------------------------------------------

!

       subroutine timestep_RK2(px,py,pz,u,v,w,fct,dfct,pm,pn,dz,dt,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny,nz
       integer(kind=4)                    ,intent(in)   :: np
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0,k0
       real(kind=8) ,dimension(nx-1,ny,nz,2),intent(in)   :: u
       real(kind=8) ,dimension(nx,ny-1,nz,2),intent(in)   :: v
       real(kind=8) ,dimension(nx,ny,nz+1,2),intent(in)   :: w
       real(kind=8) ,dimension(nx,ny,nz,2),intent(in)   :: dz
       real(kind=8) ,dimension(nx,ny),intent(in)   :: pm,pn
       real(kind=8)   ,dimension(np)      ,intent(inout):: px,py,pz
       real(kind=8)                       ,intent(in)   :: dt,fct,dfct
!      local
       real(kind=8)            :: dpx,dpy,dpz
       integer(kind=8)                  :: ip


!f2py intent(inout) px,py,pz
!f2py intent(in)   u,v,w
!f2py intent(in)   fct,dfct,dt
!f2py intent(in)   pm,pn,dz
!f2py intent(in)   npmx
!f2py intent(in)   nx,ny,nz
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   np

       do ip = 1,min(npmx,np)
        
         if (.not.(isnan(px(ip))) ) then       
           
           ! Forward Euler / predictor
           CALL advance_3d_lin(px(ip),py(ip),pz(ip),u,v,w,
     &          fct,pm,pn,dz,dt,npmx,i0,j0,k0,nx,ny,nz,np,dpx,dpy,dpz)

           ! midpoint rule / corrector
           CALL advance_3d_lin(px(ip)+0.5*dpx,py(ip)+0.5*dpy,pz(ip)+0.5*dpz,
     &          u,v,w,fct+0.5*dfct,pm,pn,dz,dt,npmx,i0,j0,k0,nx,ny,nz,np,dpx,dpy,dpz)

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
!  Runge-Kutta 4 time-stepping
!---------------------------------------------------------------------

!

       subroutine timestep_RK4(px,py,pz,u,v,w,fct,dfct,pm,pn,dz,dt,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny,nz
       integer(kind=4)                    ,intent(in)   :: np
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0,k0
       real(kind=8) ,dimension(nx-1,ny,nz,2),intent(in) :: u
       real(kind=8) ,dimension(nx,ny-1,nz,2),intent(in) :: v
       real(kind=8) ,dimension(nx,ny,nz+1,2),intent(in) :: w
       real(kind=8) ,dimension(nx,ny,nz,2),intent(in)   :: dz
       real(kind=8) ,dimension(nx,ny),intent(in)        :: pm,pn
       real(kind=8)   ,dimension(np)      ,intent(inout):: px,py,pz
       real(kind=8)                       ,intent(in)   :: dt,fct,dfct
!      local
       real(kind=8)  ,dimension(4)            :: dpx,dpy,dpz
       real(kind=8)                           :: coef
       integer(kind=8)                        :: ip

!f2py intent(inout) px,py,pz
!f2py intent(in)   u,v,w
!f2py intent(in)   fct,dfct,dt
!f2py intent(in)   pm,pn,dz
!f2py intent(in)   npmx
!f2py intent(in)   nx,ny,nz
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   np

       coef = 1/6.

       do ip = 1,min(npmx,np)
        
         if (.not.(isnan(px(ip))) ) then       

           ! Forward Euler / predictor
           CALL advance_3d_lin(px(ip),py(ip),pz(ip),u,v,w,
     &          fct,pm,pn,dz,dt,npmx,i0,j0,k0,nx,ny,nz,np,dpx(0),dpy(0),dpz(0))

           ! backward Euler / corrector
           CALL advance_3d_lin(px(ip)+0.5*dpx(0),py(ip)+0.5*dpy(0),pz(ip)+0.5*dpz(0),
     &          u,v,w,fct+0.5*dfct,pm,pn,dz,dt,npmx,i0,j0,k0,nx,ny,nz,np,dpx(1),dpy(1),dpz(1))

           ! midpoint rule / predictor
           CALL advance_3d_lin(px(ip)+0.5*dpx(1),py(ip)+0.5*dpy(1),pz(ip)+0.5*dpz(1),
     &          u,v,w,fct+0.5*dfct,pm,pn,dz,dt,npmx,i0,j0,k0,nx,ny,nz,np,dpx(2),dpy(2),dpz(2))

           ! Corrector
           CALL advance_3d_lin(px(ip)+dpx(2),py(ip)+dpy(2),pz(ip)+dpz(2),
     &          u,v,w,fct+0.5*dfct,pm,pn,dz,dt,npmx,i0,j0,k0,nx,ny,nz,np,dpx(3),dpy(3),dpz(3))


           !---------------------------------------
           ! Update position
           !---------------------------------------        

           px(ip) = px(ip) + coef * (dpx(0) + 2 * dpx(1) + 2 * dpx(2) + dpx(3))
           py(ip) = py(ip) + coef * (dpy(0) + 2 * dpy(1) + 2 * dpy(2) + dpy(3))
           pz(ip) = pz(ip) + coef * (dpz(0) + 2 * dpz(1) + 2 * dpz(2) + dpz(3))

         endif

       enddo


       end


!---------------------------------------------------------------------!
!  Compute displacement given u,v,w at particule position 
!  with linear interpolation in space and time
!---------------------------------------------------------------------
!
!

       subroutine advance_3d_lin(px,py,pz,u,v,w,fct,pm,pn,dz,dt,npmx,i0,j0,k0,nx,ny,nz,np,dpx,dpy,dpz)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny,nz
       integer(kind=4)                    ,intent(in)   :: np
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0,k0
       real(kind=8) ,dimension(nx-1,ny,nz,2),intent(in)   :: u
       real(kind=8) ,dimension(nx,ny-1,nz,2),intent(in)   :: v
       real(kind=8) ,dimension(nx,ny,nz+1,2),intent(in)   :: w
       real(kind=8) ,dimension(nx,ny,nz,2),intent(in)   :: dz
       real(kind=8) ,dimension(nx,ny),intent(in)   :: pm,pn
       real(kind=8)                       ,intent(in):: px,py,pz
       real(kind=8)                       ,intent(out):: dpx,dpy,dpz
       real(kind=8)                       ,intent(in)   :: dt,fct
!      local
       integer(kind=8)                  :: i,j,k
       integer(kind=8)                  :: i_u,j_v,k_w
       real(kind=8) ,dimension(2,2,2,2) :: wt4
       real(kind=8) ,dimension(2,2,2)   :: wt3
       real(kind=8) ,dimension(2,2)     :: wt2
       real(kind=8)                     :: fcx,fcy,fcz
       real(kind=8)                     :: fcx_u,fcy_v,fcz_w
       real(kind=8)                     :: pu,pv,pw,pdz,ppm,ppn


!f2py intent(in)   u,v,w
!f2py intent(in) px,py,pz
!f2py intent(out) dpx,dpy,dpz
!f2py intent(in)   pm,pn,dt
!f2py intent(in)   dz
!f2py intent(in)   npmx
!f2py intent(in)   nx,ny,nz
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   np


           i = max(1,min(floor(px+1)-i0,nx-1))
           j = max(1,min(floor(py+1)-j0,ny-1))
           k = max(1,min(floor(pz+1-0.5)-k0,nz-1))

           i_u = max(1,min(floor(px+1-0.5)-i0,nx-2))
           j_v = max(1,min(floor(py+1-0.5)-j0,ny-2))
           k_w = max(1,min(floor(pz+1)-k0,nz))

           !---------------------------------------
           ! 1. Linear interpolation in space and time
           !---------------------------------------        
 
           fcx = px+1 - i - i0;
           fcy = py+1 - j - j0;
           fcz = pz+1-0.5 - k - k0;

           fcx_u = px+0.5 - i - i0;
           fcy_v = py+0.5 - j - j0;
           fcz_w = pz+1 - k - k0;

           !---------------------------------------
           ! Compute velocities and level depth at particle position
           !---------------------------------------        

           CALL linear_4d(fcx_u,fcy,fcz,fct,wt4)
           pu   = sum(u(i:i+1,j:j+1,k:k+1,:)*wt4)

           !print *,'pu,dt',pu,i,fcx_u

           CALL linear_4d(fcx,fcy_v,fcz,fct,wt4)
           pv   = sum(v(i:i+1,j:j+1,k:k+1,:)*wt4)

           CALL linear_4d(fcx,fcy,fcz_w,fct,wt4)
           pw   = sum(w(i:i+1,j:j+1,k:k+1,:)*wt4)

           CALL linear_4d(fcx,fcy,fcz,fct,wt4)
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

!
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


           wt(1,1) = (1-fcy)*(1-fcx);
           wt(1,2) = (1-fcy)*   fcx;
           wt(2,1) =    fcy *(1-fcx);
           wt(2,2) =    fcy *   fcx;


       end


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  Interpolate T,S at each particle position (same than interp_3d)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_3d_ts(pvar1,pvar2,px,py,pz,var1,var2,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny,nz
       integer(kind=4)                    ,intent(in)   :: np
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
!f2py intent(in)   np

       do ip = 1,np
         if (.not.(isnan(px(ip))) ) then

           i = max(1,min(floor(px(ip)+1)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1)-j0,ny-1))
           k = max(1,min(floor(pz(ip)+1-0.5)-k0,nz-1))

           fcx = px(ip)+1 - i - i0;
           fcy = py(ip)+1 - j - j0;
           fcz = pz(ip)+0.5 - k - k0;

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


       subroutine interp_3d(pvar1,px,py,pz,var1,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny,nz
       integer(kind=4)                    ,intent(in)   :: np
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
!f2py intent(in)   np

       do ip = 1,np
         if (.not.(isnan(px(ip))) ) then

           i = max(1,min(floor(px(ip)+1)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1)-j0,ny-1))
           k = max(1,min(floor(pz(ip)+1-0.5)-k0,nz-1))

           fcx = px(ip)+1 - i - i0;
           fcy = py(ip)+1 - j - j0;
           fcz = pz(ip)+0.5 - k - k0;

           CALL linear_3d(fcx,fcy,fcz,wt3)
           pvar1(ip)   = sum(var1(i:i+1,j:j+1,k:k+1)*wt3)

         endif
       enddo
!
       end    


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Interpolate a 2D variable at each horizontal particle position
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       subroutine interp_2d(pvar1,px,py,var1,npmx,i0,j0,nx,ny,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny
       integer(kind=4)                    ,intent(in)   :: np
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
!f2py intent(in)   np


       do ip = 1,np
         if (.not.(isnan(px(ip))) ) then

           i = max(1,min(floor(px(ip)+1)-i0,nx-1))
           j = max(1,min(floor(py(ip)+1)-j0,ny-1))

           fcx = px(ip)+1 - i - i0;
           fcy = py(ip)+1 - j - j0;

           CALL linear_2d(fcx,fcy,wt)

           pvar1(ip)   = sum(var1(i:i+1,j:j+1)*wt)

         endif
       enddo
!
       end   




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Some soutines copied from the R_tools.F module
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# define CUBIC
# define DUKO_2001
# define INTERP_BELOW
# define INTERP_ABOVE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute z_r and z_w for NEW_S_COORD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# include "R_tools_fort_routines/zlevs.F"

# include "R_tools_fort_routines/zlevs_w.F"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Z interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# include "R_tools_fort_routines/sigma_to_z_intr_sfc.F"

# include "R_tools_fort_routines/sigma_to_z_intr_bot.F"








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







