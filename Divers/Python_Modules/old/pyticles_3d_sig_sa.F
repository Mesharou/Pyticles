!----------------------------------------------------------------------------------------------
! Fortran Routines for Pyticles
!
!----------------------------------------------------------------------------------------------
! 16/01/26:
!     - Add parameter 'ng' corresponding to number of ghost points in the horizontal grid
!----------------------------------------------------------------------------------------------


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Spatial interpolation routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# include "interp_3d_for_pyticles.F"



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Time-stepping routines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


       !---------------------------------------------------------------------!
       !  Forward Euler time-stepping
       !---------------------------------------------------------------------

       subroutine timestep_FE(px,py,pz,u,v,w,itim,fct,dfct,
     &          pm,pn,dz,dt,ng,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny,nz
       integer(kind=4)                    ,intent(in)   :: ng,np
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0,k0
       integer(kind=4) ,dimension(0:1)    ,intent(in)   :: itim
       real(kind=8) ,dimension(nx-1,ny,nz,2),intent(in)   :: u
       real(kind=8) ,dimension(nx,ny-1,nz,2),intent(in)   :: v
       real(kind=8) ,dimension(nx,ny,nz+1,2),intent(in)   :: w
       real(kind=8) ,dimension(nx,ny,nz,2),intent(in)   :: dz
       real(kind=8) ,dimension(nx,ny),intent(in)   :: pm,pn
       real(kind=8) ,dimension(np)        ,intent(inout):: px,py,pz
       real(kind=8)                       ,intent(in)   :: dt,fct,dfct
!      local
       real(kind=8)           :: dpx,dpy,dpz
       integer(kind=8)                  :: ip

!f2py intent(inout) px,py,pz
!f2py intent(in)   u,v,w
!f2py intent(in)   itim
!f2py intent(in)   fct,dfct,dt
!f2py intent(in)   pm,pn
!f2py intent(in)   dz
!f2py intent(in)   ng,npmx
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   nx,ny,nz
!f2py intent(in)   np

       do ip = 1,min(npmx,np)
        
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,itim,fct,pm,pn
     &         ,dz,dt,i0,j0,k0,nx,ny,nz,ng,np,dpx,dpy,dpz)

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

       subroutine timestep_RK2(px,py,pz,u,v,w,itim,fct,dfct,
     &          pm,pn,dz,dt,ng,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny,nz
       integer(kind=4)                    ,intent(in)   :: ng,np
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0,k0
       integer(kind=4) ,dimension(0:1)    ,intent(in)   :: itim
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
!f2py intent(in)   itim
!f2py intent(in)   fct,dfct,dt
!f2py intent(in)   pm,pn
!f2py intent(in)   dz
!f2py intent(in)   ng,npmx
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   nx,ny,nz
!f2py intent(in)   np


       do ip = 1,min(npmx,np)
        
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then       

           ! Forward Euler / predictor
           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,itim,
     &          fct,pm,pn,dz,dt,i0,j0,k0,nx,ny,nz,ng,np,dpx,dpy,dpz)

           if (.not.(isnan(dpx*dpy*dpz)) ) then

             ! midpoint rule / corrector
           CALL advance_3d(px(ip)+0.5*dpx,py(ip)+0.5*dpy,pz(ip)+0.5*dpz,
     &            u,v,w,itim,fct+0.5*dfct,pm,pn,dz,dt,i0,j0,k0,
     &            nx,ny,nz,ng,np,dpx,dpy,dpz)

           endif

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

       subroutine timestep_RK4(px,py,pz,u,v,w,itim,fct,dfct,pm,pn,
     &          dz,dt,ng,npmx,i0,j0,k0,nx,ny,nz,np,dpxi,dpyi,dpzi)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny,nz
       integer(kind=4)                    ,intent(in)   :: ng,np
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0,k0
       integer(kind=4) ,dimension(0:1)    ,intent(in)   :: itim
       real(kind=8) ,dimension(nx-1,ny,nz,2),intent(in) :: u
       real(kind=8) ,dimension(nx,ny-1,nz,2),intent(in) :: v
       real(kind=8) ,dimension(nx,ny,nz+1,2),intent(in) :: w
       real(kind=8) ,dimension(nx,ny,nz,2),intent(in)   :: dz
       real(kind=8) ,dimension(nx,ny),intent(in)        :: pm,pn
       real(kind=8) ,dimension(np)        ,intent(inout):: px,py,pz
       real(kind=8) ,dimension(np)        ,intent(out)  :: dpxi,dpyi,dpzi
       real(kind=8)                       ,intent(in)   :: dt,fct,dfct
!      local
       real(kind=8) ,dimension(0:3)            :: dpx,dpy,dpz
       real(kind=8)                           :: coef
       integer(kind=8)                        :: ip

!f2py intent(inout) px,py,pz
!f2py intent(in)   u,v,w
!f2py intent(in)   itim
!f2py intent(in)   fct,dfct,dt
!f2py intent(in)   pm,pn
!f2py intent(in)   dz
!f2py intent(in)   ng,npmx
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   nx,ny,nz
!f2py intent(in)   np
!f2py intent(out)  dpxi,dpyi,dpzi

       coef = 1./6.

       do ip = 1,min(npmx,np)

         if (.not.(isnan(px(ip)*py(ip)*pz(ip)))) then 

           ! Forward Euler / predictor
           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,
     &          itim,fct,pm,pn,dz,dt,i0,j0,k0,
     &          nx,ny,nz,ng,np,dpx(0),dpy(0),dpz(0))

           if (.not.(isnan(dpx(0)*dpy(0)*dpz(0)))) then 

           ! backward Euler / corrector
           CALL advance_3d(px(ip)+0.5*dpx(0),py(ip)+0.5*dpy(0),
     &          pz(ip)+0.5*dpz(0),
     &          u,v,w,itim,fct+0.5*dfct,pm,pn,dz,dt,i0,j0,k0,
     &          nx,ny,nz,ng,np,dpx(1),dpy(1),dpz(1))

            if (.not.(isnan(dpx(1)*dpy(1)*dpz(1)))) then 

           ! midpoint rule / predictor
           CALL advance_3d(px(ip)+0.5*dpx(1),py(ip)+0.5*dpy(1),
     &          pz(ip)+0.5*dpz(1),
     &          u,v,w,itim,fct+0.5*dfct,pm,pn,dz,dt,i0,j0,k0,
     &          nx,ny,nz,ng,np,dpx(2),dpy(2),dpz(2))

             if (.not.(isnan(dpx(2)*dpy(2)*dpz(2)))) then 

           ! Corrector
           CALL advance_3d(px(ip)+dpx(2),py(ip)+dpy(2),
     &          pz(ip)+dpz(2),
     &          u,v,w,itim,fct+1.*dfct,pm,pn,dz,dt,i0,j0,k0,
     &          nx,ny,nz,ng,np,dpx(3),dpy(3),dpz(3))

             endif
            endif
           endif


           !---------------------------------------
           ! Update position
           !---------------------------------------


           dpxi(ip) = coef * (dpx(0) + 2 * dpx(1) + 2 * dpx(2) + dpx(3))
           dpyi(ip) = coef * (dpy(0) + 2 * dpy(1) + 2 * dpy(2) + dpy(3))
           dpzi(ip) = coef * (dpz(0) + 2 * dpz(1) + 2 * dpz(2) + dpz(3))

           px(ip) = px(ip) + dpxi(ip)
           py(ip) = py(ip) + dpyi(ip)
           pz(ip) = pz(ip) + dpzi(ip)

          endif

       enddo


       end


        !---------------------------------------------------------------------!
        !  Adams-Bashforth 2 time-stepping
        !---------------------------------------------------------------------
        ! For Adams-methods we need to keep the previous dpx,dpy,dpz values 
        ! 

       subroutine timestep_AB2(px,py,pz,dpx,dpy,dpz,iab,
     &          u,v,w,itim,fct,dfct,pm,pn,dz,dt,ng,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none 

!      import/export
       integer(kind=4)                     ,intent(in)   :: nx,ny,nz,ng,np,npmx
       integer(kind=4)                     ,intent(in)   :: i0,j0,k0
       integer(kind=4) ,dimension(0:1)     ,intent(in)   :: itim
       real(kind=8),dimension(nx-1,ny,nz,2),intent(in)   :: u
       real(kind=8),dimension(nx,ny-1,nz,2),intent(in)   :: v
       real(kind=8),dimension(nx,ny,nz+1,2),intent(in)   :: w
       real(kind=8),dimension(nx,ny,nz,2)  ,intent(in)   :: dz
       real(kind=8),dimension(nx,ny)       ,intent(in)   :: pm,pn
       real(kind=8)   ,dimension(np)       ,intent(inout):: px,py,pz
       real(kind=8)                        ,intent(in)   :: dt,fct,dfct

!      import/export for AB scheme
       integer(kind=4),dimension(0:1)        ,intent(in):: iab
       real(kind=8) ,dimension(np,0:1)      ,intent(inout):: dpx,dpy,dpz

!      local
       real(kind=8)             :: coef
       integer(kind=8)          :: ip
       real(kind=8)             :: dpxi,dpyi,dpzi


!f2py intent(inout) px,py,pz
!f2py intent(in)   u,v,w
!f2py intent(in)   fct,dfct,dt
!f2py intent(in)   pm,pn,dz
!f2py intent(in)   ng,npmx
!f2py intent(in)   nx,ny,nz
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   itim
!f2py intent(in)   np
!f2py intent(in)   iab
!f2py intent(inout) dpx,dpy,dpz


       coef = 0.5

       do ip = 1,min(npmx,np)
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           ! Forward Euler / predictor
           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,itim,
     &          fct,pm,pn,dz,dt,i0,j0,k0,nx,ny,nz,ng,np,
     &          dpx(ip,iab(1)),dpy(ip,iab(1)),dpz(ip,iab(1)))

! 
!            if (.not.(isnan(dpx(ip,iab(0)))) ) then

            dpxi =  coef * (3. * dpx(ip,iab(1)) - dpx(ip,iab(0)))
            dpyi =  coef * (3. * dpy(ip,iab(1)) - dpy(ip,iab(0)))
            dpzi =  coef * (3. * dpz(ip,iab(1)) - dpz(ip,iab(0)))
! 
!            else
! 
!             dpx(ip,iab(1)) = dpxi
!             dpy(ip,iab(1)) = dpyi
!             dpz(ip,iab(1)) = dpzi
! 
!            endif

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
        !  Adams-Bashforth 3 time-stepping
        !---------------------------------------------------------------------
        ! For Adams-methods we need to keep the previous dpx,dpy,dpz values 
        ! 

       subroutine timestep_AB3(px,py,pz,dpx,dpy,dpz,iab,
     &          u,v,w,itim,fct,dfct,pm,pn,dz,dt,ng,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none 

!      import/export
       integer(kind=4)                     ,intent(in)   :: nx,ny,nz,ng,np,npmx
       integer(kind=4)                     ,intent(in)   :: i0,j0,k0
       integer(kind=4) ,dimension(0:1)     ,intent(in)   :: itim
       real(kind=8),dimension(nx-1,ny,nz,2),intent(in)   :: u
       real(kind=8),dimension(nx,ny-1,nz,2),intent(in)   :: v
       real(kind=8),dimension(nx,ny,nz+1,2),intent(in)   :: w
       real(kind=8),dimension(nx,ny,nz,2)  ,intent(in)   :: dz
       real(kind=8),dimension(nx,ny)       ,intent(in)   :: pm,pn
       real(kind=8)   ,dimension(np)       ,intent(inout):: px,py,pz
       real(kind=8)                        ,intent(in)   :: dt,fct,dfct

!      import/export for AB scheme
       integer(kind=4),dimension(0:2)        ,intent(inout):: iab
       real(kind=8) ,dimension(np,0:2)      ,intent(inout):: dpx,dpy,dpz

!      local
       real(kind=8)             :: coef
       integer(kind=8)          :: ip
       real(kind=8)             :: dpxi,dpyi,dpzi


!f2py intent(inout) px,py,pz
!f2py intent(in)   u,v,w
!f2py intent(in)   fct,dfct,dt
!f2py intent(in)   pm,pn,dz
!f2py intent(in)   ng,npmx
!f2py intent(in)   nx,ny,nz
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   itim
!f2py intent(in)   np
!f2py intent(inout)   iab
!f2py intent(inout) dpx,dpy,dpz


       coef = 1./12.

       do ip = 1,min(npmx,np)
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           ! Forward Euler / predictor
           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,itim,
     &          fct,pm,pn,dz,dt,i0,j0,k0,nx,ny,nz,ng,np,
     &          dpx(ip,iab(2)),dpy(ip,iab(2)),dpz(ip,iab(2)))

            dpxi =  coef * (23.*dpx(ip,iab(2)) - 16.*dpx(ip,iab(1))
     &                             + 5.*dpx(ip,iab(0)))
            dpyi =  coef * (23.*dpy(ip,iab(2)) - 16.*dpy(ip,iab(1))
     &                             + 5.*dpy(ip,iab(0)))
            dpzi =  coef * (23.*dpz(ip,iab(2)) - 16.*dpz(ip,iab(1))
     &                             + 5.*dpz(ip,iab(0)))


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
        !  Adams-Bashforth 4 time-stepping
        !---------------------------------------------------------------------
        ! For Adams-methods we need to keep the previous dpx,dpy,dpz values 
        ! 

       subroutine timestep_AB4(px,py,pz,dpx,dpy,dpz,iab,
     &          u,v,w,itim,fct,dfct,pm,pn,dz,dt,ng,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none 

!      import/export
       integer(kind=4)                     ,intent(in)   :: nx,ny,nz,ng,np,npmx
       integer(kind=4)                     ,intent(in)   :: i0,j0,k0
       integer(kind=4) ,dimension(0:1)     ,intent(in)   :: itim
       real(kind=8),dimension(nx-1,ny,nz,2),intent(in)   :: u
       real(kind=8),dimension(nx,ny-1,nz,2),intent(in)   :: v
       real(kind=8),dimension(nx,ny,nz+1,2),intent(in)   :: w
       real(kind=8),dimension(nx,ny,nz,2)  ,intent(in)   :: dz
       real(kind=8),dimension(nx,ny)       ,intent(in)   :: pm,pn
       real(kind=8)   ,dimension(np)       ,intent(inout):: px,py,pz
       real(kind=8)                        ,intent(in)   :: dt,fct,dfct

!      import/export for AB scheme
       integer(kind=4),dimension(0:3)        ,intent(inout):: iab
       real(kind=8) ,dimension(np,0:3)      ,intent(inout):: dpx,dpy,dpz

!      local
       real(kind=8)             :: coef
       integer(kind=8)          :: ip
       real(kind=8)            :: dpxi,dpyi,dpzi


!f2py intent(inout) px,py,pz
!f2py intent(in)   u,v,w
!f2py intent(in)   fct,dfct,dt
!f2py intent(in)   pm,pn,dz
!f2py intent(in)   ng,npmx
!f2py intent(in)   nx,ny,nz
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   itim
!f2py intent(in)   np
!f2py intent(inout)   iab
!f2py intent(inout) dpx,dpy,dpz

       coef = 1./24.

       do ip = 1,min(npmx,np)
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           ! Forward Euler / predictor
           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,itim,
     &          fct,pm,pn,dz,dt,i0,j0,k0,nx,ny,nz,ng,np,
     &          dpx(ip,iab(3)),dpy(ip,iab(3)),dpz(ip,iab(3)))
! 
!            if (.not.(isnan(dpx(ip,iab(0)))) ) then



           ! Adams-Bashforth 4
            dpxi =   coef * (55.*dpx(ip,iab(3)) - 59.*dpx(ip,iab(2))
     &                             + 37.*dpx(ip,iab(1)) - 9.*dpx(ip,iab(0)))
            dpyi =   coef * (55.*dpy(ip,iab(3)) - 59.*dpy(ip,iab(2))
     &                             + 37.*dpy(ip,iab(1)) - 9.*dpy(ip,iab(0)))
            dpzi =   coef * (55.*dpz(ip,iab(3)) - 59.*dpz(ip,iab(2))
     &                             + 37.*dpz(ip,iab(1)) - 9.*dpz(ip,iab(0)))

! 
!            else if (.not.(isnan(dpx(ip,iab(1)))) ) then
! 
!             if (ip.eq.1) print *,'time step 2'
!             dpxi =   1/12. * (23.*dpx(ip,iab(3)) - 16.*dpx(ip,iab(2))
!     &                              + 5.*dpx(ip,iab(1)) )
!             dpyi =   1/12. * (23.*dpy(ip,iab(3)) - 16.*dpy(ip,iab(2))
!     &                              + 5.*dpy(ip,iab(1)) )
!             dpzi =   1/12. * (23.*dpz(ip,iab(3)) - 16*dpz(ip,iab(2))
!     &                              + 5.*dpz(ip,iab(1)) )
! 
!            else if (.not.(isnan(dpx(ip,iab(2)))) ) then
! 
!             if (ip.eq.1) print *,'time step 1'
!             dpxi =   0.5 * (3. * dpx(ip,iab(3)) - dpx(ip,iab(2)))
!             dpyi =   0.5 * (3. * dpy(ip,iab(3)) - dpy(ip,iab(2)))
!             dpzi =   0.5 * (3. * dpz(ip,iab(3)) - dpz(ip,iab(2)))
! 
!            else
! 
!             if (ip.eq.1) print *,'time step 0'
!             dpxi = dpx(ip,iab(3))
!             dpyi = dpy(ip,iab(3))
!             dpzi = dpz(ip,iab(3))
! 
!            endif


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
        !  Adams-Bashforth 4 Predictor
        ! + Adams-Moulton 5 Corrector
        !---------------------------------------------------------------------
        ! For Adams-methods we need to keep the previous dpx,dpy,dpz values 
        ! 

       subroutine timestep_ABM4(px,py,pz,dpx,dpy,dpz,iab,
     &          u,v,w,itim,fct,dfct,pm,pn,dz,dt,ng,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none 

!      import/export
       integer(kind=4)                     ,intent(in)   :: nx,ny,nz,ng,np,npmx
       integer(kind=4)                     ,intent(in)   :: i0,j0,k0
       integer(kind=4) ,dimension(0:1)     ,intent(in)   :: itim
       real(kind=8),dimension(nx-1,ny,nz,2),intent(in)   :: u
       real(kind=8),dimension(nx,ny-1,nz,2),intent(in)   :: v
       real(kind=8),dimension(nx,ny,nz+1,2),intent(in)   :: w
       real(kind=8),dimension(nx,ny,nz,2)  ,intent(in)   :: dz
       real(kind=8),dimension(nx,ny)       ,intent(in)   :: pm,pn
       real(kind=8)   ,dimension(np)       ,intent(inout):: px,py,pz
       real(kind=8)                        ,intent(in)   :: dt,fct,dfct

!      import/export for AB scheme
       integer(kind=4),dimension(0:3)        ,intent(inout):: iab
       real(kind=8) ,dimension(np,0:3)      ,intent(inout):: dpx,dpy,dpz

!      local
       real(kind=8)             :: coef
       integer(kind=8)          :: ip
       real(kind=8)             :: dpxc,dpyc,dpzc
       real(kind=8)             :: dpxi,dpyi,dpzi

!f2py intent(inout) px,py,pz
!f2py intent(in)   u,v,w
!f2py intent(in)   fct,dfct,dt
!f2py intent(in)   pm,pn,dz
!f2py intent(in)   npmx
!f2py intent(in)   nx,ny,nz
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   itim
!f2py intent(in)   ng,np
!f2py intent(inout)   iab
!f2py intent(inout) dpx,dpy,dpz



       coef = 1./24.

       do ip = 1,min(npmx,np)
         if (.not.(isnan(px(ip))) ) then

           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,itim,
     &          fct,pm,pn,dz,dt,i0,j0,k0,nx,ny,nz,ng,np,
     &          dpx(ip,iab(3)),dpy(ip,iab(3)),dpz(ip,iab(3)))

           !---------------------------------------
           ! Predictor
           !---------------------------------------     
!    
!            if (.not.(isnan(dpx(ip,iab(0)))) ) then

            dpxi =   coef * (55.*dpx(ip,iab(3)) - 59.*dpx(ip,iab(2))
     &                     + 37.*dpx(ip,iab(1)) - 9.*dpx(ip,iab(0)))
            dpyi =   coef * (55.*dpy(ip,iab(3)) - 59.*dpy(ip,iab(2))
     &                     + 37.*dpy(ip,iab(1)) - 9.*dpy(ip,iab(0)))
            dpzi =   coef * (55.*dpz(ip,iab(3)) - 59.*dpz(ip,iab(2))
     &                     + 37.*dpz(ip,iab(1)) - 9.*dpz(ip,iab(0)))
! 
!            else if (.not.(isnan(dpx(ip,iab(1)))) ) then
! 
!             if (ip.eq.1) print *,'time step 2'
!             dpxi =   1/12. * (23.*dpx(ip,iab(3)) - 16.*dpx(ip,iab(2))
!     &                              + 5.*dpx(ip,iab(1)) )
!             dpyi =   1/12. * (23.*dpy(ip,iab(3)) - 16.*dpy(ip,iab(2))
!     &                              + 5.*dpy(ip,iab(1)) )
!             dpzi =   1/12.  * (23.*dpz(ip,iab(3)) - 16*dpz(ip,iab(2))
!     &                              + 5.*dpz(ip,iab(1)) )
! 
!            else if (.not.(isnan(dpx(ip,iab(2)))) ) then
! 
!             if (ip.eq.1) print *,'time step 1'
!             dpxi =   0.5 * (3. * dpx(ip,iab(3)) - dpx(ip,iab(2)))
!             dpyi =   0.5 * (3. * dpy(ip,iab(3)) - dpy(ip,iab(2)))
!             dpzi =   0.5 * (3. * dpz(ip,iab(3)) - dpz(ip,iab(2)))
! 
!            else
! 
!             if (ip.eq.1) print *,'time step 0'
!             dpxi = dpxi
!             dpyi = dpyi
!             dpzi = dpzi
! 
!            endif



           !---------------------------------------
           ! Compute velocities at new time-step
           !---------------------------------------  

           CALL advance_3d(px(ip)+dpxi,py(ip)+dpyi,pz(ip)+dpzi,
     &          u,v,w,itim,
     &          fct+dfct,pm,pn,dz,dt,i0,j0,k0,nx,ny,nz,ng,np,dpxc,dpyc,dpzc)


           !--------------------------------------- 
           ! Corrector
           !--------------------------------------- 
! 
!            if (.not.(isnan(dpx(ip,iab(0)))) ) then

            dpxi = 1./270. * ( 19.* dpxi + 251.* (  1./24. * (
     &                   9. * dpxc + 19. * dpx(ip,iab(3))
     &                 - 5. * dpx(ip,iab(2)) + dpx(ip,iab(1)) )))

            dpyi = 1./270. * ( 19.* dpyi + 251.* (  1./24. * (
     &                   9. * dpyc + 19. * dpy(ip,iab(3))
     &                 - 5. * dpy(ip,iab(2)) + dpy(ip,iab(1)) )))

            dpzi = 1./270. * ( 19.* dpzi + 251.* ( 1./24. * (
     &                   9. * dpzc + 19. * dpz(ip,iab(3))
     &                 - 5. * dpz(ip,iab(2)) + dpz(ip,iab(1)) )))
! 
! 
!            else if (.not.(isnan(dpx(ip,iab(1)))) ) then
! 
!             print *,'time step 2'
!             dpxi = 1./270. * ( 19.* dpxi + 251.* ( 1./12. * (
!     &                    5. * dpxc + 8. * dpx(ip,iab(3))
!     &                  - 1. * dpx(ip,iab(2))  )))
!             dpyi = 1./270. * ( 19.* dpyi + 251.* (  1./12. * (
!     &                    5. * dpyc + 8. * dpy(ip,iab(3))
!     &                  - 1. * dpy(ip,iab(2))  )))
!             dpzi = 1./270. * ( 19.* dpzi + 251.* (  1./12. * (
!     &                    5. * dpzc + 8. * dpz(ip,iab(3))
!     &                  - 1. * dpz(ip,iab(2))  )))
! 
!            else if (.not.(isnan(dpx(ip,iab(2)))) ) then
! 
!             print *,'time step 1'
!             dpxi = 1./270. * ( 19.* dpxi + 251.* (  1./2. * (
!     &                    1. * dpxc + 1. * dpx(ip,iab(3)) )))
!             dpyi = 1./270. * ( 19.* dpyi + 251.* (  1./2. * (
!     &                    1. * dpyc + 1. * dpy(ip,iab(3)) )))
!             dpzi = 1./270. * ( 19.* dpzi + 251.* (  1./2. * (
!     &                    1. * dpzc + 1. * dpz(ip,iab(3)) )))
! 
!            else
! 
!             print *,'time step 0'
!             pxi = 1./270. * ( 19.* pxi + 251.* ( px(ip)  +  dpxc  ))
!             pyi = 1./270. * ( 19.* pyi + 251.* ( py(ip)  +  dpyc  ))
!             pzi = 1./270. * ( 19.* pzi + 251.* ( pz(ip)  +  dpzc  ))
! 
!            endif


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
!!Omega
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# include "R_tools_fort_routines/get_omega.F"









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





