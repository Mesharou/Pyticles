!----------------------------------------------------------------------------------------------
! Fortran Routines for Pyticles
!
!----------------------------------------------------------------------------------------------
! 05/10/16:
!     - add mask as a shared array and use it in the time_step routines
!     - add the # define CHECK_MASK in pyticles_3d_sig_sa.F
!     - add subroutines interp_3d_u / interp_3d_u
! 16/01/26:
!     - Add parameter 'ng' corresponding to number of ghost points in the horizontal grid
!----------------------------------------------------------------------------------------------

# define CHECK_MASK


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
# ifdef CHECK_MASK
     &          mask,
# endif
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
# ifdef CHECK_MASK
       real(kind=8) ,dimension(nx,ny),intent(in)        :: mask
# endif
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
# ifdef CHECK_MASK
!f2py intent(in)   mask
# endif

       do ip = 1,min(npmx,np)
        
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,itim,fct,pm,pn
     &         ,dz,dt,i0,j0,k0,nx,ny,nz,ng,np,dpx,dpy,dpz)

           !---------------------------------------
           ! Update position
           !---------------------------------------
# ifdef CHECK_MASK
           call check_mask(mask,px(ip),py(ip),dpx,dpy
     &                                    ,ng,npmx,i0,j0,nx,ny)
# endif
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
# ifdef CHECK_MASK
     &          mask,
# endif
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
# ifdef CHECK_MASK
       real(kind=8) ,dimension(nx,ny),intent(in)        :: mask
# endif
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
# ifdef CHECK_MASK
!f2py intent(in)   mask
# endif

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
# ifdef CHECK_MASK
           call check_mask(mask,px(ip),py(ip),dpx,dpy
     &                                    ,ng,npmx,i0,j0,nx,ny)
# endif
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
# ifdef CHECK_MASK
     &          mask,
# endif
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
# ifdef CHECK_MASK
       real(kind=8) ,dimension(nx,ny),intent(in)        :: mask
# endif
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
# ifdef CHECK_MASK
!f2py intent(in)   mask
# endif
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

# ifdef CHECK_MASK
           !write(*,*) 'in RK4',px(ip),py(ip),dpxi(ip),dpyi(ip)
           call check_mask(mask,px(ip),py(ip),dpxi(ip),dpyi(ip)
     &                                    ,ng,npmx,i0,j0,nx,ny)
# endif

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
# ifdef CHECK_MASK
     &          mask,
# endif
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
# ifdef CHECK_MASK
       real(kind=8) ,dimension(nx,ny),intent(in)        :: mask
# endif
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
# ifdef CHECK_MASK
!f2py intent(in)   mask
# endif
!f2py intent(inout) dpx,dpy,dpz


       coef = 0.5

       do ip = 1,min(npmx,np)
         if (.not.(isnan(px(ip)*py(ip)*pz(ip))) ) then

           ! Forward Euler / predictor
           CALL advance_3d(px(ip),py(ip),pz(ip),u,v,w,itim,
     &          fct,pm,pn,dz,dt,i0,j0,k0,nx,ny,nz,ng,np,
     &          dpx(ip,iab(1)),dpy(ip,iab(1)),dpz(ip,iab(1)))



            dpxi =  coef * (3. * dpx(ip,iab(1)) - dpx(ip,iab(0)))
            dpyi =  coef * (3. * dpy(ip,iab(1)) - dpy(ip,iab(0)))
            dpzi =  coef * (3. * dpz(ip,iab(1)) - dpz(ip,iab(0)))

            !---------------------------------------
            ! Update position
            !--------------------------------------- 

# ifdef CHECK_MASK
           call check_mask(mask,px(ip),py(ip),dpxi,dpyi
     &                                    ,ng,npmx,i0,j0,nx,ny)
# endif
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
# ifdef CHECK_MASK
     &          mask,
# endif
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
# ifdef CHECK_MASK
       real(kind=8) ,dimension(nx,ny),intent(in)        :: mask
# endif
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
# ifdef CHECK_MASK
!f2py intent(in)   mask
# endif
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

# ifdef CHECK_MASK
           call check_mask(mask,px(ip),py(ip),dpxi,dpyi
     &                                    ,ng,npmx,i0,j0,nx,ny)
# endif

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
# ifdef CHECK_MASK
     &          mask,
# endif
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
# ifdef CHECK_MASK
       real(kind=8) ,dimension(nx,ny),intent(in)        :: mask
# endif
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
# ifdef CHECK_MASK
!f2py intent(in)   mask
# endif
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


            !---------------------------------------
            ! Update position
            !--------------------------------------- 

# ifdef CHECK_MASK
           call check_mask(mask,px(ip),py(ip),dpxi,dpyi
     &                                    ,ng,npmx,i0,j0,nx,ny)
# endif

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
# ifdef CHECK_MASK
     &          mask,
# endif
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
# ifdef CHECK_MASK
       real(kind=8) ,dimension(nx,ny),intent(in)        :: mask
# endif
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
# ifdef CHECK_MASK
!f2py intent(in)   mask
# endif
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


            dpxi =   coef * (55.*dpx(ip,iab(3)) - 59.*dpx(ip,iab(2))
     &                     + 37.*dpx(ip,iab(1)) - 9.*dpx(ip,iab(0)))
            dpyi =   coef * (55.*dpy(ip,iab(3)) - 59.*dpy(ip,iab(2))
     &                     + 37.*dpy(ip,iab(1)) - 9.*dpy(ip,iab(0)))
            dpzi =   coef * (55.*dpz(ip,iab(3)) - 59.*dpz(ip,iab(2))
     &                     + 37.*dpz(ip,iab(1)) - 9.*dpz(ip,iab(0)))


           !---------------------------------------
           ! Compute velocities at new time-step
           !---------------------------------------  

           CALL advance_3d(px(ip)+dpxi,py(ip)+dpyi,pz(ip)+dpzi,
     &          u,v,w,itim,
     &          fct+dfct,pm,pn,dz,dt,i0,j0,k0,nx,ny,nz,ng,np,dpxc,dpyc,dpzc)


           !--------------------------------------- 
           ! Corrector
           !--------------------------------------- 

            dpxi = 1./270. * ( 19.* dpxi + 251.* (  1./24. * (
     &                   9. * dpxc + 19. * dpx(ip,iab(3))
     &                 - 5. * dpx(ip,iab(2)) + dpx(ip,iab(1)) )))

            dpyi = 1./270. * ( 19.* dpyi + 251.* (  1./24. * (
     &                   9. * dpyc + 19. * dpy(ip,iab(3))
     &                 - 5. * dpy(ip,iab(2)) + dpy(ip,iab(1)) )))

            dpzi = 1./270. * ( 19.* dpzi + 251.* ( 1./24. * (
     &                   9. * dpzc + 19. * dpz(ip,iab(3))
     &                 - 5. * dpz(ip,iab(2)) + dpz(ip,iab(1)) )))


            !---------------------------------------
            ! Update position
            !--------------------------------------- 

# ifdef CHECK_MASK
           call check_mask(mask,px(ip),py(ip),dpxi,dpyi
     &                                    ,ng,npmx,i0,j0,nx,ny)
# endif

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
       !  Check mask routine (check if particle is overshooting into mask)
       !---------------------------------------------------------------------


       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny
       integer(kind=4)                    ,intent(in)   :: ng,npmx
       real(kind=8),dimension(nx,ny)      ,intent(in)   :: mask
       real(kind=8)                       ,intent(in)   :: px,py
       real(kind=8)                       ,intent(inout):: dpx,dpy
       integer(kind=4)                    ,intent(in)   :: i0,j0
!      local
       integer(kind=8)                                  :: i,j
       real(kind=8),dimension(2,2)                      :: wt
       real(kind=8)                                     :: fcx,fcy,pmask
       real(kind=8),dimension(4)                        :: dmask


!f2py intent(inout) dpx,dpy
!f2py intent(in)    mask,px,py
!f2py intent(in)    npmx
!f2py intent(in)    i0,j0
!f2py intent(in)    ng
!f2py intent(in)    nx,ny

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

          dmask(1) = (1-fcy)*mask(i,j)   + fcy*mask(i,j+1)
          dmask(2) = (1-fcy)*mask(i+1,j) + fcy*mask(i+1,j+1)
          dmask(3) = (1-fcx)*mask(i,j)   + fcx*mask(i+1,j)
          dmask(4) = (1-fcx)*mask(i,j+1) + fcx*mask(i+1,j+1)

          !write(*,*) 'dmask(1),dmask(2) ', dmask(1),dmask(2)
          !write(*,*) 'dmask(3),dmask(4) ', dmask(3),dmask(4)

          if ( ((dpx.gt.0).and.((dmask(1)-dmask(2)).ge.0.5)) .or.
     &         ((dpx.lt.0).and.((dmask(2)-dmask(1)).ge.0.5)) ) dpx = 0.

          if ( ((dpy.gt.0).and.((dmask(3)-dmask(4)).ge.0.5)) .or.
     &         ((dpy.lt.0).and.((dmask(4)-dmask(3)).ge.0.5)) ) dpy = 0.

          !write(*,*) 'correcting mask for ', px,py
          !write(*,*) 'corrected ', dpx,dpy

        endif

       end







       !---------------------------------------------------------------------!
       !  Runge-Kutta 4 time-stepping for 2D advection
       !---------------------------------------------------------------------

!
 

       subroutine timestep2d_RK4(px,py,u,v,itim,fct,dfct,pm,pn,
# ifdef CHECK_MASK
     &          mask,
# endif
     &          dt,ng,npmx,i0,j0,nx,ny,np,dpxi,dpyi)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny
       integer(kind=4)                    ,intent(in)   :: ng,np
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0
       integer(kind=4) ,dimension(0:1)    ,intent(in)   :: itim
       real(kind=8) ,dimension(nx-1,ny,2),intent(in)    :: u
       real(kind=8) ,dimension(nx,ny-1,2),intent(in)    :: v
       real(kind=8) ,dimension(nx,ny),intent(in)        :: pm,pn
# ifdef CHECK_MASK
       real(kind=8) ,dimension(nx,ny),intent(in)        :: mask
# endif
       real(kind=8) ,dimension(np)        ,intent(inout):: px,py
       real(kind=8) ,dimension(np)        ,intent(out)  :: dpxi,dpyi
       real(kind=8)                       ,intent(in)   :: dt,fct,dfct
!      local
       real(kind=8) ,dimension(0:3)            :: dpx,dpy
       real(kind=8)                           :: coef
       integer(kind=8)                        :: ip

!f2py intent(inout) px,py
!f2py intent(in)   u,v
!f2py intent(in)   itim
!f2py intent(in)   fct,dfct,dt
!f2py intent(in)   pm,pn
!f2py intent(in)   ng,npmx
!f2py intent(in)   i0,j0
!f2py intent(in)   nx,ny
!f2py intent(in)   np
# ifdef CHECK_MASK
!f2py intent(in)   mask
# endif
!f2py intent(out)  dpxi,dpyi

       coef = 1./6.

       do ip = 1,min(npmx,np)

         if (.not.(isnan(px(ip)*py(ip)))) then 

           ! Forward Euler / predictor
           CALL advance_2d(px(ip),py(ip),u,v,
     &          itim,fct,pm,pn,dt,i0,j0,
     &          nx,ny,ng,np,dpx(0),dpy(0))

           if (.not.(isnan(dpx(0)*dpy(0)))) then 

           ! backward Euler / corrector
           CALL advance_2d(px(ip)+0.5*dpx(0),py(ip)+0.5*dpy(0),
     &          u,v,itim,fct+0.5*dfct,pm,pn,dt,i0,j0,
     &          nx,ny,ng,np,dpx(1),dpy(1))

            if (.not.(isnan(dpx(1)*dpy(1)))) then 

           ! midpoint rule / predictor
           CALL advance_2d(px(ip)+0.5*dpx(1),py(ip)+0.5*dpy(1),
     &          u,v,itim,fct+0.5*dfct,pm,pn,dt,i0,j0,
     &          nx,ny,ng,np,dpx(2),dpy(2))

             if (.not.(isnan(dpx(2)*dpy(2)))) then 

           ! Corrector
           CALL advance_2d(px(ip)+dpx(2),py(ip)+dpy(2),
     &          u,v,itim,fct+1.*dfct,pm,pn,dt,i0,j0,
     &          nx,ny,ng,np,dpx(3),dpy(3))

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

# define CUBIC
# define DUKO_2001
# define INTERP_BELOW
# define INTERP_ABOVE

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute z_r and z_w for NEW_S_COORD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# include "R_tools_fort_routines/zlevs.F"

# include "R_tools_fort_routines/zlevs_w.F"

# include "R_tools_fort_routines/zlevs_croco_new.F"

# include "R_tools_fort_routines/zlevs_croco_new_w.F"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute z_r and z_w for OLD_S_COORD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# include "R_tools_fort_routines/zlevs_croco_old.F"

# include "R_tools_fort_routines/zlevs_croco_old_w.F"

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
!!rho1_eos
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# include "R_tools_fort_routines/rho1_eos.F"








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





