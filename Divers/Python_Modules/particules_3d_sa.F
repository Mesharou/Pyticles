
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!! 
!! cpp particules_3d_sa.F particules_3d_sa.f90 ; f2py --f90flags="-extend_source -O1 " -DF2PY_REPORT_ON_ARRAY_COPY=1 -c -m particules_3d_sa particules_3d_sa.f90
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!----------------------------------------------------------------------------------------------
       subroutine seed(px,py,pz,npmx,mask,lev0,lev1,nnx,nny,nnlev,nx,ny,np)
       implicit none 
!      import/export
       integer(kind=4)                 ,intent(in)   :: nx,ny
       integer(kind=4)                 ,intent(in)   :: nnx, nny, nnlev
       integer(kind=4)                 ,intent(in)   :: np
       integer(kind=4)                 ,intent(in)   :: lev0, lev1
       integer(kind=4),dimension(nx,ny),intent(in)   :: mask
       real(kind=8)   ,dimension(np)   ,intent(inout):: px,py,pz
       integer(kind=4)                 ,intent(inout):: npmx
!      local
       integer(kind=4) :: ip,i,j,id, k
       integer(kind=4),allocatable,dimension(:,:) :: grid

!f2py intent(inout) px,py,pz
!f2py intent(in)    mask
!f2py intent(in)    lev0,lev1
!f2py intent(inout) npmx
!f2py intent(in)    nx,ny, nnx, nny, nnlev
!f2py intent(in)    np

!      print *,'fo: ',npmx
       allocate(grid(nx,ny))
       grid = 0


       !! find empty cells
       do ip = 1,npmx  !! max particle index that is currently used
          i = floor(px(ip))
          j = floor(py(ip))
          grid(i,j) = 1
       enddo
!
!      empty cells in the mask will get a new particle
       id = 1
       do i = 1,nx,nnx
        do j = 1,ny,nny
         do k = lev0,lev1,nnlev
           if (mask(i,j).eq.1.and.grid(i,j).eq.0) then
             !! grab a particle from inbetween, if available
             do ip = id,npmx 
               if ( isnan(px(ip)) ) then
                 print * ,'found one'
                 px(ip) = i-1
                 py(ip) = j-1
                 pz(ip) = k
                 grid(i,j) = 1
                 exit
               endif
             enddo
             !! grab a particle from the end
             if (grid(i,j).eq.0.and.npmx.lt.np) then
               npmx = npmx + 1
               ip = npmx
               px(npmx) = i-1
               py(npmx) = j-1
               pz(npmx) = k
             endif
           endif
        enddo
       enddo
       enddo
       deallocate(grid)
!      print *,'fo: ',npmx

*
       end


!----------------------------------------------------------------------------------------------
       subroutine advance_3d(px,py,pz,u0,v0,w0,u1,v1,w1,fct,dx,dy,dz,dt,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny,nz
       integer(kind=4)                    ,intent(in)   :: np
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0,k0
       real(kind=8)   ,dimension(nx,ny,nz),intent(in)   :: u0,v0,w0
       real(kind=8)   ,dimension(nx,ny,nz),intent(in)   :: u1,v1,w1
       real(kind=8)   ,dimension(np)      ,intent(inout):: px,py,pz
       real(kind=8)                       ,intent(in)   :: dx,dy,dz,dt,fct
!      local
       integer(kind=8)                  :: ip,jp,kp,i,j,k
       real(kind=8)   ,dimension(2,2,2) :: wt
       real(kind=8)                     :: fcx,fcy,fcz,dxi,dyi,dzi,pu,pv,pw


!f2py intent(in)   u0,v0,w0
!f2py intent(in)   u1,v1,w1
!f2py intent(inout) px,py,pz
!f2py intent(in)   dx,dy,dz,dt
!f2py intent(in)   npmx
!f2py intent(in)   nx,ny,nz
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   np

       dxi  = 1/dx
       dyi  = 1/dy
       dzi  = 1/dz


       do ip = 1,min(npmx,np)
         if (.not.(isnan(px(ip))) ) then

           i = floor(px(ip)+1)-i0
           j = floor(py(ip)+1)-j0
           k = floor(pz(ip)+1)-k0

           if (isnan(w1(i,j,k))) then

            print *,'we have a nan at', px(ip), py(ip), pz(ip)
            pz(ip)=-1

           else

           fcx = px(ip)+1 - i - i0;
           fcy = py(ip)+1 - j - j0;
           fcz = pz(ip)+1 - k - k0;

           wt(1,1,1) = (1-fcz)*(1-fcy)*(1-fcx);
           wt(1,1,2) =    fcz *(1-fcy)*(1-fcx);
           wt(1,2,1) = (1-fcz)*   fcy *(1-fcx);
           wt(1,2,2) =    fcz *  fcy *(1-fcx);
           wt(2,1,1) = (1-fcz)*(1-fcy)*   fcx ;
           wt(2,1,2) =    fcz *(1-fcy)*   fcx ;
           wt(2,2,1) = (1-fcz)*   fcy *   fcx ;
           wt(2,2,2) =    fcz *   fcy *   fcx ;      

           if (k.ge.nz) then  
           !The particule has reach the surface 
            pu   = (1-fct)*sum(u0(i:i+1,j:j+1,k)*wt(:,:,1))+fct*sum(u1(i:i+1,j:j+1,k)*wt(:,:,1))
            pv   = (1-fct)*sum(v0(i:i+1,j:j+1,k)*wt(:,:,1))+fct*sum(v1(i:i+1,j:j+1,k)*wt(:,:,1))
            pw   = min((1-fct)*sum(w0(i:i+1,j:j+1,k)*wt(:,:,1))+fct*sum(w1(i:i+1,j:j+1,k)*wt(:,:,1)),0.)
           else
            pu   = (1-fct)*sum(u0(i:i+1,j:j+1,k:k+1)*wt)+fct*sum(u1(i:i+1,j:j+1,k:k+1)*wt)
            pv   = (1-fct)*sum(v0(i:i+1,j:j+1,k:k+1)*wt)+fct*sum(v1(i:i+1,j:j+1,k:k+1)*wt)
            pw   = (1-fct)*sum(w0(i:i+1,j:j+1,k:k+1)*wt)+fct*sum(w1(i:i+1,j:j+1,k:k+1)*wt)
           endif

           px(ip) = px(ip) + dt*pu*dxi
           py(ip) = py(ip) + dt*pv*dyi
           pz(ip) = pz(ip) + dt*pw*dzi


         endif

!           if (pz(ip)- k0>nz-1) pz(ip)=nz-1-k0
!           if ((px(ip)- i0>nx-1) .or. (px(ip)<0)) print*,'px',ip,px(ip),py(ip),pz(ip), pw,dt*pw*dzi
!           if ((py(ip)- j0>ny-1) .or. (py(ip)<0)) print*,'py',ip,px(ip),py(ip),pz(ip), pw,dt*pw*dzi
!           if ((pz(ip)- k0>nz-1) .or. (pz(ip)<0)) print*,'pz',ip,px(ip),py(ip),pz(ip), pw,dt*pw*dzi
!           print *,ip, px(ip),py(ip),pz(ip)

         endif
       enddo
!
       end
!----------------------------------------------------------------------------------------------
       subroutine advance_3d_multi(px,py,pz,u,v,w,fct,dx,dy,dz,dt,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nx,ny,nz
       integer(kind=4)                    ,intent(in)   :: np
       integer(kind=4)                    ,intent(in)   :: npmx
       integer(kind=4)                    ,intent(in)   :: i0,j0,k0
       real(kind=8) ,dimension(nx,ny,nz,2),intent(in)   :: u,v,w
       real(kind=8)   ,dimension(np)      ,intent(inout):: px,py,pz
       real(kind=8)                       ,intent(in)   :: dx,dy,dz,dt,fct
!      local
       integer(kind=8)                  :: ip,jp,kp,i,j,k
       real(kind=8) ,dimension(2,2,2,2) :: wt
       real(kind=8)                     :: fcx,fcy,fcz,dxi,dyi,dzi,pu,pv,pw


!f2py intent(in)   u,v,w
!f2py intent(inout) px,py,pz
!f2py intent(in)   dx,dy,dz,dt
!f2py intent(in)   npmx
!f2py intent(in)   nx,ny,nz
!f2py intent(in)   i0,j0,k0
!f2py intent(in)   np

       dxi  = 1/dx
       dyi  = 1/dy
       dzi  = 1/dz


       do ip = 1,min(npmx,np)
         if (.not.(isnan(px(ip))) ) then

           i = floor(px(ip)+1)-i0
           j = floor(py(ip)+1)-j0
           k = floor(pz(ip)+1)-k0

           if (isnan(w(i,j,k,2))) then

            print *,'we have a nan at', px(ip), py(ip), pz(ip)
            pz(ip)=-1

           else

           fcx = px(ip)+1 - i - i0;
           fcy = py(ip)+1 - j - j0;
           fcz = pz(ip)+1 - k - k0;

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

           if (k.ge.nz) then  
           !The particule has reach the surface 
            pu   = sum(u(i:i+1,j:j+1,k,:)*wt(:,:,1,:))
            pv   = sum(v(i:i+1,j:j+1,k,:)*wt(:,:,1,:))
            pw   = min(sum(w(i:i+1,j:j+1,k,:)*wt(:,:,1,:)),0.)
           else
            pu   = sum(u(i:i+1,j:j+1,k:k+1,:)*wt)
            pv   = sum(v(i:i+1,j:j+1,k:k+1,:)*wt)
            pw   = sum(w(i:i+1,j:j+1,k:k+1,:)*wt)
           endif

           px(ip) = px(ip) + dt*pu*dxi
           py(ip) = py(ip) + dt*pv*dyi
           pz(ip) = pz(ip) + dt*pw*dzi


         endif

!           if (pz(ip)- k0>nz-1) pz(ip)=nz-1-k0
!           if ((px(ip)- i0>nx-1) .or. (px(ip)<0)) print*,'px',ip,px(ip),py(ip),pz(ip), pw,dt*pw*dzi
!           if ((py(ip)- j0>ny-1) .or. (py(ip)<0)) print*,'py',ip,px(ip),py(ip),pz(ip), pw,dt*pw*dzi
!           if ((pz(ip)- k0>nz-1) .or. (pz(ip)<0)) print*,'pz',ip,px(ip),py(ip),pz(ip), pw,dt*pw*dzi
!           print *,ip, px(ip),py(ip),pz(ip)

         endif
       enddo
!
       end

!----------------------------------------------------------------------------------------------
       subroutine interp_3d(pvar1,pvar2,px,py,pz,var1,var2,npmx,i0,j0,k0,nx,ny,nz,np)
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
       real(kind=8)   ,dimension(2,2,2) :: wt
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
           i = floor(px(ip)+1)-i0
           j = floor(py(ip)+1)-j0
           k = floor(pz(ip)+1)-k0
           fcx = px(ip)+1 - i - i0;
           fcy = py(ip)+1 - j - j0;
           fcz = pz(ip)+1 - k - k0;
           wt(1,1,1) = (1-fcz)*(1-fcy)*(1-fcx);
           wt(1,1,2) =    fcz *(1-fcy)*(1-fcx);
           wt(1,2,1) = (1-fcz)*   fcy *(1-fcx);
           wt(1,2,2) =    fcz *  fcy *(1-fcx);
           wt(2,1,1) = (1-fcz)*(1-fcy)*   fcx ;
           wt(2,1,2) =    fcz *(1-fcy)*   fcx ;
           wt(2,2,1) = (1-fcz)*   fcy *   fcx ;
           wt(2,2,2) =    fcz *   fcy *   fcx ;
           if (k.ge.nz) then   
             pvar1(ip)   = sum(var1(i:i+1,j:j+1,k)*wt(:,:,1))
             pvar2(ip)   = sum(var2(i:i+1,j:j+1,k)*wt(:,:,1) )
           else
             pvar1(ip)   = sum(var1(i:i+1,j:j+1,k:k+1)*wt)
             pvar2(ip)   = sum(var2(i:i+1,j:j+1,k:k+1)*wt )
           endif
         endif
       enddo
!
       end    




!----------------------------------------------------------------------------------------------
       subroutine oneterp_3d(pvar1,px,py,pz,var1,npmx,i0,j0,k0,nx,ny,nz,np)
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
       real(kind=8)   ,dimension(2,2,2) :: wt
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
           i = floor(px(ip)+1)-i0
           j = floor(py(ip)+1)-j0
           k = floor(pz(ip)+1)-k0
           fcx = px(ip)+1 - i - i0;
           fcy = py(ip)+1 - j - j0;
           fcz = pz(ip)+1 - k - k0;
           wt(1,1,1) = (1-fcz)*(1-fcy)*(1-fcx);
           wt(1,1,2) =    fcz *(1-fcy)*(1-fcx);
           wt(1,2,1) = (1-fcz)*   fcy *(1-fcx);
           wt(1,2,2) =    fcz *  fcy *(1-fcx);
           wt(2,1,1) = (1-fcz)*(1-fcy)*   fcx ;
           wt(2,1,2) =    fcz *(1-fcy)*   fcx ;
           wt(2,2,1) = (1-fcz)*   fcy *   fcx ;
           wt(2,2,2) =    fcz *   fcy *   fcx ;
           if (k.ge.nz) then   
             pvar1(ip)   = sum(var1(i:i+1,j:j+1,k)*wt(:,:,1))
           else
             pvar1(ip)   = sum(var1(i:i+1,j:j+1,k:k+1)*wt)
           endif
         endif
       enddo
!
       end    


!----------------------------------------------------------------------------------------------
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
           i = floor(px(ip)+1)-i0
           j = floor(py(ip)+1)-j0
           fcx = px(ip)+1 - i - i0;
           fcy = py(ip)+1 - j - j0;
           wt(1,1) =    (1-fcy)*(1-fcx);
           wt(1,2) =       fcy *(1-fcx);
           wt(2,1) =    (1-fcy)*   fcx ;
           wt(2,2) =       fcy *   fcx ;
           pvar1(ip)   = sum(var1(i:i+1,j:j+1)*wt)
         endif
       enddo
!
       end   



!----------------------------------------------------------------------------------------------
       subroutine interp_1d(pvar1,pz,var1,npmx,nz,np)
       implicit none 
!      import/export
       integer(kind=4)                    ,intent(in)   :: nz
       integer(kind=4)                    ,intent(in)   :: np
       integer(kind=4)                    ,intent(in)   :: npmx
       real(kind=8)   ,dimension(nz),intent(in)   :: var1
       real(kind=8)   ,dimension(np)      ,intent(out):: pvar1
       real(kind=8)   ,dimension(np)      ,intent(in):: pz
!      local
       integer(kind=8)                  :: kp,k,ip
       real(kind=8)   ,dimension(2) :: wt
       real(kind=8)                     :: fcz

!f2py intent(out)   pvar1
!f2py intent(in)   var1
!f2py intent(in)   pz
!f2py intent(in)   npmx
!f2py intent(in)   nz
!f2py intent(in)   np

       do ip = 1,np
         if (.not.(isnan(pz(ip))) ) then
           k = floor(pz(ip)+1)
           fcz = pz(ip)+1 - k 
           wt(1) = (1-fcz)
           wt(2) =    fcz
           if (k.ge.nz) then   
             pvar1(ip)   = var1(k)
           else
             pvar1(ip)   = sum(var1(k:k+1)*wt)
           endif
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


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Z interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# include "R_tools_fort_routines/sigma_to_z_intr_sfc.F"

# include "R_tools_fort_routines/sigma_to_z_intr_bot.F"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute w
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# include "R_tools_fort_routines/get_wvlcty.F"






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





