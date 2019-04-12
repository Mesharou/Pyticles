# 1 "particules_3d_sa.F"
# 1 "<built-in>"
# 1 "<command-line>"
# 1 "/usr/include/stdc-predef.h" 1 3 4
# 1 "<command-line>" 2
# 1 "particules_3d_sa.F"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!! cpp particules_3d_sa.F particules_3d_sa.f90 ; f2py --f90flags="-extend_source -O1 " -DF2PY_REPORT_ON_ARRAY_COPY=1 -c -m particules_3d_sa particules_3d_sa.f90
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




!----------------------------------------------------------------------------------------------
       subroutine seed(px,py,pz,npmx,mask,lev0,lev1,nnx,nny,nnlev,nx,ny,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny
       integer(kind=4) ,intent(in) :: nnx, nny, nnlev
       integer(kind=4) ,intent(in) :: np
       integer(kind=4) ,intent(in) :: lev0, lev1
       integer(kind=4),dimension(nx,ny),intent(in) :: mask
       real(kind=8) ,dimension(np) ,intent(inout):: px,py,pz
       integer(kind=4) ,intent(inout):: npmx
! local
       integer(kind=4) :: ip,i,j,id, k
       integer(kind=4),allocatable,dimension(:,:) :: grid

!f2py intent(inout) px,py,pz
!f2py intent(in) mask
!f2py intent(in) lev0,lev1
!f2py intent(inout) npmx
!f2py intent(in) nx,ny, nnx, nny, nnlev
!f2py intent(in) np

! print *,'fo: ',npmx
       allocate(grid(nx,ny))
       grid = 0


       !! find empty cells
       do ip = 1,npmx !! max particle index that is currently used
          i = floor(px(ip))
          j = floor(py(ip))
          grid(i,j) = 1
       enddo
!
! empty cells in the mask will get a new particle
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
! print *,'fo: ',npmx

*
       end


!----------------------------------------------------------------------------------------------
       subroutine advance_3d(px,py,pz,u0,v0,w0,u1,v1,w1,fct,dx,dy,dz,dt,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz
       integer(kind=4) ,intent(in) :: np
       integer(kind=4) ,intent(in) :: npmx
       integer(kind=4) ,intent(in) :: i0,j0,k0
       real(kind=8) ,dimension(nx,ny,nz),intent(in) :: u0,v0,w0
       real(kind=8) ,dimension(nx,ny,nz),intent(in) :: u1,v1,w1
       real(kind=8) ,dimension(np) ,intent(inout):: px,py,pz
       real(kind=8) ,intent(in) :: dx,dy,dz,dt,fct
! local
       integer(kind=8) :: ip,jp,kp,i,j,k
       real(kind=8) ,dimension(2,2,2) :: wt
       real(kind=8) :: fcx,fcy,fcz,dxi,dyi,dzi,pu,pv,pw


!f2py intent(in) u0,v0,w0
!f2py intent(in) u1,v1,w1
!f2py intent(inout) px,py,pz
!f2py intent(in) dx,dy,dz,dt
!f2py intent(in) npmx
!f2py intent(in) nx,ny,nz
!f2py intent(in) i0,j0,k0
!f2py intent(in) np

       dxi = 1/dx
       dyi = 1/dy
       dzi = 1/dz


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
           wt(1,1,2) = fcz *(1-fcy)*(1-fcx);
           wt(1,2,1) = (1-fcz)* fcy *(1-fcx);
           wt(1,2,2) = fcz * fcy *(1-fcx);
           wt(2,1,1) = (1-fcz)*(1-fcy)* fcx ;
           wt(2,1,2) = fcz *(1-fcy)* fcx ;
           wt(2,2,1) = (1-fcz)* fcy * fcx ;
           wt(2,2,2) = fcz * fcy * fcx ;

           if (k.ge.nz) then
           !The particule has reach the surface
            pu = (1-fct)*sum(u0(i:i+1,j:j+1,k)*wt(:,:,1))+fct*sum(u1(i:i+1,j:j+1,k)*wt(:,:,1))
            pv = (1-fct)*sum(v0(i:i+1,j:j+1,k)*wt(:,:,1))+fct*sum(v1(i:i+1,j:j+1,k)*wt(:,:,1))
            pw = min((1-fct)*sum(w0(i:i+1,j:j+1,k)*wt(:,:,1))+fct*sum(w1(i:i+1,j:j+1,k)*wt(:,:,1)),0.)
           else
            pu = (1-fct)*sum(u0(i:i+1,j:j+1,k:k+1)*wt)+fct*sum(u1(i:i+1,j:j+1,k:k+1)*wt)
            pv = (1-fct)*sum(v0(i:i+1,j:j+1,k:k+1)*wt)+fct*sum(v1(i:i+1,j:j+1,k:k+1)*wt)
            pw = (1-fct)*sum(w0(i:i+1,j:j+1,k:k+1)*wt)+fct*sum(w1(i:i+1,j:j+1,k:k+1)*wt)
           endif

           px(ip) = px(ip) + dt*pu*dxi
           py(ip) = py(ip) + dt*pv*dyi
           pz(ip) = pz(ip) + dt*pw*dzi


         endif

! if (pz(ip)- k0>nz-1) pz(ip)=nz-1-k0
! if ((px(ip)- i0>nx-1) .or. (px(ip)<0)) print*,'px',ip,px(ip),py(ip),pz(ip), pw,dt*pw*dzi
! if ((py(ip)- j0>ny-1) .or. (py(ip)<0)) print*,'py',ip,px(ip),py(ip),pz(ip), pw,dt*pw*dzi
! if ((pz(ip)- k0>nz-1) .or. (pz(ip)<0)) print*,'pz',ip,px(ip),py(ip),pz(ip), pw,dt*pw*dzi
! print *,ip, px(ip),py(ip),pz(ip)

         endif
       enddo
!
       end
!----------------------------------------------------------------------------------------------
       subroutine advance_3d_multi(px,py,pz,u,v,w,fct,dx,dy,dz,dt,npmx,i0,j0,k0,nx,ny,nz,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nx,ny,nz
       integer(kind=4) ,intent(in) :: np
       integer(kind=4) ,intent(in) :: npmx
       integer(kind=4) ,intent(in) :: i0,j0,k0
       real(kind=8) ,dimension(nx,ny,nz,2),intent(in) :: u,v,w
       real(kind=8) ,dimension(np) ,intent(inout):: px,py,pz
       real(kind=8) ,intent(in) :: dx,dy,dz,dt,fct
! local
       integer(kind=8) :: ip,jp,kp,i,j,k
       real(kind=8) ,dimension(2,2,2,2) :: wt
       real(kind=8) :: fcx,fcy,fcz,dxi,dyi,dzi,pu,pv,pw


!f2py intent(in) u,v,w
!f2py intent(inout) px,py,pz
!f2py intent(in) dx,dy,dz,dt
!f2py intent(in) npmx
!f2py intent(in) nx,ny,nz
!f2py intent(in) i0,j0,k0
!f2py intent(in) np

       dxi = 1/dx
       dyi = 1/dy
       dzi = 1/dz


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

           if (k.ge.nz) then
           !The particule has reach the surface
            pu = sum(u(i:i+1,j:j+1,k,:)*wt(:,:,1,:))
            pv = sum(v(i:i+1,j:j+1,k,:)*wt(:,:,1,:))
            pw = min(sum(w(i:i+1,j:j+1,k,:)*wt(:,:,1,:)),0.)
           else
            pu = sum(u(i:i+1,j:j+1,k:k+1,:)*wt)
            pv = sum(v(i:i+1,j:j+1,k:k+1,:)*wt)
            pw = sum(w(i:i+1,j:j+1,k:k+1,:)*wt)
           endif

           px(ip) = px(ip) + dt*pu*dxi
           py(ip) = py(ip) + dt*pv*dyi
           pz(ip) = pz(ip) + dt*pw*dzi


         endif

! if (pz(ip)- k0>nz-1) pz(ip)=nz-1-k0
! if ((px(ip)- i0>nx-1) .or. (px(ip)<0)) print*,'px',ip,px(ip),py(ip),pz(ip), pw,dt*pw*dzi
! if ((py(ip)- j0>ny-1) .or. (py(ip)<0)) print*,'py',ip,px(ip),py(ip),pz(ip), pw,dt*pw*dzi
! if ((pz(ip)- k0>nz-1) .or. (pz(ip)<0)) print*,'pz',ip,px(ip),py(ip),pz(ip), pw,dt*pw*dzi
! print *,ip, px(ip),py(ip),pz(ip)

         endif
       enddo
!
       end

!----------------------------------------------------------------------------------------------
       subroutine interp_3d(pvar1,pvar2,px,py,pz,var1,var2,npmx,i0,j0,k0,nx,ny,nz,np)
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
       real(kind=8) ,dimension(2,2,2) :: wt
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
           i = floor(px(ip)+1)-i0
           j = floor(py(ip)+1)-j0
           k = floor(pz(ip)+1)-k0
           fcx = px(ip)+1 - i - i0;
           fcy = py(ip)+1 - j - j0;
           fcz = pz(ip)+1 - k - k0;
           wt(1,1,1) = (1-fcz)*(1-fcy)*(1-fcx);
           wt(1,1,2) = fcz *(1-fcy)*(1-fcx);
           wt(1,2,1) = (1-fcz)* fcy *(1-fcx);
           wt(1,2,2) = fcz * fcy *(1-fcx);
           wt(2,1,1) = (1-fcz)*(1-fcy)* fcx ;
           wt(2,1,2) = fcz *(1-fcy)* fcx ;
           wt(2,2,1) = (1-fcz)* fcy * fcx ;
           wt(2,2,2) = fcz * fcy * fcx ;
           if (k.ge.nz) then
             pvar1(ip) = sum(var1(i:i+1,j:j+1,k)*wt(:,:,1))
             pvar2(ip) = sum(var2(i:i+1,j:j+1,k)*wt(:,:,1) )
           else
             pvar1(ip) = sum(var1(i:i+1,j:j+1,k:k+1)*wt)
             pvar2(ip) = sum(var2(i:i+1,j:j+1,k:k+1)*wt )
           endif
         endif
       enddo
!
       end




!----------------------------------------------------------------------------------------------
       subroutine oneterp_3d(pvar1,px,py,pz,var1,npmx,i0,j0,k0,nx,ny,nz,np)
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
       real(kind=8) ,dimension(2,2,2) :: wt
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
           i = floor(px(ip)+1)-i0
           j = floor(py(ip)+1)-j0
           k = floor(pz(ip)+1)-k0
           fcx = px(ip)+1 - i - i0;
           fcy = py(ip)+1 - j - j0;
           fcz = pz(ip)+1 - k - k0;
           wt(1,1,1) = (1-fcz)*(1-fcy)*(1-fcx);
           wt(1,1,2) = fcz *(1-fcy)*(1-fcx);
           wt(1,2,1) = (1-fcz)* fcy *(1-fcx);
           wt(1,2,2) = fcz * fcy *(1-fcx);
           wt(2,1,1) = (1-fcz)*(1-fcy)* fcx ;
           wt(2,1,2) = fcz *(1-fcy)* fcx ;
           wt(2,2,1) = (1-fcz)* fcy * fcx ;
           wt(2,2,2) = fcz * fcy * fcx ;
           if (k.ge.nz) then
             pvar1(ip) = sum(var1(i:i+1,j:j+1,k)*wt(:,:,1))
           else
             pvar1(ip) = sum(var1(i:i+1,j:j+1,k:k+1)*wt)
           endif
         endif
       enddo
!
       end


!----------------------------------------------------------------------------------------------
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
           i = floor(px(ip)+1)-i0
           j = floor(py(ip)+1)-j0
           fcx = px(ip)+1 - i - i0;
           fcy = py(ip)+1 - j - j0;
           wt(1,1) = (1-fcy)*(1-fcx);
           wt(1,2) = fcy *(1-fcx);
           wt(2,1) = (1-fcy)* fcx ;
           wt(2,2) = fcy * fcx ;
           pvar1(ip) = sum(var1(i:i+1,j:j+1)*wt)
         endif
       enddo
!
       end



!----------------------------------------------------------------------------------------------
       subroutine interp_1d(pvar1,pz,var1,npmx,nz,np)
       implicit none
! import/export
       integer(kind=4) ,intent(in) :: nz
       integer(kind=4) ,intent(in) :: np
       integer(kind=4) ,intent(in) :: npmx
       real(kind=8) ,dimension(nz),intent(in) :: var1
       real(kind=8) ,dimension(np) ,intent(out):: pvar1
       real(kind=8) ,dimension(np) ,intent(in):: pz
! local
       integer(kind=8) :: kp,k,ip
       real(kind=8) ,dimension(2) :: wt
       real(kind=8) :: fcz

!f2py intent(out) pvar1
!f2py intent(in) var1
!f2py intent(in) pz
!f2py intent(in) npmx
!f2py intent(in) nz
!f2py intent(in) np

       do ip = 1,np
         if (.not.(isnan(pz(ip))) ) then
           k = floor(pz(ip)+1)
           fcz = pz(ip)+1 - k
           wt(1) = (1-fcz)
           wt(2) = fcz
           if (k.ge.nz) then
             pvar1(ip) = var1(k)
           else
             pvar1(ip) = sum(var1(k:k+1)*wt)
           endif
         endif
       enddo


       end
# 466 "particules_3d_sa.F"
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
# 480 "particules_3d_sa.F" 2


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
# 487 "particules_3d_sa.F" 2

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
# 489 "particules_3d_sa.F" 2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute w
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# 1 "R_tools_fort_routines/get_wvlcty.F" 1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute vertical velocity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!






      subroutine get_wvlcty(Lm,Mm,N,u,v, z_r,z_w,pm,pn
     & ,Wvlc)


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


      real*8 Wrk(0:Lm+1,0:N),Wvlc(0:Lm+1,0:Mm+1,N),
     & Wxi(1:Lm+1,0:Mm+1),Weta(0:Lm+1,1:Mm+1)

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
# 36 "R_tools_fort_routines/get_wvlcty.F" 2




Cf2py intent(in) Lm,Mm,N, u,v,z_r,z_w,pm,pn
Cf2py intent(out) Wvlc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      imin=0
      imax=Lm+1
      jmin=0
      jmax=Mm+1

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



      do j=jmin,jmax-1

        do i=imin,imax
          Wrk(i,0)=0.
        enddo

        do k=1,N,+1
          do i=imin,imax-1
            Wrk(i,k)=Wrk(i,k-1)-pm(i,j)*pn(i,j)*(
     & FlxU(i+1,j,k)-FlxU(i,j,k)
     & +FlxV(i,j+1,k)-FlxV(i,j,k))
          enddo
        enddo

        do i=imin+1,imax
          Wvlc(i,j,N)=+0.375*Wrk(i,N) +0.75*Wrk(i,N-1)
     & -0.125*Wrk(i,N-2)
        enddo
        do k=N-1,2,-1
          do i=imin+1,imax
            Wvlc(i,j,k)=+0.5625*(Wrk(i,k )+Wrk(i,k-1))
     & -0.0625*(Wrk(i,k+1)+Wrk(i,k-2))
          enddo
        enddo
        do i=imin+1,imax
          Wvlc(i,j, 1)= -0.125*Wrk(i,2) +0.75*Wrk(i,1)
     & +0.375*Wrk(i,0)
        enddo
      enddo
!
! Compute and add contributions due to (quasi-)horizontal motions
! along S=const surfaces by multiplying horizontal velocity
! components by slops S-coordinate surfaces:
!
      do k=1,N
        do j=jmin,jmax
          do i=imin+1,imax
            Wxi(i,j)=u(i,j,k)*(z_r(i,j,k)-z_r(i-1,j,k))
     & *(pm(i,j)+pm(i-1,j))
          enddo
        enddo
        do j=jmin+1,jmax
          do i=imin,imax
            Weta(i,j)=v(i,j,k)*(z_r(i,j,k)-z_r(i,j-1,k))
     & *(pn(i,j)+pn(i,j-1))
          enddo
        enddo
        do j=jmin+1,jmax-1
          do i=imin+1,imax-1
            Wvlc(i,j,k)=Wvlc(i,j,k)+0.25*(Wxi(i,j)+Wxi(i+1,j)
     & +Weta(i,j)+Weta(i,j+1))
          enddo
        enddo
      enddo
# 174 "R_tools_fort_routines/get_wvlcty.F"
      return
      end
# 496 "particules_3d_sa.F" 2






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
