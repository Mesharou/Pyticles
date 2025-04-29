
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      subroutine get_tracer_evolution (Lm,Mm,N,u,v, z_r,z_w
     &  ,pm,pn,dt
     &  ,t,stflx, srflx, ghat, swr_frac, W, Akt
     &  ,TXadv,TYadv,TVadv,THdiff,TVmix,TForc)    


      integer Lm,Mm,N,NT, i,j,k
     &      ,istr,iend,jstr,jend
     &      ,imin,imax,jmin,jmax
     &      ,itemp,isalt,dt,itrc

      parameter (NT=2)  
      parameter (itemp=1,isalt=2)  
      parameter (epsil=1.E-16)

      ! INPUTS
      real*8 t(0:Lm+1,0:Mm+1,N,NT)


      real*8 stflx(0:Lm+1,0:Mm+1,NT)
     &      ,srflx(0:Lm+1,0:Mm+1)
     &      ,ghat(0:Lm+1,0:Mm+1,N)

      real*8 u(1:Lm+1,0:Mm+1,N),  v(0:Lm+1,1:Mm+1,N)
     &      ,z_r(0:Lm+1,0:Mm+1,N),  z_w(0:Lm+1,0:Mm+1,0:N)
     &      ,pm(0:Lm+1,0:Mm+1), pn(0:Lm+1,0:Mm+1)
     &      ,W(0:Lm+1,0:Mm+1,0:N)
     &      ,Akt(0:Lm+1,0:Mm+1,0:N)

      ! OUTPUTS
      real*8 TXadv(0:Lm+1,0:Mm+1,N,NT)
     &      ,TYadv(0:Lm+1,0:Mm+1,N,NT)
     &      ,TVadv(0:Lm+1,0:Mm+1,N,NT)
     &      ,THdiff(0:Lm+1,0:Mm+1,N,NT)
     &      ,TVmix(0:Lm+1,0:Mm+1,N,NT)
     &      ,TForc(0:Lm+1,0:Mm+1,N,NT)


      ! LOCAL
      real*8 wrk1(0:Lm+1,0:Mm+1), wrk2(0:Lm+1,0:Mm+1)
     &      ,FX(0:Lm+1,0:Mm+1), FE(0:Lm+1,0:Mm+1)
     &      ,WORK(0:Lm+1,0:Mm+1)
     &      ,FC(0:Lm+1,0:N), DC(0:Lm+1,0:N)
     &      ,CF(0:Lm+1,0:N)
     &      ,Hz(0:Lm+1,0:Mm+1,N)
     &      ,dn_u(0:Lm+1,0:Mm+1), dm_v(0:Lm+1,0:Mm+1)
     &      ,FlxU(1:Lm+1,0:Mm+1,N),  FlxV(0:Lm+1,1:Mm+1,N)
     &      ,swr_frac(0:Lm+1,0:Mm+1,0:N)
     &      ,tnew(0:Lm+1,0:Mm+1,N,NT)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


Cf2py intent(in) Lm,Mm,N, u,v,z_r,z_w,pm,pn,dt,t,stflx, srflx, ghat, swr_frac,W, Akt
Cf2py intent(out) TXadv,TYadv,TVadv,THdiff,TVmix,TForc

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        istr=1
        iend=Lm
        jstr=1
        jend=Mm


        imin=0
        imax=Lm+1
        jmin=0
        jmax=Mm+1



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      do j=jmin,jmax
        do i=imin,imax
          do k=0,N,+1
              W(i,j,k)=W(i,j,k)/(pm(i,j)*pn(i,j))
            enddo
          enddo
        enddo


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
     &      * u(i,j,k) 
            enddo
          enddo
      enddo



      do j=jmin+1,jmax
        do i=imin,imax
            dm_v(i,j) = 2./(pm(i,j)+pm(i,j-1))
            do k=1,N,+1
              FlxV(i,j,k) = 0.5*(Hz(i,j,k)+Hz(i,j-1,k))*dm_v(i,j)
     &      * v(i,j,k) 
            enddo
          enddo
      enddo



      !!call get_swr_frac (Lm,Mm,N, Hz, swr_frac )


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! This code is valid only if the ROMS simulations has been run with UPSTREAM scheme
! So it is not valid if ADV_ISONEUTRAL was used
!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First centered scheme
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     


      do itrc=1,NT
        do k=1,N
 
# include "compute_horiz_tracer_fluxes_centered.h"
 
          do j=jstr,jend
            do i=istr,iend

              THdiff(i,j,k,itrc) = FX(i+1,j)-FX(i,j)
     &                        + FE(i,j+1)-FE(i,j)

            enddo
          enddo                 !-->  discard FX,FE
        enddo
      enddo




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! then using UPSTREAM scheme
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


      do itrc=1,NT
        do k=1,N
 
# include "compute_horiz_tracer_fluxes_upstream.h"
 
          do j=jstr,jend
            do i=istr,iend
              tnew(i,j,k,itrc)=Hz(i,j,k)*t(i,j,k,itrc)
     &                     -dt*pm(i,j)*pn(i,j)*( FX(i+1,j)-FX(i,j)
     &                                          +FE(i,j+1)-FE(i,j)
     &                                                           )
              TXadv(i,j,k,itrc) = -(FX(i+1,j)-FX(i,j))
              TYadv(i,j,k,itrc) = -(FE(i,j+1)-FE(i,j))


              THdiff(i,j,k,itrc) = THdiff(i,j,k,itrc)
     &                      -( FX(i+1,j)-FX(i,j)
     &                        + FE(i,j+1)-FE(i,j)
     &                                          )
            enddo
          enddo                 !-->  discard FX,FE
        enddo
      enddo
      


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      do j=jstr,jend
        do itrc=1,NT
 

# include "compute_vert_tracer_fluxes.h"
 
          do k=1,N            ! Apply vertical advective fluxes.
            do i=istr,iend

               tnew(i,j,k,itrc)=tnew(i,j,k,itrc)-dt*pm(i,j)
     &                             *pn(i,j)*(FC(i,k)-FC(i,k-1))

              TVadv(i,j,k,itrc) = -(FC(i,k)-FC(i,k-1))

            enddo
          enddo               !--> discard FC

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



!
! Add surface and bottom fluxes
!
#ifdef SEA_ICE_NOFLUX
          do i=istr,iend
            cff = 1.
            if (itrc.eq.itemp .and.
     &          (tnew(i,j,N,itrc) .le. -1.8)  ) cff = 0.

             tnew(i,j,N,itrc)=tnew(i,j,N,itrc)+dt*stflx(i,j,itrc)*cff

            TForc(i,j,N,itrc)= stflx(i,j,itrc)*cff/(pm(i,j)*pn(i,j))
            do k=N-1,1,-1
                TForc(i,j,k,itrc)= 0.
            enddo

          enddo
#else
!
! Add surface and bottom fluxes
!
          do i=istr,iend
            tnew(i,j,N,itrc)=tnew(i,j,N,itrc)+dt*stflx(i,j,itrc)

            TForc(i,j,N,itrc)= stflx(i,j,itrc)/(pm(i,j)*pn(i,j))
            do k=N-1,1,-1
                TForc(i,j,k,itrc)= 0.
            enddo

          enddo
#endif
  

        !write(*,*) 'sans deconner...',j,jend



# ifdef LMD_KPP
!
! Add the solar radiation flux in temperature equation. Also compute
! the nonlocal transport flux for unstable (convective) forcing
! conditions into matrix DC when using the Large et al. 1994 KPP
! scheme.
!


          if (itrc.eq.itemp) then
            do k=N-1,1,-1
              do i=istr,iend
                cff=srflx(i,j)*swr_frac(i,j,k)
#  ifdef LMD_NONLOCAL
     &                 -ghat(i,j,k)*(stflx(i,j,itemp)-srflx(i,j))
#  endif
                 tnew(i,j,k+1,itemp)=tnew(i,j,k+1,itemp) -dt*cff
                 tnew(i,j,k  ,itemp)=tnew(i,j,k  ,itemp) +dt*cff
                TForc(i,j,k+1,itrc)= TForc(i,j,k+1,itrc)
     &                 -cff/(pm(i,j)*pn(i,j))
                TForc(i,j,k,itrc)= TForc(i,j,k,itrc)
     &                 +cff/(pm(i,j)*pn(i,j))
              enddo
            enddo




c??
c?? WARNING: the above implies that light (if any) reaching all the
c?? way to the bottom is entirely absorbed within the botom-most grid
c?? box, rather that reflected/scatered back to the water column. It 
c?? is not clear, however, how to deal with this situation correctly
c??
c??         do i=istr,iend
c??           cff=srflx(i,j)*swr_frac(i,j,0)
c??#  ifdef LMD_NONLOCAL
c??  &               -ghat(i,j,0)*(stflx(i,j,itemp)-srflx(i,j))
c??#  endif
c??           t(i,j,1,itemp)=t(i,j,1,itemp) -dt*cff
c??         enddo


#  if defined LMD_NONLOCAL && defined SALINITY
          elseif (itrc.eq.isalt) then
            do k=N-1,1,-1
              do i=istr,iend
                cff=-dt*ghat(i,j,k)*stflx(i,j,isalt)
                 tnew(i,j,k+1,isalt)=tnew(i,j,k+1,isalt) -cff
                 tnew(i,j,k  ,isalt)=tnew(i,j,k  ,isalt) +cff
                TForc(i,j,k+1,itrc)= TForc(i,j,k+1,itrc)
     &                 -cff/(dt*pm(i,j)*pn(i,j))       
                TForc(i,j,k,itrc)= TForc(i,j,k,itrc)
     &                 +cff/(dt*pm(i,j)*pn(i,j))

              enddo
            enddo


#  endif
          endif
# endif


        !write(*,*) 'Forc OK',j,jend

! Perform implicit time step for vertical diffusion,
!
!   dq(k)     1     [         q(k+1)-q(k)             q(k)-q(k-1) ]
!  ------ = ----- * [ Akt(k)* ----------- - Akt(k-1)* ----------- ]
!    dt     Hz(k)   [            dz(k)                   dz(k-1)  ]
!
! where q(k) represents tracer field t(:,:,k,:,itrc). Doing so
! implies solution of a tri-diagonal system
!
!     -FC(k-1)*q_new(k-1) +[Hz(k)+FC(k-1)+FC(k)]*q_new(k)
!                       -FC(k)*q_new(k+1) = Hz(k)*q_old(k)
!
!                dt*Akt(k)
! where FC(k) = ----------- is normalized diffusivity coefficient 
!                  dz(k)
!
! defined at W-points; q_new(k) is the new-time-step (unknown) tracer
! field; q_old(k) is old-time-step tracer (known).  As long as
! vertical diffusivity Akt(k) is nonnegative, the tri-diagonal matrix
! is diagonally dominant which guarantees stability of a Gaussian
! elimination procedure, (e.g., Richtmeyer annd  Morton, 1967).
! Top and bottom boundary conditions are assumed to be no-flux,
! effectively Akt(N)=Akt(0)=0, hence FC(N)=FC(1)=0. This leads to
! equations for top and bottom grid boxes; 
!
!   -FC(N-1)*q_new(N-1) +[Hz(N)+FC(N-1)]*q_new(N) = Hz(N)*q_old(N)
!
!          [Hz(1)+FC(1)]*q_new(1) -FC(1)*q_new(2) = Hz(1)*q_old(1)
!
! The FC(N)=FC(0)=0 boundary conditions does not mean that physical
! boundary conditions are no flux: the forcing fluxes have been
! applied explicitly above.  Instead, the no-flux condition should
! be interpreted as that the implicit step merely redistributes the
! tracer concentration throughout the water column. At this moment
! the content of array t(:,:,:,itrc) has meaning of Hz*tracer.
! After the implicit step it becomes just tracer. 
!
          do k=1,N
            do i=istr,iend
               TVmix(i, j, k, itrc)=tnew(i,j,k,itrc)
            enddo
         enddo


          do i=istr,iend
            FC(i,1)=dt*Akt(i,j,1)/(z_r(i,j,2)-z_r(i,j,1))
            cff=1./(Hz(i,j,1)+FC(i,1))
            CF(i,1)=cff*FC(i,1)
            DC(i,1)=cff*tnew(i,j,1,itrc)
          enddo
          do k=2,N-1,+1
            do i=istr,iend
              FC(i,k)=dt*Akt(i,j,k)/(z_r(i,j,k+1)-z_r(i,j,k))
              cff=1./( Hz(i,j,k) +FC(i,k)+FC(i,k-1)*(1.-CF(i,k-1)) )
              CF(i,k)=cff*FC(i,k)
              DC(i,k)=cff*(tnew(i,j,k,itrc)+FC(i,k-1)*DC(i,k-1))
            enddo
          enddo

          do i=istr,iend
             tnew(i,j,N,itrc)=( tnew(i,j,N,itrc) +FC(i,N-1)
     &         *DC(i,N-1) )/(Hz(i,j,N)+FC(i,N-1)*(1.-CF(i,N-1)))
          enddo
          do k=N-1,1,-1
            do i=istr,iend
              tnew(i,j,k,itrc)=DC(i,k)+CF(i,k)*tnew(i,j,k+1,itrc)
            enddo
          enddo           !--> discard FC,CF,DC

          do k=1,N
            do i=istr,iend
              TVmix(i,j,k,itrc) = 
     &            -(TVmix(i,j,k,itrc)-tnew(i,j,k,itrc)*Hz(i,j,k))
     &                                        /(dt*pm(i,j)*pn(i,j))
            enddo
          enddo    




       enddo      ! <-- itrc 

       !write(*,*) 'the end..',j,jend


      enddo      ! <-- j


!
! Set lateral boundary conditions; nudge toward tracer climatology;
! apply land-sea mask and exchange periodic boundary conditions.
!
      do itrc=1,NT

!---------------------------------------------------------------
! Compute the tendency term of tracer diagnostics
! Divide all diagnostic terms by the cell volume
! (Hz(i,j,k,itrc)/(pm(i,j).*pn(i,j)). There after the unit
! of diagnostic terms will be: (unit of tracers)* s-1.
!
! Note: the Horizontal mixing term is computed in t3dmix
! where Trate is updated accordingly
!---------------------------------------------------------------

       do k=1,N
         do j=jstr,jend
           do i=istr,iend


              cff=pm(i,j)*pn(i,j)/Hz(i,j,k)
              TXadv(i,j,k,itrc)=TXadv(i,j,k,itrc)*cff
              TYadv(i,j,k,itrc)=TYadv(i,j,k,itrc)*cff
              TVadv(i,j,k,itrc)=TVadv(i,j,k,itrc)*cff
              TVmix(i,j,k,itrc)=TVmix(i,j,k,itrc)*cff
              THdiff(i,j,k,itrc)=THdiff(i,j,k,itrc)*cff
              TForc(i,j,k,itrc)=TForc(i,j,k,itrc)*cff


           enddo
         enddo

       !write(*,*) 'the end..',k,N

       enddo



      enddo  ! <-- itrc

!---------------------------------------------------------------

       !write(*,*) 'cooooool'


      return
      end




















