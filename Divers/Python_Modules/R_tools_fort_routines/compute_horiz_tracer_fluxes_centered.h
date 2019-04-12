! This code segment computes horizontal fluxes for tracer variables.
! Basically it interpolates tracer values from their native locations
! on C grid to horizontal velocity points. Curently three options are
! supported: 4-point symmetric fourth-order method (default); 3-point
! upstream-biased parabolic interpolation (UPSTREAM); and 4-point
! scheme where arithmetic averaging of elementary differences is
! replaced by harmonic averaging (AKIMA), resulting in mid-point
! values bounded by two nearest values at native location, regardless
! of grid-scale roughness of the interpolated field, while still
! retaining asymptotic fourth-order behavior for smooth fields.
! This code is extracted into a special module ibecause it is used
! twice, in predictor and corrector substeps for tracer variables. 
! 




# undef  UPSTREAM
# undef AKIMA




c--# define CONST_TRACERS

# ifdef UPSTREAM
#  define curv WORK
# else
#  define grad WORK
# endif

          do j=jstr,jend
            do i=imin,imax
              FX(i,j)=(t(i,j,k,itrc)-t(i-1,j,k,itrc))
            enddo
          enddo            !--> discard imin,imax
# ifndef EW_PERIODIC
          if (WESTERN_EDGE) then
            do j=jstr,jend
              FX(istr-1,j)=FX(istr,j)
            enddo
          endif
          if (EASTERN_EDGE) then
            do j=jstr,jend
              FX(iend+2,j)=FX(iend+1,j)
            enddo
          endif
# endif
          do j=jstr,jend
            do i=istr-1,iend+1
# if defined UPSTREAM
              curv(i,j)=FX(i+1,j)-FX(i,j)
# elif defined AKIMA
              cff=2.*FX(i+1,j)*FX(i,j)
              if (cff.gt.epsil) then
                grad(i,j)=cff/(FX(i+1,j)+FX(i,j))
              else
                grad(i,j)=0.
              endif
# else
              grad(i,j)=0.5*(FX(i+1,j)+FX(i,j))
# endif
            enddo
          enddo             !--> discard FX
          do j=jstr,jend
            do i=istr,iend+1
# ifdef UPSTREAM
              FX(i,j)=0.5*(t(i,j,k,itrc)+t(i-1,j,k,itrc))
     &                                                  *FlxU(i,j,k)
     &            -0.166666666666*( curv(i-1,j)*max(FlxU(i,j,k),0.)
     &                             +curv(i  ,j)*min(FlxU(i,j,k),0.))
# else
              FX(i,j)=0.5*( t(i,j,k,itrc)+t(i-1,j,k,itrc)
     &                      -0.333333333333*( grad(i,j)-grad(i-1,j))
     &                                                 )*FlxU(i,j,k)
# endif
            enddo           !--> discard curv,grad, keep FX
          enddo
 

          do j=jmin,jmax
            do i=istr,iend
              FE(i,j)=(t(i,j,k,itrc)-t(i,j-1,k,itrc))
            enddo
          enddo         !--> discard jmin,jmax

          do j=jstr-1,jend+1
            do i=istr,iend
# if defined UPSTREAM
              curv(i,j)=FE(i,j+1)-FE(i,j)
# elif defined AKIMA
              cff=2.*FE(i,j+1)*FE(i,j)
              if (cff.gt.epsil) then
                grad(i,j)=cff/(FE(i,j+1)+FE(i,j))
              else
                grad(i,j)=0.
              endif
# else
              grad(i,j)=0.5*(FE(i,j+1)+FE(i,j))
# endif
            enddo
          enddo            !--> discard FE
 
          do j=jstr,jend+1
            do i=istr,iend
# ifdef UPSTREAM
              FE(i,j)=0.5*(t(i,j,k,itrc)+t(i,j-1,k,itrc))
     &                                                  *FlxV(i,j,k)
     &            -0.166666666666*( curv(i,j-1)*max(FlxV(i,j,k),0.)
     &                             +curv(i,j  )*min(FlxV(i,j,k),0.))
# else
              FE(i,j)=0.5*( t(i,j,k,itrc)+t(i,j-1,k,itrc)
     &                       -0.333333333333*(grad(i,j)-grad(i,j-1))
     &                                                 )*FlxV(i,j,k)
# endif
            enddo
          enddo            !--> discard curv,grad, keep FE
