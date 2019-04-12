 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! ROMS ROUTINES
!!
!! copied from actual ROMS scripts
!!
!! compile with:
!! "cpp R_tools_fort.F R_tools_fort.f"
!! "f2py -DF2PY_REPORT_ON_ARRAY_COPY=1 -c -m R_tools_fort R_tools_fort.f" for python use
!!
!! print R_tools_fort.rho_eos.__doc__
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




! Included are:
!
!   subroutine rho_eos(Lm,Mm,N, T,S, z_r,z_w,rho0, 
!     &                    rho1,qp1,rho,bvf)
!   subroutine rho_grad(Lm,Mm,N, T,S, z_r,z_w,rho0,pm,pn,
!     &                    rho1,qp1,drdz,drdx,drdy)
!   subroutine sigma_to_z_intr (Lm,Mm,N, nz, z_r, z_w, rmask, var,
!     &                    z_lev, var_zlv, imin,jmin,kmin, FillValue)
!   subroutine zlevs(Lm,Mm,N, h,zeta, hc, Cs_r, Cs_w,z_r,z_w)
!    
!

# include "R_tools_fort_routines_gula/cppdefs.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Compute density anomaly (adapted from rho_eos.F)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# include "R_tools_fort_routines_gula/rho_eos_nozw.F"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# include "R_tools_fort_routines/interp_1d.F"



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute PV
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

c---# include "R_tools_fort_routines_gula/get_diagsPV_sol1.F"

c---# include "R_tools_fort_routines_gula/get_diagsPV_sol2.F"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute PV fluxes
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



# include "R_tools_fort_routines_gula/old/get_J1_sol1.F"
# include "R_tools_fort_routines_gula/old/get_J2_sol1.F"
# include "R_tools_fort_routines_gula/old/get_Jbot_sol1.F"

# include "R_tools_fort_routines_gula/get_Jbot_fromdiags.F"
# include "R_tools_fort_routines_gula/get_Jdiab_fromdiags.F"

# include "R_tools_fort_routines_gula/old/get_J1_sol2.F"
# include "R_tools_fort_routines_gula/old/get_J2_sol2.F"
# include "R_tools_fort_routines_gula/old/get_Jbot_sol2.F"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute Barotropic equation components
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# include "R_tools_fort_routines_gula/get_bot.F"

# include "R_tools_fort_routines_gula/get_bot_croco.F"

# include "R_tools_fort_routines_gula/get_bpt.F"

# include "R_tools_fort_routines_gula/get_u_prsgrd.F"

# include "R_tools_fort_routines_gula/get_vortplanet.F"

# include "R_tools_fort_routines_gula/get_vortstretch.F"

# include "R_tools_fort_routines_gula/get_intvortplanet.F"

# include "R_tools_fort_routines_gula/get_vortstretch_sol2.F"

# include "R_tools_fort_routines_gula/get_vortplantot.F"

# include "R_tools_fort_routines_gula/get_vortplantot_sol2.F"

# include "R_tools_fort_routines_gula/get_vortplantot_sol2_test.F"

# include "R_tools_fort_routines_gula/get_adv_sol2.F"



# include "R_tools_fort_routines_gula/get_uvgrid.F"

# include "R_tools_fort_routines_gula/get_adv_mix.F"

# include "R_tools_fort_routines_gula/get_adv_4th.F"


# include "R_tools_fort_routines_gula/get_fwb.F"

# include "R_tools_fort_routines_gula/get_fws.F"

# include "R_tools_fort_routines_gula/get_fdivub.F"

# include "R_tools_fort_routines_gula/get_fwdivub.F"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute Barotropic equation components for vertically averaged flow
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# include "R_tools_fort_routines_gula/get_bpt_mean.F"

# include "R_tools_fort_routines_gula/get_vortplantot_sol2_mean.F"

# include "R_tools_fort_routines_gula/get_adv_sol2_mean.F"

# include "R_tools_fort_routines_gula/get_adv_sol3_mean.F"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute stuff
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# include "R_tools_fort_routines_gula/get_rotv.F"

# include "R_tools_fort_routines_gula/get_hbbls_from_AKt.F"

# include "R_tools_fort_routines_gula/get_absvrt.F"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# include "R_tools_fort_routines_gula/get_tracer_evolution.F"

# include "R_tools_fort_routines_gula/get_uv_evolution.F"

# include "R_tools_fort_routines_gula/get_uv_evolution_old.F"

# include "R_tools_fort_routines/get_hbl.F"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# include "R_tools_fort_routines/get_swr_frac.F"

# include "R_tools_fort_routines_gula/get_ghat.F"

# include "R_tools_fort_routines/alfabeta.F"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# include "R_tools_fort_routines_gula/get_kediss.F"

# include "R_tools_fort_routines_gula/get_kediss_2d.F"

# include "R_tools_fort_routines_gula/visc3d_GP.F"

# include "R_tools_fort_routines_gula/visc3d_S.F"
# include "R_tools_fort_routines_gula/visc3d_S_baham.F"

# include "R_tools_fort_routines_gula/get_akv.F"

# include "R_tools_fort_routines_gula/lmd_vmix.F"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# include "R_tools_fort_routines_gula/solve_ttw.F"

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












