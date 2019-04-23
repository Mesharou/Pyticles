 
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

# include "R_tools_fort_routines/cppdefs.h"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Compute density anomaly (adapted from rho_eos.F)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# include "R_tools_fort_routines/rho_eos.F"

# include "R_tools_fort_routines/rho1_eos.F"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Compute Brunt-Vaissala frequency (adapted from rho_eos.F)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# include "R_tools_fort_routines/bvf_eos.F"

# include "R_tools_fort_routines/bvf_lineos.F"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Compute buoyancy (adapted from rho_eos.F)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# include "R_tools_fort_routines/get_buoy.F"


# include "R_tools_fort_routines/get_buoy1.F"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!compute neutral density gradients 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# include "R_tools_fort_routines/rho_grad.F"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!Z interpolation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# include "R_tools_fort_routines/sigma_to_z_intr_sfc.F"

# include "R_tools_fort_routines/sigma_to_z_intr_bot.F"

# include "R_tools_fort_routines/sigma_to_z_intr_bot_2d.F"

# include "R_tools_fort_routines/sigma_to_z_intr.F"

# include "R_tools_fort_routines/interp_1d.F"

# include "R_tools_fort_routines/sigma_to_z_intr_bounded.F"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute z_r and z_w for NEW_S_COORD
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

# include "R_tools_fort_routines/zlevs.F"

# include "R_tools_fort_routines/zlev.F"

# include "R_tools_fort_routines/zlevs_agrif.F"

# include "R_tools_fort_routines/zlevs_kau.F"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute Various
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# include "R_tools_fort_routines/get_rot.F"

# include "R_tools_fort_routines/get_grad.F"



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute vertical velocity
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# include "R_tools_fort_routines/get_wvlcty.F"

# include "R_tools_fort_routines/get_omega.F"


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Compute ROMS stuffs
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


# include "R_tools_fort_routines/get_swr_frac.F"

# include "R_tools_fort_routines_gula/get_ghat.F"

# include "R_tools_fort_routines/alfabeta.F"


# include "R_tools_fort_routines_gula/get_hbbl.F"

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




