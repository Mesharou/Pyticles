! Auxiliary module "compute_tile_bounds.h":
!---------- ------ ------------------------
! Compute bounds designed to cover interior points of an array
! for shared memory subdomain partitioning (tiling.)
!
! input: tile -- usually from 0 to NSUB_X*NSUB_E-1 -- indicates
!                the specified subdomain;  tile=NSUB_X*NSUB_E
!                corresponds to the whole domain of RHO points
!                treated as a single block.
! outputs: istr,iend -- starting and ending indices of subdomain
!          jstr,jend    tile in XI- and ETA-directions.
!
      integer istr,iend, jstr,jend

        istr=1
        iend=Lm
        jstr=1
        jend=Mm

 
