!
! Apply horizontal smoothing operator to wrk, while avoiding land-
! masked values which is accomplished by expressing everything in
! terms of elementary differences, subject to masking by U,V-rules.
! In the absense of masking the stencil of smoothing operator has
! the following weights, depending on coefficient settings in the
! code segment below:
!
! cff = 1/8, 1/4           cff = 1/12, 3/16         cff = 0, 1/8
!
!    1/16   1/8  1/16        1/32   1/8  1/32               1/8
!
!    1/8   1/4  1/8           1/8   3/8  1/8           1/8  1/4  1/8
!
!    1/16   1/8  1/16        1/32   1/8  1/32               1/8
!
!  2D 1-2-1-Hanning              isotropic                5-point
!   window smoother              Laplacian               Laplacian
!
! in all three cases the smoothing operator suppresses cheque-board
! mode in just one iteration; however, only the first one eliminates
! the 1D (flat-front) 2dx-modes in one iteration;   the two others
! attenuate 1D 2dx-mode by factors of 1/4 and 1/2 per iteration
! respectively.
!
      do j=jstr,jend
        do i=istr,iend
          FX(i,j)=(wrk(i,j)-wrk(i-1,j))
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          FE(i,j)=(wrk(i,j)-wrk(i,j-1))
        enddo
      enddo
      cff=1.D0/12.D0
      cff1=3.D0/16.D0
      do j=jstr,jend
        do i=istr,iend
          FE1(i,j)=FE(i,j)+cff*( FX(i+1,j)+FX(i,j-1)
     &                          -FX(i,j)-FX(i+1,j-1))
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          FX(i,j)=FX(i,j)+cff*( FE(i,j+1)+FE(i-1,j)
     &                         -FE(i,j)-FE(i-1,j+1))
        enddo
      enddo
      do j=jstr,jend
        do i=istr,iend
          wrk(i,j)=wrk(i,j)+cff1*( FX(i+1,j)-FX(i,j)
     &                          +FE1(i,j+1)-FE1(i,j))
        enddo            !--> discard FX,FE,FE1
      enddo

