!!
!! Luke McCulloch
!! Module to compute Points on a plate
!! geo.f90
!! Inputs Number of points in X, NI
!!       Number of points in Y, NJ
!!       Physical (RECTangular) Plate Dimentions, lx, ly
!! Output The Array, "Points"
!!

MODULE geo

  use precise, only : defaultp
  use input

  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp
  !REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: pts, cells, t

CONTAINS

  Subroutine line(dx, L, xo, npts, pts, ncells, cells)
    ! Subroutine to Compute Points on a line
    
    INTEGER :: i       ! dummy index
    INTEGER :: npts    ! # of points
    INTEGER :: ncells  ! # of cells
    INTEGER :: nt      ! # of timesteps

    REAL(WP), ALLOCATABLE, DIMENSION(:,:) :: pts, cells
    real(wp) :: L, xo, dx, dt
    
    dx=L/(npts-1) !fence post!

    pts(1,1)=xo
    Do i = 2,npts
       pts(1,i)=pts(1,i-1)+dx
    End Do
    
    cells(1,1)=xo+(dx/2.)-dx
    Do i = 2,ncells
       cells(1,i)=cells(1,i-1)+dx
    End Do
    write(*,*) 'cells'
    write(*,'(1D16.8)') cells

  End Subroutine line




  Subroutine plate( NI, NJ, lx, ly, dx, dy, X, Y, Points )!!! Set lx and ly to 1 for initial tests
    ! Subroutine to Compute Points on a Plate
    
    INTEGER :: NI, NJ, i, j, ipt
    REAL(wp) :: lx, ly, dx, dy, Xo, Yo!, X, Y  ! Number of X points, Number of Y Points, Length of plate in x, Length of plate in Y, step x, step y
    
    REAL(wp), DIMENSION(2,NJ*NI) :: Points    ! Matrix of Points on the Plate
    REAL(wp), DIMENSION(NI,NJ) :: X  
    REAL(wp), DIMENSION(NI,NJ) :: Y      

    dx = lx/(NI-1)
    dy = ly/(NJ-1)

    Xo=0.
    Yo=0.

 !How do you want it?  Scalars or Vectors...
    ipt=1  ! point counter
    Do j = 1,NJ
       Do i = 1,NI
          X(i,j) = Xo+REAL(dx*i)-dx
          Y(i,j) = Yo+REAL(dy*j)-dy
          Points(:,ipt)=(/ X(i,j), Y(i,j) /)  ! Vector Version - Not used
          ipt=ipt+1
       End Do
       !ipt=ipt+1
    End Do



  END SUBROUTINE plate
END MODULE geo
  
