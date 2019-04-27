!! 
!! Module to specify initial conditions in a flow
!! September 2012
!! TLM
!!
Module initial_conditions

  use precise, only : defaultp
  IMPLICIT NONE
  
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp

Contains

  Subroutine InitialConditions1D(npts, ncells, pts, cells,  HL, HRoHL, H, u, uH)
    ! Subroutine to Compute Points on a line
    
    INTEGER :: npts, ncells, i, n

    REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: pts    ! 1 x npts array of points
    REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: cells  ! 1 x npts array of cells
    REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) :: H      ! 1 x npts x nt array cell H 
    REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) :: u      ! 1 x npts x nt array cell u velocities
    REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) :: uH     ! 1 x npts x nt array cell uH 
    real(wp) :: HL, HR, HRoHL, tol


    tol=(1.E-10)
    
    HR = HL*HRoHL
    Write (*,*) '============================================================================'
    write(*,*)'                     Establishing initial conditions:'
    write(*,*)'ncells = ',ncells
   

    Do i = 1,ncells
       u(1,i,1) = 0.
       if (cells(1,i)<0.) then
          H(1,i,1) = HL
       else
          H(1,i,1) = HR
       end if
       uH(1,i,1) = u(1,i,1)*H(1,i,1)

    End Do

    do i=1,ncells
       write(*,*)'i=',i,' H = ',H(1,i,1)
    end do
    print*
    do i=1,ncells
       write(*,*)'i=',i,' uH = ',uH(1,i,1)
    end do
    print*
    do i=1,ncells
       write(*,*)'i=',i,' u = ',u(1,i,1)
    end do
    print*
   

    

  End Subroutine InitialConditions1D

End Module initial_conditions
