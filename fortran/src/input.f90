!!
!! Luke McCulloch
!! Module to compute Points on a plate
!! geo.f90
!! Inputs Number of points in X, NI
!!       Number of points in Y, NJ
!!       Physical (RECTangular) Plate Dimentions, lx, ly
!! Output The Array, "Points"
!!

MODULE input

  use precise, only : defaultp

  IMPLICIT NONE
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp



CONTAINS



  SUBROUTINE input_discrete( inputfile, outputfile, title, &
                         L, xo, npts, pts, ncells, cells, &
                         nt, H, u, uH, HL, HRoHL, H_hat, u_hat,C )

    !This is a subroutine to read in the problem geometry

     INTEGER :: npts    ! number of cell boundaries
     INTEGER :: ncells  ! number of cell centers
     INTEGER :: nt      ! number of time steps
     
     CHARACTER(len=*) :: inputfile
     CHARACTER(len=*) :: outputfile
     CHARACTER(len=*) :: title

     LOGICAL :: flexists     ! a logical variable. I .true. or .false.

     REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) :: H      ! this is the memory intensive way...
     REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) :: uH     ! this is the memory intensive way...
     REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) :: u      ! this is the memory intensive way...
     REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: H_hat  ! this is the memory intensive way...
     REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: u_hat  ! this is the memory intensive way...
     REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: pts    ! cell boundary locations
     REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: cells  ! cell center locations
     REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: k      ! thermal conductivity == heat flux
     real(wp) :: HL                                    ! heat flux is constant for now.
     real(wp) :: L                                     ! length of the domain
     real(wp) :: xo                                    ! initial x position
     real(wp) :: rho_c                                 ! density times specific heat
     real(wp) :: HRoHL                                 ! part of Robin's B.C
     real(wp) :: C

     write(*,*) 'Start Subroutine input'

     OPEN(10,file=inputfile) ! open the file
     
     ! Read the data----------------------------------------------------
     read(10,'(AAAAA)') title
     read(10,*) L, xo
     read(10,*) npts, ncells, nt
     read(10,*) HL, HRoHL
     read(10,*) C

     !!phantom nodes
     ncells=ncells+2

     !! Allocate to read in geo file in higher dimensional versions...
     ALLOCATE( pts(1,npts) )
     ALLOCATE( cells(1,ncells) )
     !ALLOCATE( h(1,ncells,nt+1))
     !ALLOCATE( uH(1,ncells,nt+1))
     !ALLOCATE( u(1,ncells,nt+1))
     ALLOCATE( h_hat(1,npts))
     ALLOCATE( u_hat(1,npts))

     close(10)
     ! Finished Reading the data------------------------------------------

     ! Write the fifi data -----------------------------------------------
     OPEN(20,file=outputfile)
     !write(20,*)
     !write(20,*) '---------------------------------------------------------'
     !write(20,*) title
     !WRITE(20, FMT='(AI8)') 'ncells = ', ncells
     !WRITE(20, FMT='(AI8)') 'HR/HL =  ', HRoHL
     !WRITE(20, FMT='(f19.10)') HRoHL
     !write(20,*) '--------------------------------------------------------'
     
     write(*,*) 'Ending Subroutine input'

  END SUBROUTINE input_discrete



!!$  SUBROUTINE output( outputfile, iter, sourcesink, cw, profilewave )
!!$  
!!$     CHARACTER(len=*) :: outputfile
!!$     !real(wp), allocatable, dimension(:,:) :: op
!!$     real(wp) :: sourcesink, cw
!!$     integer :: iter
!!$     real(wp), dimension(2,128) :: profilewave
!!$
!!$     Open(20,file= outputfile, Access= 'append', status='old')
!!$     write(20,*)
!!$     Write(20,*) ' Results From the Wigley Hull Linear Free Surface Panel Code '
!!$     write(20,*) '---------------------------------------------------------------'
!!$     write(20,*) ' (2.) Output data using the subroutine "output" in module "io" '
!!$     write(20,*)


!!$     
!!$  END SUBROUTINE output
     
!! End of Subroutine to write the data ---------------------------------
!! End of Subroutine to write the data ---------------------------------

END MODULE input
