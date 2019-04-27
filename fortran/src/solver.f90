!!
!! Luke McCulloch
!! solver.f90
!! 11-4-2012
!! Serial Approximate Riemann Solver




PROGRAM solver

  use precise, only : defaultp   ! Module to handle precision
  use constants                  ! Module for pi
  use input                      ! io module
  use geo                        ! geometry
  use initial_conditions         ! initial conditions
  use inv                        ! my implementation of the Thomas Algorithm
  use analytic                   ! exact solution

  IMPLICIT NONE  

  integer, parameter :: WP=defaultp
  real, parameter :: e=2.718281828459045
  INTEGER, DIMENSION(3), PARAMETER :: it1 = (/ 40, 80, 120 /)   ! output times

  INTEGER :: narg    ! # of command line arguments
  INTEGER :: i       ! dummy index in space
  INTEGER :: j       ! dummy index in space
  INTEGER :: n       ! dummy index in time
  INTEGER :: ex_n    ! seperate time index for exact sol.
  INTEGER :: npts    ! # of points
  INTEGER :: ncells  ! # of cells
  INTEGER :: nt      ! # of timesteps
  INTEGER :: t_1sec  ! index at approx 1 sec
  INTEGER :: waveI   ! index location of the max Left Eigenvalue

  ! Efficient Placement of these declarations...?
  CHARACTER(len=24) :: inputfile
  CHARACTER(len=24) :: outputfile
  CHARACTER(len=30) :: title
  CHARACTER(len=30) :: title2
  CHARACTER(len=30) :: scheme_type

  LOGICAL :: flexists     ! a logical variable. I .true. or .false.
  LOGICAL :: godunov      ! true for theta == integer 1 else theta is real

  REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: pts          ! 1 x npts array of points
  REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: cells        ! 1 x ncells = npts+1  array of cells

  REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: lamda        ! 2 x ncells = vector of eigenvalues
  REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: F            ! 2 x ncells = Flux vector for each cell
  REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: Q            ! 2 x ncells = storage vector for current H,uH
  REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: alpha        ! 1 x ncells = npts+1  array of cells

  REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) :: H            ! 1 x ncells array cell
  REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) :: u            ! 1 x ncells array cell 
  REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: H_hat        ! 1 x npts face centered
  REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: u_hat        ! 1 x npts face centered

  REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) :: uH           ! 1 x ncells array cell 
  REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) :: H_exact      ! 1 x ncells array cell 



  REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) :: error        ! 1 x ncells array cell temperatures
  REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: k            ! 1 x npts array of heat fluxes 
  REAL(WP), ALLOCATABLE, DIMENSION(:)     :: am           ! sub   - diagonal
  REAL(WP), ALLOCATABLE, DIMENSION(:)     :: bm           ! main  - diagonal
  REAL(WP), ALLOCATABLE, DIMENSION(:)     :: cm           ! super - diagonal
  REAL(WP), ALLOCATABLE, DIMENSION(:)     :: dm           ! implicit RHS

  
  real(wp) :: g        ! gravitational acceleration
  real(wp) :: L        ! length
  real(wp) :: xo       ! starting position
  real(wp) :: dx       ! cell width
  real(wp) :: dt       ! timestep
  real(wp) :: dt_max   ! timestep max
  real(wp) :: HL       ! thermal conductivity cnst
  real(wp) :: HRoHL    ! part of Robin's B.C
  real(wp) :: C        ! Courant Number
  real(wp) :: ev_max   ! max eigenvale abs value from analytic solution
  real(wp) :: time     ! actual time
  real(wp) :: svtm     ! saved time
  real(wp) :: tol      ! tolerance value
  real(wp) :: tout     ! print at this time
  real(wp) :: waveloc  ! x location of max left eigenvalue

  real(wp) :: store1p  ! store a variable
  real(wp) :: store2p  ! store a variable
  real(wp) :: store1n  ! store a variable
  real(wp) :: store2n  ! store a variable


  real(wp) :: dummy1   ! dummy vbl
  real(wp) :: dummy2   ! dummy vbl
  real(wp) :: dummy3   ! dummy vbl
  real(wp) :: dummy4   ! dummy vbl
  real(wp) :: dummy5   ! dummy vbl
  real(wp) :: dummy6   ! dummy vbl  


!!$  real(wp) :: rho_c    ! density times specific heat
!!$  real(wp) :: a        ! part of Robin's B.C.
!!$  real(wp) :: b        ! part of Robin's B.C. 
!!$  real(wp) :: t_s      ! part of Robin's B.C.
!!$  real(wp) :: cnst1    ! placeholder for calculations
!!$  real(wp) :: cnst2    ! placeholder for calculations
!!$  real(wp) :: delta_q  ! change
!!$  real(wp) :: a1       ! i-1 term
!!$  real(wp) :: b1       ! i   term
!!$  real(wp) :: c1       ! i+1 term
!!$  real(wp) :: theta    ! scheme variable
!!$  real(wp) :: RHS      ! explicit Physics
!!$  real(wp) :: Su       ! source used for explicit B.C. -> bk pg. 250 ex. 8.1
!!$  real(wp) :: ae       ! book shorthand
!!$  real(wp) :: aw       ! book shorthand
!!$  real(wp) :: apo      ! book shorthand
  
  real(wp) :: start
  real(wp) :: finish




  !! Command Line argument input stuff... zzz --zz---z--------------------------------------!
  write(*,*) ''
  WRITE(6,'(A)') 'Project 1, v 1.0, Sept 2012, Main Program, Luke McCulloch '
  narg = command_argument_count()
  WRITE(6, FMT='(AI3A)') 'we have ', narg, ' command line arguments '
  IF ( narg > 1 ) THEN

     CALL get_command_argument(1, inputfile)

     CALL get_command_argument(2, outputfile)
     
     !CALL get_command_argument(3, scheme_type)

     

  ELSE
     WRITE(6,'(A)') ' Input file or Output file missing!'
     WRITE(6,'(A)') ' Usage:>  ./test <inputfile> <outputfile>' ! <scheme_type> '
     STOP
  ENDIF
  INQUIRE(file=inputfile, exist=flexists)
  !! Parse the input & output files
  IF (flexists) THEN
     !! Startup:
     WRITE(6,'(AA)') ' input file,  ', inputfile
     WRITE(6,*) 'output file, ', outputfile  
     !WRITE(6,*) 'scheme type, ', scheme_type
     WRITE(6,'(A)') 'ENME 6728 Unsteady Diffusion Scheme'
     Write(6,*)     ' 8-25-2012 Implementation By Luke McCulloch'
     write(6,*) ''
     Write(6,*)     ' Soving the Shallow Water Equations in 1D'
     write(6,*) '-----------------------------------------------------'
     call input_discrete( inputfile, outputfile, title, &
                         L, xo, npts, pts, ncells, cells, &
                         nt, H, u, uH, HL, HRoHL, H_hat, u_hat,C)


     write(*,*)'HL=',HL
     write(*,*)'HR/HL=',HRoHL
     ! Parse the Scheme Type
     godunov = .True.

  ENDIF
  !! End command line argument stuff---------------------------------------------------------!

  

  write(*,'(A4,D16.4,A8,I4,A10,I4,A12,I4,A12,D16.4)') '  L= ',L,'  npts= ', npts,'ncells= ',ncells
  write(*,'(A12,I4,A12,D16.4)') ' nt= ',nt
  
  

  write(*,*) ''
  call line(dx, L, xo, npts, pts, ncells, cells)
  write(6,'(A4,D16.8)') '  dx =', dx
  write(*,*) ''
  write(*,*) '  Face Locations:'
  write(*,'(1D16.8)') pts



  If (godunov .eqv. .true.) then
     !dt_max = 0.9
     !write(*,*) 'Explicit scheme, dt max =', dt_max
     !write(*,*) 'Choose a time step less than this'
     !write(*,*) 'Please input your choice of time step size as a real number :>'
     write(*,*) 'Using Cournot Number of 0.9:'
     OPEN(10,file='ev.dat')
     read(10,'(AAAAA)') title2
     read(10,*) ev_max
     !read(*,*)  dt
     close(10)
     dt = C*dx/(ev_max+1)
     write(*,*) 'Using C =',C,', and lamda_max = ',ev_max,'dt = ',dt
     !write(*,*) 'dt=',dt
  end if

  write(*,*) 'Please input your choice of output time as a real number, in seconds :>'
  read(*,*)  tout
  nt=1!int(tout/dt)+1
  write(*,*)'nt*dt=',nt*dt,'nt=', nt
  
  ALLOCATE( h(1,ncells,nt+1))
  ALLOCATE( uH(1,ncells,nt+1))
  ALLOCATE( u(1,ncells,nt+1))
  !pause
!!$
!!$  If (implicit .eqv. .true.) then
!!$     write(*,*) 'Implicit scheme, time step is independent of mesh size'
!!$     write(*,*) 'Please input your choice of time step as a real number:>'
!!$     read(*,*)  dt
!!$  end if


  write(*,*) ''
  call InitialConditions1D(npts, ncells, pts, cells,  HL, HRoHL, H, u, uH)
  write(*,*) '  Initial H:'
  do i=1,ncells
     write(*,*)'cell loc=',cells(1,i),'initial H=', H(1,i,1),'initial U = ',u(1,i,1),'initial uH=',uH(1,i,1)
  end do
  !write(*,'(A4,D16.4)') '  k= ',k
  !stop
  ALLOCATE( lamda(2,npts) )
  ALLOCATE( F(2,npts) )
  ALLOCATE( Q(2,npts) )
  ALLOCATE( alpha(2,npts) )




  g=9.81 !! gravitational acceleration, m/s/s
  !tol = .001
!!========================================================================================
  !If (explicit .eqv. .True.) Then
  call cpu_time(start)
  Write (*,*) '------------------Begin Gudenov Shallow Water Solver------------------------'
  Write (*,*) '============================================================================'

  
  u(1,1,1)=1.0
  time=0.
  !! #Time Loop -> changed for do while!
  n=1
  do while (time<tout)
!  do n=1,nt

     !! #Cell face "pts" loop to update Roe avg vbls at the faces
     dummy5=1.
     do i=1,npts
        H_hat(1,i)=(H(1,i+1,n)+H(1,i,n))/2.
        
        u_hat(1,i)=( u(1,i+1,n)*sqrt(H(1,i+1,n)) + u(1,i,n)*sqrt(H(1,i,n)) )/&
             (sqrt(H(1,i+1,n))+sqrt(H(1,i,n)))

        lamda(1,i) = u_hat(1,i)-sqrt(g*H_hat(1,i))
        lamda(2,i) = u_hat(1,i)+sqrt(g*H_hat(1,i))


        dummy1 = 1.0/( lamda(2,i) - lamda(1,i) )
        dummy2 = H(1,i+1,n)-H(1,i,n)
        dummy3 = uH(1,i+1,n) - uH(1,i,n)
        dummy4 = u(1,i,n)*uH(1,i,n) + g*((H(1,i,n))**2)/2.

        alpha(1,i) = dummy1* ( lamda(2,i)*dummy2 - dummy3 )
        alpha(2,i) = dummy1* ( -lamda(1,i)*dummy2 + dummy3 )

        !! Define the flux at each interface
        F(1,i) = uH(1,i,n) + min(lamda(1,i),0.)*alpha(1,i) + min(lamda(2,i),0.)*alpha(2,i)
        F(2,i) = dummy4 + min(lamda(1,i),0.)*alpha(1,i)*lamda(1,i) + min(lamda(2,i),0.)*alpha(2,i)*lamda(2,i)

        dummy6 = max( abs(lamda(1,i)) , abs(lamda(2,i)) )
        dummy5 = max(dummy5,dummy6)
        !write(*,*)'lamda1= ',lamda(1,i),'lamda2= ',lamda(2,i),'','Flux1=',F(1,i),'Flux2=',F(2,i)
        !write(*,*)'alpha1= ',alpha(1,i),'alpha2= ',alpha(2,i)

     end do !End Riemann flux update
     !pause

     !! #Update cell centered values H, uH, and u
     do i=2,ncells-1
        H(1,i,n+1)  = H(1,i,n) + (dt/dx)*(F(1,i-1)-F(1,i))
        uH(1,i,n+1) = uH(1,i,n) + (dt/dx)*(F(2,i-1)-F(2,i))
        u(1,i,n+1)  = uH(1,i,n+1)/H(1,i,n+1)

        !write(*,*)'i= ',i,' H= ',H(1,i,n),'uH=',uH(1,i,n),' u= ',u(1,i,n)
     end do !End conservation update

     !Update the LHS ghost cell:
     H(1,1,n+1)  = H(1,2,n+1)
     u(1,1,n+1)  = -u(1,2,n+1)
     uH(1,1,n+1) = H(1,1,n+1)*u(1,1,n+1)
     !Update the RHS ghost cell:
     H(1,ncells,n+1)  = H(1,ncells-1,n+1)
     u(1,ncells,n+1)  = -u(1,ncells-1,n+1)
     uH(1,ncells,n+1) = u(1,ncells,n)*H(1,ncells,n)
     
     ! This is the time at which the updated values exist
     !time = time+dt
     store1p=0.
     store2p=0.
     store1n=0.
     store2n=0.
     if (n==nt) then
        t_1sec = n
        svtm   = time
        !WRITE(*,*) 't_1sec = ',t_1sec, 't=',time
        !write(*,*)'--------------------Soulution out:------------------------------'
        !do i=1,ncells
        !   write(*,*)'i= ',i,'cell = ',cells(1,i),' H= ',H(1,i,n+1),'uH=',uH(1,i,n+1),' u= ',u(1,i,n+1)
        !end do
        do i=1,npts
           store1p=max(store1p,lamda(1,i))
           store2p=max(store2p,lamda(2,i))
           store1n=max(store1n,-lamda(1,i))
           store2n=max(store2n,-lamda(2,i))
           !write(*,*)'i= ',i,'cell = ',cells(1,i),' eigen1= ', lamda(1,i),'  eigen2= ', lamda(2,i)
        end do
        store1p=max(store1p,store2p)
        store1n=-max(store1n,store2n)

        !print*

        !print*
        !write(*,*)'eigen1 max = ',store1n,'eigen2 max=',store1p
        

     end if

        !Update in place
        do i=1,ncells
           H(1,i,n)  = H(1,i,n+1)
           uH(1,i,n) = uH(1,i,n+1)
           u(1,i,n)  = u(1,i,n+1)
           
           !write(*,*)'i= ',i,' H= ',H(1,i,n),'uH=',uH(1,i,n),' u= ',u(1,i,n)
        end do !End conservation update

     dt = C*dx/(store1p)
     time = time+dt
  !end time loop
  end do
  write(*,*)'eigen1 max = ',store1n,'eigen2 max=',store1p

  do i=1,npts
     if (abs((lamda(1,i)-store1n)) < .00000000000000001) then
        !write(*,*) lamda(1,i), store1n, ((lamda(1,i)-store1n))
        waveloc = pts(1,i)
        waveI=i
        !write(*,*) waveloc, waveI
     end if
     
     !write(*,*)'i= ',i,'lamda1= ',lamda(1,i),'lamda2= ',lamda(2,i),'','Flux1=',F(1,i),'Flux2=',F(2,i)
     write(*,*)'i= ',i,'uhat= ',H_hat(1,i),'Hhat= ',H_hat(1,i)
  end do
  call cpu_time(finish)
!!========================================================================================



!!----------------------------Exact Solution----------------------------------------
  ALLOCATE( H_exact(1,ncells,3) )
  ALLOCATE( error(1,ncells,3) )



!!$  print*
!!$   write(*,*)'----------------------------------------------------'
!!$   write(*,*) 'Comparison,         Numerical    to        Exact'
!!$  do i=1,ncells
!!$     write(*,*) n, H(1,i,n+1), H_exact(1,i,n)
!!$  end do
!!$  print*
!!$  write(*,*) '  % Error'
!!$  do i=1,ncells
!!$     write(*,*) error(1,i,n)
!!$  end do

  
!!----------------------------End Exact Solution----------------------------------------

  print '("Computation Time = ",f19.10," seconds.")',finish-start
  ! Question 1 output style
  ! output Gnuplot file-------------------------Solution---------------------------
  OPEN(UNIT=1, FILE='x_pos.dat',FORM='FORMATTED',STATUS='REPLACE')!,IOSTAT=ios)
  DO i=2, ncells-1
     WRITE(UNIT=1,FMT='(F19.10)')cells(1,i)
  End do
  WRITE(UNIT=1,FMT='(" ")')
  CLOSE(UNIT=1)

  ! output Gnuplot file-------------------------Solution---------------------------
  OPEN(UNIT=2, FILE='u.dat',FORM='FORMATTED',STATUS='REPLACE')!,IOSTAT=ios)
  DO i=2, ncells-1
     WRITE(UNIT=2,FMT='(F19.10)') u(1,i,t_1sec)
  End do
  WRITE(UNIT=2,FMT='(" ")')
  CLOSE(UNIT=2)
!!$
  OPEN(UNIT=3, FILE='H.dat',FORM='FORMATTED',STATUS='REPLACE')!,IOSTAT=ios)
  DO  i=2, ncells-1
     WRITE(UNIT=3,FMT='(F19.10)')H(1,i,t_1sec)
  End do
  WRITE(UNIT=3,FMT='(" ")')
  CLOSE(UNIT=3)

  OPEN(UNIT=4, FILE='uH.dat',FORM='FORMATTED',STATUS='REPLACE')!,IOSTAT=ios)
  DO  i=2, ncells-1
     WRITE(UNIT=4,FMT='(F19.10)')uH(1,i,t_1sec)
  End do
  WRITE(UNIT=4,FMT='(" ")')
  CLOSE(UNIT=4)

  OPEN(UNIT=5, FILE='time.dat',FORM='FORMATTED',STATUS='REPLACE')!,IOSTAT=ios)
  WRITE(UNIT=5,FMT='(F19.10)')svtm
  WRITE(UNIT=5,FMT='(" ")')
  CLOSE(UNIT=5)

  OPEN(UNIT=5, FILE='LeftWaveSpeed.dat',FORM='FORMATTED',STATUS='REPLACE')!,IOSTAT=ios)
  WRITE(UNIT=5,FMT='(F19.10,F19.10,I5)')store1n, waveloc, waveI
  WRITE(UNIT=5,FMT='(" ")')
  CLOSE(UNIT=5)




  

  IF (ALLOCATED(pts))          DEALLOCATE(pts)
  IF (ALLOCATED(cells))        DEALLOCATE(cells) 
  IF (ALLOCATED(lamda))        DEALLOCATE(lamda)
  IF (ALLOCATED(F))            DEALLOCATE(F)
  IF (ALLOCATED(Q))            DEALLOCATE(Q)
  IF (ALLOCATED(alpha))        DEALLOCATE(alpha)

  IF (ALLOCATED(H))            DEALLOCATE(H) 
  IF (ALLOCATED(u))            DEALLOCATE(u) 
  IF (ALLOCATED(H_hat))        DEALLOCATE(H_hat) 
  IF (ALLOCATED(u_hat))        DEALLOCATE(u_hat) 
  IF (ALLOCATED(uH))           DEALLOCATE(uH) 
  IF (ALLOCATED(H_exact))      DEALLOCATE(H_exact)


  IF (ALLOCATED(k))            DEALLOCATE(k)
  IF (ALLOCATED(am))           DEALLOCATE(am)
  IF (ALLOCATED(bm))           DEALLOCATE(bm)
  IF (ALLOCATED(cm))           DEALLOCATE(cm)
  IF (ALLOCATED(dm))           DEALLOCATE(dm)


END PROGRAM solver


