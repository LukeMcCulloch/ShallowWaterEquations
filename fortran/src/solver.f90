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
  integer, parameter :: num_eqn = 2
  integer, parameter :: num_waves = 2
  real(wp), parameter :: e=2.718281828459045
  real(wp), parameter :: atol = 0.00000000000000001    ! a new tolerance
  INTEGER, DIMENSION(3), PARAMETER :: it1 = (/ 40, 80, 120 /)   ! output times

  INTEGER :: narg    ! # of command line arguments
  INTEGER :: i       ! dummy index in space
  INTEGER :: j       ! dummy index in space
  INTEGER :: n       ! dummy index in time
  INTEGER :: ex_n    ! seperate time index for exact sol.
  INTEGER :: npts    ! # of points 64001 interfaces (1D boundary is a point ghost nodes need to outer boundary)
  INTEGER :: ncells  ! # of cells 64002 centers (cell -> cell center)
  INTEGER :: nt      ! # of timesteps
  INTEGER :: t_1sec  ! index at approx 1 sec
  INTEGER :: waveI   ! index location of the max Left Eigenvalue

  INTEGER :: m,mw,mx ! entropy dummies

  ! Efficient Placement of these declarations...?
  CHARACTER(len=24) :: inputfile
  CHARACTER(len=24) :: outputfile
  CHARACTER(len=30) :: title
  CHARACTER(len=30) :: title2
  CHARACTER(len=30) :: scheme_type

  LOGICAL :: flexists     ! a logical variable. I .true. or .false.
  LOGICAL :: godunov      ! true for theta == integer 1 else theta is real

  integer, ALLOCATABLE, DIMENSION(:)      :: mthlim       ! LaVeque limiter selector

  REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: pts          ! 1 x npts array of points (cell boundaries)
  REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: cells        ! 1 x ncells = npts+1  array of cells

  REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: lamda        ! 2 x ncells = vector of eigenvalues
  REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: F            ! 2 x ncells = Flux vector for each cell
  REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) :: Q            ! 2 x ncells = storage vector for current H,uH
  REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: qL,qR        ! left and right q(i) on the ith cell
                                                          ! these agree with Q(i), the cell average.
  REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: alpha        ! 1 x ncells = npts+1  array of cells
  REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) :: wave         ! Riemann problem waves (num_eqn, num_wave, ncells)

  REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) :: H            ! 1 x ncells array cell - height
  REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) :: uH           ! 1 x ncells array cell - momentum
  REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) :: u            ! 1 x ncells array cell - speed
  REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: H_hat        ! 1 x npts face centered
  REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: u_hat        ! 1 x npts face centered
  REAL(WP), ALLOCATABLE, DIMENSION(:,:,:) :: H_exact      ! 1 x ncells array cell 

  REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: s            ! entropy fix
  REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: amdq         ! left-going flux-differences
  REAL(WP), ALLOCATABLE, DIMENSION(:,:)   :: apdq         ! right-going flux-differences



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
  real(wp) :: dtdx     ! dt/dx
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
  !real(wp) :: dummy2   ! dummy vbl
  !sreal(wp) :: dummy3   ! dummy vbl
  real(wp) :: dummy4   ! dummy vbl
  real(wp) :: dummy5   ! dummy vbl
  real(wp) :: dummy6   ! dummy vbl  

  real(wp) :: delta(2) !

  real(wp) :: cbar     ! sqrt(gh)

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

  !entropy fix
  real(wp) :: s0 ! entropy
  real(wp) :: s1 ! entropy
  real(wp) :: s2 ! entropy
  real(wp) :: s3 ! entropy

  real(wp) :: hr1 ! entropy
  real(wp) :: uhr1 ! entropy
  real(wp) :: hl2 ! entropy
  real(wp) :: uhl2 ! entropy
  real(wp) :: sfract ! entropy

  real(wp) :: df      ! entropy
  real(wp) :: dtdxave ! entropy
  



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
     WRITE(6,'(A)') 'ENME 6728 Unsteady Shallow Water Scheme'
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
  !write(*,'(1D16.8)') pts



  If (godunov .eqv. .true.) then
     !dt_max = 0.9
     !write(*,*) 'Explicit scheme, dt max =', dt_max
     !write(*,*) 'Choose a time step less than this'
     !write(*,*) 'Please input your choice of time step size as a real number :>'
     !write(*,*) 'Using Cournot Number of 0.9:'
     OPEN(10,file='ev.dat')
     read(10,'(AAAAA)') title2
     read(10,*) ev_max
     !read(*,*)  dt
     close(10)
     dt = C*dx/(ev_max+1)
     write(*,*) 'Using C =',C,', and lamda_max = ',ev_max,'dt = ',dt
     !write(*,*) 'dt=',dt
  end if

  !write(*,*) 'Please input your choice of output time as a real number, in seconds :>'
  !read(*,*)  tout
  tout = 1.
  !nt = !int(tout/dt)+1 
  nt=1!do time stepping in place.
  write(*,*)'nt*dt=',nt*dt,'nt=', nt
  
  ALLOCATE( h(1,ncells,nt+1))
  ALLOCATE( uH(1,ncells,nt+1))
  ALLOCATE( u(1,ncells,nt+1))
  ALLOCATE( Q(num_eqn,ncells,nt+1) )
  allocate( qL(num_eqn,ncells))
  allocate( qR(num_eqn,ncells))
  allocate( wave(num_eqn, num_waves, ncells ) )
  allocate( s(num_waves, ncells ) )
  allocate( amdq(num_eqn, ncells ) )
  allocate( apdq(num_eqn, ncells ) )
  allocate( mthlim(num_waves) )

!! limiter settings
! mthlim(mw) = 0 for no limiter
!             = 1 for minmod
!             = 2 for superbee
!             = 3 for van Leer
!             = 4 for monotonized centered
  mthlim(1) = 1
  mthlim(2) = 1


  write(*,*) ''
  call InitialConditions1D(npts, ncells, pts, cells,  HL, HRoHL, H, u, uH, Q, qL, qR)
  !call InitialConditions1D(npts, ncells, pts, cells,  HL, HRoHL, H, u, uH, Q, qL, qR)
  !write(*,*) '  Initial H:'
  !do i=1,ncells
  !   write(*,*)'cell loc=',cells(1,i),'initial H=', H(1,i,1),'initial U = ',u(1,i,1),'initial uH=',uH(1,i,1)
  !end do
  !write(*,'(A4,D16.4)') '  k= ',k
  !stop
  ALLOCATE( lamda(num_eqn,npts) )
  ALLOCATE( F(num_eqn,npts) )
  ALLOCATE( alpha(num_eqn,npts) )

  print*, 'npts   = ', npts
  print*, 'ncells = ', ncells


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
   !write(*,*) time

     !! #Cell face "pts" loop to update Roe avg vbls at the faces
     dummy5=1.

     !------------------------------------------------------
     ! Main loop of the Riemann solver.
     do i=1,npts
        ! hbar, LaVeque 15.32 
        H_hat(1,i)=(H(1,i+1,n)+H(1,i,n))/2.
        
        ! Roe Average Velocity, LaVeque 15.35
        u_hat(1,i)=( u(1,i+1,n)*sqrt(H(1,i+1,n)) + u(1,i,n)*sqrt(H(1,i,n)) )/&
             (sqrt(H(1,i+1,n))+sqrt(H(1,i,n)))
        !mixing conventions.... sorry!
        cbar = sqrt(0.5d0*g*(qr(1,i+1) + ql(1,i)))  !  0.5 * g * (u(1,i) + u(1,i-1)) where u(1) is h

        !! eigen values, LaVeque 15.36
        lamda(1,i) = u_hat(1,i)-sqrt(g*H_hat(1,i))
        lamda(2,i) = u_hat(1,i)+sqrt(g*H_hat(1,i))
        

        ! ubar = (sqrt(u_l) + sqrt(u_r)) / (sqrt(h_l) + sqrt(h_r))
        ! cbar = sqrt( 0.5 * g * (h_l + h_r))
        !
        ! ubar = u_hat
        ! cbar = sqrt(g*H_hat(1,i))



        dummy1 = 1.0/( lamda(2,i) - lamda(1,i) )
        delta(1) = H(1,i+1,n)-H(1,i,n) !delta1
        delta(2) = uH(1,i+1,n) - uH(1,i,n) !delta2
        dummy4 = u(1,i,n)*uH(1,i,n) + g*((H(1,i,n))**2)/2.

        !! alpha coefficients, LaVeque 15.39. eignevector expansion of delta(1) delta(2)
        !alpha(1,i) = dummy1* (  lamda(2,i)*delta(1) - delta(2) )
        !alpha(2,i) = dummy1* ( -lamda(1,i)*delta(1) + delta(2) )
!
        alpha(1,i) = 0.5d0*(-delta(2) + (u_hat(1,i) + cbar) * delta(1))/cbar
        alpha(2,i) = 0.5d0*( delta(2) - (u_hat(1,i) - cbar) * delta(1))/cbar



        ! Finally, compute the waves:
        wave(1,1,i) = alpha(1,i)
        wave(2,1,i) = alpha(1,i)*(u_hat(1,i) - cbar)
        s(1,i) = u_hat(1,i) - cbar
        
        wave(1,2,i) = alpha(2,i)
        wave(2,2,i) = alpha(2,i)*(u_hat(1,i) + cbar)
        s(2,i) = u_hat(1,i) + cbar
        ! needed for HR schemes and for entropy fix

        !! Define the flux at each interface
        !F(1,i) = uH(1,i,n) + min(lamda(1,i),0.)*alpha(1,i) + min(lamda(2,i),0.)*alpha(2,i)
        !F(2,i) = dummy4 + min(lamda(1,i),0.)*alpha(1,i)*lamda(1,i) + min(lamda(2,i),0.)*alpha(2,i)*lamda(2,i)




     end do !End Riemann flux update, non entropy
     !------------------------------------------------------



     !------------------------------------------------------
     !! entropy fix for transonic rarefaction (FIXME!)
      !do 200 i=1,npts
      do 200 i = 2, ncells-1 

         ! u-c in left state
         s0 = qr(2,i-1)/qr(1,i-1) - sqrt(g*qr(1,i-1))  !uh/u - sqrt(hu)

        ! check for fully supersonic case:
        if (s0 >= 0.0 .and. s(1,i) > 0.0)  then
            !all right-going
            do m=1,2
               !amdq,apdq = 0.d0 !left and right going waves
               amdq(m,i) = 0.d0
            enddo
            go to 200
        endif

        ! u-c to right of 1-wave
        hr1  = qr(1,i-1) + wave(1,1,i)
        uhr1 = qr(2,i-1) + wave(2,1,i)
        s1 =  uhr1/hr1 - dsqrt(g*hr1)
        if (s0 < 0.d0 .and. s1 > 0.d0) then
            ! transonic rarefaction in the 1-wave
            !write(*,*) 'found transonic rarefaction in the 1 wave'
            sfract = s0 * (s1-s(1,i)) / (s1-s0)
            !write(*,*) 's0 ',s0,'s1 ',s1,'s(1,i) ',s(1,i)
            !write(*,*) 'sfract ',sfract
            !write(*,*) 'time = ', time

        else if (s(1,i) < 0.d0) then
            ! 1-wave is leftgoing
            !write(*,*) '1 wave is leftgoing'
            sfract = s(1,i)
        else
            ! 1-wave is rightgoing
            !write(*,*) '1 wave is rightgoing'
            sfract = 0.d0   !# this shouldn't happen since s0 < 0
        endif

        do m=1,2
            amdq(m,i) = sfract*wave(m,1,i)
            !if (s0 < 0.d0 .and. s1 > 0.d0) then
            !   write(*,*) 'amdq(m,i) = ',amdq(m,i)
            !endif
            enddo


          
        ! -------------------------------------------------------
        ! check 2-wave:
        ! ---------------
        ! u+c in right state  (cell i)
        s3 = ql(2,i)/ql(1,i) + dsqrt(g*ql(1,i))
                      
        ! u+c to left of 2-wave
        hl2  = ql(1,i) - wave(1,2,i)
        uhl2 = ql(2,i) - wave(2,2,i)
        s2 = uhl2/hl2 + dsqrt(g*hl2)
                          
        if (s2 < 0.d0 .and. s3 > 0.d0) then
            ! transonic rarefaction in the 2-wave
            !write(*,*) 'found transonic rarefaction in the 2 wave'
            sfract = s2 * (s3-s(2,i)) / (s3-s2)
        else if (s(2,i) < 0.d0) then
            ! 2-wave is leftgoing
            !write(*,*) '2 wave is leftgoing'
            sfract = s(2,i)
        else
            ! 2-wave is rightgoing
            !write(*,*) '2 wave is rightgoing'
            go to 200
        endif
    
        do m=1,2
            amdq(m,i) = amdq(m,i) + sfract*wave(m,2,i)
            !if (s2 < 0.d0 .and. s3 > 0.d0) then
            !   write(*,*) 'amdq(m,i) = ',amdq(m,i)
            !endif
            enddo

    200 enddo


    ! compute the rightgoing flux differences:
    ! df = SUM s*wave   is the total flux difference and apdq = df - amdq

    do m=1,2
        do i=1,npts
        !do i = 2, ncells-1 
            df = 0.d0
            do mw=1,num_waves
                df = df + s(mw,i)*wave(m,mw,i)
            enddo
            apdq(m,i) = df - amdq(m,i)
         enddo
    enddo

!    !End Riemann entropy correction  
! !  ============================================

!     !! apply correction
     dtdx = dt/dx

    do i = 1, ncells-1
    !do i = 2, ncells
      ! do m = 1, num_eqn ! left boundary
      !     Q(m,i-1,n) = Q(m,i-1,n) - dtdx*amdq(m,i)
      ! end do 
      ! do m = 1, num_eqn ! right boundary
      !    Q(m,i,n) = Q(m,i,n) - dtdx*apdq(m,i)
      ! end do
      do m = 1, num_eqn ! left boundary
          Q(m,i,n) = Q(m,i,n) - dtdx*amdq(m,i)
      end do 
      do m = 1, num_eqn ! right boundary
         Q(m,i+1,n) = Q(m,i+1,n) - dtdx*apdq(m,i)
      end do
      
      !write(*,*)'i= ',i, 'amdq(m,i)=',amdq(m,i), ' apdq(m,i)= ',apdq(m,i)
   end do

!     # compute maximum wave speed:
   C = 0.d0
   do i = 2, ncells-1 
       do mw=1,num_waves
       !          # if s>0 use dtdx(i) to compute CFL,
       !          # if s<0 use dtdx(i-1) to compute CFL:
           C = dmax1(C, dtdx*s(mw,i), -dtdx*s(mw,i))
       end do
   end do



   !! call second order limiter:
!  ============================================
   call limiter(ncells,num_eqn,num_waves,1,ncells, wave,s,mthlim)

!  apply the f
!  ============================================
   !do i=1,mx+1
   do i=1,npts
   !do i = 2, ncells-1 
   !do i = 1, ncells
      do m = 1,num_eqn
          f(m,i) = 0.d0
      end do
      dtdxave = dtdx ! 0.5d0 * ( dtdx(left) + dtdx(right) )
      do mw=1,num_waves
          do m=1,num_eqn
              f(m,i) = f(m,i) + 0.5d0 * dabs(s(mw,i)) &
              * (1.d0 - dabs(s(mw,i))*dtdxave) * wave(m,mw,i)
          end do
      end do
  end do


!  # Update cell centered values H, uH, and u
!  # i.e. Q by differencing correction fluxes
!  ============================================
   do i = 2, ncells-1 
      do m = 1, num_eqn
        Q(m,i,n+1) = Q(m,i,n) + dtdx * ( F(m,i-1) - F(m,i) )
      end do
      !Q(1,i,n+1) = Q(1,i,n) + dtdx * ( F(1,i-1) - F(1,i) )
      !Q(2,i,n+1) = Q(2,i,n) + dtdx * ( F(2,i-1) - F(2,i) )
      H(1,i,n+1)  = Q(1,i,n+1)
      uH(1,i,n+1) = Q(2,i,n+1)
      u(1,i,n+1)  = uH(1,i,n+1)/H(1,i,n+1)
   end do
   



     !! #Update cell centered values H, uH, and u
   !   do i=2,ncells-1
   !      H(1,i,n+1)  = H(1,i,n) + (dt/dx)*(F(1,i-1)-F(1,i))
   !      uH(1,i,n+1) = uH(1,i,n) + (dt/dx)*(F(2,i-1)-F(2,i))
   !      u(1,i,n+1)  = uH(1,i,n+1)/H(1,i,n+1)

   !      Q(1,i,n+1) = H(1,i,n+1)
   !      Q(2,i,n+1) = uH(1,i,n+1)


   !      !write(*,*)'i= ',i,' H= ',H(1,i,n),'uH=',uH(1,i,n),' u= ',u(1,i,n)
   !   end do !End conservation update

     !Update the LHS ghost cell:
     H(1,1, n+1)  =  H(1,2,n+1)
     u(1,1, n+1)  = -u(1,2,n+1)
     uH(1,1,n+1) =   H(1,1,n+1)*u(1,1,n+1)
     Q(1,1, n+1) =   H(1,2,n+1)
     Q(2,1, n+1) =   H(1,1,n+1)*u(1,1,n+1)
     !Update the RHS ghost cell:
     H(1,ncells, n+1)  = H(1,ncells-1,n+1)
     u(1,ncells, n+1)  = -u(1,ncells-1,n+1)
     uH(1,ncells,n+1) = u(1,ncells,n)*H(1,ncells,n)
     Q(1,ncells, n+1) = H(1,ncells-1,n+1)
     Q(2,ncells, n+1) = u(1,ncells,n)*H(1,ncells,n)

     
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

        Q(1,i,n) = Q(1,i,n+1)
        Q(2,i,n) = Q(2,i,n+1) 

        qL(:,i) = Q(:,i,n+1)
        qR(:,i) = Q(:,i,n+1)
         
        !write(*,*)'i= ',i,' H= ',H(1,i,n),'uH=',uH(1,i,n),' u= ',u(1,i,n)
     end do !End conservation update

     dt = C*dx/(store1p)
     time = time+dt
  !end time loop
  end do
  write(*,*)'eigen1 max = ',store1n,'eigen2 max=',store1p

  do i=1,npts
     if (abs((lamda(1,i)-store1n)) < atol) then
        !write(*,*) lamda(1,i), store1n, ((lamda(1,i)-store1n))
        waveloc = pts(1,i)
        waveI=i
        !write(*,*) waveloc, waveI
     end if
     
     !write(*,*)'i= ',i,'lamda1= ',lamda(1,i),'lamda2= ',lamda(2,i),'','Flux1=',F(1,i),'Flux2=',F(2,i)
     !write(*,*)'i= ',i,'uhat= ',H_hat(1,i),'Hhat= ',H_hat(1,i)
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

  IF (ALLOCATED(qL))           DEALLOCATE(qL)
  IF (ALLOCATED(qR))           DEALLOCATE(qR)
  IF (ALLOCATED(wave))         DEALLOCATE(wave)
  IF (ALLOCATED(s))            DEALLOCATE(s)
  IF (ALLOCATED(amdq))         DEALLOCATE(amdq)
  IF (ALLOCATED(apdq))         DEALLOCATE(apdq)

  if (allocated(mthlim))       deallocate(mthlim)

END PROGRAM solver


