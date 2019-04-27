MODULE analytic

  USE precise, ONLY : defaultp
  use constants                  


  IMPLICIT NONE    
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp
  real, parameter, private    :: e=2.718281828459045

CONTAINS


  FUNCTION lamda_n(n,L) result(lamda)
    INTEGER  :: n
    REAL(wp) :: L,lamda
    lamda = (2.*n-1)*pi/(2.*L)
  end function lamda_n


  FUNCTION alpha_(k,rho_c) result(alpha)
    REAL(wp) :: k,rho_c,alpha
    alpha = k/rho_c
  end function alpha_
  


  subroutine t_exact(L,k_bar,rho_c,n,dt,ncells,loc,temp_exact)

    INTEGER :: i          ! summation dummy
    INTEGER :: j          ! summation dummy
    INTEGER :: n          ! timestep
    INTEGER :: ncells     ! # of cells

    real, parameter :: e=2.718281828459045

    real(wp) :: L         ! domain length
    real(wp) :: rho_c
    real(wp) :: k_bar
    real(wp) :: alpha
    real(wp) :: lamda
    real(wp) :: rational_term
    real(wp) :: dt
    real(wp) :: t

    REAL(wp), DIMENSION(ncells)   :: loc         ! x locations -> pass cells
    REAL(wp), DIMENSION(ncells)   :: temp_exact  ! solution

    alpha = alpha_(k_bar,rho_c)
    t=n*dt

    print '("1")'
    do i=1,ncells
       do j=1,30
          lamda            = lamda_n(j,L)
          rational_term    = ((-1)**(j+1)) / (-1+2.*j)
          temp_exact(i)    = temp_exact(i) + rational_term * e**(-alpha*(lamda**2)*t) *cos(lamda*loc(i))
          !print '("3")'
       end do
       
       temp_exact(i) =  (800./pi)*temp_exact(i)
    end do

  end subroutine t_exact



END MODULE analytic
