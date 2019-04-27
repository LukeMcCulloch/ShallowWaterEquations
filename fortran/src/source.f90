MODULE source

  USE precise, ONLY : defaultp


  IMPLICIT NONE    ! makes sure the compiler complains about 
                   ! undeclared variables
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp
  real, parameter, private :: e=2.718281828459045

CONTAINS


  FUNCTION source1(X,Y) result(S1)
    !
    ! Function to compute a 
    ! Scalar source term
    !

    REAL(wp) :: X, Y, S1

    S1=X*(e**Y)
    
   END FUNCTION source1


   Function sourceC(c, dt, dx) result(nu)

    REAL(wp) :: c, dt, dx, nu
     
     nu=c*dt/dx

   End Function sourceC



END MODULE source
