!!
!! Simple Fortran 90 module with program constants
!! You could have delcared pi directly in the program units that
!! need pi. However the purpose of modules is to make pieces of the
!! code reusable and easier to maintain.
!!
!! September 2012, 
!! 
MODULE constants

  use precise, only : defaultp
  IMPLICIT NONE
  
  INTEGER, PARAMETER, PRIVATE :: WP=defaultp

  ! This defines the number pi. The key word parameter makes sure, 
  ! the number pi is not overridden when the program is executed
  REAL(WP), PARAMETER :: one=1.0  
  REAL(WP), PARAMETER :: four=4.0  
  REAL(WP), PARAMETER :: pi=Four*ATAN(one)  ! using a function to define the
                                     ! parameter value is a Fortran 2003
                                     ! extension supported by g95 but
                                     !  not all other compilers
  !real(wp), parameter :: pi=3.141592653589793 ! this is an alternative form

END MODULE constants
