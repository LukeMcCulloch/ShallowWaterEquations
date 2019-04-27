MODULE precise
!
!
! Module to handle the question of precision gracefully.
! Put the following lines in your code.
!    use precise, only : defaultp
!    integer, parameter, private :: WP=defaultp
! You have to omit the private keyword outside of modules
!    use precise, only : defaultp
!    integer, parameter :: WP=defaultp
!
! Now you can declare al real variables as
!    real(wp) :: ....
! If you intend to change from single to double precision just
! change the last integer statement below from
!   defaultp=sp to
!   defaultp=dp
! recompile all program units and you are done.
!
  IMPLICIT NONE
  INTEGER, PARAMETER :: SP = KIND(1.0E+0)
  INTEGER, PARAMETER :: DP = KIND(1.0D+0)
  INTEGER, PARAMETER :: DEFAULTP = DP        ! Default working precision

END MODULE precise
