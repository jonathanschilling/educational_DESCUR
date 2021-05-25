      SUBROUTINE getcarg(narg, arg, numargs)
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER, INTENT(in)  :: narg
      INTEGER, INTENT(out) :: numargs
      CHARACTER(LEN=*), INTENT(out) :: arg
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: numchars
!-----------------------------------------------

      numargs = iargc()
      call getarg(narg, arg)

      END SUBROUTINE getcarg
