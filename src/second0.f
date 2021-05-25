      SUBROUTINE second0(stime)
      USE stel_kinds
      IMPLICIT NONE

!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), INTENT(out) :: stime
      INTEGER :: cnt, cnt_rate
!-----------------------------------------------
      CALL SYSTEM_CLOCK(count=cnt, count_rate=cnt_rate)
      IF (cnt_rate .ne. 0) THEN
         stime = REAL(cnt, rprec)/cnt_rate
      ELSE
         stime = 0
      END IF
      END SUBROUTINE second0
