      SUBROUTINE second0(stime)
      USE stel_kinds
      IMPLICIT NONE
#if defined(MPI_OPT)
      INCLUDE 'mpif.h'
#endif

!----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), INTENT(out) :: stime
      INTEGER :: cnt, cnt_rate
!-----------------------------------------------
!     USE CPU_TIME IF ACTUAL CPU USAGE DESIRED
!     CALL CPU_TIME(stime)
!     RETURN
#if defined(MPI_OPT)
      stime=MPI_Wtime()
#else
      CALL SYSTEM_CLOCK(count=cnt, count_rate=cnt_rate)
      IF (cnt_rate .ne. 0) THEN
         stime = REAL(cnt, rprec)/cnt_rate
      ELSE
         stime = 0
      END IF
#endif
      END SUBROUTINE second0
