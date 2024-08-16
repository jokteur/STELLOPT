      FUNCTION pcurr (xx)
      USE stel_kinds
      USE vmec_input, ONLY: ac, bloat
C-----------------------------------------------
      REAL(rprec), PARAMETER :: one = 1        
      INTEGER     :: i, ioff
      REAL(rprec) :: xx, pcurr, x
C-----------------------------------------------
!
!     NOTE:  AC COEFFICIENTS OBTAINED IN THREED1 FILE
!            BY MATCHING TO <JTOR> * dV/dPHI ~ SUM[x^(i+1) * ac(i)/(i+1)], i=0...UBOUND(ac)
!
      x = MIN (ABS(xx * bloat), one)
      ioff = LBOUND(ac,1)

      pcurr = 0

      DO i = UBOUND(ac,1), ioff, -1
         pcurr = x*pcurr + ac(i)/(i-ioff+1)
      END DO
      pcurr = x*pcurr

      END FUNCTION pcurr
