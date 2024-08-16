      FUNCTION piota (x)
      USE stel_kinds
      USE vmec_input, ONLY: ai
C-----------------------------------------------
      INTEGER     :: i
      REAL(rprec) :: x, piota
C-----------------------------------------------
      piota = 0

      DO i = UBOUND(ai,1), LBOUND(ai,1), -1
         piota = x*piota + ai(i)
      END DO

      END FUNCTION piota
