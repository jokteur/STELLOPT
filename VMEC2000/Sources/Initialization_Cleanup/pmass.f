      FUNCTION pmass (xx)
      USE stel_kinds
      USE vmec_input, ONLY: am, bloat, pres_scale
      USE vparams, ONLY: mu0
C-----------------------------------------------
      INTEGER     :: i
      REAL(rprec) :: xx, pmass, x
C-----------------------------------------------
!     NOTE: On entry, am is in pascals. pmass internal units are mu0*pascals (B**2 units)
      x = MIN (ABS(xx * bloat), 1._dp)

      pmass = 0

      DO i = UBOUND(am,1), LBOUND(am,1), -1
         pmass = x*pmass + am(i)
      END DO

      pmass = mu0*pres_scale*pmass

      END FUNCTION pmass

      
      FUNCTION photp (xx)
      USE stel_kinds
      USE vmec_input, ONLY: ah, bloat
C-----------------------------------------------
      INTEGER     :: i
      REAL(rprec) :: xx, photp, x
C-----------------------------------------------
!     NOTE: On entry, ah is dimensionless
!     HOT PARTICLE PRESSURE (RATIO TO ISOTROPIC PRESSURE)

      x = MIN (ABS(xx * bloat), 1._dp)

      photp = 0

      DO i = UBOUND(ah,1), LBOUND(ah,1), -1
         photp = x*photp + ah(i)
      END DO

      photp = photp

      END FUNCTION photp

      FUNCTION ptrat (xx)
      USE stel_kinds
      USE vmec_input, ONLY: at, bloat
C-----------------------------------------------
      INTEGER     :: i
      REAL(rprec) :: xx, ptrat, x
C-----------------------------------------------
!     NOTE: On entry, at is dimensionless
!     HOT PARTICLE T-perp/T-par

      x = MIN (ABS(xx * bloat), 1._dp)

      ptrat = 0

      DO i = UBOUND(at,1), LBOUND(at,1), -1
         ptrat = x*ptrat + at(i)
      END DO

      ptrat = ptrat

      END FUNCTION ptrat
