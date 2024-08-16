      MODULE stel_constants

      USE stel_kinds

!----------------------------------------------------------------------
!  Mathematical constants
!----------------------------------------------------------------------

      REAL(rprec), PARAMETER :: pi=3.14159265358979323846264338328_rprec
      REAL(rprec), PARAMETER :: pio2=pi/2
      REAL(rprec), PARAMETER :: twopi=2*pi
      REAL(rprec), PARAMETER :: sqrt2=1.41421356237309504880168872_rprec
      REAL(rprec), PARAMETER :: degree=twopi / 360
      REAL(rprec), PARAMETER :: one=1
      REAL(rprec), PARAMETER :: zero=0
 
!----------------------------------------------------------------------
!  Physical constants
!------------------------------------------------------------------

      REAL(rprec), PARAMETER :: mu0 = 2 * twopi * 1.0e-7_rprec

      END MODULE stel_constants
