      MODULE date_and_computer
      USE system_mod, ONLY: getenv
      USE safe_open_mod
      IMPLICIT NONE
      CHARACTER(LEN=3), DIMENSION(12), PARAMETER :: months =            &
        (/ 'Jan','Feb','Mar','Apr','May','Jun',                         &
           'Jul','Aug','Sep','Oct','Nov','Dec' /)
      CHARACTER(LEN=100) :: computer, os, os_release
     
      CONTAINS

      SUBROUTINE GetComputerInfo
#if defined(WIN32)
      computer = ' Window_NT'
      os       = ' MS Windows 2000'
      os_release  = ' 5.00'
#else
!     SAL - USE GNU F77 Intrisic Procedures
      CALL getenv('HOST',computer)
      CALL getenv('OSTYPE',os)
      CALL getenv('HOSTTYPE',os_release)
#endif
      END SUBROUTINE GetComputerInfo

      END MODULE date_and_computer
