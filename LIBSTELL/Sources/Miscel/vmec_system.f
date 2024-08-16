      SUBROUTINE vmec_system(cmd, ierror)
      INTEGER, OPTIONAL :: ierror
      INTEGER :: ireturn
      CHARACTER(LEN=*), INTENT(in) :: cmd

#if defined(CRAY)
      INTEGER, EXTERNAL :: ishell
      ireturn = ishell(TRIM(cmd))
#elif defined(RISC)
      CALL system(TRIM(cmd), ireturn)
#elif defined(IRIX64)
      CALL system(TRIM(cmd))
      ireturn = 0
#elif defined(LINUX) || defined(OSF1)
!      INTEGER, EXTERNAL :: system
      INTEGER :: system
      ireturn = system(TRIM(cmd))
#elif defined(WIN32) || defined(SUNOS)
      INTEGER, EXTERNAL :: system
      ireturn = system(TRIM(cmd))
#else
      INTEGER, EXTERNAL :: system
      ireturn = system(TRIM(cmd) // CHAR(0))
#endif
      IF (PRESENT(ierror)) ierror = ireturn

      END SUBROUTINE vmec_system
