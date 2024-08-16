      INTEGER FUNCTION vmec_chdir(new_path)
      IMPLICIT NONE
      CHARACTER*(*), INTENT(in) :: new_path
#if defined(CRAY)
      INTEGER :: ilen
      iLEN = 0
      CALL pxfchdir(new_path, ilen, vmec_chdir)
#else
      INTEGER, EXTERNAL :: chdir
 
#if defined(SUNOS)
      vmec_chdir = chdir(TRIM(new_path))
#else
      vmec_chdir = chdir(TRIM(new_path) // char(0))
#endif
#endif
      END FUNCTION vmec_chdir
