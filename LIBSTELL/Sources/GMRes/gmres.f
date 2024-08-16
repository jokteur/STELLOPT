      SUBROUTINE gmres (n, m, icntl, cntl, yAx, x0, b, info)
      USE stel_kinds, ONLY: rprec
      USE stel_constants, ONLY: one, zero
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: n, m,  icntl(9), info(3)
      REAL(rprec), INTENT(in)    :: b(n)
      REAL(rprec), INTENT(inout) :: x0(n)
      REAL(rprec) :: cntl(5)
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: revcom, colx, coly, colz, nbscal
      INTEGER :: irc(5), jcount
      INTEGER, PARAMETER :: matvec=1, precondLeft=2, 
     1                      precondRight=3, dotProd=4
      INTEGER :: nout, lwork
      REAL(rprec)  :: rinfo(2)
      REAL(rprec), TARGET, ALLOCATABLE :: work(:)
      REAL(rprec), POINTER :: sx(:), sy(:), sz(:)
C-----------------------------------------------
      EXTERNAL yAx
C-----------------------------------------------
!
!     EASY-TO-USE WRAPPER FOR GMRES DRIVER CALL
!
!     X0: on input, initial guess if icntl(6) == 1
!         on output, solution of Ax = b
!         NOTE: it is not overwritten UNTIL the end of this routine
!
      lwork = m**2 + m*(n+6) + 5*n + 1
      ALLOCATE (work(lwork), stat=nout)
      IF (nout .ne. 0) STOP 'Allocation error in gmres!'
      work = 0
      IF (icntl(6) .eq. 1) work(1:n) = x0
      work(n+1:2*n) = b(1:n)

*****************************************
** Reverse communication implementation
*****************************************
*
 10   CONTINUE
      CALL drive_dgmres(n,n,m,lwork,work,
     &                  irc,icntl,cntl,info,rinfo)
      revcom = irc(1)
      colx   = irc(2)
      coly   = irc(3)
      colz   = irc(4)
      nbscal = irc(5)
      sx => work(colx:);  sy => work(coly:);  sz => work(colz:)

      IF (revcom.eq.matvec) THEN
* perform the matrix vector product
*        work(colz) <-- A * work(colx)
         CALL yAx (sx, sz, n) 
         GOTO 10
*
      ELSE IF (revcom.eq.precondLeft) THEN
* perform the left preconditioning
*         work(colz) <-- M^{-1} * work(colx)
!         CALL dcopy(n,work(colx),1,work(colz),1)
         CALL dcopy(n,sx,1,sz,1)
         GOTO 10
*
      ELSE IF (revcom.eq.precondRight) THEN
* perform the right preconditioning
!         CALL dcopy(n,work(colx),1,work(colz),1)
         CALL dcopy(n,sx,1,sz,1)
         GOTO 10
*
      ELSE IF (revcom.eq.dotProd) THEN
*      perform the scalar product
*      work(colz) <-- work(colx) work(coly)
*
         CALL dgemv('C',n,nbscal,ONE,sx,n,sy,1,ZERO,sz,1)
         GOTO 10
      ENDIF

*******************************
* dump the solution to a file for debugging
*******************************
!  JDH Commented out below 2008-05-15
!      GOTO 100
!
!      nout = 11
!      open(nout,FILE='sol_dTestgmres',STATUS='unknown')
!      if (icntl(5).eq.0) then
!        write(nout,*) 'Orthogonalization : MGS'
!      elseif (icntl(5).eq.1) then
!        write(nout,*) 'Orthogonalization : IMGS'
!      elseif (icntl(5).eq.2) then
!        write(nout,*) 'Orthogonalization : CGS'
!      elseif (icntl(5).eq.3) then
!        write(nout,*) 'Orthogonalization : ICGS'
!      endif
!      write(nout,*) 'Restart : ', m
!      write(nout,*) 'info(1) = ',info(1),'  info(2) = ',info(2)
!      write(nout,*) 'rinfo(1) = ',rinfo(1),'  rinfo(2) = ',rinfo(2)
!      write(nout,*) 'Optimal workspace = ', info(3)
!      write(nout,*) 'Solution : '
!      do jcount=1,n
!        write(nout,*) work(jcount)
!      enddo
!      write(nout,*) '   '
!*
!100   continue
!*

      x0 = work(1:n)

      DEALLOCATE (work)

      END SUBROUTINE gmres
