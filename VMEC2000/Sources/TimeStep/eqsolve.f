      SUBROUTINE eqsolve(ier_flag, lscreen)
      USE vmec_main
      USE vmec_params, ONLY: ntmax, ns4, jac75_flag, norm_term_flag,
     1                       bad_jacobian_flag, more_iter_flag,
     2                       successful_term_flag
      USE precon2d, ONLY: ScratchFile, lswap2disk, ictrl_prec2d
      USE directaccess, ONLY: DeleteDAFile
      USE gmres_mod, ONLY: nfcn
      USE realspace
      USE vsvd
      USE xstuff
! Add below JDH 2010-08-03
      USE vmec_history
#if defined(SKS)
      USE parallel_include_module
      USE parallel_vmec_module, ONLY: ZeroLastNType
#endif
      USE vacmod, ONLY: nuv, nuv3
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      INTEGER :: ier_flag
      LOGICAL :: lscreen
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL(dp), PARAMETER :: p98 = 0.98_dp, p96 = 0.96_dp
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
#if defined (SKS)      
      INTEGER :: i, j, k, l, lk
#else
      INTEGER, PARAMETER :: rank=0, grank=0
#endif
      REAL(dp) :: w1, r00s, w0, wdota, r0dot, teqsolon, teqsoloff
      LOGICAL :: liter_flag, lreset_internal
C-----------------------------------------------
!
!                INDEX OF LOCAL VARIABLES
!
!        iequi   counter used to call -EQFOR- at end of run
!        ijacob  counter for number of times jacobian changes sign
!        irst    counter monitoring sign of jacobian; resets R, Z, and
!                Lambda when jacobian changes sign and decreases time step
!        signgs  sign of Jacobian : must be =1 (right-handed) or =-1 (left-handed)

!        iterj   stores position in main iteration loop (j=1,2)
!        itfsq   counter for storing FSQ into FSQT for plotting
!        ivac    counts number of free-boundary iterations
!        ndamp   number of iterations over which damping is averaged
!        meven   parity selection label for even poloidal modes of R and Z
!        modd    parity selection label for odd poloidal modes of R and
!        gc      stacked array of R, Z, Lambda Spectral force coefficients (see readin for stack order)
!        xc      stacked array of scaled R, Z, Lambda Fourier coefficients

#if defined(SKS)
      CALL second0(teqsolon)
#endif

      liter_flag = iter2 .eq. 1
 1000 CONTINUE

      itfsq = 0
      rsfac   = one
      w1      = zero
      r00s    = zero
      gphifac = zero
      grmse   = zero

!
!     COMPUTE INITIAL R, Z AND MAGNETIC FLUX PROFILES
!
   20 CONTINUE
!
!     RECOMPUTE INITIAL PROFILE, BUT WITH IMPROVED AXIS/OR RESTART 
!     FROM INITIAL PROFILE, BUT WITH A SMALLER TIME-STEP
!
      IF (irst .EQ. 2) THEN
#if defined(SKS)
        IF (PARVMEC) THEN
          CALL ZeroLastNType(pxc)
          CALL profil3d_par(pxc(1),pxc(1+irzloff),lreset_internal,
     1                       .FALSE.)
        ELSE
#endif
          xc = 0
          CALL profil3d(xc(1),xc(1+irzloff),lreset_internal,.FALSE.)
#if defined(SKS)
        END IF
#endif
        irst = 1
        IF (liter_flag) CALL restart_iter(delt0r)
      END IF
!      IF (liter_flag) CALL restart_iter(delt0r)
      liter_flag = .true.
      ier_flag = norm_term_flag

!
!
!     FORCE ITERATION LOOP
!
      iter_loop: DO WHILE (liter_flag)
!
!     ADVANCE FOURIER AMPLITUDES OF R, Z, AND LAMBDA
         CALL evolve (delt0r, ier_flag, liter_flag, lscreen)

         IF (ijacob.eq.0 .and. (ier_flag.eq.bad_jacobian_flag
     1       .or. irst.eq.4) .and. ns.ge.3) THEN
            IF (lscreen) THEN
               IF (ier_flag .eq. bad_jacobian_flag) THEN
                 IF (grank.EQ.0) 
     1           PRINT *, ' INITIAL JACOBIAN CHANGED SIGN!'
               END IF
               IF (grank.EQ.0) PRINT *,
     1            ' TRYING TO IMPROVE INITIAL MAGNETIC AXIS GUESS'
            END IF
#if defined (SKS)
            IF (PARVMEC) THEN
              CALL guess_axis_par (pr1, pz1, pru0, pzu0, lscreen)
            ELSE
#endif
              CALL guess_axis (r1, z1, ru0, zu0)
#if defined (SKS)
            ENDIF
#endif            
            lreset_internal = .true.
            ijacob = 1
            irst = 2
            GOTO 20
         ELSE IF (ier_flag.NE.norm_term_flag .AND. 
     1            ier_flag.NE.successful_term_flag) THEN
            RETURN
         ENDIF

#ifdef _ANIMEC
      w0 = wb + wpar/(gamma-one)
#else
      w0 = wb + wp/(gamma - one)
#endif

!
!     ADDITIONAL STOPPING CRITERION (set liter_flag to FALSE)
!
         IF (ijacob .eq. 25) THEN
            irst = 2
            CALL restart_iter(delt0r)
!            delt0r = p98*delt   !changed in restart
            IF (lscreen) PRINT 120, delt0r
            irst = 1
            GOTO 1000
         ELSE IF (ijacob .eq. 50) THEN
            irst = 2
            CALL restart_iter(delt0r)
!            delt0r = p96*delt
            IF (lscreen) PRINT 120, delt0r
            irst = 1
            GOTO 1000
         ELSE IF (ijacob .ge. 75) THEN
            ier_flag = jac75_flag
            liter_flag = .false.
         ELSE IF (iter2.ge.niter .and. liter_flag) THEN
            ier_flag = more_iter_flag
            liter_flag = .false.
         ENDIF

!       Store force residual, wdot for plotting
         wdota = ABS(w0 - w1)/w0
         CALL FLUSH(6)
#if defined(SKS)
         CALL MPI_Bcast(r00,1,MPI_REAL8,0,NS_COMM,MPI_ERR)
#endif
         r0dot = ABS(r00 - r00s)/r00
         r00s = r00
         w1 = w0
         IF (ivac .eq. 1) THEN
            IF (grank .EQ. 0) THEN
              IF (lscreen) PRINT 110, iter2
#if defined (SKS)
              WRITE (nthreed, 110) iter2
#endif
            END IF
            ivac = ivac + 1
         ENDIF

!     NOTE: PRINTOUT clobbers gc!
!     Increment time step and printout every nstep iterations
         IF (MOD(iter2,nstep).eq.0 .or. iter2.eq.1 .or.
     1       .not.liter_flag) CALL printout(iter2, delt0r, w0, lscreen)
         iter2 = iter2 + 1
!         iterc = iterc + 1
! JDH 2012-06-20 ^^^ iterc is a cumulative iteration counter. Used in V3FIT.
!   Never reset to 1

! JDH 2010-08-03: Call to vmec_history_store moved here from evolve.f
!  Stores fsq values and other, for later post-processing
         CALL vmec_history_store(delt0r)
         CALL flush(6)
!
!       STORE FSQ FOR PLOTTING. EVENTUALLY, STORE FOR EACH RADIAL MESH
!
         IF (MOD(iter2,niter/nstore_seq + 1).eq.0 .and. ns.eq.
     1      ns_array(multi_ns_grid)) THEN
            IF (itfsq .lt. nstore_seq) THEN
              itfsq = itfsq + 1
              fsqt(itfsq) = fsqr + fsqz
              wdot(itfsq) = MAX(wdota,c1pm13)
            END IF
         END IF
      END DO iter_loop

!SPH (021711): V3FITA - SAVE STATE FOR RESTART IF PRECONDITIONER IS ON

       IF (l_v3fit) THEN

!JDH 2011-09-14. Correct logic error.

!         IF (ictrl_prec2d .eq. 0) THEN
!            lqmr = (itype_precon .ge. 2)
!         ELSE
         IF (ictrl_prec2d .gt. 0) THEN
            CALL restart_iter(delt0r)
         END IF
      END IF

      IF (lSwap2Disk) CALL DeleteDAFile(ScratchFile)

      IF (grank .EQ. 0) THEN
	   WRITE (nthreed, 60) w0*twopi**2, wdota, r0dot
         IF (lrecon) WRITE (nthreed, 70) r00*fsqsum0/wb
         IF (nfcn .GT. 0) WRITE (nthreed, 80) nfcn
	END IF

#if defined(SKS)
      CALL second0(teqsoloff)
      eqsolve_time = eqsolve_time + (teqsoloff-teqsolon)
#endif

   60 FORMAT(/,' MHD Energy = ',1p,e13.6,3x, 'd(ln W)/dt = ',1p,e10.3,
     1         3x,'d(ln R0)/dt = ',e10.3)
   70 FORMAT(' Average radial force balance: Int[FR(m=0)]',
     1   '/Int(B**2/R) = ',1p,e12.5,' (should tend to zero)'/)
   80 FORMAT(' Function calls in GMRES: ',i5)
  110 FORMAT(/,2x,'VACUUM PRESSURE TURNED ON AT ',i4,' ITERATIONS'/)
  120 FORMAT(2x,'HAVING A CONVERGENCE PROBLEM: RESETTING DELT TO ',f8.3,
     1  /,2x,'If this does NOT resolve the problem, try changing ',
     2       '(decrease OR increase) the value of DELT')

      END SUBROUTINE eqsolve
