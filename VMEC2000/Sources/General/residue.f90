      SUBROUTINE residue (gcr, gcz, gcl)
      USE vmec_main, p5 => cp5
      USE vmec_params, ONLY: rss, zcs, rsc, zcc,                        &
                             meven, modd, ntmax, signgs
      USE realspace, ONLY: phip
      USE vsvd
      USE xstuff
      USE precon2d
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor,0:mpol1,ntmax), INTENT(inout) :: &
        gcr, gcz, gcl
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER, PARAMETER :: n0 = 0, m0 = 0, m1 = 1
      INTEGER, PARAMETER :: n3d = 0, nasym = 1
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------
      INTEGER :: nsfix, iflag, jedge, delIter
      REAL(rprec) :: r1, tnorm
!-----------------------------------------------
!
!     IMPOSE M=1 MODE CONSTRAINT TO MAKE THETA ANGLE
!     INVARIANT TO PHI-SHIFTS (AND THETA SHIFTS FOR ASYMMETRIC CASE)
!     ( ZCS = RSS, ZCC = RSC ARE THE CORRECT POLAR RELATIONS)
!
!     SYMMETRIC PERTURBATIONS (BASED ON POLAR RELATIONS):
!        RSS(n) = ZCS(n), n != 0
!     ASYMMETRIC PERTURBATIONS:
!        RSC(n) = ZCC(n), all n 
!
!     INTERNALLY:
!        XC(rss) = .5*(Rss + Zcs), XC(zcs) = .5*(Rss - Zcs)
!        XC(rsc) = .5*(Rsc + Zcc), XC(zcc) = .5*(Rsc - Zcc)
!     THIS IMPLIES THE CONSTRAINT
!        3D ONLY : GC(zcs) = 0;  
!        3D, ASYM: GC(zcc) = 0
!
      IF (lthreed) CALL constrain_m1(gcr(:,:,1,rss), gcz(:,:,1,zcs))
      IF (lasym)   CALL constrain_m1(gcr(:,:,1,rsc), gcz(:,:,1,zcc))

!     PRECONDITIONER MUST BE CALCULATED USING RAW (UNPRECONDITIONED) FORCES
      IF (ictrl_prec2d .ge. 2) RETURN

!
!     PUT FORCES INTO PHIFSAVE UNITS USED BY PRECONDITIONERS, FNORM
!
      IF (phifac .eq. zero) THEN
         STOP 'phifac = 0 in residue'
      ELSE
         tnorm = phifsave/phifac           !put all forces into phifac=phifsave units
      END IF


      IF (lrecon) THEN
!
!       MOVE R(n=0,m=0) TO SATISFY LIMITER OR AXIS POSITION
!       USE XC(NEQS2) TO STORE THIS PERTURBATION
!       TO SATISFY FORCE BALANCE AT JS=1, ADJUST PFAC IN RADFOR
!       ALSO, SCALE TOROIDAL FLUX AT EDGE TO MATCH MINOR RADIUS

        r1 = SUM(gcr(:ns,n0,m0,1))
        fsqsum0 = signgs*hs*r1/r0scale
        nsfix = 1                   !fix origin for reconstruction mode
        gcr = gcr * tnorm**2
        gcz = gcz * tnorm**2
        gcl = gcl * tnorm
        IF (iopt_raxis.gt.0 .and. iresidue.eq.2                        &
           .and. fsq.lt.fopt_axis) iresidue = 3
        IF (iresidue .lt. 3) gcr(nsfix,n0,m0,1) = zero
      ELSE
!
!     ADJUST PHIEDGE
!
         IF (imovephi .gt. 0) CALL movephi1 (gphifac)
      ENDIF
      gc(neqs1) = gphifac

!
!     COMPUTE INVARIANT RESIDUALS
!
      r1 = one/(2*r0scale)**2
      jedge = 0    
!SPH-JAH013108: MUST INCLUDE EDGE FORCE (INITIALLY) FOR V3FITA TO WORK
!ADD A V3FIT RELATED FLAG? ADD fsq criterion first
      delIter = iter2-iter1

      IF (l_v3fit) THEN
!  Coding for when run by V3FIT. Needed for correct computation
!  of partial derivatives
         IF (iter2-iter1.lt.50) jedge = 1
      ELSE
!  Coding for VMEC2000 run stand-alone
         IF (delIter.lt.50 .and.                                        &
            (fsqr+fsqz).lt.1.E-6_dp) jedge = 1
      ENDIF

      CALL getfsq (gcr, gcz, fsqr, fsqz, r1*fnorm, jedge)

      fsql = fnormL*SUM(gcl*gcl)
      fedge = r1*fnorm*SUM(gcr(ns,:,:,:)**2 + gcz(ns,:,:,:)**2)

!
!     PERFORM PRECONDITIONING AND COMPUTE RESIDUES
!
      IF (ictrl_prec2d .eq. 1) THEN

         CALL block_precond(gc)

         IF (.not.lfreeb .and. ANY(gcr(ns,:,:,:) .ne. zero))            &
            STOP 'gcr(ns) != 0 for fixed boundary in residue'
         IF (.not.lfreeb .and. ANY(gcz(ns,:,:,:) .ne. zero))            &
            STOP 'gcz(ns) != 0 for fixed boundary in residue'
         IF (ANY(gcl(:,1:,0,zsc) .ne. zero))                            &
            STOP 'gcl(m=0,n>0,sc) != 0 in residue'
         IF (lthreed) THEN
            IF (ANY(gcl(:,n0,:,zcs) .ne. zero))                         &
            STOP 'gcl(n=0,m,cs) != 0 in residue'
         END IF

         fsqr1 = SUM(gcr*gcr)
         fsqz1 = SUM(gcz*gcz)
         fsql1 = SUM(gcl*gcl)

      ELSE
!        m = 1 constraint scaling
         IF (lthreed) CALL scale_m1(gcr(:,:,1,rss), gcz(:,:,1,zcs))
         IF (lasym)   CALL scale_m1(gcr(:,:,1,rsc), gcz(:,:,1,zcc))

         iflag = 0
         CALL scalfor (gcr, arm, brm, ard, brd, crd, iflag)
         iflag = 1
         CALL scalfor (gcz, azm, bzm, azd, bzd, crd, iflag)

         CALL getfsq (gcr, gcz, fsqr1, fsqz1, fnorm1, m1)
!OLD         CALL getfsq (gcr, gcz, fsqr1, fsqz1, one, m1)

         gcl = faclam*gcl
         fsql1 = hs*SUM(gcl*gcl)

!
!     SOFT-STARTUP AVOIDS SOME "HANG-UPS" FOR LARGE RESIDUALS (04/2002)
!
!         fac = one / (one+fsqr+fsqz+fsql)
!         gcr = fac*gcr
!         gcz = fac*gcz

      ENDIF

      END SUBROUTINE residue

      SUBROUTINE constrain_m1(gcr, gcz)
      USE vmec_main, p5 => cp5 
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor), INTENT(inout) :: gcr, gcz
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      REAL, PARAMETER :: FThreshold = 1.E-6_dp
      REAL(rprec) :: temp(ns,0:ntor)
!-----------------------------------------------
!
!     COMPUTE INTERNAL gr, gz
!     NOTE: internal gz => 0 for both values of lconm1 (although gz is different)
!     FOR lconm1=T, gcr(internal) = gcr+gcz, gcz(internal) = gcr-gcz->0
!
      IF (lconm1) THEN
         temp = gcr
         gcr = osqrt2*(gcr + gcz)
         gcz = osqrt2*(temp- gcz)
      ELSE
         gcz = 0
      END IF

      IF (fsqz.lt.FThreshold) gcz = 0
 
      END SUBROUTINE constrain_m1

      SUBROUTINE scale_m1(gcr, gcz)
      USE vmec_main
      IMPLICIT NONE
!-----------------------------------------------
!   D u m m y   A r g u m e n t s
!-----------------------------------------------
      REAL(rprec), DIMENSION(ns,0:ntor), INTENT(inout) :: gcr, gcz
!-----------------------------------------------
!   L o c a l   P a r a m e t e r s
!-----------------------------------------------
      INTEGER :: n
      REAL(rprec) :: fac(ns)
!-----------------------------------------------
      IF (.not.lconm1) RETURN

      fac = osqrt2*ard(:ns,2)/(ard(:ns,2) + azd(:ns,2))
      DO n = 0, ntor
         gcr(:,n) = fac*gcr(:,n)
      END DO

      fac = osqrt2*azd(:ns,2)/(ard(:ns,2) + azd(:ns,2))
      DO n = 0, ntor
         gcz(:,n) = fac*gcz(:,n)
      END DO
  
      END SUBROUTINE scale_m1
