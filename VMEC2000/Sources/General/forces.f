      SUBROUTINE forces
      USE vmec_main, p5 => cp5
      USE realspace
      USE vforces
      IMPLICIT NONE
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(rprec), PARAMETER :: p25 = p5*p5, dphids=p25
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: l, ndim, lk, js, jsa, jsb
      REAL(rprec), DIMENSION(:), POINTER :: lv_e, lu_e, lu_o,
     1    bsqr, gvvs, guvs, guus
C-----------------------------------------------
      ndim = 1+nrzt
!
!     POINTER ALIASES
!
!      gcon => z1(:,0)
      lv_e => crmn_e; lu_e => czmn_e; lu_o => czmn_o

      bsqr => extra1(:,1);  gvvs => extra2(:,1)
      guvs => extra3(:,1);  guus => extra4(:,1)

!
!     ON ENTRY, ARMN=ZU,BRMN=ZS,AZMN=RU,BZMN=RS,LU=R*BSQ,LV = BSQ*SQRT(G)/R12
!     HERE, XS (X=Z,R) DO NOT INCLUDE DERIVATIVE OF EXPLICIT SQRT(S)
!     BSQ = |B|**2/2 + p
!     GIJ = (BsupI * BsupJ) * SQRT(G)  (I,J = U,V)
!     IT IS ESSENTIAL THAT LU,LV AT j=1 ARE ZERO INITIALLY
!
!     SOME OF THE BIGGER LOOPS WERE SPLIT TO FACILITATE CACHE
!     HITS, PIPELINING ON RISCS
!
!     FOR OPTIMIZATION ON CRAY, MUST USE COMPILER DIRECTIVES TO
!     GET VECTORIZATION OF LOOPS INVOLVING POINTERS!
!
!
!     ORIGIN OF VARIOUS TERMS
!
!     LU :  VARIATION OF DOMINANT .5*(RU-odd*Zodd - ZU-odd*Rodd) TERM
!           IN JACOBIAN
!
!     LV :  VARIATION OF R-TERM IN JACOBIAN
!
!     GVV:  VARIATION OF R**2-TERM AND Rv**2,Zv**2 IN gvv
!
!     GUU, GUV: VARIATION OF Ru, Rv, Zu, Zv IN guu, guv
!
      lu_e(1:ndim:ns) = 0; lv_e(1:ndim:ns) = 0
      guu(1:ndim:ns)  = 0; guv(1:ndim:ns)  = 0; gvv(1:ndim:ns) = 0
      armn_o(ndim)  = armn_e(ndim) *shalf(ndim)
      azmn_o(ndim)  = azmn_e(ndim) *shalf(ndim)
      brmn_o(ndim)  = brmn_e(ndim) *shalf(ndim)
      bzmn_o(ndim)  = bzmn_e(ndim) *shalf(ndim)

!$omp parallel
!$omp do
!$omp+ private(lk,jsa,jsb)
      DO lk=1,nznt
      jsa = 1 + (lk-1) * ns
      jsb = jsa -1 + ns
      guus(jsa:jsb) = guu(jsa:jsb)*shalf(jsa:jsb)
      guvs(jsa:jsb) = guv(jsa:jsb)*shalf(jsa:jsb)
      gvvs(jsa:jsb) = gvv(jsa:jsb)*shalf(jsa:jsb)

      armn_e(jsa:jsb)  = ohs*armn_e(jsa:jsb) * lu_e(jsa:jsb)
      azmn_e(jsa:jsb)  =-ohs*azmn_e(jsa:jsb) * lu_e(jsa:jsb)
      brmn_e(jsa:jsb)  = brmn_e(jsa:jsb) * lu_e(jsa:jsb)
      bzmn_e(jsa:jsb)  =-bzmn_e(jsa:jsb) * lu_e(jsa:jsb)
      bsqr(jsa:jsb)    = dphids*lu_e(jsa:jsb)/shalf(jsa:jsb)

      armn_o(jsa:jsb)  = armn_e(jsa:jsb) *shalf(jsa:jsb)
      azmn_o(jsa:jsb)  = azmn_e(jsa:jsb) *shalf(jsa:jsb)
      brmn_o(jsa:jsb)  = brmn_e(jsa:jsb) *shalf(jsa:jsb)
      bzmn_o(jsa:jsb)  = bzmn_e(jsa:jsb) *shalf(jsa:jsb)

      END DO
!$omp end do
!$omp end parallel

!     CONSTRUCT CYLINDRICAL FORCE KERNELS
!     NOTE: presg(ns+1) == 0, AND WILL BE "FILLED IN" AT EDGE
!     FOR FREE-BOUNDARY BY RBSQ
!
      DO lk=1,nznt
!DIR$ IVDEP
        DO js = 1,ns
         l = js + (lk-1) * ns
         guu(l) = p5*(guu(l) + guu(l+1))
         gvv(l) = p5*(gvv(l) + gvv(l+1))
         bsqr(l) = bsqr(l) + bsqr(l+1)
         guus(l) = p5*(guus(l) + guus(l+1))
         gvvs(l) = p5*(gvvs(l) + gvvs(l+1))
         armn_e(l) = armn_e(l+1) - armn_e(l) + p5*(lv_e(l) + lv_e(l+1))
         azmn_e(l) = azmn_e(l+1) - azmn_e(l)
         brmn_e(l) = p5*(brmn_e(l) + brmn_e(l+1))
         bzmn_e(l) = p5*(bzmn_e(l) + bzmn_e(l+1))
        END DO
      END DO
!$omp parallel
!$omp do
!$omp+ private(lk,jsa,jsb)
      DO lk=1,nznt
      jsa = 1 + (lk-1) * ns
      jsb = jsa -1 + ns
      armn_e(jsa:jsb) = armn_e(jsa:jsb) - (gvvs(jsa:jsb)*r1(jsa:jsb,1)
     1              + gvv(jsa:jsb)*r1(jsa:jsb,0))
      brmn_e(jsa:jsb) = brmn_e(jsa:jsb) + bsqr(jsa:jsb)*z1(jsa:jsb,1)
     1       -(guus(jsa:jsb)*ru(jsa:jsb,1) + guu(jsa:jsb)*ru(jsa:jsb,0))
      bzmn_e(jsa:jsb) = bzmn_e(jsa:jsb) - (bsqr(jsa:jsb)*r1(jsa:jsb,1)
     1       + guus(jsa:jsb)*zu(jsa:jsb,1) + guu(jsa:jsb)*zu(jsa:jsb,0))
      lv_e(jsa:jsb) = lv_e(jsa:jsb)*shalf(jsa:jsb)
      lu_o(jsa:jsb) = dphids*lu_e(jsa:jsb)
      END DO
!$omp end do
!$omp end parallel

      DO lk=1,nznt
!DIR$ IVDEP
        DO js = 1,ns
         l = js + (lk-1) * ns
         armn_o(l) = armn_o(l+1) - armn_o(l) - zu(l,0)*bsqr(l)
     1             + p5*(lv_e(l) + lv_e(l+1))
         azmn_o(l) = azmn_o(l+1) - azmn_o(l) + ru(l,0)*bsqr(l)
         brmn_o(l) = p5*(brmn_o(l) + brmn_o(l+1))
         bzmn_o(l) = p5*(bzmn_o(l) + bzmn_o(l+1))
         lu_o(l)   = lu_o(l) + lu_o(l+1)
        END DO
      END DO

!$omp parallel
!$omp do
!$omp+ private(lk,jsa,jsb)
      DO lk=1,nznt
      jsa = 1 + (lk-1) * ns
      jsb = jsa -1 + ns
      guu(jsa:jsb)  = guu(jsa:jsb) * sqrts(jsa:jsb)**2
      bsqr(jsa:jsb) = gvv(jsa:jsb) * sqrts(jsa:jsb)**2

      armn_o(jsa:jsb) = armn_o(jsa:jsb) - (zu(jsa:jsb,1)*lu_o(jsa:jsb)
     1      + bsqr(jsa:jsb)*r1(jsa:jsb,1) + gvvs(jsa:jsb)*r1(jsa:jsb,0))
      azmn_o(jsa:jsb) = azmn_o(jsa:jsb) + ru(jsa:jsb,1)*lu_o(jsa:jsb)
      brmn_o(jsa:jsb) = brmn_o(jsa:jsb) + z1(jsa:jsb,1)*lu_o(jsa:jsb)
     1       -(guu(jsa:jsb)*ru(jsa:jsb,1) + guus(jsa:jsb)*ru(jsa:jsb,0))
      bzmn_o(jsa:jsb) = bzmn_o(jsa:jsb) - (r1(jsa:jsb,1)*lu_o(jsa:jsb)
     1       + guu(jsa:jsb)*zu(jsa:jsb,1) + guus(jsa:jsb)*zu(jsa:jsb,0))
      END DO
!$omp end do
!$omp end parallel

      IF (lthreed) THEN
      DO lk=1,nznt
!DIR$ IVDEP
        DO js = 1,ns
         l = js + (lk-1) * ns
            guv(l)  = p5*(guv(l) + guv(l+1))
            guvs(l) = p5*(guvs(l) + guvs(l+1))
        END DO
      END DO

!$omp parallel
!$omp do
!$omp+ private(lk,jsa,jsb)
      DO lk=1,nznt
       jsa = 1 + (lk-1) * ns
       jsb = jsa -1 + ns
         brmn_e(jsa:jsb) = brmn_e(jsa:jsb) 
     1      - (guv(jsa:jsb)*rv(jsa:jsb,0) + guvs(jsa:jsb)*rv(jsa:jsb,1))
         bzmn_e(jsa:jsb) = bzmn_e(jsa:jsb) 
     1      - (guv(jsa:jsb)*zv(jsa:jsb,0) + guvs(jsa:jsb)*zv(jsa:jsb,1))
         crmn_e(jsa:jsb) = guv(jsa:jsb) *ru(jsa:jsb,0) 
     1                 + gvv(jsa:jsb) *rv(jsa:jsb,0)
     2      + gvvs(jsa:jsb)*rv(jsa:jsb,1) + guvs(jsa:jsb)*ru(jsa:jsb,1)
         czmn_e(jsa:jsb) = guv(jsa:jsb) *zu(jsa:jsb,0) 
     1                 + gvv(jsa:jsb) *zv(jsa:jsb,0)
     2      + gvvs(jsa:jsb)*zv(jsa:jsb,1) + guvs(jsa:jsb)*zu(jsa:jsb,1)
         guv(jsa:jsb) = guv(jsa:jsb) *sqrts(jsa:jsb)*sqrts(jsa:jsb)
         brmn_o(jsa:jsb) = brmn_o(jsa:jsb) 
     1      - (guvs(jsa:jsb)*rv(jsa:jsb,0) + guv(jsa:jsb)*rv(jsa:jsb,1))
         bzmn_o(jsa:jsb) = bzmn_o(jsa:jsb) 
     1      - (guvs(jsa:jsb)*zv(jsa:jsb,0) + guv(jsa:jsb)*zv(jsa:jsb,1))
         crmn_o(jsa:jsb) = guvs(jsa:jsb)*ru(jsa:jsb,0) 
     1                 + gvvs(jsa:jsb)*rv(jsa:jsb,0)
     2       + bsqr(jsa:jsb)*rv(jsa:jsb,1) + guv(jsa:jsb) *ru(jsa:jsb,1)
         czmn_o(jsa:jsb) = guvs(jsa:jsb)*zu(jsa:jsb,0) 
     1                 + gvvs(jsa:jsb)*zv(jsa:jsb,0)
     2       + bsqr(jsa:jsb)*zv(jsa:jsb,1) + guv(jsa:jsb) *zu(jsa:jsb,1)
      END DO
!$omp end do
!$omp end parallel
       ENDIF

!
!     COMPUTE CONSTRAINT FORCE KERNELS
!
!$omp parallel
!$omp do
!$omp+ private(lk,jsa,jsb)
      DO lk=1,nznt
      jsa = 1 + (lk-1) * ns
      jsb = jsa -1 + ns
      rcon(jsa:jsb,0) = (rcon(jsa:jsb,0) - rcon0(jsa:jsb))*gcon(jsa:jsb)
      zcon(jsa:jsb,0) = (zcon(jsa:jsb,0) - zcon0(jsa:jsb))*gcon(jsa:jsb)
      brmn_e(jsa:jsb) = brmn_e(jsa:jsb) + rcon(jsa:jsb,0)
      bzmn_e(jsa:jsb) = bzmn_e(jsa:jsb) + zcon(jsa:jsb,0)
      brmn_o(jsa:jsb) = brmn_o(jsa:jsb)+ rcon(jsa:jsb,0)*sqrts(jsa:jsb)
      bzmn_o(jsa:jsb) = bzmn_o(jsa:jsb)+ zcon(jsa:jsb,0)*sqrts(jsa:jsb)
      rcon(jsa:jsb,0) =  ru0(jsa:jsb) * gcon(jsa:jsb)
      zcon(jsa:jsb,0) =  zu0(jsa:jsb) * gcon(jsa:jsb)
      rcon(jsa:jsb,1) = rcon(jsa:jsb,0) * sqrts(jsa:jsb)
      zcon(jsa:jsb,1) = zcon(jsa:jsb,0) * sqrts(jsa:jsb)
      END DO
!$omp end do
!$omp end parallel
!
!     ASSIGN EDGE FORCES (JS = NS) FOR FREE BOUNDARY CALCULATION
!
      IF (ivac .ge. 1) THEN
         armn_e(ns:nrzt:ns) = armn_e(ns:nrzt:ns) 
     1                      + zu0(ns:nrzt:ns)*rbsq(1:nznt)
         armn_o(ns:nrzt:ns) = armn_o(ns:nrzt:ns) 
     1                      + zu0(ns:nrzt:ns)*rbsq(1:nznt)
         azmn_e(ns:nrzt:ns) = azmn_e(ns:nrzt:ns) 
     1                      - ru0(ns:nrzt:ns)*rbsq(1:nznt)
         azmn_o(ns:nrzt:ns) = azmn_o(ns:nrzt:ns) 
     1                      - ru0(ns:nrzt:ns)*rbsq(1:nznt)
!         fz00_edge = SUM(wint(ns:nrzt:ns)*ru0(ns:nrzt:ns)*rbsq(1:nznt))
      ENDIF

 100  CONTINUE

      END SUBROUTINE forces
