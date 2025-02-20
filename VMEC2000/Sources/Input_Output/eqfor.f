      SUBROUTINE eqfor(br, bz, bsubu, bsubv, tau, rzl_array, ier_flag)
      USE vmec_main
#ifdef _FLOW
     1              , machsq=>bcrit
#endif
      USE vmec_params
      USE realspace
      USE vforces, r12 => armn_o, bsupu => crmn_e, bsupv => czmn_e,
     1   gsqrt => azmn_o, bsq => bzmn_o, izeta => azmn_e, 
     2   brho => bzmn_e, bphi => czmn_o, curtheta => brmn_e
#ifdef _ANIMEC
     3  ,tau_an => brmn_o
#endif
      USE vacmod
      USE vsvd
      USE vspline
      USE csplinx
      USE vmec_io
      USE mgrid_mod
      USE fbal
      USE v3f_vmec_comm
      USE stel_constants, ONLY: pi

      USE parallel_include_module
#if defined(SKS)
      USE xstuff, ONLY: xc, pxc 
#else
      USE xstuff, ONLY: xc
#endif
      USE safe_open_mod
      USE timer_sub
      IMPLICIT NONE
C-----------------------------------------------
C   D u m m y   A r g u m e n t s
C-----------------------------------------------
      INTEGER :: ier_flag
      REAL(dp), DIMENSION(ns,nznt,0:1), INTENT(in) :: bsubu, bsubv
      REAL(dp), DIMENSION(nrzt), INTENT(out) :: br, bz
      REAL(dp), DIMENSION(nrzt), INTENT(out) :: tau
      REAL(dp), DIMENSION(ns,0:ntor,0:mpol1,3*ntmax),
     1  TARGET, INTENT(in) :: rzl_array
C-----------------------------------------------
C   L o c a l   P a r a m e t e r s
C-----------------------------------------------
      REAL(dp), PARAMETER :: c1p5=1.5_dp
C-----------------------------------------------
C   L o c a l   V a r i a b l e s
C-----------------------------------------------
      INTEGER :: i, icount, itheta, js, js1, l, loff,
     1   lpi, lt, n, n1, nchicur, nchiiota0, noff,
     2   nout, nsort, iv, iu, lk, nplanes
      REAL(dp), DIMENSION(:), POINTER ::
     1   rmags, zmags, rmaga, zmaga
      REAL(dp), DIMENSION(:,:,:), POINTER :: rmncc,zmnsc
      REAL(dp), DIMENSION(ns) :: phi1, chi1, jPS2
      REAL(dp) :: modb(nznt)
      REAL(dp), DIMENSION(:), ALLOCATABLE ::
     1   btor_vac, btor1, dbtor, phat, t12u, guu_1u, surf_area, 
     2   r3v, redge, rbps1u, bpol2vac, phipf_loc
      REAL(dp) :: aminr1, aminr2, aminr2in, anorm,
     1   aspectratio, betai, betstr, scaling_ratio,
     2   bminz2, bminz2in, btor, iotamax, musubi,
     3   bzcalc, bzin, chisq, chiwgt, cur0,
     4   delphid_exact, delta1, delta2, delta3, denwgt, lambda,
     5   dlogidsi, dmusubi_meas, er, es, fac, facnorm, factor, fgeo,
     6   fmax, fmin, flao, fpsi0, pavg, pitchc, pitchm,
     7   pprime, qedge, qmin1, qmin2, qmin3, qzero,
     8   raxis0, rcalc, rcen, rcenin, rgeo, rjs, 
     9   rjs1, rlao, rqmin1, rqmin2, rshaf, rshaf1, rshaf2, s11, s12,
     A   s13, s2, s3, sigr0, sigr1, sigz1, smaleli,
     B   splintx, splints, sqmin, sumbpol, sumbtot, sumbtor, sump,
     C   sump2, sump20, t1, tz, jpar_perp=0, jparPs_perp=0,
     D   tol, toroidal_flux, vnorm, vprime, wght0, xmax,
     E   xmida, xmidb, xmin, rzmax, rzmin, zxmax, zxmin, zaxis0,
     F   zmax, zmin, yr1u, yz1u, waist(2), height(2)
#ifdef _ANIMEC
     G  ,sumpar, sumper, sumpp2, sumpb2, betpar, betper, bet_s_to,
     H   bet_s_ew, bmax_t, bmin_t, betplp, betplb
#elif defined _FLOW
     I  ,betath, betarot, betato, sumrot, sumrt2, betrstr
#endif
      REAL(dp) :: d_of_kappa, tmpxc, rmssum
      INTEGER :: istat1, OFU, j, k
C-----------------------------------------------
C   E x t e r n a l   F u n c t i o n s
C-----------------------------------------------
       EXTERNAL splintx,splints
C-----------------------------------------------
      CALL second0 (teqfon)
!
!     POINTER ASSOCIATIONS
!
      rmags => rzl_array(1,:,0,rcc)
      zmags => rzl_array(1,:,0,zcs+ntmax)
      rmncc => rzl_array(:,:,:,rcc)
      zmnsc => rzl_array(:,:,:,zsc+ntmax) 
      IF (lasym) THEN
        rmaga => rzl_array(1,:,0,rcs)
        zmaga => rzl_array(1,:,0,zcc+ntmax)
      END IF

!     crmn_o => bss on half grid
      CALL bss (r12, bzmn, brmn, azmn, armn, crmn_o, bsupu, bsupv,
     1          br, bphi, bz)

!
!     STORE EDGE VALUES OF B-FIELD
!
      IF ((lfreeb .and. ivac.gt.1) .or. ledge_dump) THEN
         IF (ALLOCATED(bredge)) DEALLOCATE (bredge, bpedge, bzedge)
         ALLOCATE (bredge(2*nznt), bpedge(2*nznt), bzedge(2*nznt),
     1             stat=i)
         IF (i .ne. 0) STOP 'Error in EQFOR allocating bredge'
         DO iv = 1,nzeta
            DO iu = 1,ntheta3
               lk = iv + nzeta*(iu-1)
               n1 = ns*lk
               bredge(lk) = 1.5_dp*br(n1)   - p5*br(n1-1)
               bpedge(lk) = 1.5_dp*bphi(n1) - p5*bphi(n1-1)
               bzedge(lk) = 1.5_dp*bz(n1)   - p5*bz(n1-1)
            END DO
         END DO
      END IF

#ifdef _ANIMEC
!EVALUATE MIRROR STABILITY CRITERION; BSQ ==> MAGNETIC PRESSURE
      CALL mirror_crit(tau_an,bsq)
#elif defined _FLOW
        betath = wp/wb
        betarot = wrot/wb
        betato = betath + betarot
        write (nthreed,107) betath,betarot,betato
 107    format(" thermal beta =",1pe13.4," beta-rotation =",1pe13.4,
     &  " total beta =",1pe13.4)
#endif
!
!     NOTE: JXBFORCE ROUTINE MUST BE CALLED TO COMPUTE IZETA, JDOTB
!           ON OUTPUT, J, IZETA, JDOTB ARE IN MKS UNITS (1/MU0 FACTOR)
!
!     CAUTION: THIS CALL WILL WRITE OVER br, bz 
!

      CALL jxbforce (bsupu, bsupv, bsubu, bsubv, crmn_o, rcon, zcon, 
     1        gsqrt, bsq, curtheta, izeta, brho, sigma_an, ier_flag
#ifdef _ANIMEC
     2       ,pp1, pp2, ppar, onembc
#elif defined _FLOW
     3       ,pp1, pp2, prot, protrsq
#endif
     4              )

!
!     HALF-MESH VOLUME-AVERAGED BETA
!

      tau(1) = 0
      tau(2:nrzt) = signgs*wint(2:nrzt)*gsqrt(2:nrzt)
      DO i = 2, ns
#ifdef _ANIMEC
!BSQ contains only magnetic pressure set in call to routine mirror_crit
         s2 = SUM(bsq(i:nrzt:ns)*tau(i:nrzt:ns))/vp(i)
#elif defined _FLOW
         s2 = SUM((bsq(i:nrzt:ns)-prot(i:nrzt:ns))*tau(i:nrzt:ns))/vp(i)
#else
         s2 = SUM(bsq(i:nrzt:ns)*tau(i:nrzt:ns))/vp(i) - pres(i)
#endif
         overr(i) = SUM(tau(i:nrzt:ns)/r12(i:nrzt:ns)) / vp(i)
         beta_vol(i) = pres(i)/s2
      END DO

      betaxis = c1p5*beta_vol(2) - p5*beta_vol(3)
#ifdef _FLOW
      betaxis = betaxis * exp(p5*machsq)
#endif


      IF (grank.EQ. 0) WRITE (nthreed, 5)
 5    FORMAT(/,' NOTE:  S=normalized toroidal flux (0 - 1)',/,
     1         '        U=poloidal angle (0 - 2*pi)',/,
     1         '        V=geometric toroidal angle (0 - 2*pi)',/,
     1         '       <RADIAL FORCE> = d(Ipol)/dPHI',
     1         ' - IOTA*d(Itor)/dPHI - dp/dPHI * d(VOL)/dPHI',/,
     1         '                      = d(VOL)/dPHI*[<JSUPU>',
     1         ' - IOTA*<JSUPV> - SIGN(JAC)*dp/dPHI]',/,
     1         '       (NORMED TO SUM OF INDIVIDUAL TERMS)',//,
     1         '      S     <RADIAL    TOROIDAL      IOTA     ',
     1         ' <JSUPU>    <JSUPV>     d(VOL)/',
#ifdef _ANIMEC
     2         '   d(P_||)/    <M>     PRESF    <BSUBU>    <BSUBV>',
#else
     2         '   d(PRES)/    <M>     PRESF    <BSUBU>    <BSUBV>',
#endif
     3         '      <J.B>      <B.B>',/,
     4         '             FORCE>      FLUX                  ',
     5         '                       d(PHI) ',
     6         '    d(PHI)                             ',/,148('-'),/)

      ALLOCATE (phipf_loc(ns))

      phipf_loc(1) = twopi*signgs*(c1p5*phip(2) - p5*phip(3))
      presf(1) = c1p5*pres(2) - p5*pres(3)
      DO i = 2,ns1
         presf(i) = p5*(pres(i) + pres(i+1))
         phipf_loc(i) = p5*twopi*signgs*(phip(i) + phip(i+1))
      END DO
      presf(ns) = c1p5*pres(ns)- p5*pres(ns-1)
      phipf_loc(ns) = twopi*signgs*(c1p5*phip(ns) - p5*phip(ns1))

      phi1(1) = zero
      chi1(1) = zero
      DO i = 2, ns
         phi1(i) = phi1(i-1) + hs*phip(i)
         chi1(i) = chi1(i-1) + hs*(phip(i)*iotas(i))
      END DO

      chi = twopi*chi1

!	WRITE (36, 201) (i, phi1(i), chi1(i), iotas(i), i=1, ns)
!201   FORMAT (i4, 1p,3e14.6)

      CALL calc_fbal(bsubu, bsubv)

#ifdef _FLOW
!pmap calculated in calc_fbal contains <p(s,R)> on full mesh 
      presf = signgs*pmap
      presf(ns) = 2 * presf(ns-1)   -  presf(ns-2)
      presf(1)  = 2 * presf(2)      -  presf(3)
#endif
      bucof(1) = 0
      bvcof(1) = c1p5*bvco(2) - p5*bvco(3)
!
!     NOTE:  jcuru, jcurv on FULL radial mesh coming out of calc_fbal
!            They are local (surface-averaged) current densities (NOT integrated in s)
!            jcurX = (dV/ds)/twopi**2 <JsupX>   for X=u,v
!
      DO i = 2,ns1
         equif(i) = equif(i)*vpphi(i)/(ABS(jcurv(i)*chipf(i))
     1            + ABS(jcuru(i)*phipf(i))+ABS(presgrad(i)*vpphi(i)))
         bucof(i) = p5*(buco(i) + buco(i+1))
         bvcof(i) = p5*(bvco(i) + bvco(i+1))
      END DO

      bucof(ns) = c1p5*buco(ns) - p5*buco(ns1)
      bvcof(ns) = c1p5*bvco(ns) - p5*bvco(ns1)

      equif(1) = two*equif(2) - equif(3)
      jcuru(1) = two*jcuru(2) - jcuru(3)
      jcurv(1) = two*jcurv(2) - jcurv(3)
      presgrad(1)  = two*presgrad(2) - presgrad(3)
      presgrad(ns) = two*presgrad(ns1) - presgrad(ns1-1)
      vpphi(1)  = two*vpphi(2) - vpphi(3)
      vpphi(ns) = two*vpphi(ns1) - vpphi(ns1-1)
      equif(ns) = two*equif(ns1) - equif(ns1-1)
      jcuru(ns) = two*jcuru(ns1) - jcuru(ns1-1)
      jcurv(ns) = two*jcurv(ns1) - jcurv(ns1-1)
!     NOTE: phipf = phipf_loc/(twopi), phipf_loc ACTUAL (twopi factor) Toroidal flux derivative
!     SPH/JDH (060211): remove twopi factors from <JSUPU,V> (agree with output in JXBOUT file)
      fac = twopi*signgs
      DO js = 1, ns
         es = (js - 1)*hs
         cur0 = fac*vpphi(js)*twopi              !==dV/ds = dV/dPHI * d(PHI/ds)  (V=actual volume)
         IF(rank.EQ.0) WRITE (nthreed, 30) es, equif(js), fac*phi1(js),
     1     iotaf(js), jcuru(js)/vpphi(js)/mu0, jcurv(js)/vpphi(js)/mu0,
     2     cur0/phipf_loc(js), presgrad(js)/phipf_loc(js)/mu0,
     3     specw(js), presf(js)/mu0, bucof(js), bvcof(js), jdotb(js), 
     4     bdotb(js)
      END DO
 30   FORMAT(1p,2e10.2,2e12.4,4e11.3,0p,f7.3,1p,5e11.3)

#ifdef _FLOW
      WRITE (nthreed, 10)
 10   FORMAT(//,148('-'),/,
     1    '        S           OMEGA      TEMPERATURE',
     2    '   U-FUNCTION     <p(s,R)>   <1/2rho_mV_tor^2>',/,148('-'),/)
      pmap(1)  = c1p5*pmap(2)    - p5*pmap(3)
      anorm    = pmap(ns)           !temporarily store this value of pmap
      pmap(ns) = c1p5*pmap(ns-1) - p5*pmap(ns-2)
      pd(1)    = c1p5*pd(2)      - p5*pd(3)
      pd(ns)   = c1p5*pd(ns-1)   - p5*pd(ns-2)
      DO js = 2, ns
         es = (js - c1p5)*hs
         WRITE (nthreed, 35) es, omega(js), tpotb(js), 
     1                   0.5*machsq*rotfot(js)/(r00*r00)
     2                  ,0.5*signgs*(pmap(js)+pmap(js-1))
     3                  ,0.5*signgs*(pd(js)  +pd(js-1))
      END DO
 35   FORMAT(1p, 1p6e14.6)
      pmap(ns) = anorm              !restore original value of pmap(ns)
#endif
      DEALLOCATE (phipf_loc)

!
!     MAKE SURE WOUT FILE DOES NOT REQUIRE ANY STUFF COMPUTED BELOW....
!
      IF (ier_flag .NE. successful_term_flag) RETURN

!
!     Calculate mean (toroidally averaged) poloidal cross section area & toroidal flux
!
      anorm = twopi*hs
      vnorm = twopi*anorm
      toroidal_flux = anorm * SUM(bsupv(2:nrzt)*tau(2:nrzt))

!
!     Calculate poloidal circumference and normal surface area and aspect ratio
!     Normal is | dr/du X dr/dv | = SQRT [R**2 guu + (RuZv - RvZu)**2]
!
      ALLOCATE (guu_1u(nznt), surf_area(nznt))
      guu_1u(:nznt) = ru0(ns:nrzt:ns)*ru0(ns:nrzt:ns) +
     1   zu0(ns:nrzt:ns)*zu0(ns:nrzt:ns)
      surf_area(:nznt) = wint(ns:nrzt:ns)*SQRT(guu_1u(:nznt))
      circum_p = twopi*SUM(surf_area(:nznt))
      surf_area(:nznt) = wint(ns:nrzt:ns)*SQRT(
     1     + (r1(ns:nrzt:ns,0) + r1(ns:nrzt:ns,1))**2*guu_1u(:nznt)
     2     +((rv(ns:nrzt:ns,0) + rv(ns:nrzt:ns,1))*zu0(ns:nrzt:ns) 
     3     - (zv(ns:nrzt:ns,0) + zv(ns:nrzt:ns,1))*ru0(ns:nrzt:ns))**2)
      surf_area_p = twopi**2*SUM(surf_area(:nznt))
      DEALLOCATE (guu_1u)

      aspect = aspectratio()

!     Also, estimate mean elongation of plasma from the following relations
!     for an axisymmetric torus with elliptical cross section and semi-axes
!     a and a * kappa (kappa >= 1)
!     
!     surf_area _p = 2*pi*R * 2*pi*a ctwiddle(kappa_p)
!     volume_p    = 2*pi*R * pi*a ** 2 * kappa_p
!     cross_area _p =   pi*a ** 2 * kappa_p
!
!     The cirumference of an ellipse of semi-axes a and a * kappa_p is 
!        2 * pi * a ctwiddle(kappa_p)
!     The exact form for ctwiddle is 4 E(1 - kappa_p^2) / (2 pi), where
!      E is the complete elliptic integral of the second kind 
!     (with parameter argument m, not modulus argument k)
!
!  The coding below implements an approximate inverse of the function
!  d(kappa) = ctwiddle(kappa) / sqrt(kappa)
!  The approximate inverse is
!	   kappa = 1 + (pi^2/8) * (d^2+sqrt(d^4-1)-1)
!  Note that the variable aminor_p, for an elliptic cross section, 
!  would be a * sqrt(kappa)
!
      d_of_kappa = surf_area_p * aminor_p / ( 2 * volume_p)
      kappa_p = 1 + (pi * pi / 8) * (d_of_kappa ** 2 + SQRT                    
     &      (ABS(d_of_kappa ** 4 - 1)) -1)
      vvc_kappa_p = kappa_p ! Save result for v3fit.

      aminr1 = 2*volume_p/surf_area_p

!
!
!     OUTPUT BETAS, INDUCTANCES, SAFETY FACTORS, ETC.
!     (EXTRACTED FROM FQ-CODE, 9-10-92)
!
!     b poloidals (cylindrical estimates)
!
      rcen = p5*(router + rinner)               !geometric center
      n = 0
      n1 = n + 1
      rcenin = DOT_PRODUCT(rmncc(ns,n1,:mpol1+1:2),
     1                     mscale(:mpol1:2)*nscale(n))

      l = (mpol1+1)/2
      ALLOCATE (t12u(l))
      t12u(:l) = mscale(1:mpol1:2)*nscale(n)
      aminr2in = DOT_PRODUCT(rmncc(ns,n1,2:mpol1+1:2),t12u(:l))
      bminz2in = DOT_PRODUCT(zmnsc(ns,n1,2:mpol1+1:2),t12u(:l))
      bminz2 = DOT_PRODUCT(zmnsc(ns,n1,2:mpol1+1:2),t12u(:l))
      DEALLOCATE (t12u)
      aminr1 = SQRT(two*volume_p/(twopi*twopi*r00))  !vol av minor radius
      aminr2 = p5*(router - rinner)             !geometric plasma radius
!
!       cylindrical estimates for beta poloidal
!
#ifdef _FLOW
      sump = vnorm*SUM(prot(2:nrzt)*tau(2:nrzt))
#else
      sump = vnorm*SUM(vp(2:ns)*pres(2:ns))
#endif
      pavg = sump/volume_p
      ppeak = presf(1)
      factor = 2*pavg
!
!       delphid_exact = Integral[ (Bvac - B) * dSphi ]
!       rshaf [= RT in Eq.(12), Phys Fluids B 5 (1993) 3119]
!
!       Note: tau = |gsqrt|*wint
!
      ALLOCATE (btor_vac(nznt), btor1(nznt), dbtor(nznt), 
     1          phat(nznt), redge(nznt))
      delphid_exact = zero                         !Eq. 20 in Shafranov
      musubi = zero
      rshaf1 = zero
      rshaf2 = zero
      DO js = 2, ns
         btor_vac(:nznt) = rbtor/r12(js:nrzt:ns)
         btor1(:nznt) = r12(js:nrzt:ns)*bsupv(js:nrzt:ns)
         delphid_exact = delphid_exact + SUM(
     1     (btor_vac(:nznt)/r12(js:nrzt:ns) -
     2      bsupv(js:nrzt:ns))*tau(js:nrzt:ns))
         dbtor(:nznt) = btor1(:nznt)**2 - btor_vac(:nznt)**2
         musubi = musubi - SUM(dbtor(:nznt)*tau(js:nrzt:ns))
         phat(:nznt) = bsq(js:nrzt:ns) - p5*btor_vac(:nznt)**2
         phat(:nznt) = (phat(:nznt) - dbtor(:nznt))*tau(js:nrzt:ns)
         rshaf1 = rshaf1 + SUM(phat(:nznt))
         rshaf2 = rshaf2 + SUM(phat(:nznt)/r12(js:nrzt:ns))
      END DO

      redge(:nznt) = r1(ns:nrzt:ns,0) + r1(ns:nrzt:ns,1)
      IF (lfreeb .and. ivac.gt.1) THEN
         phat = bsqvac(:) - p5*(bsubvvac/redge(:))**2
      ELSE
         phat = c1p5*bsq(ns:nrzt:ns) - p5*bsq(ns-1:nrzt:ns)
     1        - p5*(rbtor/redge(:))**2
      END IF

      DEALLOCATE (btor_vac, btor1, dbtor)

      delphid_exact = anorm*delphid_exact
      rshaf = rshaf1/rshaf2
      fpsi0 = c1p5*bvco(2) - p5*bvco(3)
      b0 = fpsi0/r00

      IF (lrecon) THEN
      IF (apres .ne. aminr1) THEN
         WRITE (*, 50) (apres/aminr1)**2
         IF(rank.EQ.0) WRITE (nthreed, 50) (apres/aminr1)**2
      ENDIF
   50 FORMAT(/'  Multiply final PHIEDGE by ',f7.3,
     1   ' to make apres = aminr1')
!
!     Output some measured data
!
      raxis0 = SUM(raxis_cc(0:ntor))
      zaxis0 = SUM(zaxis_cs(0:ntor))
      cur0 = signgs*iotaf(1)*fpsi0/r00**2*(dkappa + one/dkappa)
      IF (lpprof.AND.rank.EQ.0) WRITE (nthreed, 60) presfac*pfac,
     1   (-100.*(r00 - rthompeak))
      IF(rank.EQ.0) WRITE (nthreed, 65) b0, cur0/mu0
   60 FORMAT(//,' Input pressure scaled by ',f6.3/,
     1   ' Pressure profile shifted ',f6.2,' cm.',
     2   ' relative to magnetic axis')
   65 FORMAT(//' B-PHI(R=Raxis,Z=0) = ',f8.3,' [Wb/M**2]',/,
     1   ' J-PHI(R=Raxis,Z=0) = iota(0)*Bt0/(mu0*R0)*(1+k**2)/k = ',1p
     2   e10.3,' [A/M**2]',2/,' Comparison of Input vs Output',8x,
     3   'Input',9x,'VMEC Output',/,1x,71('-'))
      IF(rank.EQ.0) WRITE (nthreed, 70) 1.e-6_dp*currv/mu0,
     1   1.e-6_dp*ctor/mu0, phiedge, toroidal_flux, 
     1   1.e3_dp*phidiam, 1.e3_dp*delphid, 1.e-3_dp*pthommax,
     2   1.e-3_dp*ppeak/mu0, rcenin, rcen, raxis0, r00, zaxis0, z00,
     2   rc0mse,
     3   apres, aminr1, aminr2in, aminr2, bminz2in, bminz2, rinner,
     4   router, 1.e3_dp*delphid_exact
   70 FORMAT(' Toroidal Current     =',2(10x,f10.3),'  [MA]',/,
     1   ' Edge Toroidal Flux   =',2(10x,f10.3),'  [Wb]',/,
     2   ' Diamagnetic Flux(*)  =',2(10x,f10.3),'  [mWb]',/,
     3   ' Peak Pressure        =',2(10x,f10.3),'  [KPa]',/,
     4   ' Geometric Center     =',2(10x,f10.3),'  [M]',/,
     5   ' Magnetic R Axis      =',2(10x,f10.3),'  [M]',/,
     6   ' Magnetic Z Axis      =',2(10x,f10.3),'  [M]',/,
     7   ' MSE R Axis           =',10x,f10.3,20x,'  [M]',/,
     8   ' Minor Radius (apres) =',2(10x,f10.3),'  [M]',/,
     9   ' Minor Radius (a)     =',2(10x,f10.3),'  [M]',/,
     .   ' Minor Radius (b)     =',2(10x,f10.3),'  [M]',/,
     1   ' Inboard  Midplane R  =',30x,f10.3,'  [M]',/,
     2   ' Outboard Midplane R  =',30x,f10.3,'  [M]',/,50('-')/,
     3   ' * Exact diamagnetic flux (not based on equilibrium) = ',f8.3,
     4   '  [mWb]')

!       Calculate components of total chi-squared
      total_chi_cur = (currv - ctor)**2/sigma_current**2
      nchicur = 1
      IF (iphidiam .eq. 1) total_chi_dia = (phidiam - delphid)**2/
     1   sigma_delphid**2
      nchidia = 1
      END IF   !!IF(lrecon)

      rmax_surf = MAXVAL(r1(ns:nrzt:ns,0)+r1(ns:nrzt:ns,1))
      rmin_surf = MINVAL(r1(ns:nrzt:ns,0)+r1(ns:nrzt:ns,1))
      zmax_surf = MAXVAL(ABS(z1(ns:nrzt:ns,0)+z1(ns:nrzt:ns,1)))

      DO js = 2, ns
         modb(:nznt) = SQRT(two*(bsq(js:nrzt:ns)-pres(js)))
         CALL bextrema (modb, bmin(1,js), bmax(1,js), nzeta, ntheta2)
      END DO

!
!     output geometrical, |B| quantities
!
      CALL elongation (r1, z1, waist, height)

      IF(rank.EQ.0) WRITE (nthreed, 75) bmin(1,ns), bmax(1,ns),
     1    bmin(ntheta2,ns), bmax(ntheta2,ns)
   75 FORMAT(/
     1   ' Magnetic field modulation (averaged over toroidal angle)',/,
     2   1x,71('-')/,' BMIN(u=0)             = ',f14.6/
     3   ' BMAX(u=0)             = ',f14.6/' BMIN(u=pi)            = ',
     4   f14.6/' BMAX(u=pi)            = ',f14.6/)

#ifdef _ANIMEC
      sumbtot = 2*vnorm*SUM(bsq(2:nrzt)*tau(2:nrzt))
#else
      sumbtot = 2*(vnorm*SUM(bsq(2:nrzt)*tau(2:nrzt)) - sump)
#endif
      sumbtor = vnorm*SUM(tau(2:nrzt)*(r12(2:nrzt)*bsupv(2:nrzt))**2) 
      sumbpol = sumbtot - sumbtor
      betapol = 2*sump/sumbpol
      sump20 = 2*sump
#ifdef _FLOW
      sump2  = vnorm*SUM(tau(2:nrzt)*prot(2:nrzt)*prot(2:nrzt))
      sumrot = vnorm*SUM(tau(2:nrzt)*protrsq(2:nrzt))
      sumrt2 = vnorm*SUM(tau(2:nrzt)*protrsq(2:nrzt)*protrsq(2:nrzt))
#else
      sump2 = SUM(pres(2:ns)*pres(2:ns)*vp(2:ns)*vnorm)
#endif
      betatot = sump20/sumbtot
      betator = sump20/sumbtor
#ifdef _ANIMEC
      sumpar = vnorm*SUM(tau(2:nrzt)*ppar(2:nrzt))
      sumper = vnorm*SUM(tau(2:nrzt)*pperp(2:nrzt))
      sumpp2 = vnorm*SUM(tau(2:nrzt)*pperp(2:nrzt)*pperp(2:nrzt))
      sumpb2 = vnorm*SUM(tau(2:nrzt)*ppar(2:nrzt)*ppar(2:nrzt))
#endif
      VolAvgB = SQRT(ABS(sumbtot/volume_p))
      IonLarmor = 0.0032_dp/VolAvgB
      jPS2(2:ns1) = (jpar2(2:ns1) - jdotb(2:ns1)**2/bdotb(2:ns1))
      jpar_perp = SUM(jpar2(2:ns1)*(vp(2:ns1) + vp(3:ns)))
      jparPS_perp = SUM(jPS2(2:ns1)*(vp(2:ns1) + vp(3:ns)))
      s2 = SUM(jperp2(2:ns1)*(vp(2:ns1) + vp(3:ns)))
      IF (s2 .ne. zero) THEN
         jpar_perp = jpar_perp/s2
         jparPS_perp = jparPS_perp/s2
      END IF
      IF (ntor .gt. 1) THEN
        IF(rank.EQ.0) WRITE (nthreed, 80) aspect, kappa_p, volume_p,
     1   cross_area_p,  surf_area_p, circum_p, Rmajor_p, Aminor_p,
     2   rmin_surf,  rmax_surf, zmax_surf, waist(1), height(1),
     3   waist(2), height(2)
      ELSE
        IF(rank.EQ.0) WRITE (nthreed, 80) aspect, kappa_p, volume_p,
     1   cross_area_p,  surf_area_p, circum_p, Rmajor_p, Aminor_p,
     2   rmin_surf,  rmax_surf, zmax_surf, waist(1), height(1)
      END IF
 80   FORMAT(/,' Geometric and Magnetic Quantities',/,1x,71('-')/,
     1   ' Aspect Ratio          = ',f14.6, /
     1   ' Mean Elongation       = ',f14.6, /
     1   ' Plasma Volume         = ',f14.6,' [M**3]',/
     2   ' Cross Sectional Area  = ',f14.6,' [M**2]',/
     3   ' Normal Surface Area   = ',f14.6,' [M**2]',/
     4   ' Poloidal Circumference= ',f14.6,' [M]',/
     5   ' Major Radius          = ',f14.6,' [M]',
     6   ' (from Volume and Cross Section)',/
     7   ' Minor Radius          = ',f14.6,' [M]',
     8   ' (from Cross Section)',/
     9   ' Minimum (inboard)  R  = ',f14.6,' [M]',/
     A   ' Maximum (outboard) R  = ',f14.6,' [M]',/
     A   ' Maximum height     Z  = ',f14.6,' [M]',/
     B   ' Waist (v = 0)   in R  = ',f14.6,' [M]',/
     B   ' Full Height(v = 0)    = ',f14.6,' [M]',:,/
     B   ' Waist (v = pi)  in R  = ',f14.6,' [M]',:,/
     B   ' Full Height(v = pi)   = ',f14.6,' [M]')
      IF (rank.EQ.0) WRITE (nthreed, 85) toroidal_flux,
     1       1.e-6_dp*ctor/mu0, rbtor, rbtor0, VolAvgB, IonLarmor,
     2       jpar_perp, jparPS_perp
 85   FORMAT(
     1   ' Toroidal Flux         = ',f14.6,' [Wb]',/
     1   ' Toroidal Current      = ',f14.6,' [MA]',/
     1   ' RBtor(s=1)            = ',f14.6,' [T-m]',/
     1   ' RBtor(s=0)            = ',f14.6,' [T-m]',/
     1   ' Volume Average B      = ',f14.6,' [T]',/
     2   ' Ion Larmor Radius     = ',f14.6,' [M] X Ti(keV)**0.5',/
     3   ' <J||**2>/<J-perp**2>  = ',f14.6,' (Vol. Averaged)',/
     4   ' <JPS**2>/<J-perp**2>  = ',f14.6,' (Vol. Averaged)',/)

      IF(rank.EQ.0) WRITE (nthreed, 90)
   90 FORMAT(/,71('-'),/,' MORE GEOMETRIC AND PHYSICS QUANTITIES',/,
     1    71('-'),/,' Toroidal Plane: Phi = 0',/,
     1    5x,'j',3x,'psi-psiaxis',9x,'a [M]',3x,'ellipticity',3x,
     2   'indentation',7x,'d-shape',4x,'rel. shift',6x,'<J||**2>/',4x,
     3   '<JPS**2>/',/,95x,
     4   '<J-perp**2>',3x,'<J-perp**2>'/,' -----',8(2x,12('-')))

      fac = twopi*hs*signgs
      psi(1) = zero
      ALLOCATE (r3v(ns-1))
      r3v(:ns-1) = fac*phip(2:ns)*iotas(2:ns)
      DO i = 1, ns - 1
         psi(1+i) = psi(i) + r3v(i)
      END DO
      DEALLOCATE (r3v)

!     nphi-plane, noff = 1,....,nzeta
      DO nplanes = 1, 2
         IF (nplanes .eq. 1) THEN
            noff = 1
         ELSE
            IF (nzeta .eq. 1) EXIT
            IF(rank.EQ.0) WRITE (nthreed, 95)
            noff = 1 + nzeta/2
         END IF
          
         ygeo(1) = zero
         DO js = 2, ns
            zmin =  HUGE(zmin)
            zmax = -HUGE(zmax)
            xmin =  HUGE(xmin)
            xmax = -HUGE(xmax)
            rzmax = zero

!           Theta = 0 to pi in upper half of X-Z plane
            DO icount = 1,2
               n1 = noff                                      !nphi-plane, n1 = noff,...,nzeta
               IF (icount .eq. 2)
     1         n1 = MOD(nzeta + 1 - noff,nzeta) + 1           !(twopi-v), reflected plane
               loff = js + ns*(n1-1)
               t1 = one
               IF (icount .eq. 2) t1 = -one
               DO itheta = 1,ntheta2
                  yr1u = r1(loff,0) + sqrts(js)*r1(loff,1)
                  yz1u = z1(loff,0) + sqrts(js)*z1(loff,1)
                  yz1u = t1*yz1u
                  IF (yz1u .ge. zmax) THEN
                     zmax = ABS(yz1u)
                     rzmax = yr1u
                  ELSE IF (yz1u .le. zmin) THEN
                     zmin = yz1u
                     rzmin = yr1u
                  END IF
                  IF (yr1u .ge. xmax) THEN
                     xmax = yr1u
                     zxmax = yz1u
                  ELSE IF (yr1u .le. xmin) THEN
                     xmin = yr1u
                     zxmin = yz1u
                  END IF
                  loff = loff + ns*nzeta
               END DO
            END DO


            lpi = ns*((noff-1) + nzeta*(ntheta2 - 1))
            lt  = ns*((noff-1))
            xmida = r1(js+lpi,0) + sqrts(js)*r1(js+lpi,1)
            xmidb = r1(js+lt,0)  + sqrts(js)*r1(js+lt,1)

            rgeo = p5*(xmidb + xmida)              !Geometric major radius
            ygeo(js) = p5*(xmidb - xmida)          !Geometric minor radius
   
            yinden(js) = (xmida - xmin)/(xmax - xmin) !Geometric indentation
            yellip(js) = (zmax - zmin)/(xmax - xmin)  !Geometric ellipticity
   
            ytrian(js) = (rgeo - rzmax)/(xmax - xmin) !Geometric triangularity
            yshift(js) = (r1(1+lt,0)-rgeo)/(xmax - xmin) !Geometric shift

            IF (jperp2(js) .eq. zero) jperp2(js) = EPSILON(jperp2(js))
            jpar_perp = jpar2(js)/jperp2(js)
            IF (js .lt. ns) THEN
               jparPS_perp = jPS2(js)/jperp2(js)
            ELSE
               jparPS_perp = zero
            END IF

            IF (nplanes .eq. 1) THEN
               IF(rank.EQ.0) WRITE (nthreed, 120) js, psi(js), ygeo(js),
     1        yellip(js), yinden(js), ytrian(js), yshift(js), jpar_perp, 
     2        jparPS_perp
            ELSE
               IF(rank.EQ.0) WRITE (nthreed, 120) js, psi(js), ygeo(js),
     1             yellip(js), yinden(js), ytrian(js), yshift(js)
            END IF

         END DO
      END DO

   95 FORMAT(/,71('-'),/,' Toroidal Plane: Phi = 180/Nfp',/,71('-'),/)
  120 FORMAT(1x,i5,6f14.5,1p,3e14.2)

      IF(rank.EQ.0) WRITE (nthreed, 130)
  130 FORMAT(//,' Magnetic Fields and Pressure',/,1x,71('-'))
      fac = p5/mu0
      IF(rank.EQ.0) WRITE (nthreed, 140) sump/mu0, pavg/mu0,
     1   fac*sumbpol,  fac*sumbpol/volume_p, fac*sumbtor, fac*sumbtor/
     2   volume_p, fac*sumbtot, fac*sumbtot/volume_p, c1p5*sump/mu0, 
     3   c1p5*pavg/mu0
  140 FORMAT(' Volume Integrals (Joules) and Volume ',
     1   'Averages (Pascals)',/,24x,'Integral',6x,'Average',/,
     2   ' pressure         = ',1p,2e14.6,/,' bpol**2 /(2 mu0) = ',
     3   2e14.6,/,' btor**2/(2 mu0)  = ',2e14.6,/,
     4   ' b**2/(2 mu0)     = ',2e14.6,/,' EKIN (3/2p)      = ',
     5   2e14.6,/)

      IF(rank.EQ.0) WRITE (nthreed, 800)
  800 FORMAT(/,' MAGNETIC AXIS COEFFICIENTS'/,
     1   '    n     rmag       zmag        rmag        zmag',/,
     2   '        (cos nv)   (sin nv)    (sin nv)    (cos nv)',/)
      loff = LBOUND(rmags,1)
      IF (rank. EQ. 0) THEN
      DO n = 0, ntor
         n1 = n + loff
         t1 = mscale(0)*nscale(n)
         tz = t1
         IF (.NOT.lthreed) tz = 0
         IF (lasym) THEN
             WRITE (nthreed, 820) n, t1*rmags(n1),
     1          tz*zmags(n1), tz*rmaga(n1), t1*zmaga(n1)
         ELSE
             WRITE (nthreed, 820) n, t1*rmags(n1),
     1          tz*zmags(n1)
         END IF
      END DO
      END IF
  820 FORMAT(i5,1p,4e12.4)

#ifdef _ANIMEC
      bmax_t = -HUGE(bmax_t)
      bmin_t =  HUGE(bmin_t)
      DO js=2,ns
        DO lk=1,nznt
            l=js +(lk-1)*ns
            bmax_t = MAX(bmax_t,sqrt(two*bsq(l)))
            bmin_t = MIN(bmin_t,sqrt(two*bsq(l)))
        END DO
      END DO
      betpar = sumpar / sumbtot
      betper = sumper / sumbtot
      bet_s_to = two*sqrt((two*sumpp2+sumpb2)/(3*volume_p)) /
     1           (sumbtot/volume_p)
      bet_s_ew = two*sqrt((sumpp2+sumpb2)/(two*volume_p)) /
     1           (sumbtot/volume_p)
      betplp = sumper / sumbpol
      betplb = sumpar / sumbpol
      WRITE (nthreed, 150) betatot, betapol, betator
     1 ,two*(two*betper+betpar)/3, two*betper, betpar+betper,
     2  bet_s_to,bet_s_ew, two*(two*betplp+betplb)/3, betplp+betplb

  150 FORMAT(/,' From volume averages over plasma, betas are',/,
     1   ' beta thermal total    = ',f14.6,/,
     2   ' beta thermal poloidal = ',f14.6,/,
     3   ' beta thermal toroidal = ',f14.6,/,
     4   ' beta total            = ',f14.6,/,
     5   ' beta diamagnetic      = ',f14.6,/,
     6   ' beta equal weight     = ',f14.6,/,
     7   ' beta-star total       = ',f14.6,/,
     8   ' beta-star equal weight= ',f14.6,/,
     9   ' beta poloidal total   = ',f14.6,/,
     A   ' beta poloidal (EW)    = ',f14.6,/   )
#else
      WRITE (nthreed, 152) betatot, betapol, betator
 152  FORMAT(/,' From volume averages over plasma, betas are',/,
     B   ' beta total    = ',f14.6,/,' beta poloidal = ',f14.6,/,
     C   ' beta toroidal = ',f14.6,/)
#endif

      betstr = two*SQRT(sump2/volume_p)/(sumbtot/volume_p)
#ifdef _FLOW
      betarot = two*sumrot / sumbtot
      betrstr = two*SQRT(sumrt2/volume_p)/(sumbtot/volume_p)
#endif

      IF(rank.EQ.0) WRITE (nthreed, 160) rbtor, betaxis, betstr
#ifdef _ANIMEC
     1      ,bmax_t, bmin_t, bmax_t/bmin_t,bcrit
#elif defined _FLOW
     2      ,betarot, betarot+betatot, betrstr, sqrt(machsq)
#endif
  160 FORMAT(' R * Btor-vac         = ',f14.6,' [Wb/M]',/,
#ifdef _ANIMEC
     1       ' Peak Thermal Beta    = ',f14.6,/,
     2       ' Beta-star (thermal)  = ',f14.6,/,
     3       ' B_max                = ',f14.6,' [T]',/,
     4       ' B_min                = ',f14.6,' [T]',/,
     5       ' B_max / B_min        = ',f14.6,/,
     6       ' B_crit               = ',f14.6,' [T]',/
#elif defined _FLOW
     7       ' Peak Beta            = ',f14.6,/,
     8       ' Beta-star            = ',f14.6,/
     9       ' Beta-rotation        = ',f14.6,/,
     A       ' Beta-pressure+flow   = ',f14.6,/,
     B       ' Beta-star-rotation   = ',f14.6,/,
     C       ' Mach number on axis  = ',f14.6,/
#else
     D       ' Peak Beta            = ',f14.6,/,
     E       ' Beta-star            = ',f14.6,/
#endif
     F        )


!
!
!     Shafranov surface integrals s1,s2
!     Plasma Physics vol 13, pp 757-762 (1971)
!     Also, s3 = .5*S3, defined in Lao, Nucl. Fusion 25, p.1421 (1985)
!     Note: if ctor = 0, use Int(Bsupu*Bsubu dV) for ctor*ctor/R
!     Phys. Fluids B, Vol 5 (1993) p 3121, Eq. 9a-9d
!

      ALLOCATE (rbps1u(nznt), bpol2vac(nznt))
      IF (lfreeb .and. ivac.gt.1) THEN
         bpol2vac = 2*bsqvac - bphiv*bphiv
      ELSE
         bpol2vac = 2*(c1p5*bsq(ns:nrzt:ns)   - p5*bsq(ns-1:nrzt:ns))
     1            -  ((c1p5*bsupv(ns:nrzt:ns) - p5*bsupv(ns-1:nrzt:ns))
     2            * redge)**2
      END IF

!     Compute current-like norm (factor) in Eq.(8), <a> * int(Bpol**2 * dA)
!     where <a> == 2*pi*Rs in Eq. 8 is the effective minor radius = Vol/Asurf
!     (corrects wrong description of Rs in paper, which is NOT the major radius)
!     This aminr1 = 1/2 the "correct" aminr1
      aminr1 = volume_p/surf_area_p
      factor = twopi**2*aminr1*SUM(bpol2vac*surf_area)
      factor = one/factor
      facnorm = factor*twopi**2

!     Lao's definition of normalization factor
      scaling_ratio = (mu0*curtor/circum_p)**2*volume_p
      scaling_ratio = scaling_ratio*factor

      rbps1u(:nznt) = facnorm*redge(:nznt)*phat(:nznt)
     1                *wint(ns:nznt*ns:ns)
      sigr0 = SUM(rbps1u(:nznt)*zu0(ns:nrzt:ns))
      sigr1 = SUM(rbps1u(:nznt)*zu0(ns:nrzt:ns)*redge(:nznt))
      sigz1 =-SUM(rbps1u(:nznt)*ru0(ns:nrzt:ns)*
     1                  (z1(ns:nrzt:ns,0) + z1(ns:nrzt:ns,1)))
      DEALLOCATE (redge, phat, rbps1u, bpol2vac, surf_area)

      er = sigr1 + sigz1
      rlao = volume_p/(twopi*cross_area_p)       !LAO, NUCL.FUS.25(1985)1421
      flao = rshaf/rlao
      fgeo = rshaf/rcen

      smaleli = factor*sumbpol
      vvc_smaleli = smaleli ! Save result for v3fit.
      betai   = 2*factor*sump
      musubi  = vnorm*factor*musubi
      dmusubi_meas = 2*twopi*factor*delphid*rbtor
      lambda = p5*smaleli + betai
      s11 = er - rshaf*sigr0              !Shafranov def. based on RT, Eq.(12)
      s12 = er - rcen*sigr0               !R = Rgeometric
      s13 = er - rlao*sigr0               !R = RLao
      s2  = sigr0*rshaf
      s3  = sigz1                         !1/2 S3 in Eq.(14c)
      delta1 = zero
      delta2 = one - fgeo
      delta3 = one - flao
      IF(rank.EQ.0) WRITE (nthreed, 168)
      IF(rank.EQ.0) WRITE (nthreed, 170) rshaf, rcen, rlao,
     1   scaling_ratio, s3, smaleli, musubi, betai, lambda
      IF (lrecon.AND.rank.EQ.0) WRITE (nthreed, 172) dmusubi_meas
      IF(rank.EQ.0) WRITE (nthreed, 174) delta1, delta2, delta3, 
     1   s11, s12, s13, s2, s2/fgeo, s2/flao, 
     5   musubi + s11,musubi + s12, 
     6   musubi + s13, 
     8   p5*s11 + s2, p5*s12 + s2/fgeo, p5*s13 + s2/flao, 
     A   p5*(3*betai+smaleli-musubi)/(s11+s2) - one, 
     B   p5*(3*betai+smaleli-musubi)/(s12+s2/fgeo) - one, 
     C   p5*(3*betai+smaleli-musubi)/(s13+s2/flao) - one, 
     D   p5*(betai+smaleli+musubi)/s2 - one, 
     E   p5*fgeo*(betai+smaleli+musubi)/s2 - one,
     F   p5*flao*(betai+smaleli+musubi)/s2 - one

  168 FORMAT(' Shafranov Surface Integrals',/
     1       ' Ref: S. P. Hirshman, Phys. Fluids B, 5, (1993) 3119',/,
     2       ' Note: s1 = S1/2, s2 = S2/2, where ',
     3       ' s1,s2 are the Shafranov definitions,',/,
     4       ' and s3 = S3/2, where S3 is Lao''s definition.',/,
     5       ' The quantity lsubi gives the ratio of volume poloidal',
     6     /,' field energy to the field energy estimated from the',
     7     /,' surface integral in Eq.8.',/,1x,22('-'),/)

  170 FORMAT(
     5   ' RT (Pressure-weighted)  = ',f14.6,' [M]',/,
     6   ' RG (Geometric)          = ',f14.6,' [M]',/,
     7   ' RL (Vol/2*pi*Area-Lao)  = ',f14.6,' [M]',/,
     8   ' Poloidal Field Energy',/,
     8   ' Normalization Ratio     = ',f14.6,' (Lao/Hirshman)',//,
     8   ' s3                      = ',f14.6,/,
     9   ' lsubi                   = ',f14.6,/,
     6   ' musubi                  = ',f14.6,/,
     7   ' betai                   = ',f14.6,/,
     9   ' lambda                  = ',f14.6,/)
  172 FORMAT(' musubi (diamagnetism)   = ',f14.6)
  174 FORMAT(/,32x,'R = RT',12x,'R = RG',12x,'R = RL',/, 
     1       20x,3(10x,8('-')),/,
     1   ' delta = 1 - RT/R     = ',3(f14.6,4x),/,
     2   ' s1                   = ',3(f14.6,4x),/,
     3   ' s2                   = ',3(f14.6,4x),/,
     8   ' betai (Mui + s1)     = ',3(f14.6,4x),/,
     A   ' lambda (s1/2 + s2)   = ',3(f14.6,4x),/,
     B   ' 1st Shafr''v relation = ',3(f14.6,4x),/,
     C   ' (3*Betai + Li - Mui)/[2*(s1+s2)] - 1',/,
     D   ' Radial force balance = ',3(f14.6,4x),/,
     E   ' (Betai + Li + Mui)/(2*s2) - 1',/)


      IF (lrecon) THEN
!
!     Safety Factors (q)
!
      qzero = one/ABS(iotaf(1))
      qedge = one/ABS(iotaf(ns))
      IF(rank.EQ.0) WRITE (nthreed, 180) qzero, qedge, qzero/qedge
  180 FORMAT(' Safety Factors',/,1x,14('-'),/,' q (on axis) = ',f12.4,/,
     1   ' q (at edge) = ',f12.4,/,' q(0)/qedge  = ',f12.4,/)

!
!     PRINT OUT IOTA, PRESSURE SPLINE COEFFICIENTS
!     (PRESSURE IN MKS UNITS, NWT/M**2)
!
      IF(rank.EQ.0) WRITE (nthreed, 190)
  190 FORMAT(/,' SUMMARY OF IOTA AND PRESSURE SPLINES'/,
     1   '  K   Spline Node       IOTA(K)      IOTA"(K)'/,
     2   '        SQRT(s)',/,3('-'),3(4x,10('-')))
      DO i = 1, isnodes
         IF(rank.EQ.0) WRITE (nthreed, 200) i, sknots(i), ystark(i),
     1                  y2stark(i)
      END DO
  200 FORMAT(i3,1p,3e14.3)
      IF(rank.EQ.0) WRITE (nthreed, 210)
  210 FORMAT(/,'  K   Spline Node       PRES(K)      PRES"(K)'/,
     1   '        SQRT(s)',/,3('-'),3(4x,10('-')))
      factor = pthommax
      DO i = 1, ipnodes
         IF(rank.EQ.0) WRITE (nthreed, 200) i, pknots(i),
     1      factor*ythom(i), factor* y2thom(i)
      END DO

!
!     PRINT-OUT MAGNETIC PITCH DATA
!
      nchimse = imse
      total_mse_chi = zero

      IF (imse .ne. 0) THEN

         IF(rank.EQ.0) WRITE (nthreed, 220)
  220    FORMAT(//,4x,'N     R(data)      Spline',3x,
     1      'atan[BZ/BT] ATAN[BZ/BT] Chi-Sq-Err',4x,
     2      'Sigma       BZ(in)     BZ(calc)       BT        q(in)',/,
     3      12x,'[m]',5x,'SQRT(s)-node',2x,'(in-deg)',2x,'(calc-deg)',
     4      20x,3(9x,'[T]')/,2x,3('-'),1x,2(3x,9('-')),8(2x,10('-'))/)

         msewgt = zero
         denwgt = zero
         DO n = 1, imse + 1
            nsort = isortr(n)
            pitchc = ATAN(starkcal(n))/dcon
            pitchm = ATAN(datastark(nsort))/dcon
            wght0 = ATAN(sigma_stark(nsort))/dcon
            js = indexs1(n)
            lt = indexu1(n)
            noff = ns*(lt - 1)
            rjs = r1(js+noff,0) + sqrts(js)*r1(js+noff,1)
            js1 = js + 1
            rjs1 = r1(js1+noff,0) + sqrts(js1)*r1(js1+noff,1)
            rcalc = (one - delso1(n))*rjs + delso1(n)*rjs1
            IF (delse1(n) .lt. zero) rcalc = zero
            IF (rcalc .ne. 0.) btor = fpsical(n)/rcalc
            bzin = TAN(dcon*pitchm)*btor
            bzcalc = TAN(dcon*pitchc)*btor
            chisq = zero
            IF (ABS(rcalc - router) .lt. epstan) THEN
               IF(rank.EQ.0) WRITE (nthreed, 230) n, rcalc, rsort0(n),
     1            pitchm, pitchc, wght0, one/(qcalc(n) + epstan)
            ELSE
               IF (ABS(rsort(n) - rcalc) .le. epstan) THEN
                  chisq = ((pitchc - pitchm)/wght0)**2
                  msewgt = msewgt + chisq
                  denwgt = denwgt + (pitchm/wght0)**2
               ENDIF
               IF(rank.EQ.0) WRITE (nthreed, 240) n, rcalc, rsort0(n),
     1            pitchm, pitchc, chisq, wght0, bzin, bzcalc, btor,
     2            one/(qmeas(n) + epstan)
            ENDIF
         END DO

         chiwgt = msewgt/(imse)
         msewgt = SQRT(msewgt/denwgt)
c                                             !total chi-squared for mse
         total_mse_chi = (imse)*chiwgt
         IF(rank.EQ.0) WRITE (nthreed, 250) chiwgt, msewgt

  230    FORMAT(i5,' ',1p,2e12.3,2e12.3,12x,e12.3,6x,
     1      '- Outer (phantom) Edge -',6x,e12.3)
  240    FORMAT(i5,' ',1p,2e12.3,8e12.3)
  250    FORMAT(/' MEAN CHI-SQ ERROR IN STARK DATA MATCH : ',1p,e10.3,/,
     1      ' RMS ERROR IN STARK DATA MATCH : ',e10.3,2/)

      ENDIF

!
!     PRINT-OUT PRESSURE DATA
!
      nchipres = itse
      total_pres_chi = zero

      IF (lpprof) THEN
         tswgt = zero
         denwgt = zero
         IF(rank.EQ.0) WRITE (nthreed, 300) presfac*pfac
  300    FORMAT(4x,'N     R(in)      R(calc)',4x,f6.2,
     1      ' X Pres(in)      Pres(calc)','  Chi-Sq Err       Sigma'/,2x
     2      ,3('-'),2(2x,10('-')),2(6x,12('-')),2(3x,9('-')),/)

         DO n = 1, itse
            js = indexs2(n)
            lt = indexu2(n)
            noff = ns*(lt - 1)
            rjs = r1(js+noff,0) + sqrts(js)*r1(js+noff,1)
            js1 = js + 1
            rjs1 = r1(js1+noff,0) + sqrts(js1)*r1(js1+noff,1)
            IF (delso2(n) .eq. (-one)) delso2(n) = one
            rcalc = (one - delso2(n))*rjs + delso2(n)*rjs1
            wght0 = sigma_thom(n)
            chisq = ((datathom(n)*pfac-pcalc(n)/mu0)/wght0)**2
            tswgt = tswgt + chisq
            denwgt = denwgt + (datathom(n)/wght0)**2
            IF(rank.EQ.0) WRITE (nthreed, 310) n, rthom(n), rcalc,
     1         presfac*pfac* datathom(n), pcalc(n)/mu0, chisq, wght0
         END DO

         chiwgt = tswgt/(itse)
         tswgt = SQRT(tswgt/denwgt)
         total_pres_chi = (itse)*chiwgt !total
         IF(rank.EQ.0) WRITE (nthreed, 320) chiwgt, tswgt
  310    FORMAT(i5,1p,2e12.3,2e18.3,2e12.3)
  320    FORMAT(/' MEAN CHI-SQ ERROR IN PRESSURE DATA MATCH: ',
     1      1p,e10.3,/,' RMS ERROR IN PRESSURE DATA MATCH: ',e10.3/)
      ENDIF

!
!     SUMMARIZE MAGNETICS DATA AND MATCH
!
      CALL magnetics_data

!
!     COMPUTE REAL TOROIDAL CURRENT ALONG MIDPLANE (MULTIPLY IZETA BY R/SQRT(G))
!     IN PHI = 0 PLANE
!
      CALL getcurmid (curmid, izeta, gsqrt, r12)

      DO nout = nthreed, nmac, (nmac - nthreed)
         IF (nout .eq. nmac) THEN
            IF (.not.lmac) CYCLE
            WRITE (nout, *)
     1      'FOLLOWING DATA EQUALLY SPACED IN TOROIDAL FLUX'
         END IF
         WRITE (nout, 700)
         IF (nout .eq. nthreed) WRITE (nout, 710)
         iotas(1) = two*iotas(2) - iotas(3)
         iotas(ns+1) = two*iotas(ns) - iotas(ns1)
         vp(1) = two*vp(2) - vp(3)
         vp(ns+1) = two*vp(ns) - vp(ns1)
         pres(1) = two*pres(2) - pres(3)
         pres(ns+1) = two*pres(ns) - pres(ns1)
         DO icount = 1, 2*ns - 1
            js = MOD(imid(icount) - 1,ns) + 1
            ageo(icount) = ygeo(js)
            phimid(icount) = toroidal_flux*(js - 1)/(ns1)
            psimid(icount) = psi(js)
            volpsi(icount) = vnorm*SUM(vp(2:js))
            IF (js .eq. 1) volpsi(icount) = zero
            dlogidsi = (iotas(js+1)-iotas(js))*ohs/iotaf(js)
            vprime = p5*vnorm*ohs*(vp(js)+vp(js+1))
            pprime = (pres(js+1)-pres(js))*ohs
            rgeo = rmid(icount) + ygeo(js)
            IF (icount .gt. ns) rgeo = rgeo - two*ygeo(js)
            alfa(icount) = -two*pprime*vprime*SQRT(two*volpsi(icount)/
     1         rgeo/twopi)/(iotaf(js)*phipf_loc(js))**2
            shear(icount) = -two*volpsi(icount)*dlogidsi/vprime
            WRITE (nout, 720) rmid(icount), ygeo(js), psi(js),
     1         volpsi(icount), qmid(icount), shear(icount),
     2         presmid(icount),
     3         alfa(icount), curmid(icount), datamse(icount)
         END DO
      END DO
  700 FORMAT(/3x,'RMID[M]',6x,'AMID[M]',6x,'PSI[Wb]',4x,'VOL[M**3]',9x,
     1   'q(r)',2x,'SHEAR(Spsi)',6x,'P(PASC)',6x,'ALF(P'')',2x,
     2   'JTOR(A/M**2)',1x,'ATAN(Bz/Bt)[deg]')
  710 FORMAT(10('-'),9(3x,10('-')))
  720 FORMAT(1p,e10.3,9(3x,e10.3))

!     Determine q-MIN on the s grid by MIN-splining on the SQRT-s grid
      nptsx = 2*ns - 1

      wmidx(:nptsx) = one
      tenmidx(:nptsx) = 0.1_dp
      rmidx(:nptsx) = rmid(:nptsx)
      qmidx(:nptsx) = qmid(:nptsx)

      tol = .005_dp
      sqmin = fmax(zero,one,splints,tol)
      iotamax = splints(sqmin)
      qmin3 = -99999.0_dp
      IF (iotamax .ne. zero) qmin3 = one/iotamax
      sqmin = sqmin**2

!     Determine q-MIN on a fine r-midplane grid
      tol = .005_dp
!     outboard side
      rqmin1 = fmin(rmid(ns),rmid(nptsx),splintx,tol)
      qmin1 = splintx(rqmin1)
      rqmin2 = fmin(rmid(1),rmid(ns),splintx,tol)!outboard side only
      qmin2 = splintx(rqmin2)
      IF(rank.EQ.0) WRITE (nthreed, 730) qmin1, rqmin1, qmin2, rqmin2,
     1              qmin3, sqmin
  730 FORMAT(//' MINIMUM Q :     OUTBOARD SIDE: QMIN= ',f6.3,'  AT R= ',
     1   f6.3,/,'                  INBOARD SIDE: QMIN= ',f6.3,'  AT R= '
     2   ,f6.3,/,'                    IN S SPACE: QMIN= ',f6.3,
     3   '  AT S= ',f6.3,/)

      nchiiota0 = 0
      total_chi_square = total_b_chi + total_saddle_chi + total_pres_chi
     1    + total_mse_chi + total_chi_cur + total_chi_dia + chisqerr(
     2   islope0)                        !Total CHISQ of equilibrium fit
      nchitot = nbfldn + nchiiota0 + nchisaddle + nchipres + nchimse +
     1   nchicur + nchidia
      total_chi_square_n = total_chi_square/MAX(1,nchitot)
      IF(rank.EQ.0) WRITE (nthreed, 900)
      IF(rank.EQ.0) WRITE (nthreed, 901) (bloopnames(n),nbfld(n),
     1  b_chi(n),b_chi(n)/MAX(1,nbfld(n)),n=1,nbsets)
      IF(rank.EQ.0) WRITE (nthreed, 902) nchisaddle, total_saddle_chi,
     1   total_saddle_chi/MAX(1,nchisaddle), nchipres, total_pres_chi,
     2   total_pres_chi/MAX(1,nchipres), nchimse, total_mse_chi,
     3   total_mse_chi/MAX(1,nchimse), nchicur, total_chi_cur,
     4   total_chi_cur/MAX(1,nchicur), nchidia, total_chi_dia,
     5   total_chi_dia, nchitot, total_chi_square, total_chi_square_n

  900 FORMAT(/,' CHI-SQUARED SUMMARY:',t25,'Type',t50,'Nchi',t65,'ChiSq'
     1   ,t84,'ChiSq/N'/,t25,'----',t50,'----',t65,'-----',t84,'-------'
     2   )
  901 FORMAT(t25,'B-Loops-',a,t50,i3,t60,1p,e12.4,t80,e12.4)
  902 FORMAT(t25,'Saddle',t50,i3,t60,1p,e12.4,t80,e12.4,/,t25,
     1   'Pressure',t50,i3,t60,e12.4,t80,e12.4,/,t25,'MSE',t50,i3,
     2    t60,e12.4,t80,e12.4,/,t25,'Ip',t50,i3,t60,e12.4,t80,
     3    e12.4,/,t25,'Diamagnetic Flux',t50,i3,t60,e12.4,t80,e12.4,/
     4   ,t25,'TOTAL',t50,i3,t60,e12.4,t80,e12.4)

      END IF          !!IF(LRECON)

      CALL second0 (teqfoff)
      timer(teqf) = timer(teqf) + teqfoff - teqfon
      fo_eqfor_time = timer(teqf)

      END SUBROUTINE eqfor
