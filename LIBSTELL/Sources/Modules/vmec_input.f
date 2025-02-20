      MODULE vmec_input
      USE vparams, ONLY: rprec, dp, mpol1d, ntord, ndatafmax
      USE vsvd0
      IMPLICIT NONE
!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!   For variable descriptions, see VMEC "readin.f" routine
!-----------------------------------------------
      INTEGER, PARAMETER :: mpol_default = 6
      INTEGER, PARAMETER :: ntor_default = 0
      INTEGER, PARAMETER :: ns_default   = 31
      INTEGER :: nfp, ncurr, nsin, niter, nstep, nvacskip, mpol, ntor,
     1           ntheta, nzeta, mfilter_fbdy, nfilter_fbdy,
     2           max_main_iterations, omp_num_threads
      INTEGER, DIMENSION(100) :: ns_array, niter_array
      INTEGER :: imse, isnodes, itse, ipnodes, iopt_raxis,
     1   imatch_phiedge, nflxs
      INTEGER, DIMENSION(nbsetsp) :: nbfld
      INTEGER, DIMENSION(nfloops) :: indxflx
      INTEGER, DIMENSION(nbcoilsp,nbsetsp) :: indxbfld
      REAL(rprec), DIMENSION(-ntord:ntord,0:mpol1d) ::
     1   rbs, zbc, rbc, zbs
      REAL(rprec) :: time_slice, curtor, delt, ftol, tcon0,
     1   gamma, phiedge, phidiam, sigma_current, sigma_delphid, tensi,
     2   tensp, tensi2, fpolyi, presfac, mseangle_offset, pres_offset,
     3   mseangle_offsetm, spres_ped, bloat, pres_scale, 
     4   prec2d_threshold
      REAL(rprec), DIMENSION(0:20) :: am, ai, ac
      REAL(rprec), DIMENSION(1:20) :: aphi
      CHARACTER(len=20) :: pcurr_type  !  len=12 -> len=20 J Hanson 2010-03-16
      CHARACTER(len=20) :: piota_type
      CHARACTER(len=20) :: pmass_type
      REAL(rprec), DIMENSION(ndatafmax) :: am_aux_s, am_aux_f,                 
     &   ai_aux_s, ai_aux_f, ac_aux_s, ac_aux_f

!     ANISOTROPIC AMPLITUDES: AH=PHOT/PTHERMAL, AT=TPERP/TPAR
!     bcrit: hot particle energy deposition value for |B|
      REAL(rprec), DIMENSION(0:20) :: ah, at       
      REAL(rprec)  :: bcrit
      CHARACTER(len=20) :: pt_type     ! SAL  For Ani/Flow
      CHARACTER(len=20) :: ph_type     ! SAL  For Ani/Flow,
      REAL(rprec), DIMENSION(ndatafmax) :: ah_aux_s, ah_aux_f,
     &                                     at_aux_s, at_aux_f

      REAL(rprec), DIMENSION(0:ntord) :: raxis, zaxis                !!Backwards compatibility: Obsolete
      REAL(rprec), DIMENSION(0:ntord) :: raxis_cc, raxis_cs,
     1                                   zaxis_cc, zaxis_cs
      REAL(rprec), DIMENSION(100) :: ftol_array
      REAL(rprec), DIMENSION(nigroup), TARGET :: extcur ! V3FIT needs a pointer to this.
      REAL(rprec), DIMENSION(nmse) :: mseprof
      REAL(rprec), DIMENSION(ntse) :: rthom, datathom, sigma_thom
      REAL(rprec), DIMENSION(nmse) :: rstark, datastark,
     1    sigma_stark
      REAL(rprec), DIMENSION(nfloops) :: dsiobt, sigma_flux
      REAL(rprec), DIMENSION(nbcoilsp,nbsetsp) :: bbc, sigma_b
      REAL(rprec), DIMENSION(ndatafmax) :: psa, pfa, isa, ifa
      LOGICAL :: lpofr, lmac, lfreeb, lrecon, loldout, ledge_dump,
     1           lasym, lforbal,lrfp, lmovie,  lmove_axis,
     2           lwouttxt, ldiagno,       ! J.Geiger: for txt- and diagno-output
     3           lmoreiter,               ! J.Geiger: if force residuals are not fulfilled add more iterations.
     4           lfull3d1out,             ! J.Geiger: to force full 3D1-output
     5           l_v3fit=.false.,
     6           lspectrum_dump, loptim           !!Obsolete
      LOGICAL :: lgiveup                  ! inserted M.Drevlak
      REAL(rprec) :: fgiveup              ! inserted M.Drevlak, giveup-factor for ftolv
      LOGICAL :: lbsubs                   ! J Hanson See jxbforce coding
      
      CHARACTER(len=200) :: mgrid_file
      CHARACTER(len=200) :: trip3d_file   ! SAL - TRIP3D
      CHARACTER(len=10)  :: precon_type
      CHARACTER(len=120) :: arg1
      CHARACTER(len=100) :: input_extension

      LOGICAL :: lnyquist = .TRUE.    !=false, suppress nyquist stuff; CZHU 2021.03.31

      NAMELIST /indata/ mgrid_file, time_slice, nfp, ncurr, nsin,
     1   niter, nstep, nvacskip, delt, ftol, gamma, am, ai, ac, aphi,
     1   pcurr_type, pmass_type, piota_type,
     1   am_aux_s, am_aux_f, ai_aux_s, ai_aux_f, ac_aux_s, ac_aux_f,  ! J Hanson 2010-03-16
     1   ah, at, bcrit,                                               ! WAC (anisotropic pres)
     1   ph_type, ah_aux_s, ah_aux_f,
     1   pt_type, at_aux_s, at_aux_f,
     2   rbc, zbs, rbs, zbc, spres_ped, pres_scale, raxis_cc, zaxis_cs, 
     3   raxis_cs, zaxis_cc, mpol, ntor, ntheta, nzeta, mfilter_fbdy,
     3   nfilter_fbdy, niter_array,
     4   ns_array, ftol_array, tcon0, precon_type, prec2d_threshold,
     4   curtor, sigma_current, extcur, omp_num_threads,
     5   phiedge, psa, pfa, isa, ifa, imatch_phiedge, iopt_raxis, 
     6   tensi, tensp, mseangle_offset, mseangle_offsetm, imse, 
     7   isnodes, rstark, datastark, sigma_stark, itse, ipnodes, 
     8   presfac, pres_offset, rthom, datathom, sigma_thom, phidiam, 
     9   sigma_delphid, tensi2, fpolyi, nflxs, indxflx, dsiobt, 
     A   sigma_flux, nbfld, indxbfld, bloat, raxis, zaxis,
     A   bbc, sigma_b, lpofr, lforbal, lfreeb, lmove_axis, lrecon, lmac, 
     C   lmovie,                                                           ! S Lazerson 2010
     D   lasym, ledge_dump, lspectrum_dump, loptim, lrfp,
     E   loldout, lwouttxt, ldiagno, lfull3d1out, max_main_iterations,     ! J Geiger 2010-05-04
     D   lgiveup,fgiveup,                                                  ! M.Drevlak 2012-05-10
     E   lbsubs,                                                           ! 2014-01-12 See jxbforce
     F   trip3d_file,                                                      ! SAL - TRIP3D
     G   lnyquist

      NAMELIST /mseprofile/ mseprof

      CONTAINS

      SUBROUTINE read_indata_namelist (iunit, istat)
      INTEGER, INTENT(IN) :: iunit
      INTEGER, INTENT(OUT) :: istat

!
!     INITIALIZATIONS
!
      omp_num_threads = 8
      gamma = 0
      spres_ped = 1
      mpol = mpol_default
      ntor = ntor_default
      ntheta = 0;  nzeta = 0
      ns_array = 0;  ns_array(1) = ns_default
      niter_array = -1;
      bloat = 1
      rbc = 0;  rbs = 0; zbs = 0; zbc = 0
      time_slice = 0
      nfp = 1
      ncurr = 0
      nsin = ns_default
      niter = 100
      nstep = 10
      nvacskip = 1
      delt = 1
      ftol = 1.E-10_dp
      ftol_array = 0;  ftol_array(1) = ftol
      am = 0; ai = 0; ac = 0; aphi = 0; aphi(1) = 1
      pres_scale = 1
      raxis_cc = 0; zaxis_cs = 0; raxis_cs = 0; zaxis_cc = 0;
      mfilter_fbdy = -1; nfilter_fbdy = -1
      tcon0 = 1
      precon_type = 'NONE'; prec2d_threshold = 1.E-30_dp
      curtor = 0; 
      extcur = 0;  phiedge = 1;
      mgrid_file = 'NONE'
      trip3d_file = 'NONE' ! SAL - TRIP3D
      lfreeb = .true.
      lmove_axis = .true.
      lmac = .false.
      lforbal = .false.
      lasym = .false.
      lrfp = .false.
      loldout = .false.        ! J Geiger 2010-05-04 start
      ldiagno = .false.
      lgiveup = .false.        ! inserted M.Drevlak
      fgiveup = 3.E+01_dp      ! inserted M.Drevlak
      lbsubs = .false.         ! J Hanson. See jxbforce coding
      lfull3d1out = .false.
      lmovie = .false.         ! S Lazerson for making movie files
      lmoreiter = .false.      ! default value if no max_main_iterations given.
      max_main_iterations = 1  ! to keep a presumably expected standard behavior.
#if defined(NETCDF)
      lwouttxt = .false.       ! to keep functionality as expected with netcdf
#else
      lwouttxt = .true.        ! and without netcdf
#endif

      pcurr_type = 'power_series'
      piota_type = 'power_series'
      pmass_type = 'power_series'

!     ANISTROPY PARAMETERS
      bcrit = 1
      at(0) = 1;  at(1:) = 0
      ah = 0
      ph_type = 'power_series'
      pt_type = 'power_series'
      ah_aux_s(:) = -1
      at_aux_s(:) = -1
      am_aux_s(:) = -1
      ac_aux_s(:) = -1
      ai_aux_s(:) = -1
      

!
!     BACKWARDS COMPATIBILITY
!
      raxis = 0;  zaxis = 0

      IF (iunit.eq.-327) RETURN
      
      READ (iunit, nml=indata, iostat=istat)

      IF (ALL(niter_array == -1)) niter_array = niter
      WHERE (raxis .ne. 0._dp) 
         raxis_cc = raxis
      ELSEWHERE
         raxis_cc = raxis_cc
      ENDWHERE
      WHERE (zaxis .ne. 0._dp) 
         zaxis_cs = zaxis
      ELSEWHERE
         zaxis_cs = zaxis_cs
      ENDWHERE

      raxis_cs(0) = 0; zaxis_cs(0) = 0

      IF(max_main_iterations .gt. 1) lmoreiter=.true.  !J Geiger: if more iterations are requested.

      END SUBROUTINE read_indata_namelist

      SUBROUTINE read_mse_namelist (iunit, istat)
      INTEGER :: iunit, istat

      READ (iunit, nml=mseprofile, iostat=istat)

      END SUBROUTINE read_mse_namelist
      
      SUBROUTINE write_indata_namelist (iunit, istat)
      IMPLICIT NONE
      INTEGER, INTENT(in) :: iunit
      INTEGER, INTENT(inout) :: istat
      INTEGER :: iftol,i,n,m
      INTEGER, DIMENSION(1) :: ins
      CHARACTER(LEN=*), PARAMETER :: outboo  = "(2X,A,1X,'=',1X,L1)"
      CHARACTER(LEN=*), PARAMETER :: outint  = "(2X,A,1X,'=',1X,I0)"
      CHARACTER(LEN=*), PARAMETER :: outint1 = "(2X,A,1X,'=',1X,I1.1)"
      CHARACTER(LEN=*), PARAMETER :: outint2 = "(2X,A,1X,'=',1X,I2.2)"
      CHARACTER(LEN=*), PARAMETER :: outint3 = "(2X,A,1X,'=',1X,I3.3)"
      CHARACTER(LEN=*), PARAMETER :: outint4 = "(2X,A,1X,'=',1X,I4.4)"
      CHARACTER(LEN=*), PARAMETER :: outint5 = "(2X,A,1X,'=',1X,I5.5)"
      CHARACTER(LEN=*), PARAMETER :: outint6 = "(2X,A,1X,'=',1X,I6.6)"
      CHARACTER(LEN=*), PARAMETER :: outflt="(2X,A,1X,'=',1X,ES22.12E3)"
      CHARACTER(LEN=*), PARAMETER :: outexp="(2X,A,1X,'=',1X,ES22.12E3)"
      IF (istat < 0) RETURN
      WRITE(iunit,'(A)') '!----- Runtime Parameters -----'
      WRITE(iunit,'(A)') '&INDATA'
      WRITE(iunit,outflt) 'DELT',delt
      WRITE(iunit,outint) 'NITER',niter
      WRITE(iunit,outint) 'NSTEP',nstep
      WRITE(iunit,outflt) 'TCON0',tcon0
      ins = MAXLOC(ns_array)
      WRITE(iunit,'(a,(1p,4i14))')  '  NS_ARRAY =    ',
     1     (ns_array(i), i=1,ins(1))
      iftol = 1
      DO WHILE(ftol_array(iftol).ne.0 .and. iftol.lt.100)
         iftol = iftol + 1
      END DO
      WRITE(iunit,'(a,(1p,4e14.6))')'  FTOL_ARRAY =  ',
     1     (ftol_array(i), i=1,iftol - 1)
      ins = MINLOC(niter_array)
      IF (ins(1) > 1)
     1 WRITE(iunit,'(a,(1p,4i14))') '  NITER_ARRAY = ',
     2      (niter_array(i), i=1,ins(1)-1)
      WRITE (iunit,'(2x,3a)') "PRECON_TYPE = '", TRIM(precon_type),"'" 
      WRITE (iunit,'(2x,a,1p,e14.6)') "PREC2D_THRESHOLD = ", 
     1                                prec2d_threshold
      WRITE(iunit,'(A)') '!----- Grid Parameters -----'
      WRITE(iunit,outboo) 'LASYM',lasym
      WRITE(iunit,outint4) 'NFP',nfp
      WRITE(iunit,outint4) 'MPOL',mpol
      WRITE(iunit,outint4) 'NTOR',ntor
      WRITE(iunit,outflt) 'PHIEDGE',phiedge
      WRITE(iunit,outint4) 'NTHETA',ntheta
      WRITE(iunit,outint4) 'NZETA',nzeta
      IF (lrfp) THEN
         WRITE(iunit,'(A)') '!----- RFP Parameters -----'
         WRITE(iunit,outboo) 'LRFP',lrfp
         WRITE (iunit,'(a,(1p,4ES22.12E3))') '  APHI = ',
     1                                   (aphi(n), n=1,SIZE(aphi))
      END IF
      WRITE(iunit,'(A)') '!----- Free Boundary Parameters -----'
      WRITE(iunit,outboo) 'LFREEB',lfreeb
      IF (lfreeb) THEN
         WRITE (iunit, '(2x,3a)') "MGRID_FILE = '",TRIM(mgrid_file),"'"
         DO n=1,SIZE(extcur)
            IF (extcur(n) == 0) CYCLE
            WRITE(iunit,'(2X,A,I3.3,A,ES22.12E3)')
     1      'EXTCUR(',n,') = ',extcur(n)
         END DO
         WRITE(iunit,outint4) 'NVACSKIP',nvacskip
         IF (TRIM(trip3d_file) /= 'NONE') WRITE (iunit, '(2x,3a)') 
     1             "TRIP3D_FILE = '",TRIM(trip3d_file),"'"  ! SAL - TRIP3D
      END IF
      WRITE(iunit,'(A)') '!----- Pressure Parameters -----'
      WRITE(iunit,outflt) 'GAMMA',gamma
      WRITE(iunit,outflt) 'BLOAT',bloat
      WRITE(iunit,outflt) 'SPRES_PED',spres_ped
      WRITE(iunit,outflt) 'PRES_SCALE',pres_scale
      WRITE(iunit,'(2x,3a)') "PMASS_TYPE = '",TRIM(pmass_type),"'"
      WRITE(iunit,'(a,(1p,4e22.14))')'  AM = ', (am(i-1), i=1,SIZE(am))
      i = minloc(am_aux_s(2:),DIM=1)
      IF (i > 4) THEN
         WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AM_AUX_S = ',
     1         (am_aux_s(n), n=1,i)
         WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AM_AUX_F = ', 
     1         (am_aux_f(n), n=1,i)
      END IF
      WRITE(iunit,'(A)') '!----- ANI/FLOW Parameters -----'
      WRITE(iunit,outflt) 'BCRIT',SQRT(bcrit)
      WRITE(iunit,'(2x,3a)') "PT_TYPE = '",TRIM(pt_type),"'"
      WRITE(iunit,'(a,(1p,4e22.14))')'  AT = ', (at(i-1), i=1,SIZE(at))
      i = minloc(at_aux_s(2:),DIM=1)
      IF (i > 4) THEN
         WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AT_AUX_S = ',
     1         (at_aux_s(n), n=1,i)
         WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AT_AUX_F = ', 
     1         (at_aux_f(n), n=1,i)
      END IF
      WRITE(iunit,'(2x,3a)') "PH_TYPE = '",TRIM(ph_type),"'"
      WRITE(iunit,'(a,(1p,4e22.14))')'  AH = ', (ah(i-1), i=1,SIZE(ah))
      i = minloc(ah_aux_s(2:),DIM=1)
      IF (i > 4) THEN
         WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AH_AUX_S = ',
     1         (ah_aux_s(n), n=1,i)
         WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AH_AUX_F = ', 
     1         (ah_aux_f(n), n=1,i)
      END IF
      WRITE(iunit,'(A)') '!----- Current/Iota Parameters -----'
      WRITE(iunit,outexp) 'CURTOR',curtor
      WRITE(iunit,outint) 'NCURR',ncurr
      WRITE (iunit, '(2x,3a)') "PIOTA_TYPE = '",TRIM(piota_type),"'"
      WRITE (iunit,'(a,(1p,4e22.14))') '  AI = ',(ai(n-1), n=1,SIZE(ai))
      i	= minloc(ai_aux_s(2:),DIM=1)
      IF (i > 4) THEN
         WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AI_AUX_S = ', 
     1         (ai_aux_s(n), n=1,i)
         WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AI_AUX_F = ', 
     1         (ai_aux_f(n), n=1,i)
      END IF
      WRITE (iunit, '(2x,3a)') "PCURR_TYPE = '",TRIM(pcurr_type),"'"
      WRITE (iunit,'(a,(1p,4ES22.12E3))')
     1 '  AC = ',(ac(n-1), n=1,SIZE(ac))
      i	= minloc(ac_aux_s(2:),DIM=1)
      IF (i > 4) THEN
         WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AC_AUX_S = ', 
     1         (ac_aux_s(n), n=1,i)
         WRITE (iunit,'(a,(1p,4ES22.12E3))') '  AC_AUX_F = ',
     1         (ac_aux_f(n), n=1,i)
      END IF
      WRITE(iunit,'(A)') '!----- Axis Parameters ----- '
      WRITE (iunit,'(a,(1p,4e22.14))') 
     1     '  RAXIS_CC = ',(raxis_cc(n), n=0,ntor)
      IF (lasym)
     1      WRITE (iunit,'(a,(1p,4ES22.12E3))') 
     2     '  RAXIS_CS = ',(raxis_cs(n), n=0,ntor)
      IF (lasym)
     1      WRITE (iunit,'(a,(1p,4ES22.12E3))') 
     2     '  ZAXIS_CC = ',(zaxis_cc(n), n=0,ntor)
      WRITE (iunit,'(a,(1p,4ES22.12E3))') 
     1     '  ZAXIS_CS = ',(zaxis_cs(n), n=0,ntor)
      WRITE(iunit,'(A)') '!----- Boundary Parameters -----'
      DO m = 0, mpol - 1
         DO n = -ntor, ntor
            IF ((rbc(n,m).ne.0) .or. (zbs(n,m).ne.0)) THEN
               WRITE(iunit,'(2(A,I4.3,A,I3.3,A,ES22.12E3))') 
     1         '  RBC(',n,',',m,') = ',rbc(n,m),
     2         '    ZBS(',n,',',m,') = ',zbs(n,m)
            END IF
            IF (.not. lasym) CYCLE
            IF ((rbs(n,m).ne.0) .or. (zbc(n,m).ne.0)) THEN
               WRITE(iunit,'(2(A,I4.3,A,I3.3,A,ES22.12E3))') 
     1         '  RBS(',n,',',m,') = ',rbs(n,m),
     2         '    ZBC(',n,',',m,') = ',zbc(n,m)
            END IF
         END DO
      END DO
      WRITE(iunit,'(A)') '/'
      RETURN
      END SUBROUTINE write_indata_namelist

      SUBROUTINE bcast_indata_namelist(local_master,local_comm,iflag)
      USE mpi_params
      USE mpi_inc
      IMPLICIT NONE
      INTEGER,INTENT(IN) :: local_master
      INTEGER,INTENT(INOUT) :: local_comm
      INTEGER,INTENT(INOUT) :: iflag
#if defined(MPI_OPT)
      iflag = 0
      ! Logicals
      CALL MPI_BCAST(lpofr,          1, MPI_LOGICAL, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(lmac,           1, MPI_LOGICAL, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(lfreeb,         1, MPI_LOGICAL, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(lrecon,         1, MPI_LOGICAL, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(loldout,        1, MPI_LOGICAL, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(ledge_dump,     1, MPI_LOGICAL, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(lasym,          1, MPI_LOGICAL, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(lforbal,        1, MPI_LOGICAL, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(lrfp,           1, MPI_LOGICAL, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(lmovie,         1, MPI_LOGICAL, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(lmove_axis,     1, MPI_LOGICAL, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(lwouttxt,       1, MPI_LOGICAL, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(ldiagno,        1, MPI_LOGICAL, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(lmoreiter,      1, MPI_LOGICAL, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(lfull3d1out,    1, MPI_LOGICAL, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(l_v3fit,        1, MPI_LOGICAL, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(lspectrum_dump, 1, MPI_LOGICAL, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(loptim,         1, MPI_LOGICAL, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(lgiveup,        1, MPI_LOGICAL, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(lbsubs,         1, MPI_LOGICAL, local_master, 
     1               local_comm, iflag)
      CALL MPI_BARRIER(local_comm,iflag)
      ! Integers
      CALL MPI_BCAST(nfp,                 1, MPI_INTEGER, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(ncurr,               1, MPI_INTEGER, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(nsin,                1, MPI_INTEGER, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(niter,               1, MPI_INTEGER, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(nstep,               1, MPI_INTEGER, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(nvacskip,            1, MPI_INTEGER, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(mpol,                1, MPI_INTEGER, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(ntor,                1, MPI_INTEGER, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(ntheta,              1, MPI_INTEGER, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(nzeta,               1, MPI_INTEGER, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(mfilter_fbdy,        1, MPI_INTEGER, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(nfilter_fbdy,        1, MPI_INTEGER, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(max_main_iterations, 1, MPI_INTEGER, local_master,
     1               local_comm, iflag)
      CALL MPI_BCAST(imse,                1, MPI_INTEGER, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(isnodes,             1, MPI_INTEGER, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(itse,                1, MPI_INTEGER, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(ipnodes,             1, MPI_INTEGER, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(iopt_raxis,          1, MPI_INTEGER, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(imatch_phiedge,      1, MPI_INTEGER, local_master, 
     1               local_comm, iflag)
      CALL MPI_BCAST(nflxs,               1, MPI_INTEGER, local_master, 
     1               local_comm, iflag)
      CALL MPI_BARRIER(local_comm,iflag)
      ! Integer Arrays
      CALL MPI_BCAST(ns_array,    100,              MPI_INTEGER, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(niter_array, 100,              MPI_INTEGER, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(nbfld,       nbsetsp,          MPI_INTEGER, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(indxflx,     nfloops,          MPI_INTEGER, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(indxbfld,    nbcoilsp*nbsetsp, MPI_INTEGER, 
     1               local_master, local_comm, iflag)
      CALL MPI_BARRIER(local_comm,iflag)
      ! Reals
      CALL MPI_BCAST(time_slice,       1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(curtor,           1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(delt,             1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(ftol,             1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(tcon0,            1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(gamma,            1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(phiedge,          1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(phidiam,          1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(sigma_current,    1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(sigma_delphid,    1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(tensi,            1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(tensp,            1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(tensi2,           1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(fpolyi,           1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(presfac,          1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(mseangle_offset,  1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(pres_offset,      1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(mseangle_offsetm, 1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(spres_ped,        1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(bloat,            1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(pres_scale,       1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(prec2d_threshold, 1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(bcrit,            1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(fgiveup,          1, MPI_DOUBLE_PRECISION, 
     1               local_master, local_comm, iflag)
      CALL MPI_BARRIER(local_comm,iflag)
      ! Real Arrays
      CALL MPI_BCAST(rbc,         (2*ntord+1)*(mpol1d+1),
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(zbs,         (2*ntord+1)*(mpol1d+1),
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(rbs,         (2*ntord+1)*(mpol1d+1),
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(zbc,         (2*ntord+1)*(mpol1d+1),
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(am,          21, 
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(ai,          21, 
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(ac,          21,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(aphi,        20,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(ah,          21,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(at,          21,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(am_aux_s,    ndatafmax,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(am_aux_f,    ndatafmax,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(ac_aux_s,    ndatafmax,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(ac_aux_f,    ndatafmax,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(ai_aux_s,    ndatafmax,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(ai_aux_f,    ndatafmax,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(ah_aux_s,    ndatafmax,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(ah_aux_f,    ndatafmax,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(at_aux_s,    ndatafmax,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(at_aux_f,    ndatafmax,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(raxis,       ntord+1,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(zaxis,       ntord+1,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(raxis_cc,    ntord+1,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(raxis_cs,    ntord+1,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(zaxis_cc,    ntord+1,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(zaxis_cs,    ntord+1,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(ftol_array,  100,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(extcur,      nigroup,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(mseprof,     nmse,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(rthom,       ntse,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(datathom,    ntse,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(sigma_thom,  ntse,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(rstark,      nmse,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(datastark,   nmse,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(sigma_stark, nmse,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(dsiobt,      nfloops,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(sigma_flux,  nfloops,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(bbc,         nbcoilsp*nbsetsp,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(sigma_b,     nbcoilsp*nbsetsp,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(psa,         ndatafmax,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(pfa,         ndatafmax,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(isa,         ndatafmax,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BCAST(ifa,         ndatafmax,
     1           MPI_DOUBLE_PRECISION, local_master, local_comm, iflag)
      CALL MPI_BARRIER(local_comm,iflag)
      ! Strings
      CALL MPI_BCAST(pcurr_type,      20,  MPI_CHARACTER, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(piota_type,      20,  MPI_CHARACTER, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(pmass_type,      20,  MPI_CHARACTER, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(pt_type,         20,  MPI_CHARACTER, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(ph_type,         20,  MPI_CHARACTER, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(mgrid_file,      200, MPI_CHARACTER, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(trip3d_file,     200, MPI_CHARACTER, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(precon_type,     10,  MPI_CHARACTER, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(arg1,            120, MPI_CHARACTER, 
     1               local_master, local_comm, iflag)
      CALL MPI_BCAST(input_extension, 100, MPI_CHARACTER, 
     1               local_master, local_comm, iflag)
      CALL MPI_BARRIER(local_comm,iflag)
      IF (iflag /= 0) RETURN
#endif
      iflag = 0
      RETURN
      END SUBROUTINE bcast_indata_namelist

      END MODULE vmec_input


