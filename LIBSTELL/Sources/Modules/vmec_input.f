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
     2           max_main_iterations
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

      REAL(rprec), DIMENSION(0:ntord) :: raxis, zaxis                !!Backwards compatibility: Obsolete
      REAL(rprec), DIMENSION(0:ntord) :: raxis_cc, raxis_cs,
     1                                   zaxis_cc, zaxis_cs
      REAL(rprec), DIMENSION(100) :: ftol_array
      REAL(rprec), DIMENSION(nigroup) :: extcur
      REAL(rprec), DIMENSION(nmse) :: mseprof
      REAL(rprec), DIMENSION(ntse) :: rthom, datathom, sigma_thom
      REAL(rprec), DIMENSION(nmse) :: rstark, datastark,
     1    sigma_stark
      REAL(rprec), DIMENSION(nfloops) :: dsiobt, sigma_flux
      REAL(rprec), DIMENSION(nbcoilsp,nbsetsp) :: bbc, sigma_b
      REAL(rprec), DIMENSION(ndatafmax) :: psa, pfa, isa, ifa
      LOGICAL :: lpofr, lmac, lfreeb, lrecon, loldout, ledge_dump,
     1           lasym, lforbal, lrfp, lmove_axis,
     2           lwouttxt, ldiagno,       ! J.Geiger: for txt- and diagno-output
     3           lmoreiter,               ! J.Geiger: if force residuals are not fulfilled add more iterations.
     4           lfull3d1out,             ! J.Geiger: to force full 3D1-output
     5           l_v3fit=.false.
     A         , lspectrum_dump, loptim           !!Obsolete
      LOGICAL :: lgiveup                  ! inserted M.Drevlak
      REAL(rprec) :: fgiveup              ! inserted M.Drevlak, giveup-factor for ftolv
      CHARACTER(len=200) :: mgrid_file
      CHARACTER(len=10)  :: precon_type
      CHARACTER(len=120) :: arg1
      CHARACTER(len=100) :: input_extension

      NAMELIST /indata/ mgrid_file, time_slice, nfp, ncurr, nsin,
     1   niter, nstep, nvacskip, delt, ftol, gamma, am, ai, ac, aphi,
     1   pcurr_type, pmass_type, piota_type,
     1   am_aux_s, am_aux_f, ai_aux_s, ai_aux_f, ac_aux_s, ac_aux_f,  ! J Hanson 2010-03-16
     1   ah, at, bcrit,                                               ! WAC (anisotropic pres)
     2   rbc, zbs, rbs, zbc, spres_ped, pres_scale, raxis_cc, zaxis_cs, 
     3   raxis_cs, zaxis_cc, mpol, ntor, ntheta, nzeta, mfilter_fbdy,
     3   nfilter_fbdy, niter_array,
     4   ns_array, ftol_array, tcon0, precon_type, prec2d_threshold,
     4   curtor, sigma_current, extcur,
     5   phiedge, psa, pfa, isa, ifa, imatch_phiedge, iopt_raxis, 
     6   tensi, tensp, mseangle_offset, mseangle_offsetm, imse, 
     7   isnodes, rstark, datastark, sigma_stark, itse, ipnodes, 
     8   presfac, pres_offset, rthom, datathom, sigma_thom, phidiam, 
     9   sigma_delphid, tensi2, fpolyi, nflxs, indxflx, dsiobt, 
     A   sigma_flux, nbfld, indxbfld, bloat, raxis, zaxis,
     A   bbc, sigma_b, lpofr, lforbal, lfreeb, lmove_axis, lrecon, lmac, 
     B   lasym, ledge_dump, lspectrum_dump, loptim, lrfp,
     C   loldout, lwouttxt, ldiagno, lfull3d1out, max_main_iterations,     ! J Geiger 2010-05-04
     D   lgiveup,fgiveup                                                   ! M.Drevlak 2012-05-10

      NAMELIST /mseprofile/ mseprof

      CONTAINS

      SUBROUTINE read_indata_namelist (iunit, istat)
      INTEGER :: iunit, istat

!
!     INITIALIZATIONS
!
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
      lfreeb = .true.
      lmove_axis = .true.
      lmac = .false.
      lforbal = .true.
      lasym = .false.
      lrfp = .false.
      loldout = .false.        ! J Geiger 2010-05-04 start
      ldiagno = .false.
      lgiveup = .false.        ! inserted M.Drevlak
      fgiveup = 3.E+01_dp      ! inserted M.Drevlak
      lfull3d1out = .false.
      lmoreiter = .false.      ! default value if no max_main_iterations given.
      max_main_iterations = 1  ! to keep a presumably expected standard behavior.
#if defined(NETCDF)
      lwouttxt = .false.       ! to keep functionality as expected with netcdf
#else
      lwouttxt = .true.        ! and without netcdf
#endif
                               ! J Geiger 2010-05-04 end

      pcurr_type = 'power_series'
      piota_type = 'power_series'
      pmass_type = 'power_series'

!     ANISTROPY PARAMETERS
      bcrit = 1
      at(0) = 1;  at(1:) = 0
      ah = 0

!
!     BACKWARDS COMPATIBILITY
!
      raxis = 0;  zaxis = 0
      
      READ (iunit, nml=indata, iostat=istat)

      IF (ALL(niter_array == -1)) niter_array = niter
      WHERE (raxis .ne. 0._dp) raxis_cc = raxis
      WHERE (zaxis .ne. 0._dp) zaxis_cs = zaxis
      IF(max_main_iterations .gt. 1) lmoreiter=.true.  !J Geiger: if more iterations are requested.

      END SUBROUTINE read_indata_namelist

      SUBROUTINE read_mse_namelist (iunit, istat)
      INTEGER :: iunit, istat

      READ (iunit, nml=mseprofile, iostat=istat)

      END SUBROUTINE read_mse_namelist

      END MODULE vmec_input


