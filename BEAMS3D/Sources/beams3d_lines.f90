!-----------------------------------------------------------------------
!     Module:        beams3d_lines
!     Authors:       S. Lazerson (lazerson@pppl.gov)
!     Date:          02/21/2012
!     Description:   This module contains the FIELDLINES field line
!                    variables.
!-----------------------------------------------------------------------
      MODULE beams3d_lines
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE stel_kinds, ONLY: rprec
      
!-----------------------------------------------------------------------
!     Module Variables
!          myline    Dummy index
!          nparticles    Number of Particles
!          nsteps    Number of integration steps along fieldline
!          R_lines   Radial locations along fieldline [m] (npoinc per field period)
!          Z_lines   Vertical locations along field line [m]
!          PHI_lines Toroidal locations along field line [radians]
!-----------------------------------------------------------------------
      IMPLICIT NONE
      LOGICAL  ::  ltherm
      INTEGER  ::  ns_prof1, ns_prof2, ns_prof3, ns_prof4, ns_prof5, nsh_prof4
      INTEGER  :: nparticles, nsteps, myline, mybeam, mytdex, myend, mystart_save, myend_save
      INTEGER  :: win_epower, win_ipower, win_ndot, win_dense, win_jprof, win_dist5d, win_dist5d_fida
      REAL(rprec) :: xlast,ylast,zlast ! for storing position
      REAL(rprec) :: moment, mycharge, myZ, mymass, myv_neut(3), my_end, &
                     myqm, rand_prob, cum_prob, tau, next_t, &
                     partvmax, fact_crit, fact_pa, fact_vsound, fact_kick, &
                     fact_coul, &
                     partpmax, h2_prof, h3_prof, h4_prof, h5_prof, r_h, z_h, p_h, e_h, pi_h, E_by_v
      LOGICAL, ALLOCATABLE     :: neut_lines(:,:)
      INTEGER, ALLOCATABLE     :: end_state(:)
      REAL(rprec), ALLOCATABLE :: shine_through(:), shine_port(:), GFactor(:), t_last(:)
      REAL(rprec), DIMENSION(:,:), POINTER :: ndot_prof(:,:),epower_prof(:,:), &
                                  ipower_prof(:,:),j_prof(:,:), dense_prof(:,:)
      REAL(rprec), DIMENSION(:,:,:,:,:,:), POINTER :: dist5d_prof
      REAL(rprec), DIMENSION(:,:,:,:,:), POINTER :: dist5d_fida
      REAL(rprec), ALLOCATABLE :: R_lines(:,:),Z_lines(:,:),PHI_lines(:,:),vll_lines(:,:),moment_lines(:,:),&
                                  S_lines(:,:),U_lines(:,:),B_lines(:,:), &
                                  vr_lines(:,:),vphi_lines(:,:),vz_lines(:,:)

      END MODULE beams3d_lines
