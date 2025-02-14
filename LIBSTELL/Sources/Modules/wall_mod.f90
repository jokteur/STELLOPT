!-----------------------------------------------------------------------
!     Module:        wall_mod
!     Authors:       S. Lazerson (lazerson@pppl.gov),
!                    D.J. Engels (d.j.engels@student.tue.nl)
!     Date:          May 2021
!     Description:   This module handles defining a wall as a set of
!                    triangular facets which can be used to calculate
!                    if and where a particle hits the mesh.
!                    Note that the user must deallocate face and vertex
!                    after calling wall_load_mn.
!-----------------------------------------------------------------------
      MODULE wall_mod
!-----------------------------------------------------------------------
!     Libraries
!-----------------------------------------------------------------------
      USE safe_open_mod
      IMPLICIT NONE

!-----------------------------------------------------------------------
!     Types
!         Wall_Type  The full wall, information about the vertices and the number of blocks
!         Block      A block is a part of the uniform grid, an area in space with triangles in it      
!-----------------------------------------------------------------------
      TYPE block
         INTEGER :: nfaces  ! number of faces in block
         LOGICAL :: isshared  ! whether or not shared memory active
         ! shared memory point to faces in block
         INTEGER :: win_face
         ! bounds of block in x/y/z
         DOUBLE PRECISION :: rmin(3), rmax(3)
         ! non-shared memory pointer to faces in block
         INTEGER,          DIMENSION(:),   POINTER :: face => null()
      END TYPE block

      TYPE wall_type
         ! number of blocks
         ! and the step that has to be done in the list to move one in x/y/z-direction
         integer :: nblocks, step(3)
         ! Bounds of total grid in x/y/z
         DOUBLE PRECISION :: rmin(3), rmax(3)
         ! Step of grid in space
         DOUBLE PRECISION :: stepsize
         ! Number of blocks in each direction
         INTEGER          :: br(3)
         TYPE (block), DIMENSION(:), POINTER :: blocks => null()
      END TYPE wall_type


!-----------------------------------------------------------------------
!     Module Variables
!         
!-----------------------------------------------------------------------
      LOGICAL            :: lwall_loaded
      INTEGER            :: nvertex, nface
      INTEGER, POINTER :: face(:,:)
      DOUBLE PRECISION, POINTER   :: vertex(:,:)
      CHARACTER(LEN=256) :: machine_string
      CHARACTER(LEN=256) :: date

      INTEGER, PRIVATE                         :: mystart, myend, mydelta, ik_min
      INTEGER, PRIVATE                         :: shar_rank, shar_size
      INTEGER, PRIVATE                         :: win_vertex, win_face

      DOUBLE PRECISION, PRIVATE, POINTER       :: invDenom(:)
                                                  
      DOUBLE PRECISION, PRIVATE, PARAMETER      :: zero = 0.0D+0
      DOUBLE PRECISION, PRIVATE, PARAMETER      :: one  = 1.0D+0
      DOUBLE PRECISION, PRIVATE, PARAMETER      :: epsilon = 1D-6

      LOGICAL, PRIVATE, PARAMETER               :: ldebug = .FALSE.
      LOGICAL, PRIVATE                          :: lverb  = .FALSE.

!------------------ Variables for naive approach
      
      INTEGER, POINTER :: ihit_array(:)
      INTEGER, PRIVATE                         :: win_fn, win_a0, win_v0, win_v1, &
                                                  win_dot00, win_dot01, win_dot11, &
                                                  win_d, win_ihit, win_invDenom

      DOUBLE PRECISION, PRIVATE, POINTER       :: FN(:,:), d(:)
      DOUBLE PRECISION, PRIVATE, POINTER       :: A0(:,:), V0(:,:), V1(:,:), V2(:,:),&
                                                  DOT00(:), DOT01(:), DOT02(:),&
                                                  DOT11(:), DOT12(:)

!------------------ Variables for accelerated approach

      TYPE(wall_type), PRIVATE :: wall
      TYPE(block), PRIVATE:: b

      INTEGER, PRIVATE                         :: nfaces_per_block = 50

      INTEGER, PRIVATE                         :: win_bface
      INTEGER, POINTER                         :: bface(:)
      INTEGER                                  :: nface_block
      
!-----------------------------------------------------------------------
!     Subroutines
!         wall_load_txt:   Loads triangular mesh from file
!         wall_load_mn:    Creates wall from harmonics
!         wall_load_seg:   Creates wall from segments
!         wall_dump:       Dumps triangulation data
!         wall_info:       Prints wall info
!         collide:         Calculates collision with wall
!                          Has to implementations, for double and float
!         uncount_wall_hit Reduces hit count for last location by one
!         wall_free:       Frees module memory
!     Subroutines (wall acceleration creation)
!         ACCELERATE_WALL:   Main function to create an accelerated wall from a vertex/face list
!         GET_BLOCK_SIZE :   Gets the desired block size
!         SET_NFACE      :   Sets the nfaces_per_block variable
!         CREATE_BLOCKS  :   Creates the grid of blocks
!         FILL_BLOCKS    :   Fills the blocks with vertices
!         WRITE_WALL     :   Writes the (accelerated) wall mesh to .txt
!-----------------------------------------------------------------------
!     Functions
!         get_wall_ik      Gets index of last hit
!         get_wall_area    Gets area of certain wall index
!-----------------------------------------------------------------------
      INTERFACE collide
         MODULE PROCEDURE collide_double, collide_float
      END INTERFACE

      INTERFACE free_mpi_array
         MODULE PROCEDURE free_mpi_array1d_int, free_mpi_array1d_flt, free_mpi_array1d_dbl, &
                          free_mpi_array2d_int, free_mpi_array2d_dbl
      END INTERFACE

      PRIVATE :: LINE_BOX_INTERSECTION
      PRIVATE :: mpialloc_1d_int,mpialloc_1d_dbl,mpialloc_2d_int,mpialloc_2d_dbl
      PRIVATE :: free_mpi_array1d_int, free_mpi_array1d_dbl, free_mpi_array2d_int, free_mpi_array2d_dbl
      CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Wall Constructor
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE INIT_BLOCK(this,xmin,xmax,ymin,ymax,zmin,zmax,nface_block,istat,comm,shar_comm)
      !-----------------------------------------------------------------------
      ! init_block: Initializes a single block of the accelerated uniform grid
      !-----------------------------------------------------------------------
      ! param[in]: this. Current block
      ! param[in]: xmin. Minimum x position of block
      ! param[in]: xmax. Maximum x position of block
      ! param[in]: ymin. Minimum y position of block
      ! param[in]: ymax. Maximum y position of block
      ! param[in]: zmin. Minimum z position of block
      ! param[in]: zmax. Maximum z position of block
      ! param[in]: nface_block. Number of faces in block
      ! param[in, out]: istat. Integer that shows error if != 0
      ! param[in, out]: comm. MPI communicator, handles communication between nodes
      ! param[in, out]: shar_comm. Shared MPI communicator, handles shared memory
      !-----------------------------------------------------------------------
#if defined(MPI_OPT)
         USE mpi
#endif
         TYPE(block) :: this
         DOUBLE PRECISION, INTENT(in) :: xmin, ymin, zmin, xmax, ymax, zmax
         INTEGER, INTENT(in) :: nface_block
         INTEGER, INTENT(inout) :: istat
         INTEGER, INTENT(inout), OPTIONAL :: comm, shar_comm
         LOGICAL :: lshar_comm
         INTEGER :: i

         lshar_comm = .FALSE.
#if defined(MPI_OPT)
         IF (PRESENT(comm)) THEN
            lshar_comm = .TRUE. 
            CALL MPI_Bcast(xmin,1,MPI_DOUBLE_PRECISION,0,shar_comm,istat)
            CALL MPI_Bcast(xmax,1,MPI_DOUBLE_PRECISION,0,shar_comm,istat)
            CALL MPI_Bcast(ymin,1,MPI_DOUBLE_PRECISION,0,shar_comm,istat)
            CALL MPI_Bcast(ymax,1,MPI_DOUBLE_PRECISION,0,shar_comm,istat)
            CALL MPI_Bcast(zmin,1,MPI_DOUBLE_PRECISION,0,shar_comm,istat)
            CALL MPI_Bcast(zmax,1,MPI_DOUBLE_PRECISION,0,shar_comm,istat)
         END IF
#endif
         ! Set bounds block
         this%rmin(1) = xmin
         this%rmax(1) = xmax
         this%rmin(2) = ymin
         this%rmax(2) = ymax
         this%rmin(3) = zmin
         this%rmax(3) = zmax
         this%nfaces = nface_block

         ! Keep track of size of wall
         DO i=1,3
            IF (this%rmin(i) < wall%rmin(i)) wall%rmin(i) = this%rmin(i)
            IF (this%rmax(i) > wall%rmax(i)) wall%rmax(i) = this%rmax(i)
         END DO

         IF (ldebug) WRITE(6, *) 'Block intialized', xmin, xmax, ymin, ymax, zmin, zmax, nface_block

         ! Only allocate if there is actually faces in this block
         IF (nface_block > 0) THEN
            IF (lshar_comm) THEN 
#if defined(MPI_OPT)
               CALL MPI_BARRIER(shar_comm,istat)
               IF (istat/=0) RETURN
               CALL mpialloc_1d_int(this%face, nface_block, shar_rank, 0, shar_comm, this%win_face)
               this%isshared = .TRUE.
#endif
            ELSE
               ALLOCATE(this%face(nface_block),STAT=istat)
               this%isshared = .FALSE.
            END IF

            ! Copy faces into block
            this%face(:) = bface(:)
            ! Block MPI
#if defined(MPI_OPT)
            IF (lshar_comm) CALL MPI_BARRIER(shar_comm,istat)
#endif
         ELSE
            IF (ldebug) WRITE(6, *) 'Skipped due to nface_block 0', shar_rank
         END IF
      END SUBROUTINE INIT_BLOCK
      
      SUBROUTINE wall_load_txt(filename,istat,verb,comm)
      !-----------------------------------------------------------------------
      ! wall_load_txt: Loads triangular mesh from file
      !-----------------------------------------------------------------------
      ! param[in]: filename. The file name to load in
      ! param[in, out]: istat. Integer that shows error if != 0
      ! param[in]: verb: Verbosity. True or false
      ! param[in, out]: comm. MPI communicator, handles shared memory
      !-----------------------------------------------------------------------
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER, INTENT(inout)       :: istat
      LOGICAL, INTENT(in), OPTIONAL :: verb
      INTEGER, INTENT(inout), OPTIONAL :: comm
      DOUBLE PRECISION :: xmin, ymin, zmin, xmax, ymax, zmax
      INTEGER :: iunit, ik, i, dex1, dex2, dex3
      INTEGER :: shar_comm
      LOGICAL :: shared, lwall_acc, lcomm

      IF (PRESENT(verb)) lverb = verb
      IF (lverb) WRITE(6,*) '-----  Creating wall mesh  -----'
      
      shar_rank = 0; shar_size = 1;
      lwall_loaded = .false.
      lcomm = .false.
      ! initialize MPI
#if defined(MPI_OPT)
      IF (PRESENT(comm)) THEN
         lcomm = .true.
         CALL MPI_COMM_SPLIT_TYPE(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shar_comm, istat)
         CALL MPI_COMM_RANK( shar_comm, shar_rank, istat )
         CALL MPI_COMM_SIZE( shar_comm, shar_size, istat)
      END IF
#endif
      ! open file, return if fails
      CALL safe_open(iunit,istat,TRIM(filename),'old','formatted')
      IF (istat/=0) RETURN
      lwall_acc = .false.
      ! read info
      IF (shar_rank == 0) THEN
         READ(iunit,'(A)') machine_string
         READ(iunit,'(A)') date
         READ(iunit,*) nvertex,nface
         ! Check if nvertex & nface  = 0 -> accelerated structure
         IF (nvertex == 0 .and. nface == 0) THEN
            IF (ldebug) WRITE(6,*) 'Accelerated'
            lwall_acc = .true.
            READ(iunit,*) nvertex,nface    
         END IF
      END IF
      ! Broadcast info to MPI and allocate vertex and face info
      IF (lcomm) THEN
#if defined(MPI_OPT)
         CALL MPI_Bcast(lwall_acc,1,MPI_LOGICAL,0,shar_comm,istat)
         CALL MPI_Bcast(nvertex,1,MPI_INTEGER,0,shar_comm,istat)
         CALL MPI_Bcast(nface,1,MPI_INTEGER,0,shar_comm,istat)
         CALL mpialloc_2d_dbl(vertex,nvertex,3,shar_rank,0,shar_comm,win_vertex)
         CALL mpialloc_2d_int(face,nface,3,shar_rank,0,shar_comm,win_face)
         shared = .true.
#endif
      ELSE
         ! if no MPI, allocate everything on one node
         ALLOCATE(vertex(nvertex,3),face(nface,3),STAT=istat)
         shared = .false.
      END IF

      ! read in the mesh on allocated memory
      IF (ldebug) WRITE(6, *) 'Vertex & face allocating & reading. MPI Rank: ', shar_rank
      IF (istat/=0) RETURN
      IF (shar_rank == 0) THEN
         DO ik = 1, nvertex
            READ(iunit,*) vertex(ik,1),vertex(ik,2),vertex(ik,3)
         END DO
         DO ik=1,nface
            READ(iunit,*) face(ik,1),face(ik,2),face(ik,3)
         END DO
      END IF

      ! allocate memory for information about the mesh
      IF (ldebug) WRITE(6, *) 'Pre-calculation allocation & reading. MPI Rank: ', shar_rank
#if defined(MPI_OPT)
      IF (lcomm) CALL MPI_BARRIER(shar_comm,istat)
      IF (istat/=0) RETURN
#endif
      IF (lcomm) THEN
         CALL mpialloc_2d_dbl(A0,nface,3,shar_rank,0,shar_comm,win_a0)
         CALL mpialloc_2d_dbl(V0,nface,3,shar_rank,0,shar_comm,win_v0)
         CALL mpialloc_2d_dbl(V1,nface,3,shar_rank,0,shar_comm,win_v1)
         CALL mpialloc_2d_dbl(FN,nface,3,shar_rank,0,shar_comm,win_fn)
         mydelta = CEILING(REAL(nface) / REAL(shar_size))
         mystart = 1 + shar_rank*mydelta
         myend   = mystart + mydelta
         IF (myend > nface) myend=nface
      ELSE
         ALLOCATE(A0(nface,3),V0(nface,3),V1(nface,3),&
                  FN(nface,3),STAT=istat)
         mystart = 1; myend = nface
      END IF

      IF (istat/=0) RETURN
      ! Precalculate information about mesh
      ! Calculate the face normal
      ! V  = Vertex1-Vertex0
      ! W  = Vertex2-Vertex0
      ! FN = VxW/|VxW|
      ! d  = -Vertex0.FN (. is dot product) (note weve absorbed the negative)
      DO ik = mystart, myend
         dex1 = face(ik,1)
         dex2 = face(ik,2)
         dex3 = face(ik,3)
         A0(ik,:) = vertex(dex1,:)
         V0(ik,:)  = vertex(dex3,:)-vertex(dex1,:)
         V1(ik,:)  = vertex(dex2,:)-vertex(dex1,:)
         FN(ik,1) = (V1(ik,2)*V0(ik,3))-(V1(ik,3)*V0(ik,2))
         FN(ik,2) = (V1(ik,3)*V0(ik,1))-(V1(ik,1)*V0(ik,3))
         FN(ik,3) = (V1(ik,1)*V0(ik,2))-(V1(ik,2)*V0(ik,1))
      END DO
#if defined(MPI_OPT)
      IF (lcomm) CALL MPI_BARRIER(shar_comm,istat)
#endif
      ! Check for zero area
      IF (ANY(SUM(FN*FN,DIM=2)==zero)) THEN
         istat=-327
         RETURN
      END IF
      ! allocate memory for information about mesh triangles 
      IF (lcomm) THEN
         CALL mpialloc_1d_dbl(DOT00,nface,shar_rank,0,shar_comm,win_dot00)
         CALL mpialloc_1d_dbl(DOT01,nface,shar_rank,0,shar_comm,win_dot01)
         CALL mpialloc_1d_dbl(DOT11,nface,shar_rank,0,shar_comm,win_dot11)
         CALL mpialloc_1d_dbl(invDenom,nface,shar_rank,0,shar_comm,win_invDenom)
         CALL mpialloc_1d_dbl(d,nface,shar_rank,0,shar_comm,win_d)
         CALL mpialloc_1d_int(ihit_array,nface,shar_rank,0,shar_comm,win_ihit)
         mydelta = CEILING(REAL(nface) / REAL(shar_size))
         mystart = 1 + shar_rank*mydelta
         myend   = mystart + mydelta
         IF (myend > nface) myend=nface
      ELSE
         ALLOCATE(DOT00(nface), DOT01(nface),&
                  DOT11(nface), invDenom(nface),&
                  STAT=istat)
         ALLOCATE(d(nface),STAT=istat)
         ALLOCATE(ihit_array(nface),STAT=istat)
         mystart = 1; myend = nface
      END IF

      ! if no error, calculate information about mesh triangles
      IF (istat/=0) RETURN
      DO ik = mystart, myend
         ihit_array(ik) = 0
         DOT00(ik) = V0(ik,1)*V0(ik,1) + V0(ik,2)*V0(ik,2) + V0(ik,3)*V0(ik,3)
         DOT01(ik) = V0(ik,1)*V1(ik,1) + V0(ik,2)*V1(ik,2) + V0(ik,3)*V1(ik,3)
         DOT11(ik) = V1(ik,1)*V1(ik,1) + V1(ik,2)*V1(ik,2) + V1(ik,3)*V1(ik,3)
         d(ik)     = FN(ik,1)*A0(ik,1) + FN(ik,2)*A0(ik,2) + FN(ik,3)*A0(ik,3)
         invDenom(ik) = one / (DOT00(ik)*DOT11(ik) - DOT01(ik)*DOT01(ik))
      END DO
! sync MPI
#if defined(MPI_OPT)
      IF (lcomm) CALL MPI_BARRIER(shar_comm, istat)
#endif
      ! If accelerated, read in uniform grid
      if (lwall_acc) THEN
         IF (ldebug) WRITE(6, *) 'Start accelerated reading. MPI Rank: ', shar_rank
         IF (shar_rank == 0) THEN
            READ(iunit, *) wall%nblocks, wall%step(1), wall%step(2), wall%step(3)
            READ(iunit, *) wall%stepsize, wall%br(1), wall%br(2), wall%br(3)
         END IF
         ! If MPI, broadcast this info again
         IF (lcomm) THEN
#if defined(MPI_OPT)
            CALL MPI_Bcast(wall%nblocks,1,MPI_INTEGER,0,shar_comm,istat)
            CALL MPI_Bcast(wall%step,3,MPI_INTEGER,0,shar_comm,istat)
            CALL MPI_Bcast(wall%stepsize,1,MPI_DOUBLE_PRECISION,0,shar_comm,istat)
            CALL MPI_Bcast(wall%br,3,MPI_INTEGER,0,shar_comm,istat)
#endif
         END IF
         ! Allocate room for all the blocks
         ALLOCATE(wall%blocks(wall%nblocks), STAT=istat)
         IF (istat/=0) RETURN

         ! This is the maximum wall size possible
         DO i=1,3 
            wall%rmin(i) = 1D+20
            wall%rmax(i) = -1D+20
         END DO

         ! Start looping over all blocks
         IF (ldebug) WRITE(6, *) 'Start blocks loop. Rank and wall info: ', shar_rank, wall%nblocks, wall%br         
         DO ik=1, wall%nblocks
            IF (ldebug) WRITE(6, *) 'In loop', shar_rank, ik
            ! Read in the bounds of the block and the number of faces
            IF (shar_rank == 0) THEN
               READ(iunit, *) xmin,xmax,ymin,ymax,zmin,zmax
               READ(iunit, *) nface_block
               IF (ldebug) WRITE(6, *) 'Block reading allocation: ', ik, nface_block
            END IF

#if defined(MPI_OPT)
            IF (lcomm) CALL MPI_Bcast(nface_block,1,MPI_INTEGER,0,shar_comm,istat)
#endif
            ! Only allocate and read faces if there are faces to read for this block
            IF (nface_block > 0) THEN
               IF (lcomm) THEN
                  CALL mpialloc_1d_int(bface,nface_block,shar_rank,0,shar_comm,win_bface)
               ELSE
                  ! if no MPI, allocate everything on one node
                  ALLOCATE(bface(nface_block),STAT=istat)
               END IF

               ! Read all the faces
               IF (istat/=0) RETURN
               IF (shar_rank == 0) THEN
                  IF (ldebug) WRITE(6, *) 'Block reading'
                  DO i=1, nface_block
                     READ(iunit,*) bface(i)
                  END DO
               END IF
            END IF

#if defined(MPI_OPT)
            IF (lcomm) CALL MPI_BARRIER(shar_comm,istat)
#endif

            ! Initialize the block
            IF (ldebug .and. shar_rank == 0) WRITE(6, *) 'Init block: ', ik
            CALL INIT_BLOCK(wall%blocks(ik),xmin,xmax,ymin,ymax,zmin,zmax,nface_block,istat,comm,shar_comm)

            IF (ldebug .and. shar_rank == 0) WRITE(6, *) 'Init block done: ', ik
            ! Also only deallocate if nface > 0
            IF (nface_block > 0) CALL free_mpi_array(win_bface, bface, shared)
         END DO
      ELSE
         ! Else create an accelerated wall manually
         IF (lcomm) THEN
            CALL ACCELERATE_WALL(istat, comm, shar_comm)
         ELSE
            CALL ACCELERATE_WALL(istat)
         END IF
#if defined(MPI_OPT)
         IF (PRESENT(comm)) CALL MPI_BARRIER(shar_comm, istat)
#endif
         IF (istat/=0) RETURN
      END IF
      ! close file
      CLOSE(iunit)
      ! sync MPI
      IF (lcomm) THEN
#if defined(MPI_OPT)
         CALL MPI_BARRIER(shar_comm, istat)
         CALL MPI_COMM_FREE(shar_comm, istat)
#endif
      END IF
      IF (ldebug) WRITE(6, *) 'Done reading wall from txt: ', shar_rank
      ! set wall as loaded and return
      lwall_loaded = .true.
      RETURN
      END SUBROUTINE wall_load_txt

      SUBROUTINE wall_load_mn(Rmn,Zmn,xm,xn,mn,nu,nv,verb,comm,Rmn2,Zmn2)
      !-----------------------------------------------------------------------
      ! wall_load_mn: Creates wall from harmonics
      !-----------------------------------------------------------------------
      ! param[in]: Rmn. Harmonics in R direction (cos)
      ! param[in]: Zmn. Harmonics in Z direction (sin)
      ! param[in]: xm. Poloidal Harmonic Array
      ! param[in]: xn. Toroidal Harmonic Array
      ! param[in]: mn. Total number of modes
      ! param[in]: nu. Total poloidal gridpoints
      ! param[in]: nv. Total toroidal gridpoints
      ! param[in]: verb: Verbosity. True or false
      ! param[in, out]: comm. MPI communicator, handles shared memory
      ! param[in]: Rmn2. Harmonics in R direction (sin)
      ! param[in]: Zmn2. Harmonics in Z direction (cos)
      !-----------------------------------------------------------------------
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: Rmn(mn), Zmn(mn), xm(mn), xn(mn)
      INTEGER, INTENT(in) :: mn, nu, nv
      LOGICAL, INTENT(in), OPTIONAL :: verb
      INTEGER, INTENT(inout), OPTIONAL :: comm
      DOUBLE PRECISION, INTENT(in), OPTIONAL :: Rmn2(mn), Zmn2(mn)
      INTEGER :: u, v, i, j, istat, dex1, dex2, dex3, ik, nv2
      INTEGER :: shar_comm
      LOGICAL :: shared, lcomm
      DOUBLE PRECISION :: pi2, th, zt, pi
      DOUBLE PRECISION, ALLOCATABLE :: r_temp(:,:),z_temp(:,:),x_temp(:,:),y_temp(:,:)


      lcomm = .FALSE.
      IF (PRESENT(comm)) lcomm = .TRUE.
      IF (PRESENT(verb)) lverb = verb
      IF (lverb) WRITE(6,*) '-----  Creating wall mesh  -----'

      ! create info usually read from file manually
      machine_string = '          HARMONICS'
      date = '      TODAY'
      pi2 = 8.0D+00 * ATAN(one)
      pi  = 4.0E+00 * ATAN(one)
      nv2 = nv/2
      shar_rank = 0; shar_size = 1;
      ! initialize MPI
      IF (lcomm) THEN
#if defined(MPI_OPT)
         CALL MPI_COMM_SPLIT_TYPE(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shar_comm, istat)
         CALL MPI_COMM_RANK( shar_comm, shar_rank, istat )
         CALL MPI_COMM_SIZE( shar_comm, shar_size, istat)
#endif
      END IF

      ! If shared memory rank not changed above, create temporary locations from harmonic calculations
      IF (shar_rank==0)THEN
         ALLOCATE(r_temp(nu,nv),z_temp(nu,nv),x_temp(nu,nv),y_temp(nu,nv))
         r_temp(:,:) = 0
         z_temp(:,:) = 0
         DO u = 1, nu
            DO v = 1, nv
               DO i = 1, mn
                  th = pi2*DBLE(u-1)/DBLE(nu)
                  zt = pi2*DBLE(v-1)/DBLE(nv)
                  r_temp(u,v) = r_temp(u,v) + Rmn(i)*DCOS(xm(i)*th+xn(i)*zt)
                  z_temp(u,v) = z_temp(u,v) + Zmn(i)*DSIN(xm(i)*th+xn(i)*zt)
               END DO
            END DO
         END DO
         IF (PRESENT(Rmn2)) THEN
            DO u = 1, nu
               DO v = 1, nv
                  DO i = 1, mn
                     th = pi2*DBLE(u-1)/DBLE(nu)
                     zt = pi2*DBLE(v-1)/DBLE(nv)
                     r_temp(u,v) = r_temp(u,v) + Rmn2(i)*DSIN(xm(i)*th+xn(i)*zt)
                     z_temp(u,v) = z_temp(u,v) + Zmn2(i)*DCOS(xm(i)*th+xn(i)*zt)
                  END DO
               END DO
            END DO
         END IF
         DO v = 1, nv
            zt = pi2*DBLE(v-1)/DBLE(nv)
            x_temp(:,v) = r_temp(:,v) * DCOS(zt)
            y_temp(:,v) = r_temp(:,v) * DSIN(zt)
         END DO
      END IF

      nvertex = nu*nv
      nface   = 2*nu*nv
      istat = 0
      ! allocate shared memory with of mesh
      ! calculate which part each node has to calculate using mydelta
      IF (lcomm) THEN
         CALL mpialloc_2d_dbl(vertex,nvertex,3,shar_rank,0,shar_comm,win_vertex)
         CALL mpialloc_2d_int(face,nface,3,shar_rank,0,shar_comm,win_face)
         shared = .true.
         mydelta = CEILING(REAL(nface) / REAL(shar_size))
         mystart = 1 + shar_rank*mydelta
         myend   = mystart + mydelta
         IF (myend > nface) myend=nface
      ELSE
         ! if no MPI, allocate everything on one node
         ALLOCATE(vertex(nvertex,3),face(nface,3),STAT=istat)
         shared = .false.
         mystart = 1; myend=nface
      END IF

      i = 1  ! Tracks vertex index
      j = 1 ! Tracks face index
      ! Do further calculations to create mesh if shared memory rank is zero
      IF (shar_rank==0)THEN
         DO v = 1, nv-1
            DO u = 1, nu-1
               vertex(i,1) = x_temp(u,v)
               vertex(i,2) = y_temp(u,v)
               vertex(i,3) = z_temp(u,v)
               face(j,1) = i
               face(j,2) = i + nu
               face(j,3) = i + 1
               j = j + 1
               face(j,1) = i + nu
               face(j,2) = i + nu + 1
               face(j,3) = i + 1
               j = j + 1
               i=i+1
            END DO
            u = nu
            vertex(i,1) = x_temp(u,v)
            vertex(i,2) = y_temp(u,v)
            vertex(i,3) = z_temp(u,v)
            face(j,1) = i
            face(j,2) = i+nu
            face(j,3) = i-nu+1
            j = j + 1
            face(j,1) = i + nu
            face(j,2) = i + 1
            face(j,3) = i-nu+1
            j = j + 1
            i=i+1
         END DO
         v = nv
         DO u = 1, nu - 1
            vertex(i,1) = x_temp(u,v)
            vertex(i,2) = y_temp(u,v)
            vertex(i,3) = z_temp(u,v)
            face(j,1) = i
            face(j,2) = i - nu*(nv-1)
            face(j,3) = i + 1
            j = j + 1
            face(j,1) = i - nu*(nv-1)
            face(j,2) = i - nu*(nv-1) + 1
            face(j,3) = i + 1
            j = j + 1
            i=i+1
         END DO
         u = nu
         vertex(i,1) = x_temp(u,v)
         vertex(i,2) = y_temp(u,v)
         vertex(i,3) = z_temp(u,v)
         face(j,1) = i
         face(j,2) = nu
         face(j,3) = i - nu + 1
         j = j + 1
         face(j,1) = nu
         face(j,2) = 1
         face(j,3) = i - nu + 1
         j = j + 1
         i=i+1
         ! remove temperary information
         DEALLOCATE(r_temp,z_temp,x_temp,y_temp)
      END IF
      ! if using MPI, wait here
#if defined(MPI_OPT)
      IF (lcomm) CALL MPI_BARRIER(shar_comm,istat)
#endif
      ! allocate memory for information about the mesh
      IF (lcomm) THEN
         CALL mpialloc_2d_dbl(A0,nface,3,shar_rank,0,shar_comm,win_a0)
         CALL mpialloc_2d_dbl(V0,nface,3,shar_rank,0,shar_comm,win_v0)
         CALL mpialloc_2d_dbl(V1,nface,3,shar_rank,0,shar_comm,win_v1)
         CALL mpialloc_2d_dbl(FN,nface,3,shar_rank,0,shar_comm,win_fn)
         mydelta = CEILING(REAL(nface) / REAL(shar_size))
         mystart = 1 + shar_rank*mydelta
         myend   = mystart + mydelta
         IF (myend > nface) myend=nface
      ELSE
         ALLOCATE(A0(nface,3),V0(nface,3),V1(nface,3),&
                  FN(nface,3),STAT=istat)
         mystart = 1; myend = nface
      END IF

      IF (istat/=0) RETURN
      ! Precalculate information about mesh
      ! Calculate the face normal
      ! W  = Vertex1-Vertex0
      ! W  = Vertex2-Vertex0
      ! FN = VxW/|VxW|
      ! d  = -Vertex0.FN (. is dot product) (not weve dropped the minus in our formulation)
      DO ik = mystart, myend
         dex1 = face(ik,1)
         dex2 = face(ik,2)
         dex3 = face(ik,3)
         A0(ik,:) = vertex(dex1,:)
         V0(ik,:)  = vertex(dex3,:)-vertex(dex1,:)
         V1(ik,:)  = vertex(dex2,:)-vertex(dex1,:)
         FN(ik,1) = (V1(ik,2)*V0(ik,3))-(V1(ik,3)*V0(ik,2))
         FN(ik,2) = (V1(ik,3)*V0(ik,1))-(V1(ik,1)*V0(ik,3))
         FN(ik,3) = (V1(ik,1)*V0(ik,2))-(V1(ik,2)*V0(ik,1))
      END DO
      ! allocate memory for information about mesh triangles 
      IF (lcomm) THEN
         CALL MPI_BARRIER(shar_comm, istat)
         CALL mpialloc_1d_dbl(DOT00,nface,shar_rank,0,shar_comm,win_dot00)
         CALL mpialloc_1d_dbl(DOT01,nface,shar_rank,0,shar_comm,win_dot01)
         CALL mpialloc_1d_dbl(DOT11,nface,shar_rank,0,shar_comm,win_dot11)
         CALL mpialloc_1d_dbl(invDenom,nface,shar_rank,0,shar_comm,win_invDenom)
         CALL mpialloc_1d_dbl(d,nface,shar_rank,0,shar_comm,win_d)
         CALL mpialloc_1d_int(ihit_array,nface,shar_rank,0,shar_comm,win_ihit)
         mydelta = CEILING(REAL(nface) / REAL(shar_size))
         mystart = 1 + shar_rank*mydelta
         myend   = mystart + mydelta
         IF (myend > nface) myend=nface
      ELSE
         ALLOCATE(DOT00(nface), DOT01(nface),&
                  DOT11(nface), invDenom(nface),&
                  STAT=istat)
         ALLOCATE(d(nface),STAT=istat)
         ALLOCATE(ihit_array(nface),STAT=istat)
         mystart = 1; myend = nface
      END IF

      ! if no error, calculate information about mesh triangles
      IF (istat/=0) RETURN
      DO ik = mystart, myend
         ihit_array(ik) = 0
         DOT00(ik) = V0(ik,1)*V0(ik,1) + V0(ik,2)*V0(ik,2) + V0(ik,3)*V0(ik,3)
         DOT01(ik) = V0(ik,1)*V1(ik,1) + V0(ik,2)*V1(ik,2) + V0(ik,3)*V1(ik,3)
         DOT11(ik) = V1(ik,1)*V1(ik,1) + V1(ik,2)*V1(ik,2) + V1(ik,3)*V1(ik,3)
         d(ik)     = FN(ik,1)*A0(ik,1) + FN(ik,2)*A0(ik,2) + FN(ik,3)*A0(ik,3)
         invDenom(ik) = one / (DOT00(ik)*DOT11(ik) - DOT01(ik)*DOT01(ik))
      END DO
      ! Set very low number of faces/block for meshes that won't be saved to .dat anyway
      CALL SET_NFACE(10)
      IF (lcomm) THEN 
#if defined(MPI_OPT)
         CALL MPI_BARRIER(comm,istat)
#endif
         ! Create an accelerated wall manually
         CALL ACCELERATE_WALL(istat, comm, shar_comm)
      ELSE
         CALL ACCELERATE_WALL(istat)
      END IF

      ! sync MPI
      IF (lcomm) THEN
#if defined(MPI_OPT)
         CALL MPI_BARRIER(shar_comm, istat)
         CALL MPI_COMM_FREE(shar_comm, istat)
#endif
      END IF
      ! set wall as loaded and return
      lwall_loaded = .true.
      RETURN
      END SUBROUTINE wall_load_mn

      SUBROUTINE wall_load_seg(npts,rseg,zseg,nphi,istat,verb,comm)
      !-----------------------------------------------------------------------
      ! wall_load_seg: Creates wall from segments
      !-----------------------------------------------------------------------
      ! param[in]: npts. Number of points in r and Z direction
      ! param[in]: rseg. Segments in r direction (with npts number of points)
      ! param[in]: zseg. Segmetns in z direction (with npts number of points)
      ! param[in]: nphi. Number of points in phi direction
      ! param[in, out]: istat. Integer that shows error if != 0
      ! param[in]: verb: Verbosity. True or false
      ! param[in, out]: comm. MPI communicator, handles shared memory
      !-----------------------------------------------------------------------
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      INTEGER, INTENT(in) :: npts, nphi
      DOUBLE PRECISION, INTENT(in) :: rseg(npts)
      DOUBLE PRECISION, INTENT(in) :: Zseg(npts)
      INTEGER, INTENT(inout)       :: istat
      LOGICAL, INTENT(in), OPTIONAL :: verb
      INTEGER, INTENT(inout), OPTIONAL :: comm
      INTEGER :: shar_comm
      LOGICAL :: shared, lcomm
      INTEGER :: nseg, ij, ik, il, im
      DOUBLE PRECISION :: dphi
      INTEGER :: dex1, dex2, dex3


      lcomm = .FALSE.
      IF (PRESENT(comm)) lcomm = .TRUE.
      IF (PRESENT(verb)) lverb = verb
      IF (lverb) WRITE(6,*) '-----  Creating wall mesh  -----'
      
      shar_rank = 0; shar_size = 1;
      dphi = 8.0D+00 * ATAN(one)/nphi
      ! initialize MPI
      IF (lcomm) THEN
#if defined(MPI_OPT)
         CALL MPI_COMM_SPLIT_TYPE(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, shar_comm, istat)
         CALL MPI_COMM_RANK( shar_comm, shar_rank, istat )
         CALL MPI_COMM_SIZE( shar_comm, shar_size, istat)
#endif
      END IF
      ! create info usually read from file manually
      machine_string = '          SEGMENTS'
      date = '      TODAY'
      nseg = npts-1
      nvertex = nseg * 4 * nphi
      nface   = nseg * 2 * nphi
      ! allocate shared memory with of mesh
      ! calculate which part each node has to calculate using mydelta
      IF (lcomm) THEN
         CALL mpialloc_2d_dbl(vertex,nvertex,3,shar_rank,0,shar_comm,win_vertex)
         CALL mpialloc_2d_int(face,nface,3,shar_rank,0,shar_comm,win_face)
         shared = .true.
         mydelta = CEILING(REAL(nface) / REAL(shar_size))
         mystart = 1 + shar_rank*mydelta
         myend   = mystart + mydelta
         IF (myend > nface) myend=nface
      ELSE
         ALLOCATE(vertex(nvertex,3),face(nface,3),STAT=istat)
         shared = .false.
         mystart = 1; myend=nface
      END IF

      ! if no error, and shared rank is zero, create triangles from input info
      IF (istat/=0) RETURN
      IF (shar_rank == 0) THEN
         il = 1; im = 1
         DO ik = 1, nphi
            DO ij = 1, nseg
               ! create the 4 points
               vertex(il,1) = rseg(ij)*cos(dphi*(ik-1))
               vertex(il,2) = rseg(ij)*sin(dphi*(ik-1))
               vertex(il,3) = zseg(ij)
               vertex(il+1,1) = rseg(ij+1)*cos(dphi*(ik-1))
               vertex(il+1,2) = rseg(ij+1)*sin(dphi*(ik-1))
               vertex(il+1,3) = zseg(ij+1)
               vertex(il+2,1) = rseg(ij)*cos(dphi*ik)
               vertex(il+2,2) = rseg(ij)*sin(dphi*ik)
               vertex(il+2,3) = zseg(ij)
               vertex(il+3,1) = rseg(ij+1)*cos(dphi*ik)
               vertex(il+3,2) = rseg(ij+1)*sin(dphi*ik)
               vertex(il+3,3) = zseg(ij+1)
               ! Create the 2 triangles
               face(im,1) = il
               face(im,2) = il+2
               face(im,3) = il+1
               face(im+1,1) = il+2
               face(im+1,2) = il+3
               face(im+1,3) = il+1
               ! Adjust index
               im = im + 2
               il = il + 4
            END DO 
         END DO
      END IF
      ! if using MPI, wait here
      ! allocate memory for information about the mesh
#if defined(MPI_OPT)
      IF (lcomm) CALL MPI_BARRIER(comm,istat)
#endif
      IF (istat/=0) RETURN
      IF (lcomm) THEN
         CALL mpialloc_2d_dbl(A0,nface,3,shar_rank,0,shar_comm,win_a0)
         CALL mpialloc_2d_dbl(V0,nface,3,shar_rank,0,shar_comm,win_v0)
         CALL mpialloc_2d_dbl(V1,nface,3,shar_rank,0,shar_comm,win_v1)
         CALL mpialloc_2d_dbl(FN,nface,3,shar_rank,0,shar_comm,win_fn)
         mydelta = CEILING(REAL(nface) / REAL(shar_size))
         mystart = 1 + shar_rank*mydelta
         myend   = mystart + mydelta
         IF (myend > nface) myend=nface
      ELSE
         ALLOCATE(A0(nface,3),V0(nface,3),V1(nface,3),&
                  FN(nface,3),STAT=istat)
         mystart = 1; myend = nface
      END IF

      IF (istat/=0) RETURN
      ! Precalculate information about mesh
      ! Calculate the face normal
      ! V  = Vertex1-Vertex0
      ! W  = Vertex2-Vertex0
      ! FN = VxW/|VxW|
      ! d  = -Vertex0.FN (. is dot product) (note weve absorbed the negative)
      DO ik = mystart, myend
         dex1 = face(ik,1)
         dex2 = face(ik,2)
         dex3 = face(ik,3)
         A0(ik,:) = vertex(dex1,:)
         V0(ik,:)  = vertex(dex3,:)-vertex(dex1,:)
         V1(ik,:)  = vertex(dex2,:)-vertex(dex1,:)
         FN(ik,1) = (V1(ik,2)*V0(ik,3))-(V1(ik,3)*V0(ik,2))
         FN(ik,2) = (V1(ik,3)*V0(ik,1))-(V1(ik,1)*V0(ik,3))
         FN(ik,3) = (V1(ik,1)*V0(ik,2))-(V1(ik,2)*V0(ik,1))
      END DO
#if defined(MPI_OPT)
      IF (lcomm) CALL MPI_BARRIER(comm,istat)
#endif
      ! Check for zero area
      IF (ANY(SUM(FN*FN,DIM=2)==zero)) THEN
         istat=-327
         RETURN
      END IF
      ! allocate memory for information about mesh triangles 
      IF (lcomm) THEN
         CALL mpialloc_1d_dbl(DOT00,nface,shar_rank,0,shar_comm,win_dot00)
         CALL mpialloc_1d_dbl(DOT01,nface,shar_rank,0,shar_comm,win_dot01)
         CALL mpialloc_1d_dbl(DOT11,nface,shar_rank,0,shar_comm,win_dot11)
         CALL mpialloc_1d_dbl(invDenom,nface,shar_rank,0,shar_comm,win_invDenom)
         CALL mpialloc_1d_dbl(d,nface,shar_rank,0,shar_comm,win_d)
         CALL mpialloc_1d_int(ihit_array,nface,shar_rank,0,shar_comm,win_ihit)
         mydelta = CEILING(REAL(nface) / REAL(shar_size))
         mystart = 1 + shar_rank*mydelta
         myend   = mystart + mydelta
         IF (myend > nface) myend=nface
      ELSE
         ALLOCATE(DOT00(nface), DOT01(nface),&
                  DOT11(nface), invDenom(nface),&
                  STAT=istat)
         ALLOCATE(d(nface),STAT=istat)
         ALLOCATE(ihit_array(nface),STAT=istat)
         mystart = 1; myend = nface
      END IF

      ! if no error, calculate information about mesh triangles
      IF (istat/=0) RETURN
      DO ik = mystart, myend
         ihit_array(ik) = 0
         DOT00(ik) = V0(ik,1)*V0(ik,1) + V0(ik,2)*V0(ik,2) + V0(ik,3)*V0(ik,3)
         DOT01(ik) = V0(ik,1)*V1(ik,1) + V0(ik,2)*V1(ik,2) + V0(ik,3)*V1(ik,3)
         DOT11(ik) = V1(ik,1)*V1(ik,1) + V1(ik,2)*V1(ik,2) + V1(ik,3)*V1(ik,3)
         d(ik)     = FN(ik,1)*A0(ik,1) + FN(ik,2)*A0(ik,2) + FN(ik,3)*A0(ik,3)
         invDenom(ik) = one / (DOT00(ik)*DOT11(ik) - DOT01(ik)*DOT01(ik))
      END DO
      ! Set very low number of faces/block for meshes that won't be saved to .dat anyway
      CALL SET_NFACE(10)
      IF (lcomm) THEN 
#if defined(MPI_OPT)
         CALL MPI_BARRIER(comm,istat)
#endif
         ! Create an accelerated wall manually
         CALL ACCELERATE_WALL(istat, comm, shar_comm)
      ELSE
         CALL ACCELERATE_WALL(istat)
      END IF

      ! sync MPI
      IF (lcomm) THEN
#if defined(MPI_OPT)
         CALL MPI_BARRIER(shar_comm, istat)
         CALL MPI_COMM_FREE(shar_comm, istat)
#endif
      END IF
      ! set wall as loaded and return
      lwall_loaded = .true.
      RETURN
      END SUBROUTINE wall_load_seg

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Wall info
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE wall_dump(filename,istat)
      !-----------------------------------------------------------------------
      ! wall_dump: Dumps triangulation data
      !-----------------------------------------------------------------------
      ! param[in]: filename. The file name to load in
      ! param[in, out]: istat. Integer that shows error if != 0
      !-----------------------------------------------------------------------
      CHARACTER(LEN=*), INTENT(in) :: filename
      INTEGER, INTENT(inout)       :: istat
      INTEGER :: iunit, ik
      ! open file
      CALL safe_open(iunit,istat,'wall_dump.'//TRIM(filename),'unknown','formatted')
      ! for every face, output info
      DO ik = 1, nface
         WRITE(iunit,'(13(ES20.10))')  A0(ik,1), A0(ik,2), A0(ik,3), &
                                       FN(ik,1), FN(ik,2), FN(ik,3),&
                                       V0(ik,1), V0(ik,2), V0(ik,3),&
                                       V1(ik,1), V1(ik,2), V1(ik,3), &
                                       d(ik)
      END DO
      WRITE(iunit,*) '#  V0(ik,1), V0(ik,2), V0(ik,3), FN(ik,1), FN(ik,2), FN(ik,3),',&
                        'V(ik,1), V(ik,2), V(ik,3), W(ik,1), W(ik,2), W(ik,3), d(ik)'
      ! close file
      CLOSE(iunit)
      RETURN
      END SUBROUTINE wall_dump

      
      SUBROUTINE wall_info(iunit)
      !-----------------------------------------------------------------------
      ! wall_info: Prints wall info
      !-----------------------------------------------------------------------
      ! param[in]: iunit. Location to print information about wall
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER, INTENT(in)    :: iunit
      INTEGER                :: i, c, max_c
      DOUBLE PRECISION       :: mean_c
      WRITE(iunit,'(A)')         ' -----  Vessel Information  -----'
      WRITE(iunit,'(3X,A,A)')    'Wall Name : ',TRIM(machine_string(10:))
      WRITE(iunit,'(3X,A,A)')    'Date      : ',TRIM(date(6:))
      WRITE(iunit,'(3X,A,I7)')   'Faces     : ',nface
      c = 0; max_c = 0
      IF (wall%nblocks > 0) THEN
         WRITE(iunit,'(3X,A,I7)')    'Blocks    : ',wall%nblocks
         DO i=1,wall%nblocks
            IF (wall%blocks(i)%nfaces > max_c) max_c = wall%blocks(i)%nfaces
            c = c + wall%blocks(i)%nfaces
         END DO
         mean_c = c / wall%nblocks
         WRITE(iunit,'(3X,A,F9.2)')    'Mean faces per block: ', mean_c
         WRITE(iunit,'(3X,A,I7)')   'Highest faces per block: ', max_c   
      END IF
          
      RETURN
      END SUBROUTINE wall_info

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Wall collide
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE collide_float(x0,y0,z0,x1,y1,z1,xw,yw,zw,lhit)
      !-----------------------------------------------------------------------
      ! collide_float: Implementation of collide for floating point values
      !-----------------------------------------------------------------------
      ! param[in]: x0. x-location of first point of the line segment to check
      ! param[in]: y0. y-location of first point of the line segment to check
      ! param[in]: z0. z-location of first point of the line segment to check
      ! param[in]: x1. x-location of second point of the line segment to check
      ! param[in]: y1. y-location of second point of the line segment to check
      ! param[in]: z1. z-location of second point of the line segment to check
      ! param[out]: xw. x-location of hit (if hit has been found)
      ! param[out]: yw. y-location of hit (if hit has been found)
      ! param[out]: zw. z-location of hit (if hit has been found)
      ! param[out]: lhit. Logical that shows if hit has been found or not
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      REAL, INTENT(in) :: x0, y0, z0, x1, y1, z1
      REAL, INTENT(out) :: xw, yw, zw
      LOGICAL, INTENT(out) :: lhit
      DOUBLE PRECISION :: x0d, y0d, z0d, x1d, y1d, z1d
      DOUBLE PRECISION :: xwd, ywd, zwd
      LOGICAL          :: lhit2
      ! function simply converts from floating point to double to help compiler
      xw=zero; yw=zero; zw=zero; lhit=.FALSE.
      x0d=x0; y0d=y0; z0d=z0
      x1d=x1; y1d=y1; z1d=z1
      CALL collide_double(x0d,y0d,z0d,x1d,y1d,z1d,xwd,ywd,zwd,lhit2)
      xw=xwd; yw=ywd; zw=zwd; lhit=lhit2
      RETURN
      END SUBROUTINE collide_float

      SUBROUTINE collide_double(x0,y0,z0,x1,y1,z1,xw,yw,zw,lhit)
      !-----------------------------------------------------------------------
      ! collide_double: Implementation of collide for double precision values
      !-----------------------------------------------------------------------
      ! param[in]: x0. x-location of first point of the line segment to check
      ! param[in]: y0. y-location of first point of the line segment to check
      ! param[in]: z0. z-location of first point of the line segment to check
      ! param[in]: x1. x-location of second point of the line segment to check
      ! param[in]: y1. y-location of second point of the line segment to check
      ! param[in]: z1. z-location of second point of the line segment to check
      ! param[out]: xw. x-location of hit (if hit has been found)
      ! param[out]: yw. y-location of hit (if hit has been found)
      ! param[out]: zw. z-location of hit (if hit has been found)
      ! param[out]: lhit. Logical that shows if hit has been found or not
      !-----------------------------------------------------------------------
      IMPLICIT NONE
      ! Line positions
      DOUBLE PRECISION, INTENT(in) :: x0, y0, z0, x1, y1, z1
      ! Hit positions and logical if hit was found
      DOUBLE PRECISION, INTENT(out) :: xw, yw, zw
      LOGICAL, INTENT(out) :: lhit
      ! Loop integers, block integers
      INTEGER :: ik, i, k1, k2, b_found
      ! In which block the ray is in x/y/z
      INTEGER :: br(3)
      ! Line direction/position vectors
      DOUBLE PRECISION :: dr(3), r0(3)
      ! Hit calculation
      DOUBLE PRECISION :: V2x, V2y, V2z, DOT02l, DOT12l, tloc, tmin, alphal, betal
      ! Whether or not a step is in negative or positive direction in x/y/z
      INTEGER :: step(3)
      ! Time until block step
      DOUBLE PRECISION :: tDelta(3), tcomp(3)
      ! Whether or not out of grid in x/y/z
      LOGICAL :: outlow(3), outhigh(3)
      xw=zero; yw=zero; zw=zero; lhit=.FALSE.
      ik_min = zero
      tmin = one + epsilon
      ! Define DR
      dr(1) = x1-x0
      dr(2) = y1-y0
      dr(3) = z1-z0
      ! initialize
      r0(1) = x0
      r0(2) = y0
      r0(3) = z0
      DO i=1,3
         outlow(i) = .false.
         outhigh(i) = .false.
         tcomp(i) = zero
         ! check direction line to determine in which direction to step through the blocks
         IF (dr(i) < zero) THEN
            step(i) = -1
         ELSE
            step(i) = 1
         END IF
      END DO

      k1 = 1; k2 = wall%nblocks
      b_found = -1

      ! Do integer division to find which block the line is in
      DO i=1,3
         br(i) = INT((r0(i) - wall%rmin(i)) / wall%stepsize) + 1
      END DO
                  
      ! Check if outside the grid anywhere
      IF (ANY(br < 1) .or. ANY(br > wall%br)) THEN

         ! Check closest block in grid
         ! Also keep track which dimension is out of bounds
         DO i=1,3
            IF (br(i) < 1) THEN 
               br(i) = 1
               outlow(i) = .TRUE.
            ELSE IF (br(i) > wall%br(i)) THEN
               br(i) = wall%br(i)
               outhigh(i) = .TRUE.
            END IF
         END DO
         
         ! Find correct block
         b_found = br(3) + (br(2) - 1) * wall%step(2) + (br(1) - 1) * wall%step(1)
         b = wall%blocks(b_found)

         ! Check where grid is entered in nearest block
         ! Check if it went into the nearest block within the current distance for each dimension
         ! Only check for dimension that were out of bounds
         DO i=1,3
            IF (outlow(i) .and. step(i) .eq. 1) THEN
               tDelta(i) = (b%rmin(i) - r0(i)) / dr(i)
               IF (tDelta(i) < tmin) THEN
                  tcomp(i) = tDelta(i) + epsilon
               END IF
            END IF

            IF (outhigh(i) .and. step(i) .eq. -1) THEN
               tDelta(i) = (b%rmax(i) - r0(i)) / dr(i)
               IF (tDelta(i) < tmin) THEN
                  tcomp(i) = tDelta(i) + epsilon
               END IF
            END IF
         END DO
         
         ! If it did not enter any dimension in time, set block to -1
         IF (ANY((tcomp == zero) .and. (outlow .or. outhigh))) THEN
            b_found = -1                         
         END IF
      ELSE
         ! If not outside grid, set b_found
         b_found = br(3) + (br(2) - 1) * wall%step(2) + (br(1) - 1) * wall%step(1)
      END IF

      ! Traverse blocks
      DO WHILE (.true.)
         ! Reset hit info
         xw=zero; yw=zero; zw=zero; lhit=.FALSE.
         ik_min = zero
         tmin = one + epsilon
         ! If outside block, exit
         IF (b_found > wall%nblocks .or. b_found < 1) EXIT
         b = wall%blocks(b_found)

         k1 = 1; k2 = b%nfaces
         ! Check every triangle
         ! Based on Badouel's algorithm
         ! Source: https://graphics.stanford.edu/courses/cs348b-98/gg/intersect.html
         DO i = k1,k2
            ! get ik by reading face number in block
            ik = b%face(i)
            ! calculate whether or not this line segment ever hits the plane of the triangle
            alphal = FN(ik,1)*dr(1) + FN(ik,2)*dr(2) + FN(ik,3)*dr(3)
            betal = FN(ik,1)*r0(1) + FN(ik,2)*r0(2) + FN(ik,3)*r0(3)
            ! tloc indicated when hit. If hit between r0 and r1, tloc between 0 and 1
            tloc = (d(ik)-betal)/alphal
            IF (tloc > one) CYCLE
            IF (tloc <= zero) CYCLE
            ! If the line segment hits the plane of the triangle
            ! calculate if it actually hits on the triangle
            V2x = x0 + tloc*dr(1) - A0(ik,1)
            V2y = y0 + tloc*dr(2) - A0(ik,2)
            V2z = z0 + tloc*dr(3) - A0(ik,3)
            DOT02l = V0(ik,1)*V2x + V0(ik,2)*V2y + V0(ik,3)*V2z
            DOT12l = V1(ik,1)*V2x + V1(ik,2)*V2y + V1(ik,3)*V2z
            alphal = (DOT11(ik)*DOT02l-DOT01(ik)*DOT12l)*invDenom(ik)
            betal  = (DOT00(ik)*DOT12l-DOT01(ik)*DOT02l)*invDenom(ik)
            ! In that case, these should be false
            IF ((alphal < -epsilon) .or. (betal < -epsilon) .or. (alphal+betal > one + epsilon)) CYCLE
            ! else check if this was the closest hit, and then store
            IF (tloc < tmin) THEN
               ik_min = ik
               tmin = tloc
            END IF
         END DO
         
         ! Check when the line will leave the block
         DO i=1,3
            IF (step(i) .eq. 1) THEN
               tDelta(i) = (b%rmax(i) - r0(i)) / dr(i)
            ELSE
               tDelta(i) = (b%rmin(i) - r0(i)) / dr(i)
            END IF
         END DO

         ! Check if leaves block before hits. 
         ! Compare with tcomp to make sure that ray always moves in positive time direction 
         IF (ANY(tDelta < tmin .and. tDelta > tcomp)) THEN
            ! If true, find which direction is smallest time until leave block
            IF (tDelta(1) < tDelta(2)) THEN
               IF (tDelta(1) < tDelta(3)) THEN
                  i = 1
               ELSE 
                  i = 3
               END IF
            ELSE
               IF (tDelta(2) < tDelta(3)) THEN
                  i = 2
               ELSE
                  i = 3
               END IF
            END IF
            ! Having found exit direction, step in that direction
            br(i) = br(i) + step(i)
            ! If leave grid, exit
            IF (br(i) > wall%br(i) .or. br(i) < 1) EXIT    
            b_found = b_found + wall%step(i) * step(i)
            tcomp(i) = tDelta(i) + epsilon
         ELSE 
            ! If none smaller than hit time, exit loop. Hit found.
            EXIT
         END IF
      END DO

      ! if any index stored, hit was found, calculate location and increment ihit_array
      IF (ik_min > zero) THEN
         lhit = .TRUE.
         xw   = x0 + tmin*dr(1)
         yw   = y0 + tmin*dr(2)
         zw   = z0 + tmin*dr(3)
         ihit_array(ik_min) = ihit_array(ik_min) + 1
      END IF
      RETURN
      END SUBROUTINE collide_double

      SUBROUTINE uncount_wall_hit
      !-----------------------------------------------------------------------
      ! uncount_wall_hit: Reduces ihit_array at last found hit location with one
      !-----------------------------------------------------------------------
         IMPLICIT NONE
         ihit_array(ik_min) = ihit_array(ik_min) - 1
      END SUBROUTINE

      INTEGER FUNCTION get_wall_ik()
      !-----------------------------------------------------------------------
      ! get_wall_ik: Gets index of last hit location
      !-----------------------------------------------------------------------
      ! return[integer]: get_wall_ik. Last hit index 
      !-----------------------------------------------------------------------
         IMPLICIT NONE
         get_wall_ik = ik_min
         RETURN
      END FUNCTION

      DOUBLE PRECISION FUNCTION get_wall_area(ik)
      !-----------------------------------------------------------------------
      ! get_wall_ik: Gets index of last hit location
      !-----------------------------------------------------------------------
      ! param[in]: ik. Index of wall location to check
      ! return[double]: get_wall_area. Area of wall location checked 
      !-----------------------------------------------------------------------
         IMPLICIT NONE
         INTEGER, INTENT(in) :: ik
         get_wall_area = 0.5*SQRT(SUM(FN(ik,:)*FN(ik,:)))
         RETURN
      END FUNCTION

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Wall acceleration creation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE ACCELERATE_WALL(istat, comm, shar_comm)
      !-----------------------------------------------------------------------
      ! ACCELERATE_WALL: Accelerates an unaccelerated wall
      !-----------------------------------------------------------------------
      ! param[in, out]: istat. Integer that shows error if != 0
      ! param[in, out]: comm. MPI communicator, handles communication between nodes
      ! param[in, out]: shar_comm. Shared MPI communicator, handles shared memory
      !-----------------------------------------------------------------------
         INTEGER, OPTIONAL :: istat
         INTEGER, INTENT(inout), OPTIONAL :: comm, shar_comm
   
         DOUBLE PRECISION :: size
         
         IF (lverb) WRITE(6, *) '-----  Creating accelerated wall  -----'
         IF (lverb) CALL WALL_INFO(6)
   
         CALL GET_BLOCK_SIZE(size)

         IF (lverb) WRITE(6, *) 'Creating blocks. Chosen size:', size
   
         CALL CREATE_BLOCKS(size, istat)
   
#if defined(MPI_OPT)
         IF (PRESENT(comm)) CALL MPI_BARRIER(shar_comm,istat)
#endif

         IF (lverb) WRITE(6, *) 'Filling blocks. Number of blocks: ', wall%nblocks
   
#if defined(MPI_OPT)
         IF (PRESENT(comm)) THEN 
            CALL FILL_BLOCKS(istat, comm, shar_comm)
            CALL MPI_BARRIER(shar_comm,istat)
         ELSE
#endif
            CALL FILL_BLOCKS(istat)
#if defined(MPI_OPT)
         END IF
#endif
   
         IF (lverb) WRITE(6, *) '-----  Done creating accelerated wall  -----'
      END SUBROUTINE ACCELERATE_WALL
   
      SUBROUTINE GET_BLOCK_SIZE(size)
      !-----------------------------------------------------------------------
      ! GET_BLOCK_SIZE: Determines the block size of the newly accelerated wall
      !-----------------------------------------------------------------------
      ! param[in, out]: size. Size of the blocks
      !-----------------------------------------------------------------------
         DOUBLE PRECISION, INTENT(out) :: size
   
         DOUBLE PRECISION :: wall_size(3), square_wall_size
         INTEGER :: i, nblocks_desired
   
         square_wall_size = MAXVAL(MAXVAL(vertex, DIM=1) - MINVAL(vertex, DIM=1)) + 0.1
         wall_size = MAXVAL(vertex, DIM=1) - MINVAL(vertex, DIM=1)
         nblocks_desired = nface / nfaces_per_block + 1
   
         DO i=1,100000
            size = square_wall_size / i
            IF (CEILING(wall_size(1)/size) * CEILING(wall_size(2)/size) * CEILING(wall_size(3)/size)> nblocks_desired) EXIT
         END DO      
      END SUBROUTINE GET_BLOCK_SIZE
   
      SUBROUTINE SET_NFACE(new_nfaces_per_block)
      !-----------------------------------------------------------------------
      ! SET_NFACE: Sets the number of faces desired per block to a new value
      !-----------------------------------------------------------------------
      ! param[in]: new_nfaces_per_block. New value desired
      !-----------------------------------------------------------------------
      
         INTEGER, INTENT(in) :: new_nfaces_per_block
   
         nfaces_per_block = new_nfaces_per_block
   
      END SUBROUTINE
   
      SUBROUTINE CREATE_BLOCKS(size, istat)
      !-----------------------------------------------------------------------
      ! CREATE_BLOCKS: Creates a grid of blocks
      !-----------------------------------------------------------------------
      ! param[in]: size. Size of the blocks
      ! param[in, out]: istat. Integer that shows error if != 0
      !-----------------------------------------------------------------------
         DOUBLE PRECISION, INTENT(in) :: size
         INTEGER, OPTIONAL :: istat
   
         DOUBLE PRECISION :: wall_size(3), rmin(3), buffer
         DOUBLE PRECISION :: tmp
         INTEGER :: i, j, xi, yi, zi
         INTEGER :: nblocks(3)
   
         DOUBLE PRECISION, POINTER :: xs(:), ys(:), zs(:)
   
         ! Find size of mesh
         rmin = MINVAL(vertex, DIM=1)
         wall_size = MAXVAL(vertex, DIM=1) - rmin

         ! For each dimension, find bounds of blocks in a very ugly way
         nblocks(1) = INT(wall_size(1)/size) + 1
         nblocks(2) = INT(wall_size(2)/size) + 1
         nblocks(3) = INT(wall_size(3)/size) + 1
         ALLOCATE(xs(nblocks(1)),STAT=istat)
         ALLOCATE(ys(nblocks(2)),STAT=istat)
         ALLOCATE(zs(nblocks(3)),STAT=istat)
         ! XS
         buffer = nblocks(1) * size - wall_size(1)
         DO j = 1, nblocks(1)
            xs(j) = rmin(1) - buffer / 2 + (j - 1) * size
         END DO
         ! YS
         buffer = nblocks(2) * size - wall_size(2)
         DO j = 1, nblocks(2)
            ys(j) = rmin(2) - buffer / 2 + (j - 1) * size
         END DO
         ! ZS
         buffer = nblocks(3) * size - wall_size(3)
         DO j = 1, nblocks(3)
            zs(j) = rmin(3) - buffer / 2 + (j - 1) * size
         END DO

         ! Set info about wall
         wall%rmin(1) = xs(1)
         wall%rmin(2) = ys(1)
         wall%rmin(3) = zs(1)

         wall%rmax(1) = xs(nblocks(1))
         wall%rmax(2) = ys(nblocks(2))
         wall%rmax(3) = zs(nblocks(3))
   
         wall%nblocks = nblocks(1) * nblocks(2) * nblocks(3)
         wall%stepsize = size
         wall%br = nblocks
         wall%step(3) = 1
         wall%step(2) = nblocks(3)
         wall%step(1) = nblocks(2) * nblocks(3)
   
         ! Set bound of each block
         ALLOCATE(wall%blocks(wall%nblocks),STAT=istat)

         DO i = 1, wall%nblocks
            zi = MOD(i-1,nblocks(3))+1
            yi = MOD(i-1,nblocks(2)*nblocks(3))
            yi = FLOOR(REAL(yi) / REAL(nblocks(3)))+1
            xi = CEILING(REAL(i) / REAL(nblocks(2)*nblocks(3)))
            wall%blocks(i)%rmin(1) = xs(xi)
            wall%blocks(i)%rmin(2) = ys(yi)
            wall%blocks(i)%rmin(3) = zs(zi)
                  
            wall%blocks(i)%rmax(1) = xs(xi) + size
            wall%blocks(i)%rmax(2) = ys(yi) + size
            wall%blocks(i)%rmax(3) = zs(zi) + size
                  
            wall%blocks(i)%nfaces = 0

         END DO

         DEALLOCATE(xs)
         DEALLOCATE(ys)
         DEALLOCATE(zs)
      END SUBROUTINE CREATE_BLOCKS
   
      SUBROUTINE FILL_BLOCKS(istat, comm, shar_comm)
      !-----------------------------------------------------------------------
      ! FILL_BLOCKS: Fills the grid of blocks with the correct faces
      !-----------------------------------------------------------------------
      ! param[in, out]: istat. Integer that shows error if != 0
      ! param[in, out]: comm. MPI communicator, handles communication between nodes
      ! param[in, out]: shar_comm. Shared MPI communicator, handles shared memory
      !-----------------------------------------------------------------------
         INTEGER, OPTIONAL :: istat
         INTEGER, INTENT(inout), OPTIONAL :: comm, shar_comm
         ! Whether or not shared memory is used
         LOGICAL :: shared, lcomm
         ! Integer for shared memory
         INTEGER :: win_counter_arr, win_mask_face

         ! Integer array that counters number of faces in each block 
         INTEGER, POINTER :: counter_arr(:)
         ! Integer array that determines how to split blocks over threads. Only used with MPI
         INTEGER, POINTER :: mysplit(:)
         ! Block bounds
         DOUBLE PRECISION :: rmin(3), rmax(3)
         ! Masks if vertex in block or face already in block
         LOGICAL, POINTER :: mask_face(:,:), mask(:)
         ! Loop integers and counter
         INTEGER :: i, j, k, counter

         ! By Sam
         INTEGER :: win_A1, win_A2
         DOUBLE PRECISION, POINTER :: A1(:,:), A2(:,:)
         
         ! Comm logit
         lcomm = .FALSE.
#if defined(MPI_OPT)
         IF (PRESENT(comm)) lcomm = .TRUE.
#endif

         ! Split blocks between threads
         ! Also allocate counter_arr
         IF (lcomm) THEN
            shared = .TRUE.
            mydelta = CEILING(REAL(wall%nblocks) / REAL(shar_size))
            mystart = 1 + shar_rank*mydelta
            myend   = mystart + mydelta
            IF (myend > wall%nblocks) myend=wall%nblocks
            CALL mpialloc_1d_int(counter_arr, wall%nblocks, shar_rank, 0, shar_comm, win_counter_arr)
            CALL mpialloc_2d_dbl(A1, nface, 3, shar_rank, 0, shar_comm, win_A1)
            CALL mpialloc_2d_dbl(A2, nface, 3, shar_rank, 0, shar_comm, win_A2)
         ELSE
            shared = .FALSE.
            mystart = 1; myend = wall%nblocks
            ALLOCATE(counter_arr(wall%nblocks))
            ALLOCATE(A1(nface,3),A2(nface,3))
         END IF

         counter_arr = 0
         ALLOCATE(mask_face(nface,3), STAT=istat)
         ALLOCATE(mask(nface), STAT=istat)

         ! Helper for vertices per face
         IF (shar_rank == 0) THEN
            A1 = V0+A0
            A2 = V1+A0
         END IF

#if defined(MPI_OPT)
         IF (lcomm) CALL MPI_BARRIER(shar_comm,istat)
#endif

         ! Find how many vertices in each block
         IF (shar_rank == 0 .and. ldebug) WRITE(6, *) 'Filling blocks: Finding the number of vertices in each block'
         DO i=mystart, myend
            ! Define bounds block
            rmin = wall%blocks(i)%rmin - epsilon
            rmax = wall%blocks(i)%rmax + epsilon

            ! Extend the domain to include the neighbors (SAL)
            !rmin = rmin - (rmax - rmin)*0.25
            !rmax = rmax + (rmax - rmin)*0.25

            ! Check vertices in block
            mask_face(:,:) = .FALSE.


            !mask_face(:,1) = (A0(:,1) < rmax(1) .and. A0(:,1) >= rmin(1) &
            !   .and. A0(:,2) < rmax(2) .and. A0(:,2) >= rmin(2) &
            !   .and. A0(:,3) < rmax(3) .and. A0(:,3) >= rmin(3))
            !mask_face(:,2) = (A1(:,1) < rmax(1) .and. A1(:,1) >= rmin(1) &
            !   .and. A1(:,2) < rmax(2) .and. A1(:,2) >= rmin(2) &
            !   .and. A1(:,3) < rmax(3) .and. A1(:,3) >= rmin(3))
            !mask_face(:,3) = (A2(:,1) < rmax(1) .and. A2(:,1) >= rmin(1) &
            !   .and. A2(:,2) < rmax(2) .and. A2(:,2) >= rmin(2) &
            !   .and. A2(:,3) < rmax(3) .and. A2(:,3) >= rmin(3))

            ! Could probably vectorize
            DO j = 1,nface
               CALL LINE_BOX_INTERSECTION(A0(j,1),A0(j,2),A0(j,3), &
                  A1(j,1),A1(j,2),A1(j,3),rmin(1),rmin(2),rmin(3),&
                  rmax(1),rmax(2),rmax(3),mask_face(j,1))
               CALL LINE_BOX_INTERSECTION(A1(j,1),A1(j,2),A1(j,3), &
                  A2(j,1),A2(j,2),A2(j,3),rmin(1),rmin(2),rmin(3),&
                  rmax(1),rmax(2),rmax(3),mask_face(j,2))
               CALL LINE_BOX_INTERSECTION(A2(j,1),A2(j,2),A2(j,3), &
                  A0(j,1),A0(j,2),A0(j,3),rmin(1),rmin(2),rmin(3),&
                  rmax(1),rmax(2),rmax(3),mask_face(j,3))
            END DO

            ! Count found faces
            mask = ANY(mask_face,2)
            counter_arr(i) = COUNT(mask)
         END DO

#if defined(MPI_OPT)
         IF (lcomm) CALL MPI_BARRIER(shar_comm,istat)
#endif

         ! Allocate face arrays for each block. Each thread executes this
         ! Also sets size of block for each thread
         IF (shar_rank == 0 .and. ldebug) WRITE(6, *) 'Filling blocks: Allocating each block'
         DO i=1, wall%nblocks
            IF (counter_arr(i) > 0) THEN
               IF (lcomm) THEN 
                  CALL mpialloc_1d_int(wall%blocks(i)%face, counter_arr(i), shar_rank, 0, shar_comm, wall%blocks(i)%win_face)
               ELSE
                  ALLOCATE(wall%blocks(i)%face(counter_arr(i)), STAT=istat)
               END IF
            END IF
            wall%blocks(i)%nfaces = counter_arr(i)
            wall%blocks(i)%isshared = shared
         END DO

         ! Efficiently divide blocks over threads
         IF (shar_rank == 0 .and. ldebug) WRITE(6, *) 'Filling blocks: Dividing work'
         IF (lcomm) THEN
            ALLOCATE(mysplit(shar_size + 1))
            DO j=0,100
               mydelta = CEILING(SUM(counter_arr) / REAL(shar_size) / REAL(100 + j) * 100D+0)
               mysplit = 0
               counter = 0
               k=2
               DO i=1, wall%nblocks
                  counter = counter + counter_arr(i)
                  IF (counter > mydelta) THEN
                     mysplit(k) = i
                     k = k + 1
                     IF (k > shar_size + 1) k = shar_size + 1
                     counter = 0
                  END IF
               END DO
               mysplit(k) = wall%nblocks
               IF (k .eq. shar_size + 1) EXIT
            END DO
            mystart = 1 + mysplit(shar_rank + 1)
            myend   = mysplit(shar_rank + 2)
            IF (mystart > wall%nblocks .OR. myend .eq. 0) myend=mystart
            DEALLOCATE(mysplit)
         ELSE
            mystart = 1; myend = wall%nblocks
         END IF
#if defined(MPI_OPT)
         IF (lcomm) CALL MPI_BARRIER(shar_comm,istat)
#endif

         ! Actually assign faces to blocks
         IF (shar_rank == 0 .and. ldebug) WRITE(6, *) 'Filling blocks: Putting faces in each block'
         DO i=mystart, myend
            ! Skip if for some reason out of bounds
            IF (i > wall%nblocks) EXIT
            ! Only add if actually faces in the block
            IF (counter_arr(i) > 0) THEN
               mask_face = .FALSE.
               ! Define bounds block
               rmin = wall%blocks(i)%rmin - epsilon
               rmax = wall%blocks(i)%rmax + epsilon

               ! Extend the domain to include the neighbors (SAL)
               !rmin = rmin - (rmax - rmin)*0.25
               !rmax = rmax + (rmax - rmin)*0.25
      
               ! Check which vertices are in block         
               !mask_face(:,1) = (A0(:,1) < rmax(1) .and. A0(:,1) >= rmin(1) &
               !   .and. A0(:,2) < rmax(2) .and. A0(:,2) >= rmin(2) &
               !   .and. A0(:,3) < rmax(3) .and. A0(:,3) >= rmin(3))
               !mask_face(:,2) = (A1(:,1) < rmax(1) .and. A1(:,1) >= rmin(1) &
               !   .and. A1(:,2) < rmax(2) .and. A1(:,2) >= rmin(2) &
               !   .and. A1(:,3) < rmax(3) .and. A1(:,3) >= rmin(3))
               !mask_face(:,3) = (A2(:,1) < rmax(1) .and. A2(:,1) >= rmin(1) &
               !   .and. A2(:,2) < rmax(2) .and. A2(:,2) >= rmin(2) &
               !   .and. A2(:,3) < rmax(3) .and. A2(:,3) >= rmin(3))

               ! Could probably vectorize
               DO j = 1,nface
                  CALL LINE_BOX_INTERSECTION(A0(j,1),A0(j,2),A0(j,3), &
                     A1(j,1),A1(j,2),A1(j,3),rmin(1),rmin(2),rmin(3),&
                     rmax(1),rmax(2),rmax(3),mask_face(j,1))
                  CALL LINE_BOX_INTERSECTION(A1(j,1),A1(j,2),A1(j,3), &
                     A2(j,1),A2(j,2),A2(j,3),rmin(1),rmin(2),rmin(3),&
                     rmax(1),rmax(2),rmax(3),mask_face(j,2))
                  CALL LINE_BOX_INTERSECTION(A2(j,1),A2(j,2),A2(j,3), &
                     A0(j,1),A0(j,2),A0(j,3),rmin(1),rmin(2),rmin(3),&
                     rmax(1),rmax(2),rmax(3),mask_face(j,3))
               END DO

               ! Mask
               mask = ANY(mask_face,2)
               
               ! Add faces to face list of wall block
               counter = 0
               DO j=1,nface
                  IF (mask(j)) THEN
                     counter = counter + 1
                     wall%blocks(i)%face(counter) = j
                  END IF
                  ! Break if all faces have been found
                  IF (counter == counter_arr(i)) EXIT
               END DO
            END IF           
         END DO
         
         ! Cleanup
#if defined(MPI_OPT)
         IF (lcomm) CALL MPI_BARRIER(shar_comm,istat)
#endif
         IF (shar_rank == 0 .and. ldebug) WRITE(6, *) 'Filling blocks: Starting cleanup fill blocks'
         DEALLOCATE(mask_face, mask)
         CALL free_mpi_array(win_counter_arr, counter_arr, shared)
         CALL free_mpi_array(win_A1, A1, shared)
         CALL free_mpi_array(win_A2, A2, shared)
         END SUBROUTINE FILL_BLOCKS
   
      SUBROUTINE WRITE_WALL(filename, istat)
      !-----------------------------------------------------------------------
      ! WRITE_WALL: Saves a mesh to file
      !-----------------------------------------------------------------------
      ! param[in]: filename. The filename loaded in. Will add "_acc" to end
      ! param[in, out]: istat. Integer that shows error if != 0
      !-----------------------------------------------------------------------      
         CHARACTER(LEN=*), INTENT(in) :: filename
         INTEGER, OPTIONAL :: istat
   
         INTEGER            :: ppos, i
         CHARACTER(LEN=256) :: filename_res
   
         ppos = scan(trim(filename),".", BACK= .true.)
         filename_res = filename(1:ppos - 1) // "_acc.dat"
   
         OPEN(10, file=filename_res, STATUS="UNKNOWN", ACTION="WRITE")
   
         WRITE(10, *) 'MACHINE: ', TRIM(machine_string(10:))
         WRITE(10, *) 'DATE: ', TRIM(date(6:))
         WRITE(10, "(I12, I12)") 0, 0
         WRITE(10, "(I12, I12)") nvertex, nface
         WRITE(10, "(E20.10, E20.10, E20.10)") TRANSPOSE(vertex)
         WRITE(10, "(I12, I12, I12)") TRANSPOSE(face)

         IF (ASSOCIATED(wall%blocks)) THEN
   
            WRITE(10, "(I12, I12, I12, I12)") wall%nblocks, wall%step
            WRITE(10, "(E20.10, I12, I12, I12)") wall%stepsize, wall%br
      
            DO i=1,wall%nblocks
               WRITE(10, "(6(E20.10))") wall%blocks(i)%rmin(1), wall%blocks(i)%rmax(1), wall%blocks(i)%rmin(2), &
                                       wall%blocks(i)%rmax(2), wall%blocks(i)%rmin(3), wall%blocks(i)%rmax(3)
               WRITE(10, "(I12)") wall%blocks(i)%nfaces
               IF (wall%blocks(i)%nfaces > 0) WRITE(10, "(I12)") wall%blocks(i)%face
            END DO
            IF (lverb) WRITE(6, *) 'Accelerated wall written to mesh file. Location: ', filename_res
         ELSE
            IF (lverb) WRITE(6, *) 'Wall written to mesh file. Location: ', filename_res
         END IF
   
         CLOSE(10)        
      END SUBROUTINE WRITE_WALL

      SUBROUTINE LINE_BOX_INTERSECTION(x1i, y1i, z1i, x2i, y2i, z2i, xmin, ymin, zmin, xmax, ymax, zmax, intersects)
      !-----------------------------------------------------------------------
      ! LINE_BOX_INTERSECTION: Determines if line intersects a 3D box
      !-----------------------------------------------------------------------
      ! param[in]: x1,y1,z1,x2,y2,z2 Line Endpoints
      ! param[in]: xmin,ymin,zmin,xmax,ymax,zmax Box bounds
      ! param[out]: intersects Does Line Intersect
      !----------------------------------------------------------------------- 
      IMPLICIT NONE
      DOUBLE PRECISION, INTENT(in) :: x1i, y1i, z1i, x2i, y2i, z2i, xmin, ymin, zmin, xmax, ymax, zmax
      LOGICAL, INTENT(out) :: intersects
      DOUBLE PRECISION :: tmin, tmax, x1, y1, z1, x2, y2, z2, dx, dy, dz, t1, t2

      intersects = .false.

      ! Force x2 > x1
      x1 = MIN(x1i,x2i)
      y1 = MIN(y1i,y2i)
      z1 = MIN(z1i,z2i)
      x2 = MAX(x1i,x2i)
      y2 = MAX(y1i,y2i)
      z2 = MAX(z1i,z2i)

      ! Helpers
      dx = one/(x2-x1)
      dy = one/(y2-y1)
      dz = one/(z2-z1)

      ! Calculate the minimum and maximum values of t for each axis
      tmin = zero
      tmax = 1.0D20
      t1 = (xmin - x1) * dx
      t2 = (xmax - x1) * dx
      tmin = max(tmin, min(min(t1, t2), tmax))
      tmax = min(tmax, max(max(t1, t2), tmin))
      t1 = (ymin - y1) * dy
      t2 = (ymax - y1) * dy
      tmin = max(tmin, min(min(t1, t2), tmax))
      tmax = min(tmax, max(max(t1, t2), tmin))
      t1 = (zmin - z1) * dz
      t2 = (zmax - z1) * dz
      tmin = max(tmin, min(min(t1, t2), tmax))
      tmax = min(tmax, max(max(t1, t2), tmin))

      ! Check if there is an intersection

      intersects = (tmin < tmax)

      RETURN

      END SUBROUTINE LINE_BOX_INTERSECTION

   

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Wall Destructors
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE BLOCK_DESTROY(this)
      !-----------------------------------------------------------------------
      ! BLOCK_DESTROY: Destroys a single block
      !-----------------------------------------------------------------------
      ! param[in, out]: this. The block to destroy 
      !-----------------------------------------------------------------------
         IMPLICIT NONE
         TYPE(block), INTENT(inout) :: this
         IF (this%nfaces > 0) CALL free_mpi_array(this%win_face, this%face, this%isshared)
         this%isshared = .FALSE.
         this%nfaces = 0
         RETURN
         END SUBROUTINE BLOCK_DESTROY

      SUBROUTINE wall_free(istat,shared_comm)
      !-----------------------------------------------------------------------
      ! wall_free: Removes wall from memory
      !-----------------------------------------------------------------------
      ! param[in, out]: istat. Integer that shows error if != 0
      ! param[in, out]: shared_comm. MPI communicator, handles shared memory
      !-----------------------------------------------------------------------
#if defined(MPI_OPT)
      USE mpi
#endif
      IMPLICIT NONE
      INTEGER, INTENT(inout) :: istat
      INTEGER, INTENT(inout), OPTIONAL :: shared_comm
      INTEGER :: i
      ! Check if have to delete accelerated wall or the normal one
      IF (ASSOCIATED(wall%blocks)) THEN
         DO i = 1, wall%nblocks
            CALL BLOCK_DESTROY(wall%blocks(i))
         END DO
         DEALLOCATE(wall%blocks)
         wall%nblocks = -1
      END IF
      IF (PRESENT(shared_comm)) THEN
#if defined(MPI_OPT)
         CALL MPI_WIN_FENCE(0,win_vertex,istat)
         CALL MPI_WIN_FREE(win_vertex,istat)
         CALL MPI_WIN_FENCE(0,win_face,istat)
         CALL MPI_WIN_FREE(win_face,istat)
         CALL MPI_WIN_FENCE(0,win_fn,istat)
         CALL MPI_WIN_FREE(win_fn,istat)
         CALL MPI_WIN_FENCE(0,win_a0,istat)
         CALL MPI_WIN_FREE(win_a0,istat)
         CALL MPI_WIN_FENCE(0,win_v0,istat)
         CALL MPI_WIN_FREE(win_v0,istat)
         CALL MPI_WIN_FENCE(0,win_v1,istat)
         CALL MPI_WIN_FREE(win_v1,istat)
         CALL MPI_WIN_FENCE(0,win_dot00,istat)
         CALL MPI_WIN_FREE(win_dot00,istat)
         CALL MPI_WIN_FENCE(0,win_dot01,istat)
         CALL MPI_WIN_FREE(win_dot01,istat)
         CALL MPI_WIN_FENCE(0,win_dot11,istat)
         CALL MPI_WIN_FREE(win_dot11,istat)
         CALL MPI_WIN_FENCE(0,win_d,istat)
         CALL MPI_WIN_FREE(win_d,istat)
         CALL MPI_WIN_FENCE(0,win_invdenom,istat)
         CALL MPI_WIN_FREE(win_invdenom,istat)
         CALL MPI_WIN_FENCE(0,win_ihit,istat)
         CALL MPI_WIN_FREE(win_ihit,istat)    
         IF (ASSOCIATED(vertex)) NULLIFY(vertex)
         IF (ASSOCIATED(face)) NULLIFY(face)
         IF (ASSOCIATED(FN)) NULLIFY(FN)
         IF (ASSOCIATED(A0)) NULLIFY(A0)
         IF (ASSOCIATED(V0)) NULLIFY(V0)
         IF (ASSOCIATED(V1)) NULLIFY(V1)
         IF (ASSOCIATED(DOT00)) NULLIFY(DOT00)
         IF (ASSOCIATED(DOT01)) NULLIFY(DOT01)
         IF (ASSOCIATED(DOT11)) NULLIFY(DOT11)
         IF (ASSOCIATED(d)) NULLIFY(d)
         IF (ASSOCIATED(ihit_array)) NULLIFY(ihit_array)
      ELSE
#endif
         IF (ASSOCIATED(FN)) DEALLOCATE(FN)
         IF (ASSOCIATED(A0)) DEALLOCATE(A0)
         IF (ASSOCIATED(V0)) DEALLOCATE(V0)
         IF (ASSOCIATED(V1)) DEALLOCATE(V1)
         IF (ASSOCIATED(DOT00)) DEALLOCATE(DOT00)
         IF (ASSOCIATED(DOT01)) DEALLOCATE(DOT01)
         IF (ASSOCIATED(DOT11)) DEALLOCATE(DOT11)
         IF (ASSOCIATED(invDenom)) DEALLOCATE(invDenom)
         IF (ASSOCIATED(d)) DEALLOCATE(d)
         IF (ASSOCIATED(vertex)) DEALLOCATE(vertex)
         IF (ASSOCIATED(face)) DEALLOCATE(face)
         IF (ASSOCIATED(ihit_array)) DEALLOCATE(ihit_array)
      END IF
      
      machine_string=''
      date=''
      nface = -1
      nvertex = -1

      lwall_loaded = .false.
      RETURN
      END SUBROUTINE wall_free

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Memory Allocation Subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE mpialloc_1d_int(array,n1,subid,mymaster,share_comm,win)
      !-----------------------------------------------------------------------
      ! mpialloc_1d_int: Allocated a 1D integer array to shared memory
      ! Taken from LIBSTELL/Sources/Modules/mpi_sharemem.f90
      ! Included here to reduce dependencies
      !-----------------------------------------------------------------------
      ! Libraries
#if defined(MPI_OPT)
      USE mpi
#endif
      USE ISO_C_BINDING
      IMPLICIT NONE
      ! Arguments
      INTEGER, POINTER, INTENT(inout) :: array(:)
      INTEGER, INTENT(in) :: n1
      INTEGER, INTENT(in) :: subid
      INTEGER, INTENT(in) :: mymaster
      INTEGER, INTENT(inout) :: share_comm
      INTEGER, INTENT(inout) :: win
      ! Variables
      INTEGER :: disp_unit, ier
      INTEGER :: array_shape(1)
#if defined(MPI_OPT)
      INTEGER(KIND=MPI_ADDRESS_KIND) :: window_size
#endif
      TYPE(C_PTR) :: baseptr
      ! Initialization
      ier = 0
      array_shape(1) = n1
      disp_unit = 1
#if defined(MPI_OPT)
      window_size = 0_MPI_ADDRESS_KIND
      IF (subid == mymaster) window_size = INT(n1,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      CALL MPI_WIN_ALLOCATE_SHARED(window_size, disp_unit, MPI_INFO_NULL, share_comm, baseptr, win ,ier)
      IF (subid /= mymaster) CALL MPI_WIN_SHARED_QUERY(win, 0, window_size, disp_unit, baseptr, ier)
      CALL C_F_POINTER(baseptr, array, array_shape)
#endif
      RETURN
      END SUBROUTINE mpialloc_1d_int

      SUBROUTINE mpialloc_1d_dbl(array,n1,subid,mymaster,share_comm,win)
      !-----------------------------------------------------------------------
      ! mpialloc_1d_int: Allocated a 1D double array to shared memory
      ! Taken from LIBSTELL/Sources/Modules/mpi_sharemem.f90
      ! Included here to reduce dependencies
      !-----------------------------------------------------------------------
      ! Libraries
#if defined(MPI_OPT)
      USE mpi
#endif
      USE ISO_C_BINDING
      IMPLICIT NONE
      ! Arguments
      DOUBLE PRECISION, POINTER, INTENT(inout) :: array(:)
      INTEGER, INTENT(in) :: n1
      INTEGER, INTENT(in) :: subid
      INTEGER, INTENT(in) :: mymaster
      INTEGER, INTENT(inout) :: share_comm
      INTEGER, INTENT(inout) :: win
      ! Variables
      INTEGER :: disp_unit, ier
      INTEGER :: array_shape(1)
#if defined(MPI_OPT)
      INTEGER(KIND=MPI_ADDRESS_KIND) :: window_size
#endif
      TYPE(C_PTR) :: baseptr
      ! Initialization
      ier = 0
      array_shape(1) = n1
      disp_unit = 1
#if defined(MPI_OPT)
      window_size = 0_MPI_ADDRESS_KIND
      IF (subid == mymaster) window_size = INT(n1,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      CALL MPI_WIN_ALLOCATE_SHARED(window_size, disp_unit, MPI_INFO_NULL, share_comm, baseptr, win ,ier)
      IF (subid /= mymaster) CALL MPI_WIN_SHARED_QUERY(win, 0, window_size, disp_unit, baseptr, ier)
      CALL C_F_POINTER(baseptr, array, array_shape)
#endif
      RETURN
      END SUBROUTINE mpialloc_1d_dbl

      SUBROUTINE mpialloc_2d_int(array,n1,n2,subid,mymaster,share_comm,win)
      !-----------------------------------------------------------------------
      ! mpialloc_1d_int: Allocated a 2D integer array to shared memory
      ! Taken from LIBSTELL/Sources/Modules/mpi_sharemem.f90
      ! Included here to reduce dependencies
      !-----------------------------------------------------------------------
      ! Libraries
#if defined(MPI_OPT)
      USE mpi
#endif
      USE ISO_C_BINDING
      IMPLICIT NONE
      ! Arguments
      INTEGER, POINTER, INTENT(inout) :: array(:,:)
      INTEGER, INTENT(in) :: n1
      INTEGER, INTENT(in) :: n2
      INTEGER, INTENT(in) :: subid
      INTEGER, INTENT(in) :: mymaster
      INTEGER, INTENT(inout) :: share_comm
      INTEGER, INTENT(inout) :: win
      ! Variables
      INTEGER :: disp_unit, ier
      INTEGER :: array_shape(2)
#if defined(MPI_OPT)
      INTEGER(KIND=MPI_ADDRESS_KIND) :: window_size
#endif
      TYPE(C_PTR) :: baseptr
      ! Initialization
      ier = 0
      array_shape(1) = n1
      array_shape(2) = n2
      disp_unit = 1
#if defined(MPI_OPT)
      window_size = 0_MPI_ADDRESS_KIND
      IF (subid == mymaster) window_size = INT(n1*n2,MPI_ADDRESS_KIND)*4_MPI_ADDRESS_KIND
      CALL MPI_WIN_ALLOCATE_SHARED(window_size, disp_unit, MPI_INFO_NULL, share_comm, baseptr, win ,ier)
      IF (subid /= mymaster) CALL MPI_WIN_SHARED_QUERY(win, 0, window_size, disp_unit, baseptr, ier)
      CALL C_F_POINTER(baseptr, array, array_shape)
#endif
      RETURN
      END SUBROUTINE mpialloc_2d_int

      SUBROUTINE mpialloc_2d_dbl(array,n1,n2,subid,mymaster,share_comm,win)
      !-----------------------------------------------------------------------
      ! mpialloc_1d_int: Allocated a 2D double array to shared memory
      ! Taken from LIBSTELL/Sources/Modules/mpi_sharemem.f90
      ! Included here to reduce dependencies
      !-----------------------------------------------------------------------
      ! Libraries
#if defined(MPI_OPT)
      USE mpi
#endif
      USE ISO_C_BINDING
      IMPLICIT NONE
      ! Arguments
      DOUBLE PRECISION, POINTER, INTENT(inout) :: array(:,:)
      INTEGER, INTENT(in) :: n1
      INTEGER, INTENT(in) :: n2
      INTEGER, INTENT(in) :: subid
      INTEGER, INTENT(in) :: mymaster
      INTEGER, INTENT(inout) :: share_comm
      INTEGER, INTENT(inout) :: win
      ! Variables
      INTEGER :: disp_unit, ier
      INTEGER :: array_shape(2)
#if defined(MPI_OPT)
      INTEGER(KIND=MPI_ADDRESS_KIND) :: window_size
#endif
      TYPE(C_PTR) :: baseptr
      ! Initialization
      ier = 0
      array_shape(1) = n1
      array_shape(2) = n2
      disp_unit = 1
#if defined(MPI_OPT)
      window_size = 0_MPI_ADDRESS_KIND
      IF (subid == mymaster) window_size = INT(n1*n2,MPI_ADDRESS_KIND)*8_MPI_ADDRESS_KIND
      CALL MPI_WIN_ALLOCATE_SHARED(window_size, disp_unit, MPI_INFO_NULL, share_comm, baseptr, win ,ier)
      IF (subid /= mymaster) CALL MPI_WIN_SHARED_QUERY(win, 0, window_size, disp_unit, baseptr, ier)
      CALL C_F_POINTER(baseptr, array, array_shape)
#endif
      RETURN
      END SUBROUTINE mpialloc_2d_dbl

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!    Memory Freeing Subroutines
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      SUBROUTINE free_mpi_array1d_int(win_local,array_local,isshared)
         IMPLICIT NONE
         LOGICAL, INTENT(in) :: isshared
         INTEGER, INTENT(inout) :: win_local
         INTEGER, POINTER, INTENT(inout) :: array_local(:)
         INTEGER :: istat
         istat=0
#if defined(MPI_OPT)
         IF (isshared) THEN
            CALL MPI_WIN_FENCE(0, win_local,istat)
            CALL MPI_WIN_FREE(win_local,istat)
            IF (ASSOCIATED(array_local)) NULLIFY(array_local)
         ELSE
#endif
            IF (ASSOCIATED(array_local)) DEALLOCATE(array_local)
#if defined(MPI_OPT)
         ENDIF
#endif
         RETURN
         END SUBROUTINE free_mpi_array1d_int
   
         SUBROUTINE free_mpi_array1d_flt(win_local,array_local,isshared)
         IMPLICIT NONE
         LOGICAL, INTENT(in) :: isshared
         INTEGER, INTENT(inout) :: win_local
         REAL, POINTER, INTENT(inout) :: array_local(:)
         INTEGER :: istat
         istat=0
#if defined(MPI_OPT)
         IF (isshared) THEN
            CALL MPI_WIN_FENCE(0, win_local,istat)
            CALL MPI_WIN_FREE(win_local,istat)
            IF (ASSOCIATED(array_local)) NULLIFY(array_local)
         ELSE
#endif
            IF (ASSOCIATED(array_local)) DEALLOCATE(array_local)
#if defined(MPI_OPT)
         ENDIF
#endif
         RETURN
         END SUBROUTINE free_mpi_array1d_flt
   
         SUBROUTINE free_mpi_array1d_dbl(win_local,array_local,isshared)
         IMPLICIT NONE
         LOGICAL, INTENT(in) :: isshared
         INTEGER, INTENT(inout) :: win_local
         DOUBLE PRECISION, POINTER, INTENT(inout) :: array_local(:)
         INTEGER :: istat
         istat=0
#if defined(MPI_OPT)
         IF (isshared) THEN
            CALL MPI_WIN_FENCE(0, win_local,istat)
            CALL MPI_WIN_FREE(win_local,istat)
            IF (ASSOCIATED(array_local)) NULLIFY(array_local)
         ELSE
#endif
            IF (ASSOCIATED(array_local)) DEALLOCATE(array_local)
#if defined(MPI_OPT)
         ENDIF
#endif
         RETURN
         END SUBROUTINE free_mpi_array1d_dbl
   
         SUBROUTINE free_mpi_array2d_int(win_local,array_local,isshared)
         IMPLICIT NONE
         LOGICAL, INTENT(in) :: isshared
         INTEGER, INTENT(inout) :: win_local
         INTEGER, POINTER, INTENT(inout) :: array_local(:,:)
         INTEGER :: istat
         istat=0
#if defined(MPI_OPT)
         IF (isshared) THEN
            CALL MPI_WIN_FENCE(0, win_local,istat)
            CALL MPI_WIN_FREE(win_local,istat)
            IF (ASSOCIATED(array_local)) NULLIFY(array_local)
         ELSE
#endif
            IF (ASSOCIATED(array_local)) DEALLOCATE(array_local)
#if defined(MPI_OPT)
         ENDIF
#endif
         RETURN
         END SUBROUTINE free_mpi_array2d_int
   
         SUBROUTINE free_mpi_array2d_dbl(win_local,array_local,isshared)
         IMPLICIT NONE
         LOGICAL, INTENT(in) :: isshared
         INTEGER, INTENT(inout) :: win_local
         DOUBLE PRECISION, POINTER, INTENT(inout) :: array_local(:,:)
         INTEGER :: istat
         istat=0
#if defined(MPI_OPT)
         IF (isshared) THEN
            CALL MPI_WIN_FENCE(0, win_local,istat)
            CALL MPI_WIN_FREE(win_local,istat)
            IF (ASSOCIATED(array_local)) NULLIFY(array_local)
         ELSE
#endif
            IF (ASSOCIATED(array_local)) DEALLOCATE(array_local)
#if defined(MPI_OPT)
         ENDIF
#endif
         RETURN
         END SUBROUTINE free_mpi_array2d_dbl
!-----------------------------------------------------------------------
!     End Module
!-----------------------------------------------------------------------
      END MODULE wall_mod
