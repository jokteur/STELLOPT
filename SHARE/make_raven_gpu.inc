######################################################################
#            Modules to be used
#            Load using following commands before running compile
#            Also set correct machine name
######################################################################
# module purge
# module load git
# module load nvhpcsdk/21
# module load impi/2021.2
# module load mkl
# module load hdf5-mpi/1.10.7
# module load netcdf-mpi/4.7.0
# module load intel/19.1.3
# module load anaconda

# export MACHINE="raven_gpu"

#######################################################################
#            Define Basic Utilities
#######################################################################
  SHELL = /bin/sh
  PWD1 = `pwd`
  MYHOME = $(HOME)/bin
  PRECOMP:= fpp -C -P -DLINUX
  COMPILE = nvfortran
  COMPILE_FREE = nvfortran
  LINK    = ifort -Bstatic -o
  LINK_AR = ar -ruv
  LINK_C  = icc -shared -Wl,-z-defs

#######################################################################
#            Define Compiler Flags
#######################################################################
  FLAGS_R = -I$(MKL_HOME)/include -O2 -Kieee  \
            -Mnodefaultunit -mp=gpu -Minfo=mp -gpu=cuda11.2,cc80
  FLAGS_D = -I$(MKL_HOME)/include -O2 -Kieee  \
            -Mnodefaultunit -g -traceback -mp=gpu -Minfo=mp -gpu=cuda11.2,cc80
  LIBS    =  ${MKLROOT}/lib/intel64/libmkl_blas95_lp64.a \
             ${MKLROOT}/lib/intel64/libmkl_lapack95_lp64.a \
             ${MKLROOT}/lib/intel64/libmkl_scalapack_lp64.a \
             -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_lp64.a \
             ${MKLROOT}/lib/intel64/libmkl_sequential.a \
             ${MKLROOT}/lib/intel64/libmkl_core.a \
             ${MKLROOT}/lib/intel64/libmkl_blacs_intelmpi_lp64.a \
             -Wl,--end-group -lpthread -lm -ldl
#######################################################################
#            MPI Options
#######################################################################
  LMPI    = T
  MPI_COMPILE = mpipgfortran
  MPI_COMPILE_FREE = mpipgfortran
  MPI_COMPILE_C = mpipgfortran
  MPI_LINK = mpipgfortran
  MPI_RUN  = srun
  MPI_RUN_OPTS = --nodes=1 --ntasks-per-node=72 --time=0:30:00 -p express
  MPI_RUN_OPTS_SM   = --nodes=1 --ntasks-per-node=72 --time=0:30:00 -p express
  MPI_RUN_OPTS_MD   = --nodes=8 --ntasks-per-node=72 --time=0:30:00 -p express
  MPI_RUN_OPTS_LG   = --nodes=32 --ntasks-per-node=72 --time=0:30:00 -p express

#######################################################################
#            NAG Options
#######################################################################
  LNAG = F
  NAG_LIB = $(NAGFLIB) $(NAGMKLLIB)

#######################################################################
#            NETCDF Options
#######################################################################
  LNETCDF = T
  NETCDF_INC = $(shell nc-config --fflags)
  NETCDF_LIB = $(shell nc-config --flibs)

#######################################################################
#            HDF5 Options
#######################################################################
  LHDF5 = T
  PRECOMP+= -DHDF5_PAR
  HDF5_INC = -I$(HDF5_HOME)/include 
  HDF5_LIB = -Wl,-rpath,$(HDF5_HOME)/lib -L$(HDF5_HOME)/lib \
             -lhdf5 -lhdf5_fortran -lhdf5_hl -lhdf5hl_fortran 

#######################################################################
#            FFTW3 Options
#######################################################################
  LFFTW3 = F
  FFTW3_INC = -I$(FFTW_HOME)/include
  FFTW3_LIB = -L$(FFTW_HOME)/lib -lfftw3 -lfftw3f \
              -Wl,-rpath,$(FFTW_HOME)/lib

#######################################################################
#             PGPLOT Options
#######################################################################
  LPGPLOT = F
  PGPLOT_INC =
  PGPLOT_LIB = 

#######################################################################
#             SILO Options
#######################################################################
  LSILO = F
  SILO_INC = -I$(SILO_DIR)/include
  SILO_LIB = -L$(SILO_DIR)/lib -lsilo

#######################################################################
#            DKES/NEO Options
#######################################################################
  LDKES = T
  LNEO  = T

#######################################################################
#            GENE Options
#######################################################################
  LGENE = F
  GENE_DIR = $(GENE_PATH)/bin/obj_cobra
  FUTILS_DIR = $(GENE_PATH)/external/cobra/futils-gene-mod/src
  LIB_GENE = libgene.a
  LIB_FUTILS = libfutils.a
  GENE_INC = -I$(GENE_DIR) \
             -I$(FUTILS_DIR)
  GENE_LIB = $(GENE_DIR)/$(LIB_GENE) \
             $(FUTILS_DIR)/$(LIB_FUTILS)

#######################################################################
#            COILOPT++ Options
#######################################################################
  LCOILOPT = F
  COILOPT_INC = -I$(COILOPT_PATH)
  COILOPTPP_DIR = $(COILOPT_PATH)
  LIB_COILOPTPP = libcoilopt++.a
  COILOPT_LIB = $(COILOPT_PATH)/$(LIB_COILOPTPP) 

#######################################################################
#            TERPSICHORE Options
#######################################################################
  LTERPSICHORE= F
  TERPSICHORE_INC = -I$(TERPSICHORE_PATH)
  TERPSICHORE_DIR = $(TERPSICHORE_PATH)
  LIB_TERPSICHORE = libterpsichore.a
  TERPSICHORE_LIB = $(TERPSICHORE_DIR)/$(LIB_TERPSICHORE)

#######################################################################
#            TRAVIS Options
#######################################################################
  LTRAVIS= F
  TRAVIS_DIR = $(TRAVIS_PATH)
  LIB_TRAVIS = libtravis64_sopt.a
  LIB_MCONF  = libmconf64.a
  TRAVIS_LIB = $(TRAVIS_DIR)/lib/$(LIB_TRAVIS) \
               $(TRAVIS_DIR)/magconf/lib/$(LIB_MCONF) -lstdc++

#######################################################################
#            SFINCS Options
#######################################################################
  LSFINCS = F
  SFINCS_INC = -I$(SFINCS_DIR)
  LIB_SFINCS = libsfincs.a

  SFINCS_LIB = $(SFINCS_DIR)/$(LIB_SFINCS)

#######################################################################
#            REGCOIL Options
#######################################################################
  LREGCOIL= F
  REGCOIL_DIR = $(REGCOIL_PATH)
  REGCOIL_INC = -I$(REGCOIL_DIR) 
  LIB_REGCOIL = libregcoil.a
  REGCOIL_LIB = $(REGCOIL_DIR)/$(LIB_REGCOIL) -fopenmp

#######################################################################
#            Available Energy Options
#######################################################################
  LAEOPT= F
  AEOPT_DIR = $(AEOPT_PATH)
  AEOPT_INC = -I$(AEOPT_DIR) 
  LIB_AEOPT = libtrapAE.a
  AEOPT_LIB = $(AEOPT_PATH)/$(LIB_AEOPT)

#######################################################################
#            LIBSTELL Shared Options
#######################################################################
  LIB_SHARE = $(LIBS) $(NETCDF_LIB) $(HDF5_LIB) 

