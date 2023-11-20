cmake_minimum_required(VERSION 3.22)

project(STELLOPT)

enable_language(Fortran)


# Third party libraries
find_package(HDF5 COMPONENTS Fortran)
find_package(MPI COMPONENTS Fortran)
set(NETCDF_F90 "YES")
set(NETCDF_F77 "YES")
include(cmake-modules/FindNetCDF.cmake)

add_compile_definitions(MPI_OPT NETCDF LHDF5)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fallow-argument-mismatch -Wno-argument-mismatch -fPIC -O2 -fexternal-blas -fbacktrace -fcheck=all,no-array-temps")
    # -fcheck=all,no-array-temps 
endif()

add_subdirectory(LIBSTELL)