cmake_minimum_required(VERSION 3.16)

project(STELLOPT)

enable_language(Fortran)


set(FLOW OFF CACHE BOOL "Activate VMEC flow version")
set(ANIMEC OFF CACHE BOOL "Activate ANIMEC version")

if (${FLOW})
    if (${ANIMEC})
        message(FATAL_ERROR "Cannot activate both FLOW and ANIMEC versions")
    endif()
endif()

if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
    # require at least gcc 9.3
    if (CMAKE_CXX_COMPILER_VERSION VERSION_LESS 9.3)
        message(FATAL_ERROR "GCC version must be at least 9.3!")
    endif()
endif()

# Third party libraries
find_package(HDF5 COMPONENTS Fortran)
find_package(MPI COMPONENTS Fortran)
set(NETCDF_F90 "YES")
set(NETCDF_F77 "YES")
include(cmake-modules/FindNetCDF.cmake)
find_package(LAPACK REQUIRED)

add_compile_definitions(NETCDF LHDF5 LINUX)

if (${FLOW})
    add_compile_definitions(_FLOW)
elseif(${ANIMEC})
    add_compile_definitions(_ANIMEC)
endif()

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -cpp -fallow-argument-mismatch -Wno-argument-mismatch -fPIC -O2 -fexternal-blas -fbacktrace -fcheck=all,no-array-temps")
    set(Fortran_FreeFORM_FLAGS "-ffree-form -ffree-line-length-none")
    set(Fortran_FixedFORM_FLAGS "-ffixed-form -ffixed-line-length-none")
endif()
#TODO: Intel, MSVC, Clang

add_subdirectory(LIBSTELL)
add_subdirectory(VMEC2000)