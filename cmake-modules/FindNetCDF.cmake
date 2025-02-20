# - Find NetCDF
# Find the native NetCDF includes and library
#
#  NETCDF_INCLUDES    - where to find netcdf.h, etc
#  NETCDF_LIBRARIES   - Link these libraries when using NetCDF
#  NETCDF_FOUND       - True if NetCDF found including required interfaces (see below)
#
# Your package can require certain interfaces to be FOUND by setting these
#
#  NETCDF_CXX         - require the C++ interface and link the C++ library
#  NETCDF_F77         - require the F77 interface and link the fortran library
#  NETCDF_F90         - require the F90 interface and link the fortran library
#
# The following are not for general use and are included in
# NETCDF_LIBRARIES if the corresponding option above is set.
#
#  NETCDF_LIBRARIES_C    - Just the C interface
#  NETCDF_LIBRARIES_CXX  - C++ interface, if available
#  NETCDF_LIBRARIES_F77  - Fortran 77 interface, if available
#  NETCDF_LIBRARIES_F90  - Fortran 90 interface, if available
#
# Normal usage would be:
#  set (NETCDF_F90 "YES")
#  find_package (NetCDF REQUIRED)
#  target_link_libraries (uses_f90_interface ${NETCDF_LIBRARIES})
#  target_link_libraries (only_uses_c_interface ${NETCDF_LIBRARIES_C})

if (NETCDF_INCLUDES AND NETCDF_LIBRARIES)
  # Already in cache, be silent
  set (NETCDF_FIND_QUIETLY TRUE)
endif (NETCDF_INCLUDES AND NETCDF_LIBRARIES)

# Netcdf offers nc-config and nf-config, which is usefull to find the libraries
find_program(NC_CONFIG nc-config)
find_program(NF_CONFIG nf-config)
if (NC_CONFIG)
  exec_program(${NC_CONFIG} ARGS --includedir OUTPUT_VARIABLE NETCDF_INC_)
  exec_program(${NC_CONFIG} ARGS --libdir OUTPUT_VARIABLE NETCDF_LIB_)
endif()
if (NF_CONFIG)
  exec_program(${NF_CONFIG} ARGS --includedir OUTPUT_VARIABLE NETCDFF_INC_)
  exec_program(${NF_CONFIG} ARGS --flibs OUTPUT_VARIABLE NF_LIBS_)
  # Extract the netcdf library path from the output of nf-config
  string(REGEX MATCH "-L(.*) -lnetcdff" TMP_OUT ${NF_LIBS_})
  string(FIND ${TMP_OUT} " -lnetcdff" TMP_STRING_POS)
  math(EXPR TMP_STRING_POS "${TMP_STRING_POS}-2")
  string(SUBSTRING ${TMP_OUT} 2 ${TMP_STRING_POS} NETCDFF_LIB_)
endif()

find_path (NETCDF_INCLUDES netcdf.h
  HINTS ${NETCDF_INC_} ENV NETCDF_DIR ENV NETCDF_INC ENV NETCDF_INCLUDES ENV NETCDF_HOME)
find_library (NETCDF_LIBRARIES_C NAMES netcdf 
  HINTS ${NETCDF_LIB_} ENV NETCDF_LIB ENV NETCDF_HOME)

mark_as_advanced(NETCDF_LIBRARIES_C)

set (NetCDF_has_interfaces "YES") # will be set to NO if we're missing any interfaces
set (NetCDF_libs "${NETCDF_LIBRARIES_C}")

get_filename_component (NetCDF_lib_dirs "${NETCDF_LIBRARIES_C}" PATH)


macro (NetCDF_check_interface lang header libs)  
  if (NETCDF_${lang})
    find_path (NETCDF_INCLUDES_${lang} NAMES ${header}
      HINTS ${NETCDFF_INC_} ENV NETCDFF_HOME ENV NETCDFF_INC ENV NETCDFF_INCLUDE NO_DEFAULT_PATH)
    find_library (NETCDF_LIBRARIES_${lang} NAMES ${libs}
      HINTS ${NETCDFF_LIB_} ENV NETCDFF_HOME ENV NETCDFF_LIB NO_DEFAULT_PATH)

    mark_as_advanced (NETCDF_INCLUDES_${lang} NETCDF_LIBRARIES_${lang})

    if (NETCDF_INCLUDES_${lang} AND NETCDF_LIBRARIES_${lang})
      list (INSERT NetCDF_libs 0 ${NETCDF_LIBRARIES_${lang}}) # prepend so that -lnetcdf is last
    else (NETCDF_INCLUDES_${lang} AND NETCDF_LIBRARIES_${lang})
      set (NetCDF_has_interfaces "NO")
      message (STATUS "Failed to find NetCDF interface for ${lang}")
    endif (NETCDF_INCLUDES_${lang} AND NETCDF_LIBRARIES_${lang})
  endif (NETCDF_${lang})
endmacro (NetCDF_check_interface)

# # NetCDF_check_interface (CXX netcdf.h netcdf)
NetCDF_check_interface (F77 netcdf.inc  netcdff)
NetCDF_check_interface (F90 netcdf.mod  netcdff)

set (NETCDF_LIBRARIES "${NetCDF_libs}" CACHE STRING "All NetCDF libraries required for interface level")

# handle the QUIETLY and REQUIRED arguments and set NETCDF_FOUND to TRUE if
# all listed variables are TRUE
include (FindPackageHandleStandardArgs)
find_package_handle_standard_args (NetCDF DEFAULT_MSG NETCDF_LIBRARIES NETCDF_INCLUDES NetCDF_has_interfaces)

set(NETCDF_LIBRARIES "${NETCDF_LIBRARIES};${NETCDF_LIBRARIES_F77};${NETCDF_LIBRARIES_F90}")
set(NETCDF_INCLUDES "${NETCDF_INCLUDES};${NETCDF_INCLUDES_F77};${NETCDF_INCLUDES_F90}")

mark_as_advanced (NETCDF_LIBRARIES NETCDF_INCLUDES)
