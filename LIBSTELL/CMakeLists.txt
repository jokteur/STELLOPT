file(GLOB libstell_src RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    Sources/Coils/*.f
    Sources/FFTpack/*.f
    Sources/GMRes/*.f
    Sources/Lsoe/*.f
    Sources/Miscel/*.f
    Sources/Modules/*.f
    Sources/NCLASS/*.f
    Sources/NumerMeth/*.f
    Sources/Optimization/*.f
    Sources/Pspline/*.f
)
file(GLOB libstell_freeform RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    Sources/Ezcdf/*.f90
    Sources/FFTpack/*.f90
    Sources/Miscel/*.f90
    Sources/Modules/*.f90
    Sources/Optimization/*.f90
    Sources/Pspline/*.f90
)

set_source_files_properties(${libstell_src} PROPERTIES COMPILE_FLAGS ${Fortran_FixedFORM_FLAGS})
set_source_files_properties(${libstell_freeform} PROPERTIES COMPILE_FLAGS ${Fortran_FreeFORM_FLAGS})

add_library(libstell STATIC ${libstell_src} ${libstell_freeform})
target_link_libraries(libstell PRIVATE ${NETCDF_LIBRARIES} MPI::MPI_Fortran ${HDF5_LIBRARIES})
target_include_directories(libstell PUBLIC ${NETCDF_INCLUDES} MPI::MPI_Fortran ${HDF5_INCLUDE_DIRS})

get_target_property(LIB_DIR libstell BINARY_DIR)
set(LIBSTELL_MOD_DIR ${LIB_DIR} PARENT_SCOPE)