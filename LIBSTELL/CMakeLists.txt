file(GLOB source_list RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
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
file(GLOB source_list_freeform RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    Sources/Ezcdf/*.f90
    Sources/FFTpack/*.f90
    Sources/Miscel/*.f90
    Sources/Modules/*.f90
    Sources/Optimization/*.f90
    Sources/Pspline/*.f90
)

set_source_files_properties(${source_list} PROPERTIES COMPILE_FLAGS "-ffixed-line-length-none")
set_source_files_properties(${source_list_freeform} PROPERTIES COMPILE_FLAGS "-ffree-form -ffree-line-length-none")

add_library(libstell STATIC ${source_list} ${source_list_freeform})
target_link_libraries(libstell PRIVATE ${NETCDF_LIBRARIES} MPI::MPI_Fortran ${HDF5_LIBRARIES})
target_include_directories(libstell PUBLIC ${NETCDF_INCLUDES} MPI::MPI_Fortran ${HDF5_INCLUDE_DIRS})