file(GLOB source_list RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    Sources/Coils/*.f
    Sources/Ezcdf/*.f90
    Sources/FFTpack/*.f
    Sources/FFTpack/*.f90
    Sources/GMRes/*.f
    Sources/Lsoe/*.f
    Sources/Miscel/*.f
    Sources/Miscel/*.f90
    Sources/Modules/*.f90
    Sources/Modules/*.f
    Sources/NCLASS/*.f
    Sources/NumerMeth/*.f
    Sources/Optimization/*.f
    Sources/Optimization/*.f90
    Sources/Pspline/*.f
    Sources/Pspline/*.f90
)

message(a ${MPI_INCLUDE_PATH})

add_library(libstell STATIC ${source_list})
target_link_libraries(libstell PRIVATE ${NETCDF_LIBRARIES} MPI::MPI_Fortran)
target_include_directories(libstell PUBLIC ${NETCDF_INCLUDES} MPI::MPI_Fortran)