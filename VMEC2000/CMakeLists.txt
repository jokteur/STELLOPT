file(GLOB source_list RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    Sources/General/*.f
    Sources/Hessian/*.f
    Sources/Initialization_Cleanup/*.f
    Sources/Input_Output/*.f
    Sources/NESTOR_vacuum/*.f
    Sources/Reconstruction/*.f
    Sources/Splines/*.f
    Sources/TimeStep/*.f
)
file(GLOB source_list_freeform RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    Sources/General/*.f90
    Sources/Hessian/*.f90
    Sources/Initialization_Cleanup/*.f90
    Sources/Input_Output/*.f90
    Sources/NESTOR_vacuum/*.f90
    Sources/Reconstruction/*.f90
    Sources/Splines/*.f90
    Sources/TimeStep/*.f90
)

set_source_files_properties(${source_list} PROPERTIES COMPILE_FLAGS ${Fortran_FixedFORM_FLAGS})
set_source_files_properties(${source_list_freeform} PROPERTIES COMPILE_FLAGS ${Fortran_FreeFORM_FLAGS})

set(VMEC_NAME "vmec_standard")
if (${ANIMEC})
    set(VMEC_NAME "vmec_animec")
elseif (${FLOW})
    set(VMEC_NAME "vmec_flow")
endif()
    

message(STATUS "source_list: ${LIBSTELL_MOD_DIR}")
add_executable(${VMEC_NAME} ${source_list} ${source_list_freeform})
target_link_libraries(${VMEC_NAME} PRIVATE ${NETCDF_LIBRARIES} MPI::MPI_Fortran ${HDF5_LIBRARIES} libstell ${LAPACK_LIBRARIES})
target_include_directories(${VMEC_NAME} PUBLIC ${NETCDF_INCLUDES} ${MPI_Fortran_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS} ${LIBSTELL_MOD_DIR})