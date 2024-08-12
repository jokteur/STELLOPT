file(GLOB source_list RELATIVE ${CMAKE_CURRENT_SOURCE_DIR}
    Sources/*.f
)

set_source_files_properties(${source_list} PROPERTIES COMPILE_FLAGS ${Fortran_FixedFORM_FLAGS})

set(VMEC_NAME "vmec_standard")
if (${ANIMEC})
    set(VMEC_NAME "vmec_animec")
elseif (${FLOW})
    set(VMEC_NAME "vmec_flow")
endif()
    

# message(STATUS "source_list: ${LIBSTELL_MOD_DIR}")
add_executable(makegrid ${source_list})
target_link_libraries(makegrid PRIVATE ${NETCDF_LIBRARIES} MPI::MPI_Fortran ${HDF5_LIBRARIES} libstell ${LAPACK_LIBRARIES})
target_include_directories(makegrid PUBLIC ${NETCDF_INCLUDES} ${MPI_Fortran_INCLUDE_DIRS} ${HDF5_INCLUDE_DIRS} ${LIBSTELL_MOD_DIR})
install(TARGETS makegrid)