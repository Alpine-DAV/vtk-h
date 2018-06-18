#==============================================================================
# See License.txt
#==============================================================================
add_subdirectory(utils)
add_subdirectory(thirdparty_builtin)

set(vtkh_core_headers
  DataSet.hpp
  Error.hpp
  vtkh.hpp
  )

set(vtkh_core_sources
  DataSet.cpp
  vtkh.cpp
  )


set(vtkh_thirdparty_libs ${VTKm_LIBRARIES} diy vtkh_utils)

if(ENABLE_OPENMP)
  list(APPEND vtkh_thirdparty_libs openmp)
endif()

if(CUDA_FOUND)
  list(APPEND vtkh_thirdparty_libs cuda)
endif()


# blt vtk dep, includes all of vtkm
blt_add_library(
  NAME vtkh_core
  SOURCES ${vtkh_core_sources}
  HEADERS ${vtkh_core_headers}
  DEPENDS_ON ${vtkh_thirdparty_libs}
  )

# Install libraries
install(TARGETS vtkh_core
  EXPORT ${VTKh_EXPORT_NAME}
  ARCHIVE DESTINATION ${VTKh_INSTALL_LIB_DIR}
  LIBRARY DESTINATION ${VTKh_INSTALL_LIB_DIR}
  RUNTIME DESTINATION ${VTKh_INSTALL_BIN_DIR}
  )

# Install headers
install(FILES ${vtkh_core_headers}
  DESTINATION ${VTKh_INSTALL_INCLUDE_DIR}/vtkh)

if(MPI_FOUND)
    # # necessary b/c of how find_cuda works
    if(CUDA_FOUND)
         set(CUDA_NVCC_FLAGS_PREV ${CUDA_NVCC_FLAGS})
         set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -D VTKH_PARALLEL")
    endif()

    blt_add_library(
      NAME vtkh_core_mpi
      SOURCES ${vtkh_core_sources}
      HEADERS ${vtkh_core_headers}
      DEPENDS_ON ${vtkh_thirdparty_libs} mpi vtkh_utils_mpi
      )

    if(CUDA_FOUND)
         SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS_PREV})
    else()
         blt_add_target_compile_flags(TO vtkh_core_mpi FLAGS "-D VTKH_PARALLEL")
    endif()

  # Install libraries
  install(TARGETS vtkh_core_mpi
    EXPORT ${VTKh_EXPORT_NAME}
    ARCHIVE DESTINATION ${VTKh_INSTALL_LIB_DIR}
    LIBRARY DESTINATION ${VTKh_INSTALL_LIB_DIR}
    RUNTIME DESTINATION ${VTKh_INSTALL_BIN_DIR}
    )
endif()

add_subdirectory(filters)
add_subdirectory(rendering)
