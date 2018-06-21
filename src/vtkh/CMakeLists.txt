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


set(vtkh_thirdparty_libs vtkm diy vtkh_utils)


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

if(ENABLE_OPENMP)
  if(CUDA_FOUND)
    blt_add_target_compile_flags(TO vtkh_utils FLAGS "-Xcompiler ${OpenMP_CXX_FLAGS}")
  else()
    blt_add_target_compile_flags(TO vtkh_utils FLAGS "${OpenMP_CXX_FLAGS}")
  endif()
endif()


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

  blt_add_library(
      NAME vtkh_core_mpi
      SOURCES ${vtkh_core_sources}
      HEADERS ${vtkh_core_headers}
      DEPENDS_ON ${vtkh_thirdparty_libs} mpi vtkh_utils_mpi
      )

  if(ENABLE_OPENMP)
      if(CUDA_FOUND)
          blt_add_target_compile_flags(TO vtkh_utils FLAGS "-Xcompiler ${OpenMP_CXX_FLAGS}")
      else()
          blt_add_target_compile_flags(TO vtkh_utils FLAGS "${OpenMP_CXX_FLAGS}")
      endif()
  endif()

  blt_add_target_compile_flags(TO vtkh_core_mpi FLAGS "-D VTKH_PARALLEL")

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
