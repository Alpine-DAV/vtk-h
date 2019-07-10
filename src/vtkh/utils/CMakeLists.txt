#==============================================================================
# See License.txt
#==============================================================================
set(vtkh_utils_headers
  Logger.hpp
  Mutex.hpp
  PNGEncoder.hpp
  StatisticsDB.hpp
  ThreadSafeContainer.hpp
  vtkm_array_utils.hpp
  vtkm_dataset_info.hpp
  )

set(vtkh_utils_sources
  Logger.cpp
  PNGEncoder.cpp
  StatisticsDB.cpp
  vtkm_dataset_info.cpp
  )

set(vtkh_utils_thirdparty_libs vtkm vtkh_lodepng)

if(CUDA_FOUND)
  list(APPEND vtkh_utils_thirdparty_libs cuda)
endif()

if(ENABLE_OPENMP)
  list(APPEND vtkh_utils_thirdparty_libs openmp)
endif()

blt_add_library(
  NAME vtkh_utils
  SOURCES ${vtkh_utils_sources}
  HEADERS ${vtkh_utils_headers}
  DEPENDS_ON ${vtkh_utils_thirdparty_libs}
  )

if(CUDA_FOUND AND NOT ENABLE_SHARED_LIBS)
  set_target_properties(vtkh_utils PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()

if(ENABLE_OPENMP)
    blt_add_target_compile_flags(TO vtkh_utils FLAGS "-DVTKH_USE_OPENMP")
    list(APPEND vtkh_utils_thirdparty_libs openmp)
endif()

# Install libraries
install(TARGETS vtkh_utils
  EXPORT ${VTKh_EXPORT_NAME}
  ARCHIVE DESTINATION ${VTKh_INSTALL_LIB_DIR}
  LIBRARY DESTINATION ${VTKh_INSTALL_LIB_DIR}
  RUNTIME DESTINATION ${VTKh_INSTALL_BIN_DIR}
  )

# Install headers
install(FILES ${vtkh_utils_headers}
        DESTINATION ${VTKh_INSTALL_INCLUDE_DIR}/vtkh/utils)

if (MPI_FOUND)

    blt_add_library(
      NAME vtkh_utils_mpi
      SOURCES ${vtkh_utils_sources}
      HEADERS ${vtkh_utils_headers}
      DEPENDS_ON ${vtkh_utils_thirdparty_libs} mpi
      )

if(CUDA_FOUND AND NOT ENABLE_SHARED_LIBS)
  set_target_properties(vtkh_utils_mpi PROPERTIES CUDA_RESOLVE_DEVICE_SYMBOLS ON)
endif()

    blt_add_target_compile_flags(TO vtkh_utils_mpi FLAGS "-DVTKH_PARALLEL")

    # Install libraries
    install(TARGETS vtkh_utils_mpi
      EXPORT ${VTKh_EXPORT_NAME}
      ARCHIVE DESTINATION ${VTKh_INSTALL_LIB_DIR}
      LIBRARY DESTINATION ${VTKh_INSTALL_LIB_DIR}
      RUNTIME DESTINATION ${VTKh_INSTALL_BIN_DIR}
      )
endif()
