#==============================================================================
# See License.txt
#==============================================================================
add_subdirectory(thirdparty_builtin)

##################################
# Hide sybols for all vtkh targets
##################################
set(CMAKE_CXX_VISIBILITY_PRESET hidden)

add_subdirectory(utils)

set(vtkh_core_headers
  DataSet.hpp
  Error.hpp
  Logger.hpp
  Timer.hpp
  StatisticsDB.hpp
  vtkh.hpp
  vtkh_exports.h
  )

set(vtkh_core_sources
  DataSet.cpp
  Logger.cpp
  Timer.cpp
  StatisticsDB.cpp
  vtkh.cpp
  )


set(vtkh_core_deps vtkm diy vtkh_utils)


if(CUDA_FOUND)
  list(APPEND vtkh_core_deps cuda)
endif()

if(ENABLE_OPENMP)
    list(APPEND vtkm_core_deps openmp)
endif()

if (ENABLE_SERIAL)
    # blt vtk dep, includes all of vtkm
    blt_add_library(
      NAME vtkh_core
      SOURCES ${vtkh_core_sources}
      HEADERS ${vtkh_core_headers}
      DEPENDS_ON ${vtkh_core_deps}
      )

    vtkm_add_target_information(vtkh_core DEVICE_SOURCES ${vtkh_core_sources})
    target_compile_definitions(vtkh_core PRIVATE VTKH_COMPILING_FLAG)
    if(ENABLE_LOGGING)
      target_compile_definitions(vtkh_core PRIVATE VTKH_ENABLE_LOGGING)
    endif()

    # Install libraries
    install(TARGETS vtkh_core
      EXPORT ${VTKh_EXPORT_NAME}
      ARCHIVE DESTINATION ${VTKh_INSTALL_LIB_DIR}
      LIBRARY DESTINATION ${VTKh_INSTALL_LIB_DIR}
      RUNTIME DESTINATION ${VTKh_INSTALL_BIN_DIR}
      )
endif()

# Install headers
install(FILES ${vtkh_core_headers}
  DESTINATION ${VTKh_INSTALL_INCLUDE_DIR}/vtkh)

if(MPI_FOUND)

  set(vtkh_core_mpi_deps vtkm diy vtkh_utils_mpi)

  if(CUDA_FOUND)
    list(APPEND vtkh_core_mpi_deps cuda)
  endif()


  blt_add_library(
      NAME vtkh_core_mpi
      SOURCES ${vtkh_core_sources}
      HEADERS ${vtkh_core_headers}
      DEPENDS_ON ${vtkh_core_mpi_deps} mpi
      )

  vtkm_add_target_information(vtkh_core_mpi DEVICE_SOURCES ${vtkh_core_sources})

  target_compile_definitions(vtkh_core_mpi PRIVATE VTKH_PARALLEL)
  target_compile_definitions(vtkh_core_mpi PRIVATE VTKH_COMPILING_FLAG)
  if(ENABLE_LOGGING)
    target_compile_definitions(vtkh_core_mpi PRIVATE VTKH_ENABLE_LOGGING)
  endif()

  # Install libraries
  install(TARGETS vtkh_core_mpi
    EXPORT ${VTKh_EXPORT_NAME}
    ARCHIVE DESTINATION ${VTKh_INSTALL_LIB_DIR}
    LIBRARY DESTINATION ${VTKh_INSTALL_LIB_DIR}
    RUNTIME DESTINATION ${VTKh_INSTALL_BIN_DIR}
    )
endif()

add_subdirectory(vtkm_filters)
add_subdirectory(filters)
add_subdirectory(rendering)
