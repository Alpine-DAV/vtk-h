#==============================================================================
# See License.txt
#==============================================================================

set(Headers
  Filter.hpp
  Clip.hpp
  MarchingCubes.hpp
  Threshold.hpp
  WaveletCompressor.hpp
  )

set(Sources
  Clip.cpp
  MarchingCubes.cpp
  Threshold.cpp
  WaveletCompressor.cpp
  )

add_library(vtkh_filters ${Sources} ${Headers})

target_include_directories(vtkh_filters
  PUBLIC
    $<BUILD_INTERFACE:${VTKh_SOURCE_DIR}>
    $<BUILD_INTERFACE:${VTKh_BINARY_DIR}>
    $<INSTALL_INTERFACE:${VTKh_INSTALL_INCLUDE_DIR}>
  )

target_link_libraries(vtkh_filters
  PUBLIC vtkh_core
  PRIVATE vtkh_dependencies_private)

# Install libraries
install(TARGETS vtkh_filters
  EXPORT ${VTKh_EXPORT_NAME}
  ARCHIVE DESTINATION ${VTKh_INSTALL_LIB_DIR}
  LIBRARY DESTINATION ${VTKh_INSTALL_LIB_DIR}
  RUNTIME DESTINATION ${VTKh_INSTALL_BIN_DIR}
  )

# Install headers
install(FILES ${Headers}
  DESTINATION ${VTKh_INSTALL_INCLUDE_DIR}/vtkh/filters)

if (ENABLE_MPI)
  add_library(vtkh_filters_par ${Sources} ${Headers})

  target_include_directories(vtkh_filters_par
    PUBLIC
      $<BUILD_INTERFACE:${VTKh_SOURCE_DIR}>
      $<BUILD_INTERFACE:${VTKh_BINARY_DIR}>
      $<INSTALL_INTERFACE:${VTKh_INSTALL_INCLUDE_DIR}>
    )
  target_link_libraries(vtkh_filters_par
    PUBLIC vtkh_core_par
    PRIVATE vtkh_dependencies_private)


  # Install libraries
  install(TARGETS vtkh_filters_par
    EXPORT ${VTKh_EXPORT_NAME}
    ARCHIVE DESTINATION ${VTKh_INSTALL_LIB_DIR}
    LIBRARY DESTINATION ${VTKh_INSTALL_LIB_DIR}
    RUNTIME DESTINATION ${VTKh_INSTALL_BIN_DIR}
    )
endif()
