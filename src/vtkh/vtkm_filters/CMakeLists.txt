#==============================================================================
# See License.txt
#==============================================================================

set(vtkm_filters_headers
  vtkmCellAverage.hpp
  vtkmClip.hpp
  vtkmCleanGrid.hpp
  vtkmClipWithField.hpp
  vtkmExtractStructured.hpp
  vtkmLagrangian.hpp
  vtkmMarchingCubes.hpp
  vtkmPointAverage.hpp
  vtkmThreshold.hpp
  vtkmVectorMagnitude.hpp
  )

set(vtkm_filters_sources
  vtkmCellAverage.cpp
  vtkmClip.cpp
  vtkmCleanGrid.cpp
  vtkmClipWithField.cpp
  vtkmExtractStructured.cpp
  vtkmLagrangian.cpp
  vtkmMarchingCubes.cpp
  vtkmPointAverage.cpp
  vtkmThreshold.cpp
  vtkmVectorMagnitude.cpp
  )

set(vtkm_filters_deps vtkm)

if(CUDA_FOUND)
    # triggers cuda compile
    list(APPEND vtkm_filters_deps cuda)
    # defend against vtkm code cuda can't reseason about
    if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS 10.0.0)
        set(particle_cuda_src "${CMAKE_CURRENT_SOURCE_DIR}/vtkmLagrangian.cpp")
        set_source_files_properties(${particle_cuda_src} PROPERTIES COMPILE_FLAGS "-Xptxas --opt-level=0")
    endif()
endif()


blt_add_library(
  NAME vtkm_compiled_filters
  SOURCES ${vtkm_filters_sources}
  HEADERS ${vtkm_filters_headers}
  DEPENDS_ON ${vtkm_filters_deps}
  )

if(ENABLE_OPENMP)
  if(CUDA_FOUND)
    blt_add_target_compile_flags(TO vtkm_compiled_filters FLAGS "-Xcompiler ${OpenMP_CXX_FLAGS}")
  else()
    blt_add_target_compile_flags(TO vtkm_compiled_filters FLAGS "${OpenMP_CXX_FLAGS}")
  endif()
endif()

# Install libraries
install(TARGETS vtkm_compiled_filters
  EXPORT ${VTKh_EXPORT_NAME}
  ARCHIVE DESTINATION ${VTKh_INSTALL_LIB_DIR}
  LIBRARY DESTINATION ${VTKh_INSTALL_LIB_DIR}
  RUNTIME DESTINATION ${VTKh_INSTALL_BIN_DIR}
  )

# Install headers
install(FILES ${vtkm_filters_headers}
  DESTINATION ${VTKh_INSTALL_INCLUDE_DIR}/vtkh/vtkm_filters)
