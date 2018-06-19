#==============================================================================
# See License.txt
#==============================================================================

set(vtkh_filters_headers
  Filter.hpp
  CellAverage.hpp
  CleanGrid.hpp
  Clip.hpp
  ClipField.hpp
  IsoVolume.hpp
  NoOp.hpp
  MarchingCubes.hpp
  PointAverage.hpp  
  Recenter.hpp
  Threshold.hpp
  Slice.hpp
  )

set(vtkh_filters_sources
  Filter.cpp
  CellAverage.cpp
  CleanGrid.cpp
  Clip.cpp
  ClipField.cpp
  IsoVolume.cpp
  NoOp.cpp
  MarchingCubes.cpp
  PointAverage.cpp
  Recenter.cpp
  Threshold.cpp
  Slice.cpp
  )

set(vtkh_filters_deps vtkh_core  )

if(CUDA_FOUND)
    # triggers cuda compile
    list(APPEND vtkh_filters_deps cuda)
endif()


blt_add_library(
  NAME vtkh_filters
  SOURCES ${vtkh_filters_sources}
  HEADERS ${vtkh_filters_headers}
  DEPENDS_ON ${vtkh_filters_deps}
  )

if(ENABLE_OPENMP)
  if(CUDA_FOUND)
    blt_add_target_compile_flags(TO vtkh_utils FLAGS "-Xcompiler ${OpenMP_CXX_FLAGS}")
  else()
    add_target_compile_flags(TO vtkh_utils FLAGS "${OpenMP_CXX_FLAGS}")
  endif()
endif()

# Install libraries
install(TARGETS vtkh_filters
  EXPORT ${VTKh_EXPORT_NAME}
  ARCHIVE DESTINATION ${VTKh_INSTALL_LIB_DIR}
  LIBRARY DESTINATION ${VTKh_INSTALL_LIB_DIR}
  RUNTIME DESTINATION ${VTKh_INSTALL_BIN_DIR}
  )

# Install headers
install(FILES ${vtkh_filters_headers}
  DESTINATION ${VTKh_INSTALL_INCLUDE_DIR}/vtkh/filters)

if (MPI_FOUND)

    set(vtkh_filters_mpi_deps vtkh_core_mpi vtkh_utils vtkh_utils_mpi)

    # necessary b/c of how find_cuda works
    if(CUDA_FOUND)
         # triggers cuda compile
         list(APPEND vtkh_filters_mpi_deps cuda)
    endif()
    
    blt_add_library(
      NAME vtkh_filters_mpi
      SOURCES ${vtkh_filters_sources}
      HEADERS ${vtkh_filters_headers}
      DEPENDS_ON ${vtkh_filters_mpi_deps}
      )
  
    if(ENABLE_OPENMP)
          if(CUDA_FOUND)
              blt_add_target_compile_flags(TO vtkh_utils FLAGS "-Xcompiler ${OpenMP_CXX_FLAGS}")
          else()
              add_target_compile_flags(TO vtkh_utils FLAGS "${OpenMP_CXX_FLAGS}")
          endif()
    endif()

    blt_add_target_compile_flags(TO vtkh_filters_mpi FLAGS "-D VTKH_PARALLEL")


    # Install libraries
    install(TARGETS vtkh_filters_mpi
            EXPORT ${VTKh_EXPORT_NAME}
            ARCHIVE DESTINATION ${VTKh_INSTALL_LIB_DIR}
            LIBRARY DESTINATION ${VTKh_INSTALL_LIB_DIR}
            RUNTIME DESTINATION ${VTKh_INSTALL_BIN_DIR})

endif()
