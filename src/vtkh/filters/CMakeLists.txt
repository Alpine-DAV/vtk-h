#==============================================================================
# See License.txt
#==============================================================================

set(vtkh_filters_headers
  Filter.hpp
  CellAverage.hpp
  CleanGrid.hpp
  Clip.hpp
  ClipField.hpp
  CompositeVector.hpp
  Gradient.hpp
  GhostStripper.hpp
  HistSampling.hpp
  Histogram.hpp
  Log.hpp
  ContourTree.hpp
  IsoVolume.hpp
  NoOp.hpp
  Lagrangian.hpp
  LagrangianInterpolation.hpp
  MarchingCubes.hpp
  Particle.hpp
  ParticleAdvection.hpp
  Integrator.hpp
  PointAverage.hpp
  Recenter.hpp
  Threshold.hpp
  Statistics.hpp
  Slice.hpp
  VectorComponent.hpp
  VectorMagnitude.hpp
  communication/BoundsMap.hpp
  communication/Communicator.hpp
  communication/MemStream.h
  )

set(vtkh_comm_filters_headers
  communication/BoundsMap.hpp
  communication/Communicator.hpp
  communication/MemStream.h
  )

set(vtkh_filters_sources
  Filter.cpp
  CellAverage.cpp
  CleanGrid.cpp
  Clip.cpp
  ClipField.cpp
  CompositeVector.cpp
  Gradient.cpp
  GhostStripper.cpp
  HistSampling.cpp
  Histogram.cpp
  Log.cpp
  ContourTree.cpp
  IsoVolume.cpp
  NoOp.cpp
  Lagrangian.cpp
  LagrangianInterpolation.cpp
  MarchingCubes.cpp
  ParticleAdvection.cpp
  PointAverage.cpp
  Recenter.cpp
  Threshold.cpp
  Slice.cpp
  Statistics.cpp
  VectorComponent.cpp
  VectorMagnitude.cpp
  communication/MemStream.cpp
  )

set(vtkh_comm_filters_sources
  communication/MemStream.cpp
  )

if (ENABLE_SERIAL)
    set(vtkh_filters_deps vtkh_core vtkh_utils vtkm_compiled_filters)

    if(CUDA_FOUND)
        # triggers cuda compile
        list(APPEND vtkh_filters_deps cuda)
    endif()

    if(ENABLE_OPENMP)
        list(APPEND vtkh_filters_deps openmp)
    endif()

    blt_add_library(
      NAME vtkh_filters
      SOURCES ${vtkh_filters_sources} ${vtkh_comm_filters_sources}
      HEADERS ${vtkh_filters_headers} ${vtkh_comm_filters_headers}
      DEPENDS_ON ${vtkh_filters_deps}
      )

    target_compile_definitions(vtkh_filters PRIVATE VTKH_COMPILING_FLAG)

    if(ENABLE_OPENMP)
        target_compile_definitions(vtkh_filters PRIVATE VTKH_USE_OPENMP)
    endif()

    if(ENABLE_LOGGING)
      target_compile_definitions(vtkh_filters PRIVATE VTKH_ENABLE_LOGGING)
    endif()

    if(ENABLE_CUDA)
        set_target_properties(vtkh_filters PROPERTIES LINKER_LANGUAGE CUDA)
    endif()

    vtkm_add_target_information(vtkh_filters DEVICE_SOURCES ${vtkh_filter_sources})

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

    # Install headers
    install(FILES ${vtkh_comm_filters_headers}
      DESTINATION ${VTKh_INSTALL_INCLUDE_DIR}/vtkh/filters/communication)
endif()

if (MPI_FOUND)

  set(vtkh_filters_mpi_deps vtkh_core_mpi vtkh_utils_mpi vtkm_compiled_filters)
  set(vtkh_mpi_headers communication/Messenger.hpp
                       communication/ParticleMessenger.hpp
                       communication/RayMessenger.hpp)

    set(vtkh_mpi_sources communication/Messenger.cpp
                         communication/ParticleMessenger.cpp
                         communication/RayMessenger.cpp)

    # necessary b/c of how find_cuda works
    if(CUDA_FOUND)
         # triggers cuda compile
         list(APPEND vtkh_filters_mpi_deps cuda)
    endif()

    if(ENABLE_OPENMP)
         list(APPEND vtkh_filters_mpi_deps openmp)
    endif()

    blt_add_library(
      NAME vtkh_filters_mpi
      SOURCES ${vtkh_filters_sources} ${vtkh_comm_filters_sources} ${vtkh_mpi_sources}
      HEADERS ${vtkh_filters_headers} ${vtkh_comm_filters_headers} ${vtkh_mpi_headers}
      DEPENDS_ON ${vtkh_filters_mpi_deps}
      )

    target_compile_definitions(vtkh_filters_mpi PRIVATE VTKH_PARALLEL)
    target_compile_definitions(vtkh_filters_mpi PRIVATE VTKH_COMPILING_FLAG)

    if(ENABLE_OPENMP)
          target_compile_definitions(vtkh_filters_mpi PRIVATE VTKH_USE_OPENMP)
    endif()

    if(ENABLE_LOGGING)
      target_compile_definitions(vtkh_filters_mpi PRIVATE VTKH_ENABLE_LOGGING)
    endif()

    vtkm_add_target_information(vtkh_filters_mpi DEVICE_SOURCES ${vtkh_filter_sources})

    if(ENABLE_CUDA)
        set_target_properties(vtkh_filters_mpi PROPERTIES LINKER_LANGUAGE CUDA)
    endif()


    # Install libraries
    install(TARGETS vtkh_filters_mpi
            EXPORT ${VTKh_EXPORT_NAME}
            ARCHIVE DESTINATION ${VTKh_INSTALL_LIB_DIR}
            LIBRARY DESTINATION ${VTKh_INSTALL_LIB_DIR}
            RUNTIME DESTINATION ${VTKh_INSTALL_BIN_DIR})

endif()
