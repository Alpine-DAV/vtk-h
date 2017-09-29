#==============================================================================
# See License.txt
#==============================================================================
set(Headers
  Image.hpp
  ImageCompositor.hpp
  MeshRenderer.hpp
  RayTracer.hpp
  Render.hpp
  Renderer.hpp
  VolumeRenderer.hpp
  compositing/Compositor.hpp
  )

set(Sources
  Image.cpp
  MeshRenderer.cpp
  RayTracer.cpp
  Render.cpp
  Renderer.cpp
  VolumeRenderer.cpp
  compositing/Compositor.cpp
  )


set(vtkh_rendering_deps vtkh_core vtkh_utils )

if(VTKm_CUDA_FOUND)
    # triggers cuda compile
    list(APPEND vtkh_rendering_deps cuda)
endif()

blt_add_library(
  NAME vtkh_rendering
  SOURCES ${Sources}
  HEADERS ${Headers}
  DEPENDS_ON ${vtkh_rendering_deps}
  )

# target_link_libraries(vtkh_rendering
#   PRIVATE vtkh_dependencies_private)
#
# target_include_directories(vtkh_rendering
#   PUBLIC
#     $<BUILD_INTERFACE:${VTKh_SOURCE_DIR}>
#     $<BUILD_INTERFACE:${VTKh_BINARY_DIR}>
#
#     $<INSTALL_INTERFACE:${VTKh_INSTALL_INCLUDE_DIR}>  # <prefix>/include
#   )

########
# TODO: We should have a macro for install targets + headers
########


# Install libraries
install(TARGETS vtkh_rendering
  EXPORT ${VTKh_EXPORT_NAME}
  ARCHIVE DESTINATION ${VTKh_INSTALL_LIB_DIR}
  LIBRARY DESTINATION ${VTKh_INSTALL_LIB_DIR}
  RUNTIME DESTINATION ${VTKh_INSTALL_BIN_DIR}
  )

# Install headers
install(FILES ${Headers}
  DESTINATION ${VTKh_INSTALL_INCLUDE_DIR}/vtkh/rendering)


#------------------------------------------------------------------------------
# Handle parallel library
#------------------------------------------------------------------------------
set(ParallelHeaders
  compositing/DirectSendCompositor.hpp
  compositing/DIYCompositor.hpp
  compositing/RadixKCompositor.hpp
  compositing/vtkh_diy_collect.hpp
  compositing/vtkh_diy_image_block.hpp
  compositing/vtkh_diy_utils.hpp
  )

set(ParallelSources
  compositing/DIYCompositor.cpp
  compositing/DirectSendCompositor.cpp
  compositing/RadixKCompositor.cpp
  )

if (ENABLE_MPI)
    
    set(vtkh_rendering_par_deps vtkh_core_par vtkh_utils )

    # necessary b/c of how find_cuda works
    if(VTKm_CUDA_FOUND)
         set(CUDA_NVCC_FLAGS_PREV ${CUDA_NVCC_FLAGS})
         set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -D PARALLEL")
         # triggers cuda compile
         list(APPEND vtkh_rendering_par_deps cuda)
    endif()
    
    blt_add_library(
        NAME vtkh_rendering_par
        SOURCES ${Sources} ${ParallelSources}
        HEADERS ${Headers} ${ParallelHeaders}
        DEPENDS_ON ${vtkh_rendering_par_deps}
    )

    if(VTKm_CUDA_FOUND)
         SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS_PREVIOUS})
    else()
         blt_add_target_compile_flags(TO vtkh_rendering_par FLAGS "-D PARALLEL")
    endif()


  # target_link_libraries(vtkh_rendering_par
  #   PRIVATE vtkh_dependencies_private)
  # target_include_directories(vtkh_rendering_par
  #   PUBLIC
  #     $<BUILD_INTERFACE:${VTKh_SOURCE_DIR}>
  #     $<BUILD_INTERFACE:${VTKh_BINARY_DIR}>
  #
  #     $<INSTALL_INTERFACE:${VTKh_INSTALL_INCLUDE_DIR}>  # <prefix>/include
  #   )

  ########
  # TODO: We should have a macro for install targets + headers
  ########


  # Install libraries
  install(TARGETS vtkh_rendering_par
    EXPORT ${VTKh_EXPORT_NAME}
    ARCHIVE DESTINATION ${VTKh_INSTALL_LIB_DIR}
    LIBRARY DESTINATION ${VTKh_INSTALL_LIB_DIR}
    RUNTIME DESTINATION ${VTKh_INSTALL_BIN_DIR}
    )

  # Install headers
  install(FILES ${ParallelHeaders}
    DESTINATION ${VTKh_INSTALL_INCLUDE_DIR}/vtkh/rendering)

endif()
