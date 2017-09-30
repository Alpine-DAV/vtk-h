#==============================================================================
# See License.txt
#==============================================================================
add_subdirectory(utils)
add_subdirectory(thirdparty_builtin)

set(Headers
  DataSet.hpp
  Error.hpp
  vtkh.hpp
  )

set(Sources
  DataSet.cpp
  vtkh.cpp
  )


set(vtkh_thirdparty_libs vtkm diy)

if(TBB_FOUND)
  list(APPEND vtkh_thirdparty_libs tbb)
endif()

if(VTKM_FOUND)
  list(APPEND vtkh_thirdparty_libs vtkm)
endif()

if(CUDA_FOUND)
  list(APPEND vtkh_thirdparty_libs cuda)
endif()


#
# add_library(vtkh_dependencies_private INTERFACE)
# if (ENABLE_TBB)
#   target_link_libraries(vtkh_dependencies_private INTERFACE tbb)
# endif()


# blt vtk dep, includes all of vtkm
blt_add_library(
  NAME vtkh_core
  SOURCES ${Sources}
  HEADERS ${Headers}
  DEPENDS_ON ${vtkh_thirdparty_libs}
  )

# target_link_libraries(vtkh_core PRIVATE vtkh_dependencies_private)
#
# target_include_directories(vtkh_core
#   PUBLIC
#     $<BUILD_INTERFACE:${VTKh_SOURCE_DIR}>
#     $<BUILD_INTERFACE:${VTKh_BINARY_DIR}>
#
#     $<INSTALL_INTERFACE:${VTKh_INSTALL_INCLUDE_DIR}>  # <prefix>/include
#   )
#

########
# TODO: We should have a macro for install targets + headers
########

# Install libraries
install(TARGETS vtkh_core
  EXPORT ${VTKh_EXPORT_NAME}
  ARCHIVE DESTINATION ${VTKh_INSTALL_LIB_DIR}
  LIBRARY DESTINATION ${VTKh_INSTALL_LIB_DIR}
  RUNTIME DESTINATION ${VTKh_INSTALL_BIN_DIR}
  )

# Install headers
install(FILES ${Headers}
  DESTINATION ${VTKh_INSTALL_INCLUDE_DIR}/vtkh)

if(MPI_FOUND)
    
    # # necessary b/c of how find_cuda works
    if(CUDA_FOUND)
         set(CUDA_NVCC_FLAGS_PREV ${CUDA_NVCC_FLAGS})
         set(CUDA_NVCC_FLAGS "${CUDA_NVCC_FLAGS} -D PARALLEL")
    endif()

    blt_add_library(
      NAME vtkh_core_par
      SOURCES ${Sources}
      HEADERS ${Headers}
      DEPENDS_ON ${vtkh_thirdparty_libs} mpi
      )

    if(CUDA_FOUND)
         SET(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS_PREV})
    else()
         blt_add_target_compile_flags(TO vtkh_core_par FLAGS "-D PARALLEL")
    endif()
    

    
  # target_link_libraries(vtkh_core_par PRIVATE vtkh_dependencies_private)
  #
  # target_compile_definitions(vtkh_core_par
  #   PUBLIC PARALLEL # rename to VTKh_IS_PARALLEL or something like that.
  #   )
  #
  # target_include_directories(vtkh_core_par
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
  install(TARGETS vtkh_core_par
    EXPORT ${VTKh_EXPORT_NAME}
    ARCHIVE DESTINATION ${VTKh_INSTALL_LIB_DIR}
    LIBRARY DESTINATION ${VTKh_INSTALL_LIB_DIR}
    RUNTIME DESTINATION ${VTKh_INSTALL_BIN_DIR}
    )
endif()

add_subdirectory(filters)
add_subdirectory(rendering)
