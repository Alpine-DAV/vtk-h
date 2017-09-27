#==============================================================================
# See License.txt
#==============================================================================
set(Headers
  PNGEncoder.hpp
  vtkm_array_utils.hpp
  vtkm_dataset_info.hpp
  )

set(Sources
  PNGEncoder.cpp
  )


blt_add_library(
  NAME vtkh_utils
  SOURCES ${Sources}
  HEADERS ${Headers}
  DEPENDS_ON vtkh_lodepng
  )

# target_link_libraries(vtkh_utils PRIVATE vtkh_lodepng)
# target_include_directories(vtkh_utils
#   PRIVATE
#     $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
#   INTERFACE
#     $<BUILD_INTERFACE:${VTKh_SOURCE_DIR}>
#     $<BUILD_INTERFACE:${VTKh_BINARY_DIR}>
#     $<INSTALL_INTERFACE:${VTKh_INSTALL_INCLUDE_DIR}>  # <prefix>/include
#   )
#


########
# TODO: We should have a macro for install targets + headers
########

# Install libraries
install(TARGETS vtkh_utils
  EXPORT ${VTKh_EXPORT_NAME}
  ARCHIVE DESTINATION ${VTKh_INSTALL_LIB_DIR}
  LIBRARY DESTINATION ${VTKh_INSTALL_LIB_DIR}
  RUNTIME DESTINATION ${VTKh_INSTALL_BIN_DIR}
  )

# Install headers
install(FILES ${Headers}
  DESTINATION ${VTKh_INSTALL_INCLUDE_DIR}/vtkh/utils)
