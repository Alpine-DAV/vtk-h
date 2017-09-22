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

add_library(vtkh_utils ${Sources})
target_link_libraries(vtkh_utils PRIVATE lodepng)


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
