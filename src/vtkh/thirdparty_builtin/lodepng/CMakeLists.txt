#==============================================================================
# See License.txt
#==============================================================================
set(Headers
  lodepng.h
  )

set(Sources
  lodepng.cpp
  )

add_library(lodepng ${Sources})
target_include_directories(lodepng INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
  $<INSTALL_INTERFACE:vtkh/thirdparty_builtin/lodepng>  # <prefix>/vtkh/thirdparty_builtin/lodepng
  )

# install target for serial ascent lib
install(TARGETS lodepng
  EXPORT ${VTKh_EXPORT_NAME}
  ARCHIVE DESTINATION ${VTKh_INSTALL_LIB_DIR}
  LIBRARY DESTINATION ${VTKh_INSTALL_LIB_DIR}
  RUNTIME DESTINATION ${VTKh_INSTALL_BIN_DIR}
  )

# Install headers
install(FILES ${Headers}
  DESTINATION ${VTKh_INSTALL_INCLUDE_DIR}/vtkh/thirdparty_builtin/lodepng
  )

# Install others files.
install(FILES License.txt README.md
  DESTINATION ${VTKh_INSTALL_SHARED_RESOURCES_DIR}/thirdparty_builtin/lodepng
  )
