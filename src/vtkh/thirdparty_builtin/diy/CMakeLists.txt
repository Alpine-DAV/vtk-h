#==============================================================================
# See License.txt
#==============================================================================
add_library(diy INTERFACE)
target_include_directories(diy INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:vtkh/thirdparty_builtin>  # <prefix>/vtkh/thirdparty_builtin/
  )

install(TARGETS diy
  EXPORT ${VTKh_EXPORT_NAME})

# Install headers
install(DIRECTORY include/diy
  DESTINATION ${VTKh_INSTALL_INCLUDE_DIR}/vtkh/thirdparty_builtin)

# Install other files.
install(FILES LEGAL.txt LICENSE.txt
  DESTINATION ${VTKh_INSTALL_SHARED_RESOURCES_DIR}/thirdparty_builtin/diy
  )
