#==============================================================================
# See License.txt
#==============================================================================

########
# TODO: blt_register library support interface libs 
########

add_library(diy INTERFACE)
target_include_directories(diy INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:${VTKh_INSTALL_INCLUDE_DIR}/vtkh/thirdparty_builtin>  # <prefix>/vtkh/thirdparty_builtin/
  )


########
# TODO: We should have a macro for install targets + headers
########


install(TARGETS diy
  EXPORT ${VTKh_EXPORT_NAME})

# Install headers
install(DIRECTORY include/diy
  DESTINATION ${VTKh_INSTALL_INCLUDE_DIR}/vtkh/thirdparty_builtin)

# Install other files.
install(FILES LEGAL.txt LICENSE.txt
  DESTINATION ${VTKh_INSTALL_SHARED_RESOURCES_DIR}/thirdparty_builtin/diy
  )
