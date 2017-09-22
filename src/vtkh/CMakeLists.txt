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


blt_add_library(
  NAME vtkh
  SOURCES ${Sources}
  HEADERS ${Headers}
  DEPENDS_ON vtkm vtkh_utils
  )

target_include_directories(vtkh
  PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}

  # not convinced these should be exposed, but leaving them for now.
  INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    $<INSTALL_INTERFACE:vtkh>  # <prefix>/vtkh
  )


# Install libraries
install(TARGETS vtkh
  EXPORT ${VTKh_EXPORT_NAME}
  ARCHIVE DESTINATION ${VTKh_INSTALL_LIB_DIR}
  LIBRARY DESTINATION ${VTKh_INSTALL_LIB_DIR}
  RUNTIME DESTINATION ${VTKh_INSTALL_BIN_DIR}
  )

# Install headers
install(FILES ${Headers}
  DESTINATION ${VTKh_INSTALL_INCLUDE_DIR}/vtkh)


#add_subdirectory(filters)
#add_subdirectory(rendering)
