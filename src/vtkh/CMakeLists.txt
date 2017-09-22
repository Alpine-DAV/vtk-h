#==============================================================================
# See License.txt
#==============================================================================


set(Headers
  vtkh.hpp
  DataSet.hpp
  Error.hpp)

set(Sources
  vtkh.cpp
  DataSet.cpp)


blt_add_library(
  NAME vtkh
  SOURCES ${Sources}
  HEADERS ${Headers}
  DEPENDS_ON vtkm
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

#add_subdirectory(utils)
#add_subdirectory(filters)
#add_subdirectory(rendering)
