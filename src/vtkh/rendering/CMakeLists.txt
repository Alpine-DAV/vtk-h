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

add_library(vtkh_rendering ${Sources} ${Headers})
target_link_libraries(vtkh_rendering PUBLIC vtkh_core vtkm_rendering)

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
