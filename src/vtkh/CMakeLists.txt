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
  DEPENDS_ON vtkm_cont)

#add_subdirectory(utils)
#add_subdirectory(filters)
#add_subdirectory(rendering)
