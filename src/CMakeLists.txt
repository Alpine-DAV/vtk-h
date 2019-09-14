#==============================================================================
# See License.txt
#==============================================================================

################################
# VTK-h
################################

################################
# cmake policy selections
################################
# avoid warnings for project commands w/o VERSION
cmake_policy(SET CMP0048 NEW)

project(VTKh VERSION "0.4.0")

set(VTKh_EXPORT_NAME "VTKhTargets")
if (NOT DEFINED VTKh_INSTALL_INCLUDE_DIR)
  set(VTKh_INSTALL_INCLUDE_DIR "include")
endif()
if (NOT DEFINED VTKh_INSTALL_CONFIG_DIR)
  set(VTKh_INSTALL_CONFIG_DIR "lib")
endif()
if (NOT DEFINED VTKh_INSTALL_LIB_DIR)
  set(VTKh_INSTALL_LIB_DIR "lib")
endif()
if (NOT DEFINED VTKh_INSTALL_BIN_DIR)
  set(VTKh_INSTALL_BIN_DIR "bin")
endif()
if (NOT DEFINED VTKh_INSTALL_SHARED_RESOURCES_DIR)
  set(VTKh_INSTALL_SHARED_RESOURCES_DIR "share/vtkh")
endif()
if (NOT DEFINED VTKh_INSTALL_CMAKE_MODULE_DIR)
  set(VTKh_INSTALL_CMAKE_MODULE_DIR "${VTKh_INSTALL_SHARED_RESOURCES_DIR}/cmake")
endif()

set(VTKh_CMAKE_MODULE_PATH ${VTKh_SOURCE_DIR}/cmake)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${VTKh_CMAKE_MODULE_PATH})

################################
# Some CMake module includes.
################################
include(CMakeDependentOption)

# Set default build types
include(VTKhBuildType)

################################
# Build Options
################################
option(BUILD_SHARED_LIBS  "Build shared libraries"    ON)
option(ENABLE_TESTS       "Build tests"               ON)
option(ENABLE_CUDA        "Enable CUDA Support"       OFF)
option(ENABLE_MPI         "Build MPI Support"         ON)
option(ENABLE_LOGGING     "Generate log files "       OFF)
option(ENABLE_SERIAL      "Build serial (non-MPI) libraries"    ON)

if(NOT ENABLE_SERIAL AND NOT ENABLE_MPI)
  message(FATAL_ERROR "No libraries are built. "
    "Please set ENABLE_SERIAL, ENABLE_MPI or both to ON")
endif()


if(ENABLE_CUDA AND BUILD_SHARED_LIBS)
  message(FATAL_ERROR "Static libraries are required when building with CUDA")
endif()
set (CUDA_SEPARABLE_COMPILATION ON CACHE BOOL "" )
# enforce the version after the ENABLE Cuda option
if(ENABLE_CUDA)
  cmake_minimum_required(VERSION 3.14)
else()
  cmake_minimum_required(VERSION 3.9)
endif()
################################
# Set BLT Options
################################

# don't use BLT's all warnings feature
set(ENABLE_ALL_WARNINGS OFF CACHE INTERNAL "")

################################
# Init BLT
################################
if (NOT BLT_CXX_STD)
  set(BLT_CXX_STD "c++11" CACHE STRING "")
endif()

# if BLT_SOURCE_DIR is not set - use "blt" as default
if(NOT BLT_SOURCE_DIR)
    set(BLT_SOURCE_DIR "blt")
endif()

# init blt using BLT_SOURCE_DIR
include(${BLT_SOURCE_DIR}/SetupBLT.cmake)

################################
# Basic CMake Setup
################################
include(cmake/CMakeBasics.cmake)

################################
# Setup 3rd Party Libs
################################
include(cmake/Setup3rdParty.cmake)

################################
# Setup tests helpers
################################
include(cmake/SetupTests.cmake)

################################
# Setup project wide includes
################################
include(cmake/SetupIncludes.cmake)

#------------------------------------------------------------------------------
# Determine source version.
include(VTKhDetermineVersion)
find_package(Git QUIET)
file(STRINGS version.txt version_txt)
extract_version_components("${version_txt}" "VTKh")
determine_version(${VTKh_SOURCE_DIR} "${GIT_EXECUTABLE}" "VTKh")
#------------------------------------------------------------------------------

add_subdirectory(vtkh)

#------------------------------------------------------------------------------
# Add a interface target that makes it easier to depend on VTKh
if (ENABLE_SERIAL)
  add_library(vtkh INTERFACE)
  target_link_libraries(vtkh INTERFACE
                        vtkh_core
                        vtkh_utils
                        vtkh_filters
                        vtkh_rendering)

  # Install libraries
  install(TARGETS vtkh EXPORT ${VTKh_EXPORT_NAME})
endif()

if (ENABLE_MPI)
  add_library(vtkh_mpi INTERFACE)
  target_link_libraries(vtkh_mpi INTERFACE
                        vtkh_core_mpi
                        vtkh_utils_mpi
                        vtkh_filters_mpi
                        vtkh_rendering_mpi)

  # Install libraries
  install(TARGETS vtkh_mpi EXPORT ${VTKh_EXPORT_NAME})
endif()


#------------------------------------------------------------------------------

if(ENABLE_TESTS)
  add_subdirectory(tests)
endif()


################################
# Create CMake importable
# exports for all of our targets
################################
export(EXPORT ${VTKh_EXPORT_NAME}
  FILE ${CMAKE_BINARY_DIR}/${VTKh_INSTALL_CONFIG_DIR}/VTKhTargets.cmake
  )

install(EXPORT ${VTKh_EXPORT_NAME}
  DESTINATION ${VTKh_INSTALL_LIB_DIR}
  FILE VTKhTargets.cmake
  )

include(CMakePackageConfigHelpers)
configure_package_config_file(
  ${VTKh_SOURCE_DIR}/cmake/VTKhConfig.cmake.in
  ${VTKh_BINARY_DIR}/${VTKh_INSTALL_CONFIG_DIR}/VTKhConfig.cmake
  INSTALL_DESTINATION ${VTKh_INSTALL_CONFIG_DIR}
  PATH_VARS
    VTKh_INSTALL_INCLUDE_DIR
    VTKh_INSTALL_CONFIG_DIR
    VTKh_INSTALL_LIB_DIR
    VTKh_INSTALL_BIN_DIR
    VTKh_INSTALL_CMAKE_MODULE_DIR
    VTKh_INSTALL_SHARED_RESOURCES_DIR
    )
write_basic_package_version_file(
  ${VTKh_BINARY_DIR}/${VTKh_INSTALL_CONFIG_DIR}/VTKhConfigVersion.cmake
  VERSION ${VTKh_VERSION}
  COMPATIBILITY ExactVersion)

install(
  FILES
    ${VTKh_BINARY_DIR}/${VTKh_INSTALL_CONFIG_DIR}/VTKhConfig.cmake
    ${VTKh_BINARY_DIR}/${VTKh_INSTALL_CONFIG_DIR}/VTKhConfigVersion.cmake
  DESTINATION ${VTKh_INSTALL_CONFIG_DIR}
  )
