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

project(VTKh VERSION "0.7.0")

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
option(ENABLE_KOKKOS      "Build KOKKOS Support"      OFF)
option(ENABLE_LOGGING     "Generate log files "       OFF)
option(ENABLE_SERIAL      "Build serial (non-MPI) libraries" ON)
option(ENABLE_FILTER_CONTOUR_TREE "Build contour tree support" OFF)

if(NOT ENABLE_SERIAL AND NOT ENABLE_MPI)
  message(FATAL_ERROR "No libraries are built. "
    "Please set ENABLE_SERIAL, ENABLE_MPI or both to ON")
endif()

if(ENABLE_CUDA)
    # ask politely a few times for SEPARABLE COMPILATION
    set (CMAKE_CUDA_SEPARABLE_COMPILATION ON CACHE BOOL "" )
    set (CUDA_SEPARABLE_COMPILATION ON CACHE BOOL "" )
endif()

# enforce the version after the ENABLE Cuda option
if(ENABLE_CUDA)
  cmake_minimum_required(VERSION 3.14)
else()
  cmake_minimum_required(VERSION 3.9)
endif()

################################
# cmake policy selection
################################
# cmake 3.18+ use CMP0104 OLD to manually specify cuda flags
#https://cmake.org/cmake/help/latest/policy/CMP0104.html
if(POLICY CMP0104)
    cmake_policy(SET CMP0104 OLD)
endif()

################################
# Init BLT
################################
# This also includes 
# our BLT defaults
include(cmake/SetupBLT.cmake)

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
# include(VTKhDetermineVersion)
# find_package(Git QUIET)
# file(STRINGS version.txt version_txt)
# extract_version_components("${version_txt}" "VTKh")
# determine_version(${VTKh_SOURCE_DIR} "${GIT_EXECUTABLE}" "VTKh")
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
# Add examples
################################
if(ENABLE_EXAMPLES)
    add_subdirectory(examples)
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
