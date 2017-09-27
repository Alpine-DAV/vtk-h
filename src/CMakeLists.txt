#==============================================================================
# See License.txt
#==============================================================================

cmake_minimum_required(VERSION 3.2)

################################
# cmake policy selections
################################
# avoid default mac osx rpath settings for cmake 3.0
cmake_policy(SET CMP0042 OLD)
# avoid warnings for project commands w/o VERSION
cmake_policy(SET CMP0048 NEW)


################################
# VTK-h
################################
project(VTKh)

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
option(ENABLE_MPI         "Build MPI Support"         ON)

################################
# Find VTKm and setup vars based on VTKm's capabilities.
################################
#find_package(VTKm REQUIRED OPTIONAL_COMPONENTS Rendering Serial CUDA TBB)

# These options that affect BLT. We set them up based on what
# capabilities were found in VTKm
cmake_dependent_option(ENABLE_CUDA "Build CUDA Support" ON "VTKm_CUDA_FOUND" OFF)
cmake_dependent_option(ENABLE_TBB "Build TBB Support" ON "VTKm_TBB_FOUND" OFF)

################################
# Set BLT Options
################################
# forward the build shared libs option to blt
set(ENABLE_SHARED_LIBS ${BUILD_SHARED_LIBS} CACHE INTERNAL "")

# don't use BLT's all warnings feature
set(ENABLE_ALL_WARNINGS OFF CACHE INTERNAL "")

################################
# Init BLT
################################
include(blt/SetupBLT.cmake)

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
add_library(vtkh INTERFACE)
target_link_libraries(vtkh INTERFACE vtkh_core vtkh_filters vtkh_rendering)

# Install libraries
install(TARGETS vtkh EXPORT ${VTKh_EXPORT_NAME})

if (ENABLE_MPI)
  add_library(vtkh_par INTERFACE)
  target_link_libraries(vtkh_par INTERFACE vtkh_core_par vtkh_filters_par
    vtkh_rendering_par)
  target_compile_definitions(vtkh_par
    INTERFACE PARALLEL # rename to VTKh_IS_PARALLEL or something like that.
    )

  # Install libraries
  install(TARGETS vtkh_par EXPORT ${VTKh_EXPORT_NAME})
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
