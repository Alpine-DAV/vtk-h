###############################################################################
# Copyright (c) 2015-2017, Lawrence Livermore National Security, LLC.
#
# Produced at the Lawrence Livermore National Laboratory
#
# LLNL-CODE-716457
#
# All rights reserved.
#
# This file is part of Ascent.
#
# For details, see: http://software.llnl.gov/ascent/.
#
# Please also read ascent/LICENSE
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the disclaimer below.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the disclaimer (as noted below) in the
#   documentation and/or other materials provided with the distribution.
#
# * Neither the name of the LLNS/LLNL nor the names of its contributors may
#   be used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY,
# LLC, THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY
# DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
# STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
# IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
###############################################################################

###############################################################################
# Setup VTKm
###############################################################################

if(NOT VTKM_DIR)
    MESSAGE(FATAL_ERROR "VTKm support needs explicit VTKM_DIR")
endif()

MESSAGE(STATUS "Looking for VTKm using VTKM_DIR = ${VTKM_DIR}")

# use VTKM_DIR to setup the options that cmake's find VTKm needs
file(GLOB VTKm_DIR "${VTKM_DIR}/lib/cmake/vtkm-*")
if(NOT VTKm_DIR)
    MESSAGE(FATAL_ERROR "Failed to find VTKm at VTKM_DIR=${VTKM_DIR}/lib/cmake/vtk-*")
endif()

find_package(VTKm REQUIRED QUIET)

if(ENABLE_CUDA AND NOT VTKm_ENABLE_CUDA)
   message(FATAL_ERROR "VTK-h CUDA support requires VTK-m with CUDA support (ENABLE_CUDA == TRUE, however VTKm_ENABLE_CUDA == FALSE")
endif()


set(VTKM_FOUND TRUE)

set(VTKM_TARGETS vtkm_cont vtkm_filter vtkm_rendering)

if(ENABLE_CUDA)
    # we need to inject the vtkm cuda flags into CMAKE_CUDA_FLAGS
    vtkm_get_cuda_flags(_fetch_vtkm_cuda_flags)
    set(CMAKE_CUDA_FLAGS  "${CMAKE_CUDA_FLAGS} ${_fetch_vtkm_cuda_flags}")
    unset(_fetch_vtkm_cuda_flags)
endif()

# VTKM does not seem to propogate includes it exposes to us, so we have to work
# around this.
file(GLOB LCL_DIR "${VTKM_DIR}/include/vtkm-*/vtkm/thirdparty/lcl/vtkmlcl/")
include_directories("${LCL_DIR}")

blt_register_library(NAME vtkm
                     LIBRARIES ${VTKM_TARGETS}
                     )
