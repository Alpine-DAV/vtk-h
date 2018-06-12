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

################################
# Basic TPL Tests
################################

list(APPEND VTKM_TESTS t_vtkm_smoke)
set(MPI_TESTS t_diy_smoke)

################################
# Add our Main Unit Tests
################################
message(STATUS "Adding thirdparty lib unit tests")
foreach(TEST ${TPL_TESTS})
    add_cpp_test(TEST ${TEST})
endforeach()

################################
# Add our VTKm Unit Tests
################################
foreach(TEST ${VTKM_TESTS})
  add_cpp_test(TEST ${TEST}_Serial SOURCE ${TEST} DEPENDS_ON vtkm)
  target_compile_definitions(${TEST}_Serial PRIVATE "VTKM_DEVICE_ADAPTER=VTKM_DEVICE_ADAPTER_SERIAL")
  if (ENABLE_OPENMP)
    add_cpp_test(TEST ${TEST}_OPENMP SOURCE ${TEST} DEPENDS_ON vtkm)
    target_compile_definitions(${TEST}_OPENMP PRIVATE "VTKM_DEVICE_ADAPTER=VTKM_DEVICE_ADAPTER_OPENMP")
  endif()
  if (ENABLE_CUDA)
    add_cuda_test(TEST ${TEST}_CUDA SOURCE ${TEST} DEPENDS_ON vtkm)
    target_compile_definitions(${TEST}_CUDA PRIVATE "VTKM_DEVICE_ADAPTER=VTKM_DEVICE_ADAPTER_CUDA")
  endif()
endforeach()

################################
# Add our Optional Unit Tests
################################
if(MPI_FOUND)
  message(STATUS "MPI enabled: Adding related unit tests")
  foreach(TEST ${MPI_TESTS})
    # this uses 2 procs
    add_cpp_mpi_test(TEST ${TEST} NUM_PROCS 2)
    target_link_libraries(${TEST} PRIVATE diy)
  endforeach()
else()
  message(STATUS "MPI disabled: Skipping related tests")
endif()


