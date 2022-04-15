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
# Setup Kokkos
###############################################################################

if(ENABLE_KOKKOS)
    if(KOKKOS_DIR)
        MESSAGE("Provided KOKKOS_DIR:${KOKKOS_DIR}")

        # check for both lib64 and lib
        if(EXISTS ${KOKKOS_DIR}/lib64/cmake/Kokkos/)
            set(KOKKOS_CMAKE_CONFIG_DIR ${KOKKOS_DIR}/lib64/cmake/Kokkos/)
        endif()

        if(EXISTS ${KOKKOS_DIR}/lib/cmake/Kokkos/)
            set(KOKKOS_CMAKE_CONFIG_DIR ${KOKKOS_DIR}/lib/cmake/Kokkos/)
        endif()

        if(NOT EXISTS ${KOKKOS_CMAKE_CONFIG_DIR}/KokkosConfig.cmake)
            MESSAGE(FATAL_ERROR "Could not find Kokkos CMake include file (${KOKKOS_CMAKE_CONFIG_DIR}/KokkosConfig.cmake)")
        endif()

    ###############################################################################
    # Import Kokkos CMake targets
    ###############################################################################
    find_package(Kokkos REQUIRED
                 NO_DEFAULT_PATH
                 COMPONENTS separable_compilation
                 PATHS ${KOKKOS_CMAKE_CONFIG_DIR})
    else()
	 MESSAGE(FATAL_ERROR "Kokkos support needs explicit KOKKOS_DIR")
    endif()	 
endif()
