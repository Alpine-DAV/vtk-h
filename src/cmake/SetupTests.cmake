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

##------------------------------------------------------------------------------
## - Builds and adds a test that uses gtest
##
## add_cpp_test( TEST test [SOURCE filename] DEPENDS_ON dep1 dep2... )
##------------------------------------------------------------------------------
function(add_cpp_test)

    set(options)
    set(singleValueArgs TEST SOURCE)
    set(multiValueArgs DEPENDS_ON)

    # parse our arguments
    cmake_parse_arguments(arg
                         "${options}"
                         "${singleValueArgs}"
                         "${multiValueArgs}" ${ARGN} )

    message(STATUS " [*] Adding Unit Test: ${arg_TEST}")

    if (NOT arg_SOURCE)
      set(arg_SOURCE ${arg_TEST})
    endif()

    blt_add_executable( NAME ${arg_TEST}
                        SOURCES ${arg_SOURCE}.cpp ${fortran_driver_source}
                        OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}
                        DEPENDS_ON "${arg_DEPENDS_ON}" gtest)

    blt_add_test( NAME ${arg_TEST}
                  COMMAND ${arg_TEST}
                    )

endfunction()


##------------------------------------------------------------------------------
## - Builds and adds a test that uses gtest
##
## add_cuda_test( TEST test DEPENDS_ON dep1 dep2... )
##------------------------------------------------------------------------------
function(add_cuda_test)

    set(options)
    set(singleValueArgs TEST)
    set(multiValueArgs DEPENDS_ON)

    # parse our arguments
    cmake_parse_arguments(arg
                         "${options}" 
                         "${singleValueArgs}" 
                         "${multiValueArgs}" ${ARGN} )

    message(STATUS " [*] Adding CUDA Unit Test: ${arg_TEST}")

    blt_add_executable( NAME ${arg_TEST}
                        SOURCES ${arg_TEST}.cpp ${fortran_driver_source}
                        OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}
                        DEPENDS_ON "${arg_DEPENDS_ON}" gtest cuda)

    blt_add_test( NAME ${arg_TEST}
                  COMMAND ${arg_TEST}
                    )

endfunction()


##------------------------------------------------------------------------------
## - Builds and adds a test that uses gtest and mpi
##
## add_cpp_mpi_test( TEST test NUM_PROCS 2 DEPENDS_ON dep1 dep2... )
##------------------------------------------------------------------------------
function(add_cpp_mpi_test)

    set(options)
    set(singleValueArgs TEST NUM_PROCS)
    set(multiValueArgs DEPENDS_ON)

    # parse our arguments
    cmake_parse_arguments(arg
                         "${options}" 
                         "${singleValueArgs}" 
                         "${multiValueArgs}" ${ARGN} )

    message(STATUS " [*] Adding Unit Test: ${arg_TEST}")


    blt_add_executable( NAME ${arg_TEST}
                        SOURCES ${arg_TEST}.cpp ${fortran_driver_source}
                        OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}
                        DEPENDS_ON "${arg_DEPENDS_ON}" gtest mpi)

    blt_add_test( NAME ${arg_TEST}
                  COMMAND ${arg_TEST}
                  NUM_MPI_TASKS ${arg_NUM_PROCS})

endfunction()


##------------------------------------------------------------------------------
## - Adds a python based unit test
##
## add_python_test( TEST test)
##------------------------------------------------------------------------------
function(add_python_test TEST)
    message(STATUS " [*] Adding Python-based Unit Test: ${TEST}")
    add_test(NAME ${TEST} COMMAND 
             ${PYTHON_EXECUTABLE} -B -m unittest -v ${TEST})
    # make sure python can pick up the modules we built
    set(PYTHON_TEST_PATH "${CMAKE_BINARY_DIR}/python-modules/:${CMAKE_CURRENT_SOURCE_DIR}")
    if(EXTRA_PYTHON_MODULE_DIRS)
        set(PYTHON_TEST_PATH "${EXTRA_PYTHON_MODULE_DIRS}:${PYTHON_TEST_PATH}")
    endif()
    set_property(TEST ${TEST} PROPERTY ENVIRONMENT  "PYTHONPATH=${PYTHON_TEST_PATH}")
endfunction(add_python_test)


##------------------------------------------------------------------------------
## - Adds a fortran based unit test
##
## add_fortran_test( TEST test DEPENDS_ON dep1 dep2... )
##------------------------------------------------------------------------------
macro(add_fortran_test)
    set(options)
    set(singleValueArgs TEST)
    set(multiValueArgs DEPENDS_ON)

    # parse our arguments
    cmake_parse_arguments(arg
                         "${options}" 
                         "${singleValueArgs}" 
                         "${multiValueArgs}" ${ARGN} )

    message(STATUS " [*] Adding Fortran Unit Test: ${arg_TEST}")
    blt_add_executable( NAME ${arg_TEST}
                        SOURCES ${arg_TEST}.f
                        OUTPUT_DIR ${CMAKE_CURRENT_BINARY_DIR}
                        DEPENDS_ON fruit "${arg_DEPENDS_ON}")

    blt_add_test( NAME ${arg_TEST}
                  COMMAND ${arg_TEST})

endmacro(add_fortran_test)
