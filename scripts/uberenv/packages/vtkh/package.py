# Copyright 2013-2019 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack import *

import sys
import os
import socket


import llnl.util.tty as tty
from os import environ as env


def cmake_cache_entry(name, value, vtype=None):
    """
    Helper that creates CMake cache entry strings used in
    'host-config' files.
    """
    if vtype is None:
        if value == "ON" or value == "OFF":
            vtype = "BOOL"
        else:
            vtype = "PATH"
    return 'set({0} "{1}" CACHE {2} "")\n\n'.format(name, value, vtype)


class Vtkh(Package,CudaPackage):
    """VTK-h is a toolkit of scientific visualization algorithms for emerging
    processor architectures. VTK-h brings together several projects like VTK-m
    and DIY2 to provide a toolkit with hybrid parallel capabilities."""

    homepage = "https://github.com/Alpine-DAV/vtk-h"
    git      = "https://github.com/Alpine-DAV/vtk-h.git"
    maintainers = ['cyrush']


    version('ascent_ver', commit='b8504989173666e009d9ceb37a690b8284e42146', submodules=True, preferred=True)
    version('develop', branch='develop', submodules=True)
    version('0.1.0', branch='develop', tag='v0.1.0', submodules=True)

    variant("shared", default=True, description="Build vtk-h as shared libs")
    variant("mpi", default=True, description="build mpi support")
    variant("serial", default=True, description="build serial (non-mpi) libraries")
    variant("tbb", default=False, description="build tbb support")
    variant("cuda", default=False, description="build cuda support")
    variant("openmp", default=(sys.platform != 'darwin'),
            description="build openmp support")

    depends_on("cmake@3.14.1:3.14.5")

    depends_on("mpi", when="+mpi")
    depends_on("intel-tbb", when="@0.1.0+tbb")
    depends_on("cuda", when="+cuda")

    depends_on("vtkm@master~tbb+openmp", when="@develop+openmp")
    depends_on("vtkm@master~tbb~openmp", when="@develop~openmp")

    depends_on("vtkm@master+cuda~tbb+openmp", when="@develop+cuda+openmp")
    depends_on("vtkm@master+cuda~tbb~openmp", when="@develop+cuda~openmp")

    depends_on("vtkm@master~tbb+openmp~shared", when="@develop+openmp~shared")
    depends_on("vtkm@master~tbb~openmp~shared", when="@develop~openmp~shared")

    depends_on("vtkm@master+cuda~tbb+openmp~shared", when="@develop+cuda+openmp~shared")
    depends_on("vtkm@master+cuda~tbb~openmp~shared", when="@develop+cuda~openmp~shared")

    depends_on("vtkm@ascent_ver~tbb+openmp", when="@ascent_ver+openmp")
    depends_on("vtkm@ascent_ver~tbb~openmp", when="@ascent_ver~openmp")

    depends_on("vtkm@ascent_ver+cuda~tbb+openmp", when="@ascent_ver+cuda+openmp")
    depends_on("vtkm@ascent_ver+cuda~tbb~openmp", when="@ascent_ver+cuda~openmp")

    depends_on("vtkm@ascent_ver~tbb+openmp~shared", when="@ascent_ver+openmp~shared")
    depends_on("vtkm@ascent_ver~tbb~openmp~shared", when="@ascent_ver~openmp~shared")

    depends_on("vtkm@ascent_ver+cuda~tbb+openmp~shared", when="@ascent_ver+cuda+openmp~shared")
    depends_on("vtkm@ascent_ver+cuda~tbb~openmp~shared", when="@ascent_ver+cuda~openmp~shared")

    # secretly change build to static when building with cuda to bypass spack variant
    # forwarding crazyness
    #conflicts('+cuda', when='+shared', msg='vtk-h must be built statically (~shared) when cuda is enabled')

    def install(self, spec, prefix):
        with working_dir('spack-build', create=True):
            cmake_args = ["../src",
                          "-DVTKM_DIR={0}".format(spec["vtkm"].prefix),
                          "-DENABLE_TESTS=OFF",
                          "-DBUILD_TESTING=OFF"]

            # shared vs static libs logic
            # force static when building with cuda
            if "+cuda" in spec:
              cmake_args.append('-DBUILD_SHARED_LIBS=OFF')
            else:
              if "+shared" in spec:
                  cmake_args.append('-DBUILD_SHARED_LIBS=ON')
              else:
                  cmake_args.append('-DBUILD_SHARED_LIBS=OFF')

            # mpi support
            if "+mpi" in spec:
                mpicc = spec['mpi'].mpicc
                mpicxx = spec['mpi'].mpicxx
                cmake_args.extend(["-DMPI_C_COMPILER={0}".format(mpicc),
                                   "-DMPI_CXX_COMPILER={0}".format(mpicxx)])
                mpiexe_bin = join_path(spec['mpi'].prefix.bin, 'mpiexec')
                if os.path.isfile(mpiexe_bin):
                    cmake_args.append("-DMPIEXEC={0}".format(mpiexe_bin))

            # openmp support
            if "+openmp" in spec:
                cmake_args.append("-DENABLE_OPENMP=ON")

            # cuda support
            if "+cuda" in spec:
                cmake_args.append("-DVTKm_ENABLE_CUDA:BOOL=ON")
                cmake_args.append("-DENABLE_CUDA:BOOL=ON")
                cmake_args.append("-DCMAKE_CUDA_HOST_COMPILER={0}".format(env["SPACK_CXX"]))
                #if 'cuda_arch' in spec.variants:
                #    cuda_arch = spec.variants['cuda_arch'].value[0]
                #    vtkm_cuda_arch = "native"
                #    arch_map = {"75":"turing", "70":"volta",
                #                "62":"pascal", "61":"pascal", "60":"pascal",
                #                "53":"maxwell", "52":"maxwell", "50":"maxwell",
                #                "35":"kepler", "32":"kepler", "30":"kepler"}
                #    if cuda_arch in arch_map:
                #      vtkm_cuda_arch = arch_map[cuda_arch]
                #    cmake_args.append(
                #        '-DVTKm_CUDA_Architecture={0}'.format(vtkm_cuda_arch))
                #else:
                #    # this fix is necessary if compiling platform has cuda, but
                #    # no devices (this's common for front end nodes on hpc clus
                #    # ters)
                #    # we choose kepler as a lowest common denominator
                #    cmake_args.append("-DVTKm_CUDA_Architecture=native")
            else:
                cmake_args.append("-DVTKm_ENABLE_CUDA:BOOL=OFF")
                cmake_args.append("-DENABLE_CUDA:BOOL=OFF")
            # use release, instead of release with debug symbols b/c vtkh libs
            # can overwhelm compilers with too many symbols
            for arg in std_cmake_args:
                if arg.count("CMAKE_BUILD_TYPE") == 0:
                    cmake_args.extend(std_cmake_args)
            cmake_args.append("-DCMAKE_BUILD_TYPE=Release")
            cmake(*cmake_args)
            make()
            make("install")

            host_cfg_fname = self.create_host_config(spec,
                                                     prefix)

            install(host_cfg_fname, prefix)

    def create_host_config(self, spec, prefix, py_site_pkgs_dir=None):
        """
        This method creates a 'host-config' file that specifies
        all of the options used to configure and build vtkh.
        """

        #######################
        # Compiler Info
        #######################
        c_compiler = env["SPACK_CC"]
        cpp_compiler = env["SPACK_CXX"]

        #######################################################################
        # By directly fetching the names of the actual compilers we appear
        # to doing something evil here, but this is necessary to create a
        # 'host config' file that works outside of the spack install env.
        #######################################################################

        sys_type = spec.architecture
        # if on llnl systems, we can use the SYS_TYPE
        if "SYS_TYPE" in env:
            sys_type = env["SYS_TYPE"]

        ##############################################
        # Find and record what CMake is used
        ##############################################

        cmake_exe = spec['cmake'].command.path

        host_cfg_fname = "%s-%s-%s-vtkh.cmake" % (socket.gethostname(),
                                                  sys_type,
                                                  spec.compiler)

        cfg = open(host_cfg_fname, "w")
        cfg.write("##################################\n")
        cfg.write("# spack generated host-config\n")
        cfg.write("##################################\n")
        cfg.write("# {0}-{1}\n".format(sys_type, spec.compiler))
        cfg.write("##################################\n\n")

        # Include path to cmake for reference
        cfg.write("# cmake from spack \n")
        cfg.write("# cmake executable path: %s\n\n" % cmake_exe)

        #######################
        # Compiler Settings
        #######################

        cfg.write("#######\n")
        cfg.write("# using %s compiler spec\n" % spec.compiler)
        cfg.write("#######\n\n")
        cfg.write("# c compiler used by spack\n")
        cfg.write(cmake_cache_entry("CMAKE_C_COMPILER", c_compiler))
        cfg.write("# cpp compiler used by spack\n")
        cfg.write(cmake_cache_entry("CMAKE_CXX_COMPILER", cpp_compiler))

        # shared vs static libs
        if "+shared" in spec:
            cfg.write(cmake_cache_entry("BUILD_SHARED_LIBS", "ON"))
        else:
            cfg.write(cmake_cache_entry("BUILD_SHARED_LIBS", "OFF"))

        #######################################################################
        # Core Dependencies
        #######################################################################

        #######################
        # VTK-h (and deps)
        #######################

        cfg.write("# vtk-m support \n")

        if "+openmp" in spec:
            cfg.write("# enable openmp support\n")
            cfg.write(cmake_cache_entry("ENABLE_OPENMP", "ON"))

        cfg.write("# vtk-m from spack\n")
        cfg.write(cmake_cache_entry("VTKM_DIR", spec['vtkm'].prefix))

        #######################################################################
        # Optional Dependencies
        #######################################################################

        #######################
        # Serial
        #######################

        if "+serial" in spec:
            cfg.write(cmake_cache_entry("ENABLE_SERIAL", "ON"))
        else:
            cfg.write(cmake_cache_entry("ENABLE_SERIAL", "OFF"))

        #######################
        # MPI
        #######################

        cfg.write("# MPI Support\n")

        if "+mpi" in spec:
            mpicc_path = spec['mpi'].mpicc
            mpicxx_path = spec['mpi'].mpicxx
            mpifc_path = spec['mpi'].mpifc
            # if we are using compiler wrappers on cray systems
            # use those for mpi wrappers, b/c  spec['mpi'].mpicxx
            # etc make return the spack compiler wrappers
            # which can trip up mpi detection in CMake 3.14
            if cpp_compiler == "CC":
                mpicc_path = "cc"
                mpicxx_path = "CC"
                mpifc_path = "ftn"
            cfg.write(cmake_cache_entry("ENABLE_MPI", "ON"))
            cfg.write(cmake_cache_entry("MPI_C_COMPILER", mpicc_path ))
            cfg.write(cmake_cache_entry("MPI_CXX_COMPILER", mpicxx_path ))
            cfg.write(cmake_cache_entry("MPI_Fortran_COMPILER", mpifc_path ))
            mpiexe_bin = join_path(spec['mpi'].prefix.bin, 'mpiexec')
            if os.path.isfile(mpiexe_bin):
                # starting with cmake 3.10, FindMPI expects MPIEXEC_EXECUTABLE
                # vs the older versions which expect MPIEXEC
                if self.spec["cmake"].satisfies('@3.10:'):
                    cfg.write(cmake_cache_entry("MPIEXEC_EXECUTABLE",
                                                mpiexe_bin))
                else:
                    cfg.write(cmake_cache_entry("MPIEXEC",
                                                mpiexe_bin))
        else:
            cfg.write(cmake_cache_entry("ENABLE_MPI", "OFF"))

        #######################
        # CUDA
        #######################

        cfg.write("# CUDA Support\n")

        if "+cuda" in spec:
            cfg.write(cmake_cache_entry("ENABLE_CUDA", "ON"))
            #cfg.write(cmake_cache_entry("VTKm_ENABLE_CUDA","ON"))
            #cfg.write(cmake_cache_entry("CMAKE_CUDA_HOST_COMPILER",''.format(env["SPACK_CXX"])))
            #if 'cuda_arch' in spec.variants:
            #    cuda_arch = spec.variants['cuda_arch'].value[0]
            #    vtkm_cuda_arch = "native"
            #    arch_map = {"75":"turing", "70":"volta",
            #                "62":"pascal", "61":"pascal", "60":"pascal",
            #                "53":"maxwell", "52":"maxwell", "50":"maxwell",
            #                "35":"kepler", "32":"kepler", "30":"kepler"}
            #    if cuda_arch in arch_map:
            #      vtkm_cuda_arch = arch_map[cuda_arch]
            #    cfg.write(cmake_cache_entry('VTKm_CUDA_Architecture','{0}'.format(vtkm_cuda_arch)))
        else:
            cfg.write(cmake_cache_entry("ENABLE_CUDA", "OFF"))

        cfg.write("##################################\n")
        cfg.write("# end spack generated host-config\n")
        cfg.write("##################################\n")
        cfg.close()

        host_cfg_fname = os.path.abspath(host_cfg_fname)
        tty.info("spack generated conduit host-config file: " + host_cfg_fname)
        return host_cfg_fname
