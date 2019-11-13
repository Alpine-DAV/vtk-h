# Copyright 2013-2019 Lawrence Livermore National Security, LLC and other
# Spack Project Developers. See the top-level COPYRIGHT file for details.
#
# SPDX-License-Identifier: (Apache-2.0 OR MIT)


from spack import *
import os


class Vtkm(CMakePackage, CudaPackage):
    """VTK-m is a toolkit of scientific visualization algorithms for emerging
    processor architectures. VTK-m supports the fine-grained concurrency for
    data analysis and visualization algorithms required to drive extreme scale
    computing by providing abstract models for data and execution that can be
    applied to a variety of algorithms across many different processor
    architectures."""

    homepage = "https://m.vtk.org/"
    url      = "https://gitlab.kitware.com/api/v4/projects/vtk%2Fvtk-m/repository/archive.tar.gz?sha=v1.3.0"
    git      = "https://gitlab.kitware.com/vtk/vtk-m.git"

    # overrider master to pin
    # version used for ascent
    version('master', tag='v1.5.0', preferred=True)
    #version('master', branch='master')
    version('1.3.0', "d9f6e274dec2ea01273cccaba356d23ca88c5a25")
    version('1.2.0', "3295fed86012226c107e1f2605ca7cc583586b63")
    version('1.1.0', "6aab1c0885f6ffaaffcf07930873d0df")

    # use release, instead of release with debug symbols b/c vtkm libs
    # can overwhelm compilers with too many symbols
    variant('build_type', default='Release', description='CMake build type',
            values=('Debug', 'Release', 'RelWithDebInfo', 'MinSizeRel'))
    variant("shared", default=True, description="build shared libs")
    variant("cuda", default=False, description="build cuda support")
    variant("doubleprecision", default=True,
            description='enable double precision')
    variant("logging", default=False, description="build logging support")
    variant("mpi", default=False, description="build mpi support")
    variant("openmp", default=False, description="build openmp support")
    variant("rendering", default=True, description="build rendering support")
    variant("tbb", default=True, description="build TBB support")
    variant("64bitids", default=False,
            description="enable 64 bits ids")

    depends_on("cmake")

    depends_on("tbb", when="+tbb")
    depends_on("cuda", when="+cuda")
    depends_on("mpi", when="+mpi")

    # secretly change build to static when building with cuda to bypass spack variant
    # forwarding crazyness
    #conflicts('+cuda', when='+shared', msg='vtk-m must be built statically (~shared) when cuda is enabled')

    def cmake_args(self):
        spec = self.spec
        options = []
        with working_dir('spack-build', create=True):
            options = ["../",
                       "-DVTKm_ENABLE_TESTING:BOOL=OFF"]
            # shared vs static libs logic
            # force building statically with cuda
            if "+cuda" in spec:
                options.append('-DBUILD_SHARED_LIBS=OFF')
            else:
              if "+shared" in spec:
                  options.append('-DBUILD_SHARED_LIBS=ON')
              else:
                  options.append('-DBUILD_SHARED_LIBS=OFF')
            # cuda support
            if "+cuda" in spec:
                options.append("-DVTKm_ENABLE_CUDA:BOOL=ON")
                options.append("-DCMAKE_CUDA_HOST_COMPILER={0}".format(env["SPACK_CXX"]))
                if 'cuda_arch' in spec.variants:
                    cuda_arch = spec.variants['cuda_arch'].value[0]
                    vtkm_cuda_arch = "native"
                    arch_map = {"75":"turing", "70":"volta",
                                "62":"pascal", "61":"pascal", "60":"pascal",
                                "53":"maxwell", "52":"maxwell", "50":"maxwell",
                                "35":"kepler", "32":"kepler", "30":"kepler"}
                    if cuda_arch in arch_map:
                      vtkm_cuda_arch = arch_map[cuda_arch]
                    options.append(
                        '-DVTKm_CUDA_Architecture={0}'.format(vtkm_cuda_arch))
            else:
                options.append("-DVTKm_ENABLE_CUDA:BOOL=OFF")

            # double precision
            if "+doubleprecision" in spec:
                options.append("-DVTKm_USE_DOUBLE_PRECISION:BOOL=ON")
            else:
                options.append("-DVTKm_USE_DOUBLE_PRECISION:BOOL=OFF")

            # logging support
            if "+logging" in spec:
                if spec.satisfies('@:1.2.0') and \
                        spec['vtkm'].version.string != 'master' and \
                        spec['vtkm'].version.string != 'ascent_ver' :
                    raise InstallError('logging is not supported for\
                            vtkm version lower than 1.3')
                options.append("-DVTKm_ENABLE_LOGGING:BOOL=ON")
            else:
                options.append("-DVTKm_ENABLE_LOGGING:BOOL=OFF")

            # mpi support
            if "+mpi" in spec:
                if spec.satisfies('@:1.2.0') and \
                        spec['vtkm'].version.string != 'master' and \
                        spec['vtkm'].version.string != 'ascent_ver':
                    raise InstallError('mpi is not supported for\
                            vtkm version lower than 1.3')
                options.append("-DVTKm_ENABLE_MPI:BOOL=ON")
            else:
                options.append("-DVTKm_ENABLE_MPI:BOOL=OFF")

            # openmp support
            if "+openmp" in spec:
                # openmp is added since version 1.3.0
                if spec.satisfies('@:1.2.0') and \
                        spec['vtkm'].version.string != 'master' and \
                        spec['vtkm'].version.string != 'ascent_ver':
                    raise InstallError('OpenMP is not supported for\
                            vtkm version lower than 1.3')
                options.append("-DVTKm_ENABLE_OPENMP:BOOL=ON")
            else:
                options.append("-DVTKm_ENABLE_OPENMP:BOOL=OFF")

            # rendering support
            if "+rendering" in spec:
                options.append("-DVTKm_ENABLE_RENDERING:BOOL=ON")
            else:
                options.append("-DVTKm_ENABLE_RENDERING:BOOL=OFF")

            # tbb support
            if "+tbb" in spec:
                # vtk-m detectes tbb via TBB_ROOT env var
                os.environ["TBB_ROOT"] = spec["tbb"].prefix
                options.append("-DVTKm_ENABLE_TBB:BOOL=ON")
            else:
                options.append("-DVTKm_ENABLE_TBB:BOOL=OFF")

            # 64 bit ids
            if "+64bitids" in spec:
                options.append("-DVTKm_USE_64BIT_IDS:BOOL=ON")
                print("64 bit ids enabled")
            else:
                options.append("-DVTKm_USE_64BIT_IDS:BOOL=OFF")
            return options
