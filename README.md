# VTK-h is now developed as part of Ascent

VTK-h is now developed as part of Ascent (https://github.com/Alpine-DAV/ascent).

It was moved as part of a great amalgamation to simplify ongoing development and support of Ascent.


VTK-h
=====

VTK-h is a toolkit of scientific visualization algorithms for emerging processor architectures. VTK-h
bring together several projects like VTK-m and DIY2 to provide a toolkit with hybrid parallel capabilities.

VTK-m Version:
==============
VTK-h use a specific version of VTK-m.
For the current version see:

https://github.com/Alpine-DAV/vtk-h/blob/develop/hashes.txt

Additionally, VTK-m has several build settings and we test with a specific combination.
When building VTK-m, be sure to set:

- `VTKm_USE_64BIT_IDS:BOOL=OFF`
- `VTKm_USE_DOUBLE_PRECISION:BOOL=ON`
- `VTKm_USE_DEFAULT_TYPES_FOR_ASCENT:BOOL=ON`
- `VTKm_NO_DEPRECATED_VIRTUAL:BOOL=ON`


Source Repo
=================

VTK-h's source is hosted on GitHub:

https://github.com/Alpine-DAV/vtk-h

License
===========

VTK-h is released under [3-clause BSD license](/LICENSE).
