#ifndef VTK_H_HPP
#define VTK_H_HPP

#include <string>

#ifdef PARALLEL
#include <mpi.h>
#endif

namespace vtkh
{

  std::string AboutVTKH();
  bool IsSerialEnabled();
  bool IsOpenMPEnabled();
  bool IsCUDAEnabled();
  void ForceSerial();
  void ForceOpenMP();
  void ForceCUDA();
  void ResetDevices();
  int GetMPIRank();
  int GetMPISize();
#ifdef PARALLEL
  void SetMPIComm(MPI_Comm mpi_comm);
  MPI_Comm GetMPIComm();
#endif

}
#endif
