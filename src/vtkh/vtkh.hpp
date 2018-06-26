#ifndef VTK_H_HPP
#define VTK_H_HPP

#include <string>

namespace vtkh
{

  std::string AboutVTKH();
  bool        IsSerialEnabled();
  bool        IsOpenMPEnabled();
  bool        IsCUDAEnabled();
  bool        IsMPIEnabled();
  
  void        ForceSerial();
  void        ForceOpenMP();
  void        ForceCUDA();
  void        ResetDevices();

  int         GetMPIRank();
  int         GetMPISize();
  
  void        SetMPICommHandle(int mpi_comm_id);
  int         GetMPICommHandle();
}
#endif
