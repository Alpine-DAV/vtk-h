#ifndef VTK_H_HPP
#define VTK_H_HPP

#include <string>

namespace vtkh
{

  std::string AboutVTKH();
  // is backend support compiled in
  bool        IsSerialAvailible();
  bool        IsOpenMPAvailible();
  bool        IsCUDAAvailible();
  // is backend enabled (e.g., ForceX)
  bool        IsSerialEnabled();
  bool        IsOpenMPEnabled();
  bool        IsCUDAEnabled();
  
  bool        IsMPIEnabled();
  
  int         CUDADeviceCount();
  void        SelectCUDADevice(int device_index);
      
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
