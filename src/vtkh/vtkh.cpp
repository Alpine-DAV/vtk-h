#include "vtkh.hpp"
#include "Error.hpp"

#include <vtkm/cont/RuntimeDeviceInformation.h>
#include <vtkm/cont/RuntimeDeviceTracker.h>
#include <vtkm/cont/DeviceAdapterListTag.h>
#include <sstream>

#ifdef VTKH_PARALLEL
#include <mpi.h>
#include <vtkh/utils/vtkh_mpi_utils.hpp>
#endif

namespace vtkh
{

static int g_mpi_comm_id = -1;

#ifdef VTKH_PARALLEL // mpi case


//---------------------------------------------------------------------------//
void
CheckCommHandle()
{
  if(g_mpi_comm_id == -1)
  {
    std::stringstream msg;
    msg<<"VTK-h internal error. There is no valid MPI comm available. ";
    msg<<"It is likely that VTKH.SetMPICommHandle(int) was not called.";
    throw Error(msg.str());
  }
}

//---------------------------------------------------------------------------//
void 
SetMPICommHandle(int mpi_comm_id)
{
  g_mpi_comm_id = mpi_comm_id;
}

//---------------------------------------------------------------------------//
int
GetMPICommHandle()
{
  CheckCommHandle();
  return g_mpi_comm_id;
}

//---------------------------------------------------------------------------//
int 
GetMPIRank()
{
  CheckCommHandle();
  int rank;
  MPI_Comm comm = GetMPIComm(); 
  MPI_Comm_rank(comm, &rank);
  return rank;
}

//---------------------------------------------------------------------------//
int 
GetMPISize()
{
  CheckCommHandle();
  int size;
  MPI_Comm comm = GetMPIComm(); 
  MPI_Comm_size(comm, &size);
  return size;
}

#else // non-mpi case

//---------------------------------------------------------------------------//
void
CheckCommHandle()
{
  std::stringstream msg;
  msg<<"VTK-h internal error. Trying to access MPI comm in non-mpi vtkh lib";
  msg<<"Are you using the right library (vtkh vs vtkh_mpi)?";
  throw Error(msg.str());
}

//---------------------------------------------------------------------------//
void 
SetMPICommHandle(int mpi_comm_id)
{
  std::stringstream msg;
  msg<<"VTK-h internal error. Trying to set MPI comm handle in non-mpi vtkh lib";
  msg<<"Are you using the right library (vtkh vs vtkh_mpi)?";
  throw Error(msg.str());
}

//---------------------------------------------------------------------------//
int 
GetMPICommHandle()
{
  std::stringstream msg;
  msg<<"VTK-h internal error. Trying to get MPI comm handle in non-mpi vtkh lib";
  msg<<"Are you using the right library (vtkh vs vtkh_mpi)?";
  throw Error(msg.str());
  return g_mpi_comm_id;
}

//---------------------------------------------------------------------------//
int 
GetMPIRank()
{
  return 1;
}

//---------------------------------------------------------------------------//
int 
GetMPISize()
{
  return 1;
}
#endif

//---------------------------------------------------------------------------//
bool
IsSerialEnabled()
{
  vtkm::cont::RuntimeDeviceInformation<vtkm::cont::DeviceAdapterTagSerial> serial;
  return serial.Exists();
}

//---------------------------------------------------------------------------//
bool
IsOpenMPEnabled()
{
  vtkm::cont::RuntimeDeviceInformation<vtkm::cont::DeviceAdapterTagOpenMP> omp;
  return omp.Exists();
}

//---------------------------------------------------------------------------//
bool
IsCUDAEnabled()
{
  vtkm::cont::RuntimeDeviceInformation<vtkm::cont::DeviceAdapterTagCuda> cuda;
  return cuda.Exists();
}

//---------------------------------------------------------------------------//
void
ForceSerial()
{
  vtkm::cont::RuntimeDeviceTracker global_tracker;
  global_tracker = vtkm::cont::GetGlobalRuntimeDeviceTracker();
  global_tracker.ForceDevice(vtkm::cont::DeviceAdapterTagSerial());
}

//---------------------------------------------------------------------------//
void
ForceOpenMP()
{
  vtkm::cont::RuntimeDeviceTracker global_tracker;
  global_tracker = vtkm::cont::GetGlobalRuntimeDeviceTracker();
  global_tracker.ForceDevice(vtkm::cont::DeviceAdapterTagOpenMP());
}

//---------------------------------------------------------------------------//
void
ForceCUDA()
{
  vtkm::cont::RuntimeDeviceTracker global_tracker;
  global_tracker = vtkm::cont::GetGlobalRuntimeDeviceTracker();
  global_tracker.ForceDevice(vtkm::cont::DeviceAdapterTagCuda());
}

//---------------------------------------------------------------------------//
void
ResetDevices()
{
  vtkm::cont::RuntimeDeviceTracker global_tracker;
  global_tracker = vtkm::cont::GetGlobalRuntimeDeviceTracker();
  global_tracker.Reset();
}

//---------------------------------------------------------------------------//
std::string
AboutVTKH()
{
  std::stringstream msg;
  msg<<"---------------- VTK-h -------------------\n";
#ifdef VTKH_PARALLEL
  int version, subversion;
  MPI_Get_version(&version, &subversion);
  msg<<"MPI version: "<<version<<"."<<subversion<<"\n";
#else
  msg<<"MPI version: n/a\n";
#endif
  msg<<"VTK-m adapters: ";

  if(IsCUDAEnabled())
  {
    msg<<"Cuda ";
  }

  if(IsOpenMPEnabled())
  {
    msg<<"OpenMP ";
  }

  if(IsSerialEnabled())
  {
    msg<<"Serial ";
  }
  msg<<"\n";
 msg<<"------------------------------------------\n";
  return msg.str();
}

}
