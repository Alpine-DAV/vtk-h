#include "vtkh.hpp"
#include "Error.hpp"
#include <vtkh/Logger.hpp>

#include <vtkm/cont/RuntimeDeviceInformation.h>
#include <vtkm/cont/RuntimeDeviceTracker.h>
#include <vtkm/cont/DeviceAdapterListTag.h>

#ifdef VTKM_CUDA
#include <cuda.h>
#endif

#include <sstream>

#ifdef VTKH_PARALLEL
#include <mpi.h>
#endif

namespace vtkh
{

static int g_mpi_comm_id = -1;


//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
#ifdef VTKH_PARALLEL // mpi case
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

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
#ifdef VTKH_ENABLE_LOGGING
  DataLogger::GetInstance()->SetRank(GetMPIRank());
#endif
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
  int rank;
  MPI_Comm comm = MPI_Comm_f2c(GetMPICommHandle());
  MPI_Comm_rank(comm, &rank);
  return rank;
}

//---------------------------------------------------------------------------//
int
GetMPISize()
{
  int size;
  MPI_Comm comm = MPI_Comm_f2c(GetMPICommHandle());
  MPI_Comm_size(comm, &size);
  return size;
}

//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//
#else // non-mpi case
//---------------------------------------------------------------------------//
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
void
CheckCommHandle()
{
  std::stringstream msg;
  msg<<"VTK-h internal error. Trying to access MPI comm in non-mpi vtkh lib.";
  msg<<"Are you using the right library (vtkh vs vtkh_mpi)?";
  throw Error(msg.str());
}

//---------------------------------------------------------------------------//
void
SetMPICommHandle(int mpi_comm_id)
{
  std::stringstream msg;
  msg<<"VTK-h internal error. Trying to set MPI comm handle in non-mpi vtkh lib.";
  msg<<"Are you using the right library (vtkh vs vtkh_mpi)?";
  throw Error(msg.str());
}

//---------------------------------------------------------------------------//
int
GetMPICommHandle()
{
  std::stringstream msg;
  msg<<"VTK-h internal error. Trying to get MPI comm handle in non-mpi vtkh lib.";
  msg<<"Are you using the right library (vtkh vs vtkh_mpi)?";
  throw Error(msg.str());
  return g_mpi_comm_id;
}

//---------------------------------------------------------------------------//
int
GetMPIRank()
{
  return 0;
}

//---------------------------------------------------------------------------//
int
GetMPISize()
{
  return 1;
}
//---------------------------------------------------------------------------//
#endif
//---------------------------------------------------------------------------//

//---------------------------------------------------------------------------//
bool
IsMPIEnabled()
{
#ifdef VTKH_PARALLEL
  return true;
#else
  return false;
#endif
}

std::string GetCurrentDevice()
{
  std::string device = "serial";
  // use the same prefered ordering as vtkm
  if(IsCUDAEnabled())
  {
    device = "cuda";
  }
  else if(IsOpenMPEnabled())
  {
    device = "openmp";
  }

  return device;
}

//---------------------------------------------------------------------------//
bool
IsSerialAvailable()
{
  vtkm::cont::RuntimeDeviceInformation info;
  return info.Exists(vtkm::cont::DeviceAdapterTagSerial());
}


//---------------------------------------------------------------------------//
bool
IsOpenMPAvailable()
{
  vtkm::cont::RuntimeDeviceInformation info;
  return info.Exists(vtkm::cont::DeviceAdapterTagOpenMP());
}

//---------------------------------------------------------------------------//
bool
IsCUDAAvailable()
{
  vtkm::cont::RuntimeDeviceInformation info;
  return info.Exists(vtkm::cont::DeviceAdapterTagCuda());
}

//---------------------------------------------------------------------------//
bool
IsSerialEnabled()
{
  vtkm::cont::RuntimeDeviceTracker &device_tracker
    = vtkm::cont::GetRuntimeDeviceTracker();
  return device_tracker.CanRunOn(vtkm::cont::DeviceAdapterTagSerial());
}


//---------------------------------------------------------------------------//
bool
IsOpenMPEnabled()
{
  vtkm::cont::RuntimeDeviceTracker &device_tracker
    = vtkm::cont::GetRuntimeDeviceTracker();
  return device_tracker.CanRunOn(vtkm::cont::DeviceAdapterTagOpenMP());
}

//---------------------------------------------------------------------------//
bool
IsCUDAEnabled()
{
  vtkm::cont::RuntimeDeviceTracker &device_tracker
    = vtkm::cont::GetRuntimeDeviceTracker();
  return device_tracker.CanRunOn(vtkm::cont::DeviceAdapterTagCuda());
}

//---------------------------------------------------------------------------//
int
CUDADeviceCount()
{
    int device_count = 0;
#ifdef VTKM_CUDA
    cudaError_t res = cudaGetDeviceCount(&device_count);
    if(res != cudaSuccess)
    {
        std::stringstream msg;
        msg << "Failed to get CUDA device count" << std::endl
            << "CUDA Error Message: "
            << cudaGetErrorString(res);
        throw Error(msg.str());
    }
    return device_count;

#else
    throw Error("Cannot fetch CUDA device count: VTK-m lacks CUDA support");
#endif
    return device_count;
}

//---------------------------------------------------------------------------//
void
SelectCUDADevice(int device_index)
{
#ifdef VTKM_CUDA
    int device_count = CUDADeviceCount();
    // make sure index is ok
    if(device_index >= 0 && device_index < device_count)
    {
        cudaError_t res = cudaSetDevice(device_index);
        if(res != cudaSuccess)
        {
            std::stringstream msg;
            msg << "Failed to set CUDA device (device index = "
                << device_index << ")" << std::endl
                << "CUDA Error Message: "
                << cudaGetErrorString(res);
            throw Error(msg.str());
        }
    }
    else
    {
        std::stringstream msg;
        msg << "Invalid CUDA device index: "
            << device_index
            << " (number of devices = "
            << device_index << ")";
        throw Error(msg.str());
    }
#else
    throw Error("Cannot set CUDA device: VTK-m lacks CUDA support");
#endif

}

//---------------------------------------------------------------------------//
void
ForceSerial()
{
  vtkm::cont::RuntimeDeviceTracker &device_tracker
    = vtkm::cont::GetRuntimeDeviceTracker();
  device_tracker.ForceDevice(vtkm::cont::DeviceAdapterTagSerial());
}

//---------------------------------------------------------------------------//
void
ForceOpenMP()
{
  vtkm::cont::RuntimeDeviceTracker &device_tracker
    = vtkm::cont::GetRuntimeDeviceTracker();
  device_tracker.ForceDevice(vtkm::cont::DeviceAdapterTagOpenMP());
}

//---------------------------------------------------------------------------//
void
ForceCUDA()
{
  vtkm::cont::RuntimeDeviceTracker &device_tracker
    = vtkm::cont::GetRuntimeDeviceTracker();
  device_tracker.ForceDevice(vtkm::cont::DeviceAdapterTagCuda());
}

//---------------------------------------------------------------------------//
void
ResetDevices()
{
  vtkm::cont::RuntimeDeviceTracker &device_tracker
    = vtkm::cont::GetRuntimeDeviceTracker();
  device_tracker.Reset();
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

  if(IsCUDAAvailable())
  {
    msg<<"Cuda (";
    if(IsCUDAEnabled())
    {
      msg<<"enabled) ";
    }
    else
    {
      msg<<"disabled) ";
    }

  }

  if(IsOpenMPAvailable())
  {
    msg<<"OpenMP (";
    if(IsOpenMPEnabled())
    {
      msg<<"enabled) ";
    }
    else
    {
      msg<<"disabled) ";
    }
  }

  if(IsSerialAvailable())
  {
    msg<<"Serial (";
    if(IsSerialEnabled())
    {
      msg<<"enabled) ";
    }
    else
    {
      msg<<"disabled) ";
    }
  }
  msg<<"\n";
 msg<<"------------------------------------------\n";
  return msg.str();
}

}
