#ifndef VTK_H_VTKM_STREAMLINE_HPP
#define VTK_H_VTKM_STREAMLINE_HPP

/*
#ifdef VTKH_PARALLEL
#include <mpi.h>
#include <vtkh/vtkh.hpp>
#include <diy/mpi.hpp>
#endif
*/

#include <vtkm/Particle.h>
#include <vtkm/cont/PartitionedDataSet.h>

namespace vtkh
{

class vtkmStreamline
{
public:
  vtkm::cont::PartitionedDataSet Run(vtkm::cont::PartitionedDataSet &input,
                                     std::string field_name,
                                     double step_size,
                                     int numSteps,
                                     const std::vector<vtkm::Particle>& seeds);
};
}
#endif
