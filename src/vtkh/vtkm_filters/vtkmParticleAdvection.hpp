#ifndef VTK_H_VTKM_PARTICLE_ADVECTION_HPP
#define VTK_H_VTKM_PARTICLE_ADVECTION_HPP

#include <vtkm/Particle.h>
#include <vtkm/cont/PartitionedDataSet.h>

namespace vtkh
{

class vtkmParticleAdvection
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
