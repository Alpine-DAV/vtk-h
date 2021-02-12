#include "vtkmParticleAdvection.hpp"

#include <vtkm/filter/ParticleAdvection.h>

namespace vtkh
{
vtkm::cont::PartitionedDataSet
vtkmParticleAdvection::Run(vtkm::cont::PartitionedDataSet &input,
                          std::string field_name,
                          double step_size,
                          int numSteps,
                          const std::vector<vtkm::Particle>& seeds)
{
#ifdef VTKH_BYPASS_VTKM_BIH
  return vtkm::cont::PartitionedDataSet();
#else
  auto seedsAH = vtkm::cont::make_ArrayHandle(seeds, vtkm::CopyFlag::Off);

  vtkm::filter::ParticleAdvection particleAdvectionFilter;
  particleAdvectionFilter.SetStepSize(step_size);
  particleAdvectionFilter.SetActiveField(field_name);
  particleAdvectionFilter.SetSeeds(seedsAH);
  particleAdvectionFilter.SetNumberOfSteps(numSteps);
  auto output = particleAdvectionFilter.Execute(input);
  return output;
#endif
}

} // namespace vtkh
