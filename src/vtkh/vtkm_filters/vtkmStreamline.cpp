#include "vtkmStreamline.hpp"

#include <vtkm/filter/Streamline.h>

namespace vtkh
{
vtkm::cont::PartitionedDataSet
vtkmStreamline::Run(vtkm::cont::PartitionedDataSet &input,
                    std::string field_name,
                    double step_size,
                    int numSteps,
                    const std::vector<vtkm::Particle>& seeds)
{
#ifdef VTKH_BYPASS_VTKM_BIH
  return vtkm::cont::PartitionedDataSet();
#else
  auto seedsAH = vtkm::cont::make_ArrayHandle(seeds, vtkm::CopyFlag::Off);

  vtkm::filter::Streamline streamlineFilter;
  streamlineFilter.SetStepSize(step_size);
  streamlineFilter.SetActiveField(field_name);
  streamlineFilter.SetSeeds(seedsAH);
  streamlineFilter.SetNumberOfSteps(numSteps);
  auto output = streamlineFilter.Execute(input);
  return output;
#endif
}

} // namespace vtkh
