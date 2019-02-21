#include "vtkmThreshold.hpp"
#include <vtkm/filter/Threshold.h>
namespace vtkh
{
vtkm::cont::DataSet
vtkmThreshold::DoIt(vtkm::cont::DataSet &input,
                    std::string field_name,
                    double min_value,
                    double max_value,
                    vtkm::filter::FieldSelection map_fields)
{
  vtkm::filter::Threshold thresholder;
  thresholder.SetUpperThreshold(max_value);
  thresholder.SetLowerThreshold(min_value);
  thresholder.SetActiveField(field_name);
  thresholder.SetFieldsToPass(map_fields);
  auto output = thresholder.Execute(input);
  vtkh::StripPermutation(output);
  return output;
}

} // namespace vtkh
