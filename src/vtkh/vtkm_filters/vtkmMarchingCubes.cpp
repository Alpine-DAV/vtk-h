#include "vtkmMarchingCubes.hpp"
#include <vtkm/filter/MarchingCubes.h>

namespace vtkh
{
vtkm::cont::DataSet
vtkmMarchingCubes::Run(vtkm::cont::DataSet &input,
                       std::string field_name,
                       std::vector<double> iso_values,
                       vtkm::filter::FieldSelection map_fields)
{
  vtkm::filter::MarchingCubes marcher;

  marcher.SetFieldsToPass(map_fields);
  marcher.SetIsoValues(iso_values);
  marcher.SetMergeDuplicatePoints(false);
  marcher.SetActiveField(field_name);

  auto output = marcher.Execute(input);
  return output;
}

} // namespace vtkh
