#include "vtkmCleanGrid.hpp"
#include <vtkm/filter/CleanGrid.h>

namespace vtkh
{
vtkm::cont::DataSet
vtkmCleanGrid::Run(vtkm::cont::DataSet &input,
                   vtkm::filter::FieldSelection map_fields)
{
  vtkm::filter::CleanGrid cleaner;
  cleaner.SetFieldsToPass(map_fields);
  auto output = cleaner.Execute(input);
  return output;
}

} // namespace vtkh
