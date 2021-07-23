#include "vtkmPointTransform.hpp"
#include <vtkm/filter/PointTransform.h>

namespace vtkh
{
vtkm::cont::DataSet
vtkmPointTransform::Run(vtkm::cont::DataSet &input,
                        vtkm::Matrix<double,4,4> &transform,
                        vtkm::filter::FieldSelection map_fields)
{
  vtkm::filter::PointTransform trans;

  trans.SetChangeCoordinateSystem(true);
  trans.SetFieldsToPass(map_fields);
  vtkm::Matrix<vtkm::FloatDefault,4,4> default_mat;
  for(int i = 0; i < 4; i++)
  {
    for(int j = 0; j < 4; j++)
    {
      default_mat[i][j] = static_cast<vtkm::FloatDefault>(transform[i][j]);
    }
  }
  trans.SetTransform(default_mat);

  auto output = trans.Execute(input);
  return output;
}

} // namespace vtkh
