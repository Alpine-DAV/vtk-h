#include "vtkmClip.hpp"
#include <vtkm/filter/ClipWithImplicitFunction.h>

namespace vtkh
{
vtkm::cont::DataSet
vtkmClip::Run(vtkm::cont::DataSet &input,
                       vtkm::cont::ImplicitFunctionHandle func,
                       vtkm::Id activeCellSet,
                       bool invert,
                       vtkm::filter::FieldSelection map_fields)
{
  vtkm::filter::ClipWithImplicitFunction clipper;

  clipper.SetImplicitFunction(func);
  clipper.SetActiveCellSetIndex(activeCellSet);
  clipper.SetInvertClip(invert);
  clipper.SetFieldsToPass(map_fields);

  auto output = clipper.Execute(input);
  return output;
}

} // namespace vtkh
