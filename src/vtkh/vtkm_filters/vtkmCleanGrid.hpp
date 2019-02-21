#ifndef VTK_H_VTKM_CLEAN_GRID_HPP
#define VTK_H_VTKM_CLEAN_GRID_HPP

#include <vtkm/cont/DataSet.h>
#include <vtkm/filter/FieldSelection.h>

namespace vtkh
{

class vtkmCleanGrid
{
public:
  vtkm::cont::DataSet Run(vtkm::cont::DataSet &input,
                          vtkm::filter::FieldSelection map_fields);
};
}
#endif
