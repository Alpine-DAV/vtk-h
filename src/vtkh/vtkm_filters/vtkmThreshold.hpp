#ifndef VTK_H_VTKM_THRESHOLD_HPP
#define VTK_H_VTKM_THRESHOLD_HPP

#include <vtkm/cont/DataSet.h>
#include <vtkm/filter/FieldSelection.h>
#include <vtkh/utils/vtkm_permutation_removal.hpp>

namespace vtkh
{

class vtkmThreshold
{
public:
  vtkm::cont::DataSet DoIt(vtkm::cont::DataSet &input,
                           std::string field_name,
                           double min_value,
                           double max_value,
                           vtkm::filter::FieldSelection map_fields);
};

}
#endif
