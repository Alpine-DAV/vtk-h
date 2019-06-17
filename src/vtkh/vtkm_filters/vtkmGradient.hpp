#ifndef VTK_H_VTKM_GRADIENT_HPP
#define VTK_H_VTKM_GRADIENT_HPP

#include <vtkm/cont/DataSet.h>
#include <vtkm/filter/FieldSelection.h>

namespace vtkh
{

struct GradientParameters
{
  bool use_point_gradient = true;
  bool compute_divergence = false;
  bool compute_vorticity  = false;
  bool compute_qcriterion = false;

  std::string output_name = "gradient";
  std::string divergence_name = "divergence";
  std::string vorticity_name = "vorticity";
  std::string qcriterion_name = "qcriterion";

};

class vtkmGradient
{
public:
  vtkm::cont::DataSet Run(vtkm::cont::DataSet &input,
                          std::string field_name,
                          GradientParameters params,
                          vtkm::filter::FieldSelection map_fields);
};
}
#endif
