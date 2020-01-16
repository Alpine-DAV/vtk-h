#ifndef VTK_H_LAGRANGIAN_INTERPOLATION_HPP
#define VTK_H_LAGRANGIAN_INTERPOLATION_HPP

#include <vtkh/vtkh_exports.h>
#include <vtkh/vtkh.hpp>
#include <vtkh/filters/Filter.hpp>
#include <vtkh/DataSet.hpp>

namespace vtkh
{

class VTKH_API LagrangianInterpolation : public Filter
{
public:
  LagrangianInterpolation();
  virtual ~LagrangianInterpolation();
  std::string GetName() const override;
  void SetField(const std::string &field_name);
//  void SetInputPath(const std::string &input_path);
  void SetWriteFrequency(const int &write_frequency);

private:
  bool AllMessagesReceived(bool *a, int num_ranks);
  bool BoundsCheck(float x, float y, float z, double *BB);
  bool BoundsCheck(float x, float y, float z, double xmin, double xmax, double ymin, double ymax, double zmin, double zmax);
  std::vector<int> GetNeighborRankList(vtkm::Id num_ranks, vtkm::Id rank, double *bbox_list);

protected:
  void PreExecute() override;
  void PostExecute() override;
  void DoExecute() override;

  std::string m_field_name;
//  std::string m_input_path;
  int m_write_frequency;
};

} //namespace vtkh
#endif
