#ifndef VTK_H_NO_OP_HPP
#define VTK_H_NO_OP_HPP

#include <vtkh/vtkh.hpp>
#include <vtkh/filters/Filter.hpp>
#include <vtkh/DataSet.hpp>

namespace vtkh
{

class Histogram : public Filter
{
public:
  Histogram();
  virtual ~Histogram();
  std::string GetName() const override;
  void SetField(const std::string &field_name);

protected:
  void PreExecute() override;
  void PostExecute() override;
  void DoExecute() override;

  std::string m_field_name;
  int m_num_bins;
};

} //namespace vtkh
#endif
