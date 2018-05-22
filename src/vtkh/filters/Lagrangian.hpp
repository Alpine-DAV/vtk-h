#ifndef VTK_H_NO_OP_HPP
#define VTK_H_NO_OP_HPP

#include <vtkh/vtkh.hpp>
#include <vtkh/filters/Filter.hpp>
#include <vtkh/DataSet.hpp>

namespace vtkh
{

class Lagrangian : public Filter
{
public:
  Lagrangian(); 
  virtual ~Lagrangian(); 
  std::string GetName() const override;
  void SetField(const std::string &field_name);

protected:
  void PreExecute() override;
  void PostExecute() override;
  void DoExecute() override;

  std::string m_field_name;
};

} //namespace vtkh
#endif
