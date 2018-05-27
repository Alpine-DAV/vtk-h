#ifndef VTK_H_LAGRANGIAN_HPP
#define VTK_H_LAGRANGIAN_HPP

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
  void SetStepSize(const double &step_size);
  void SetWriteFrequency(const int &write_frequency);

protected:
  void PreExecute() override;
  void PostExecute() override;
  void DoExecute() override;

  std::string m_field_name;
	double m_step_size;
	double m_write_frequency;
};

} //namespace vtkh
#endif
