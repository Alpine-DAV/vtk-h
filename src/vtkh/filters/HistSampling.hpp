#ifndef VTK_H_HIST_SAMPLING_HPP
#define VTK_H_HIST_SAMPLING_HPP

#include <vtkh/vtkh.hpp>
#include <vtkh/filters/Filter.hpp>
#include <vtkh/DataSet.hpp>

namespace vtkh
{

class HistSampling : public Filter
{
public:
  HistSampling(); 
  virtual ~HistSampling(); 
  std::string GetName() const override;
  void SetField(const std::string &field_name);
  std::string GetField() const;
protected:
  void PreExecute() override;
  void PostExecute() override;
  void DoExecute() override;

  std::string m_field_name;
};

} //namespace vtkh
#endif
