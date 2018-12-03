#ifndef VTK_H_LOG_HPP
#define VTK_H_LOG_HPP

#include <vtkh/vtkh.hpp>
#include <vtkh/filters/Filter.hpp>
#include <vtkh/DataSet.hpp>
#include <vtkm/Range.h>

namespace vtkh
{

class Log: public Filter
{
public:
  Log(); 
  virtual ~Log(); 
  std::string GetName() const override; 
  void SetField(const std::string &field_name);
  void SetResultField(const std::string &field_name);

  std::string GetField() const;
  std::string GetResultField() const;
protected:
  void PreExecute() override;
  void PostExecute() override;
  void DoExecute() override;
  std::string m_field_name;
  std::string m_result_name;

};

} //namespace vtkh
#endif
