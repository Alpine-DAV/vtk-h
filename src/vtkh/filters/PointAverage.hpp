#ifndef VTK_H_POINT_AVERAGE_HPP
#define VTK_H_POINT_AVERAGE_HPP

#include <vtkh/vtkh.hpp>
#include <vtkh/filters/Filter.hpp>
#include <vtkh/DataSet.hpp>

namespace vtkh
{

class PointAverage : public Filter
{
public:
  PointAverage(); 
  virtual ~PointAverage(); 
  std::string GetName() const override;
  void SetLevels(const int &levels);
  void SetField(const std::string &field_name);
  void SetOutputField(const std::string &field_name);

protected:
  void PreExecute() override;
  void PostExecute() override;
  void DoExecute() override;

  std::string m_field_name, m_output_field_name;
};

} //namespace vtkh
#endif
