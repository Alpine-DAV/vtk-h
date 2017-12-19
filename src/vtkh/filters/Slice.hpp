#ifndef VTK_H_SLICE_HPP
#define VTK_H_SLICE_HPP

#include <vtkh/vtkh.hpp>
#include <vtkh/filters/Filter.hpp>
#include <vtkh/DataSet.hpp>


namespace vtkh
{

class Slice : public Filter
{
public:
  Slice(); 
  virtual ~Slice(); 
  void SetPlane(vtkm::Vec<vtkm::Float32,3> point, vtkm::Vec<vtkm::Float32,3> normal);
protected:
  void PreExecute() override;
  void PostExecute() override;
  void DoExecute() override;
  vtkm::Vec<vtkm::Float32,3> m_point;
  vtkm::Vec<vtkm::Float32,3> m_normal;
};

} //namespace vtkh
#endif
