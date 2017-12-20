#ifndef VTK_H_CLEAN_GRID_HPP
#define VTK_H_CLEAN_GRID_HPP

#include <vtkh/vtkh.hpp>
#include <vtkh/filters/Filter.hpp>
#include <vtkh/DataSet.hpp>


namespace vtkh
{

class CleanGrid : public Filter
{
public:
  CleanGrid(); 
  virtual ~CleanGrid(); 
protected:
  void PreExecute() override;
  void PostExecute() override;
  void DoExecute() override;
};

} //namespace vtkh
#endif
