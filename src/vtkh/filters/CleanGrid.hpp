#ifndef VTK_H_CLEAN_GRID_HPP
#define VTK_H_CLEAN_GRID_HPP

#include <vtkh/vtkh_exports.h>
#include <vtkh/vtkh.hpp>
#include <vtkh/filters/Filter.hpp>
#include <vtkh/DataSet.hpp>


namespace vtkh
{

class VTKH_API CleanGrid : public Filter
{
public:
  CleanGrid();
  virtual ~CleanGrid();
  std::string GetName() const override;
protected:
  void PreExecute() override;
  void PostExecute() override;
  void DoExecute() override;
};

} //namespace vtkh
#endif
