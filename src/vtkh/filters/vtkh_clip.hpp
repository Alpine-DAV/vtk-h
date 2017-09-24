#ifndef VTK_H_CLIP_HPP
#define VTK_H_CLIP_HPP

#include "vtkh.hpp"
#include "vtkh_filter.hpp"
#include "vtkh_data_set.hpp"

#include <memory>

namespace vtkh
{

class Clip: public Filter
{
public:
  Clip(); 
  virtual ~Clip(); 
  void SetBoxClip(const vtkm::Bounds &clipping_bounds);
  void SetSphereClip(const double center[3], const double radius);
  void SetPlaneClip(const double origin[3], const double normal[3]);
  void SetCellSetIndex(vtkm::Id index);
  void SetCellSet(const std::string &cell_set);
protected:
  void PreExecute() override;
  void PostExecute() override;
  void DoExecute() override;

  struct InternalsType;
  std::shared_ptr<InternalsType> m_internals;
  std::string m_cell_set;
};

} //namespace vtkh
#endif
