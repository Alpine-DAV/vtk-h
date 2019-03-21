#ifndef VTK_H_ISO_VOLUME_HPP
#define VTK_H_ISO_VOLUME_HPP

#include <vtkh/filters/Filter.hpp>
#include <memory>

namespace vtkh
{

class IsoVolume: public Filter
{
public:
  IsoVolume();
  virtual ~IsoVolume();
  std::string GetName() const override;
  void SetRange(const vtkm::Range range);
  void SetField(const std::string field_name);
protected:
  void PreExecute() override;
  void PostExecute() override;
  void DoExecute() override;

  vtkm::Range m_range;
  std::string m_field_name;
};

} //namespace vtkh
#endif
