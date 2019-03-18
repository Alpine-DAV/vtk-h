#ifndef VTK_H_PATH_TRACE_HPP
#define VTK_H_PATH_TRACE_HPP

#include <vtkh/vtkh.hpp>
#include <vtkh/filters/Filter.hpp>
#include <vtkh/DataSet.hpp>

#include <vtkh/filters/communication/BoundsMap.hpp>
#include <vtkh/filters/communication/SpatialQuery.hpp>

#include <vtkm/rendering/raytracing/Ray.h>
#include <vtkm/rendering/Camera.h>

namespace vtkh
{

class PathTrace : public Filter
{
public:
  using vtkmRay = vtkm::rendering::raytracing::Ray<vtkm::Float32>;
  using vtkmCamera = vtkm::rendering::Camera;

  PathTrace();
  virtual ~PathTrace();
  std::string GetName() const override;
  void SetField(const std::string &field_name);
  void SetCamera(const vtkmCamera &camera);

protected:
  void PreExecute() override;
  void PostExecute() override;
  void DoExecute() override;

  std::string m_field_name;
  vtkmCamera m_camera;

  BoundsMap m_bounds_map;
  SpatialQuery m_spatial_query;


  void CreateRays(vtkmRay &rays);
};

} //namespace vtkh
#endif
