#ifndef VTK_H_PATH_TRACE_HPP
#define VTK_H_PATH_TRACE_HPP

#include <vtkh/vtkh.hpp>
#include <vtkh/filters/Filter.hpp>
#include <vtkh/DataSet.hpp>

#include <vtkh/filters/communication/BoundsMap.hpp>
#include <vtkh/filters/communication/SpatialQuery.hpp>
#include <vtkh/filters/communication/Ray.hpp>

#include <vtkm/rendering/raytracing/Ray.h>
#include <vtkm/rendering/Camera.h>

#ifdef VTKH_PARALLEL
#include <vtkh/filters/communication/RayMessenger.hpp>
#endif

namespace vtkh
{

class PathTrace : public Filter
{
public:
  using vtkmRay = vtkm::rendering::raytracing::Ray<vtkm::Float32>;
  using vtkmCamera = vtkm::rendering::Camera;
  using IncomingQueue = std::map<int, std::vector<Ray>>;
  using OutgoingQueue = std::map<int, std::vector<Ray>>;
  using DomainQueue = std::map<int, std::vector<Ray>>;
  using DomainRays = std::map<int, vtkmRay>;

  PathTrace();
  virtual ~PathTrace();
  std::string GetName() const override;
  void SetField(const std::string &field_name);
  void SetCamera(const vtkmCamera &camera);

protected:
  void PreExecute() override;
  void PostExecute() override;
  void DoExecute() override;

  void PackOutgoing(vtkmRay &rays);               // pack outgoing rays into out q
  void RouteOutgoing();                           // route out q to destinations
  void RouteIncoming(std::vector<Ray> &in_rays);  // route incoming rays to local domains
  void Recv();

  std::string m_field_name;
  vtkmCamera m_camera;

  BoundsMap m_bounds_map;
  SpatialQuery m_spatial_query;
  vtkm::Bounds m_bounds;

  OutgoingQueue m_out_q;
  IncomingQueue m_in_q;
  DomainRays    m_dom_rays;

  int m_rank;
  int m_procs;

#ifdef VTKH_PARALLEL
  RayMessenger m_messenger;
#endif

  void CreateRays(vtkmRay &rays);
};

} //namespace vtkh
#endif
