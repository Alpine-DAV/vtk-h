#ifndef VTK_H_SPATIAL_QUERY_HPP
#define VTK_H_SPATIAL_QUERY_HPP

#include <vtkh/vtkh.hpp>
#include "BoundsMap.hpp"

#include <vtkm/rendering/raytracing/Ray.h>

namespace vtkh
{

class SpatialQuery
{
public:
  SpatialQuery();
  SpatialQuery(const BoundsMap &bounds_map);

  void SetBoundsMap(const BoundsMap &bounds_map);
  //void IntersectRays(vtkm::rendering::raytracing::Ray<vtkm::Float32> &rays,
  //                   vtkm::cont::ArrayHandle<vtkm::Int32> &candidates,
  //                   vtkm::cont::ArrayHandle<vtkm::Int32> &offsets,
  //                   vtkm::cont::ArrayHandle<vtkm::Int32> &counts);
  void IntersectRays(vtkm::rendering::raytracing::Ray<vtkm::Float32> &rays);
protected:
  BoundsMap m_bounds_map;
  bool m_is_built;
  void Build();

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32,4>> m_inner_nodes;
  vtkm::cont::ArrayHandle<vtkm::Id> m_leaf_nodes;
  std::vector<vtkm::cont::ArrayHandle<vtkm::Float32>> m_aabb_handles;

};

} // namespace vtkh
#endif
