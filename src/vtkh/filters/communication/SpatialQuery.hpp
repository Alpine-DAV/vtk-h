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

  struct Result
  {
    vtkm::cont::ArrayHandle<vtkm::Int32> m_candidates;
    vtkm::cont::ArrayHandle<vtkm::Int32> m_offsets;
    vtkm::cont::ArrayHandle<vtkm::Int32> m_counts;
    vtkm::cont::ArrayHandle<vtkm::Int32> m_domain_map;

    Result(vtkm::cont::ArrayHandle<vtkm::Int32> &candidates,
           vtkm::cont::ArrayHandle<vtkm::Int32> &offsets,
           vtkm::cont::ArrayHandle<vtkm::Int32> &counts,
           vtkm::cont::ArrayHandle<vtkm::Int32> &domain_map)
      : m_candidates(candidates),
        m_offsets(offsets),
        m_counts(counts),
        m_domain_map(domain_map)
    {}

    void GetDomains(const int index, std::vector<int> &domains)
    {
      domains.clear();
      auto pcounts = m_counts.GetPortalConstControl();
      auto poffsets = m_offsets.GetPortalConstControl();
      auto pcandidates = m_candidates.GetPortalConstControl();
      auto pmap = m_domain_map.GetPortalConstControl();
      const int count = pcounts.Get(index);
      const int offset = poffsets.Get(index);
      domains.resize(count);
      for(int i = 0; i < count; ++i)
      {
        const int d_index = pcandidates.Get(offset + i);
        const int dom = pmap.Get(d_index);
        domains[i] = dom;
      }

    }
  };

  Result IntersectRays(vtkm::rendering::raytracing::Ray<vtkm::Float32> &rays);

  //void IntersectRays(vtkm::rendering::raytracing::Ray<vtkm::Float32> &rays);
protected:
  BoundsMap m_bounds_map;
  bool m_is_built;
  void Build();

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32,4>> m_inner_nodes;
  vtkm::cont::ArrayHandle<vtkm::Id> m_leaf_nodes;
  std::vector<vtkm::cont::ArrayHandle<vtkm::Float32>> m_aabb_handles;

  vtkm::cont::ArrayHandle<vtkm::Int32> m_domain_map;
};

} // namespace vtkh
#endif
