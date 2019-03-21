
#include "SpatialQuery.hpp"
#include "SpatialQueryDetail.hpp"

#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/rendering/raytracing/BoundingVolumeHierarchy.h>

namespace vtkh
{

SpatialQuery::SpatialQuery()
 : m_is_built(false)
{}

SpatialQuery::SpatialQuery(const BoundsMap &bounds_map)
 : m_bounds_map(bounds_map),
   m_is_built(false)
{
}

void SpatialQuery::SetBoundsMap(const BoundsMap &bounds_map)
{
  m_bounds_map = bounds_map;
  m_is_built = false;
}

void SpatialQuery::Build()
{
  vtkm::rendering::raytracing::AABBs aabbs;
  int size = m_bounds_map.bm.size();
  if(size == 0) std::cout<<"Don't build with size 0\n";
  aabbs.xmins.Allocate(size);
  aabbs.ymins.Allocate(size);
  aabbs.zmins.Allocate(size);
  aabbs.xmaxs.Allocate(size);
  aabbs.ymaxs.Allocate(size);
  aabbs.zmaxs.Allocate(size);

  auto x_mn = aabbs.xmins.GetPortalControl();
  auto y_mn = aabbs.ymins.GetPortalControl();
  auto z_mn = aabbs.zmins.GetPortalControl();
  auto x_mx = aabbs.xmaxs.GetPortalControl();
  auto y_mx = aabbs.ymaxs.GetPortalControl();
  auto z_mx = aabbs.zmaxs.GetPortalControl();

  int i = 0;
  m_domain_map.Allocate(size);
  auto domain_portal = m_domain_map.GetPortalControl();
  for(auto it = m_bounds_map.bm.begin(); it != m_bounds_map.bm.end(); ++it)
  {
    const vtkm::Bounds &bounds = it->second;
    x_mn.Set(i, bounds.X.Min);
    y_mn.Set(i, bounds.Y.Min);
    z_mn.Set(i, bounds.Z.Min);

    x_mx.Set(i, bounds.X.Max);
    y_mx.Set(i, bounds.Y.Max);
    z_mx.Set(i, bounds.Z.Max);
    // keep track of the aabb index to the domain id so we
    // can map back later
    domain_portal.Set(i, it->first);

    ++i;
  }

  m_aabb_handles.resize(6);
  // Copy the arrays since constructing the bvh will alter the
  // contents of the aabbs we hand it
  vtkm::cont::ArrayCopy(aabbs.xmins, m_aabb_handles[0]);
  vtkm::cont::ArrayCopy(aabbs.ymins, m_aabb_handles[1]);
  vtkm::cont::ArrayCopy(aabbs.zmins, m_aabb_handles[2]);

  vtkm::cont::ArrayCopy(aabbs.xmaxs, m_aabb_handles[3]);
  vtkm::cont::ArrayCopy(aabbs.ymaxs, m_aabb_handles[4]);
  vtkm::cont::ArrayCopy(aabbs.zmaxs, m_aabb_handles[5]);

  vtkm::rendering::raytracing::LinearBVH builder;
  builder.SetData(aabbs);
  builder.Construct();
  m_inner_nodes = builder.FlatBVH;
  m_leaf_nodes = builder.Leafs;
  m_is_built = true;
}

SpatialQuery::Result
SpatialQuery::IntersectRays(vtkm::rendering::raytracing::Ray<vtkm::Float32> &rays)
{
  if(!m_is_built)
  {
    Build();
  }

  detail::AABBData aabb_data(m_aabb_handles[0],
                             m_aabb_handles[1],
                             m_aabb_handles[2],
                             m_aabb_handles[3],
                             m_aabb_handles[4],
                             m_aabb_handles[5]);

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32,2>> min_max;
  min_max.Allocate(rays.NumRays);

  std::cout<<"done\n";
  vtkm::worklet::DispatcherMapField<detail::ClosestCandidates>(detail::ClosestCandidates())
    .Invoke(rays.Dir, rays.Origin, rays.MinDistance, rays.MaxDistance, min_max, m_inner_nodes, m_leaf_nodes, aabb_data );

  vtkm::cont::ArrayHandle<vtkm::Id> counts;
  counts.Allocate(rays.NumRays);

  vtkm::worklet::DispatcherMapField<detail::CountCandidates>(detail::CountCandidates())
    .Invoke(rays.Dir, rays.Origin, min_max, counts, m_inner_nodes, m_leaf_nodes, aabb_data);

  vtkm::cont::ArrayHandle<vtkm::Id> offsets;
  vtkm::Id sum = vtkm::cont::Algorithm::ScanExclusive(counts, offsets);

  std::cout<<"Sum "<<sum<<"\n";

  vtkm::cont::ArrayHandle<vtkm::Id> candidates;
  candidates.Allocate(sum);

  vtkm::worklet::DispatcherMapField<detail::GetCandidates>(detail::GetCandidates())
    .Invoke(rays.Dir, rays.Origin, min_max, counts, m_inner_nodes, m_leaf_nodes, aabb_data, offsets, candidates);

  //auto c = candidates.GetPortalControl();
  //auto cc = counts.GetPortalControl();
  //for(int i = 0; i < cc.Get(0); ++i)
  //{
  //  std::cout<<"Candidate "<<c.Get(i)<<"\n";
  //}
  return Result(candidates, offsets, counts, m_domain_map);

}

} // namespace vtkh

