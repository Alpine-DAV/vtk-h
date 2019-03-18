//-----------------------------------------------------------------------------
///
/// file: t_vtk-h_dataset.cpp
///
//-----------------------------------------------------------------------------

#include "gtest/gtest.h"

#include <vtkh/vtkh.hpp>
#include <vtkh/DataSet.hpp>
#include <vtkh/filters/communication/SpatialQuery.hpp>
#include <vtkm/rendering/raytracing/Ray.h>
#include "t_test_utils.hpp"

#include <iostream>

void add_ray(vtkm::rendering::raytracing::Ray<vtkm::Float32> &rays,
             const vtkm::Id &index,
             const vtkm::Vec<vtkm::Float32,3> &origin,
             const vtkm::Vec<vtkm::Float32,3> &dir,
             const vtkm::Float32 &min_distance,
             const vtkm::Float32 &max_distance)
{
  auto o_portal = rays.Origin.GetPortalControl();
  auto d_portal = rays.Dir.GetPortalControl();
  auto min_portal = rays.MinDistance.GetPortalControl();
  auto max_portal = rays.MaxDistance.GetPortalControl();
  o_portal.Set(index, origin);
  d_portal.Set(index, dir);
  min_portal.Set(index, min_distance);
  max_portal.Set(index, max_distance);

}


//----------------------------------------------------------------------------
TEST(vtkh_spatial_bounds, vtkh_basic)
{

  vtkm::Bounds b1,b2;

  vtkm::Vec<vtkm::Float32, 3> p1(0.f, 0.f, 0.f);
  vtkm::Vec<vtkm::Float32, 3> p2(1.f, 1.f, 1.f);
  b1.Include(p1);
  b1.Include(p2);

  vtkm::Vec<vtkm::Float32, 3> p3(1.f, 0.f, 0.f);
  vtkm::Vec<vtkm::Float32, 3> p4(2.f, 1.f, 1.f);

  b2.Include(p3);
  b2.Include(p4);

  vtkh::BoundsMap bmap;
  bmap.AddBlock(0,b1);
  bmap.AddBlock(1,b2);

  vtkh::SpatialQuery squery(bmap);

  vtkm::rendering::raytracing::Ray<vtkm::Float32> rays;
  rays.Resize(1);

  vtkm::Vec<vtkm::Float32,3> orig(-.1f, .5f, .5f);
  vtkm::Vec<vtkm::Float32,3> dir(1.f, 0.f, 0.f);
  vtkm::Float32 d_min = 0.f, d_max = 4.f;

  add_ray(rays, 0, orig, dir, d_min, d_max);

  squery.IntersectRays(rays);

}

//----------------------------------------------------------------------------
TEST(vtkh_spatial_bounds, vtkh_overlap)
{

  vtkm::Bounds b1,b2;

  vtkm::Vec<vtkm::Float32, 3> p1(0.f, 0.f, 0.f);
  vtkm::Vec<vtkm::Float32, 3> p2(1.f, 1.f, 1.f);
  b1.Include(p1);
  b1.Include(p2);

  vtkm::Vec<vtkm::Float32, 3> p3(0.5f, 0.f, 0.f);
  vtkm::Vec<vtkm::Float32, 3> p4(2.f, 1.f, 1.f);

  b2.Include(p3);
  b2.Include(p4);

  vtkh::BoundsMap bmap;
  bmap.AddBlock(0,b1);
  bmap.AddBlock(1,b2);

  vtkh::SpatialQuery squery(bmap);

  vtkm::rendering::raytracing::Ray<vtkm::Float32> rays;
  rays.Resize(1);

  vtkm::Vec<vtkm::Float32,3> orig(.75f, .5f, .5f);
  vtkm::Vec<vtkm::Float32,3> dir(1.f, 0.f, 0.f);
  vtkm::Float32 d_min = 0.f, d_max = 4.f;

  add_ray(rays, 0, orig, dir, d_min, d_max);

  squery.IntersectRays(rays);

}

//----------------------------------------------------------------------------
TEST(vtkh_spatial_bounds, vtkh_edge_case)
{
  vtkh::BoundsMap bmap;
  vtkm::Id3 dims(2,2,2);
  int blockid = 0;

  for(vtkm::Id z = 0; z < dims[2]; ++z)
    for(vtkm::Id y = 0; y < dims[1]; ++y)
      for(vtkm::Id x = 0; x < dims[0]; ++x)
      {
        //constexpr float eps = .000001f;
        constexpr float eps = .00f;
        vtkm::Bounds bounds;
        vtkm::Vec<vtkm::Float32, 3> point;
        point[0] = vtkm::Float32(x) - eps;
        point[1] = vtkm::Float32(y) - eps;
        point[2] = vtkm::Float32(z) - eps;
        std::cout<<point<<" ";
        bounds.Include(point);
        point[0] = vtkm::Float32(x) + 1.f + eps;
        point[1] = vtkm::Float32(y) + 1.f + eps;
        point[2] = vtkm::Float32(z) + 1.f + eps;
        bounds.Include(point);
        std::cout<<point<<"\n";
        bmap.AddBlock(blockid,bounds);
        blockid++;
      }

  vtkh::SpatialQuery squery(bmap);

  vtkm::rendering::raytracing::Ray<vtkm::Float32> rays;
  rays.Resize(1);

  vtkm::Vec<vtkm::Float32,3> orig(1.f, 1.f, 0.0f);
  vtkm::Vec<vtkm::Float32,3> dir(0.f, 1.f, 0.f);
  vtkm::Float32 d_min = 0.f, d_max = 4.f;

  add_ray(rays, 0, orig, dir, d_min, d_max);

  squery.IntersectRays(rays);

}

