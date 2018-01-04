//-----------------------------------------------------------------------------
///
/// file: t_vtk-h_slice.cpp
///
//-----------------------------------------------------------------------------

#include "gtest/gtest.h"

#include <vtkh/vtkh.hpp>
#include <vtkh/DataSet.hpp>
#include <vtkh/filters/Slice.hpp>
#include <vtkh/rendering/RayTracer.hpp>
#include <vtkh/rendering/Scene.hpp>
#include "t_test_utils.hpp"

#include <iostream>


TEST(vtkh_slice, vtkh_slice)
{
  vtkh::DataSet data_set;
 
  const int base_size = 32;
  const int num_blocks = 1; 
  
  for(int i = 0; i < num_blocks; ++i)
  {
    data_set.AddDomain(CreateTestData(i, num_blocks, base_size), i);
  }

  vtkh::Slice slicer;

  vtkm::Vec<vtkm::Float32,3> normal(.5f,.5f,.5f);
  vtkm::Vec<vtkm::Float32,3> point(16.f,16.f,16.f);
  slicer.SetPlane(point, normal);
  slicer.SetInput(&data_set);
  slicer.Update();
  vtkh::DataSet *slice  = slicer.GetOutput();

  vtkm::Bounds bounds = slice->GetGlobalBounds();
  float bg_color[4] = { 0.f, 0.f, 0.f, 1.f};
  vtkm::rendering::Camera camera;
  camera.ResetToBounds(bounds);
  vtkh::Render render = vtkh::MakeRender<vtkh::RayTracer>(512, 
                                                          512, 
                                                          camera, 
                                                          *slice, 
                                                          "slice",
                                                           bg_color);  
  vtkh::RayTracer tracer;
  tracer.SetInput(slice);
  tracer.SetField("cell_data"); 

  vtkh::Scene scene;
  scene.AddRenderer(&tracer);
  scene.AddRender(render);
  scene.Render();

  delete slice; 
}

TEST(vtkh_slice, vtkh_mulit_slice)
{
  vtkh::DataSet data_set;
 
  const int base_size = 32;
  const int num_blocks = 1; 
  
  for(int i = 0; i < num_blocks; ++i)
  {
    data_set.AddDomain(CreateTestData(i, num_blocks, base_size), i);
  }

  vtkm::Bounds bounds;
  vtkh::Scene scene;
  std::vector<vtkm::Id> domain_ids;
  
  // add the first slice
  vtkh::Slice slicer1;

  vtkm::Vec<vtkm::Float32,3> normal1(1.0f,1.f,.0f);
  vtkm::Vec<vtkm::Float32,3> point1(16.f,16.f,16.f);
  slicer1.SetPlane(point1, normal1);
  slicer1.SetInput(&data_set);
  slicer1.Update();
  vtkh::DataSet *slice1  = slicer1.GetOutput();
  domain_ids = slice1->GetDomainIds();
  bounds.Include(slice1->GetGlobalBounds());

  vtkh::RayTracer tracer1;
  tracer1.SetInput(slice1);
  tracer1.SetField("cell_data"); 
  scene.AddRenderer(&tracer1);

  // add the second slice
  vtkh::Slice slicer2;

  vtkm::Vec<vtkm::Float32,3> normal2(1.f,.0f,.0f);
  vtkm::Vec<vtkm::Float32,3> point2(16.f,16.f,16.f);
  slicer2.SetPlane(point2, normal2);
  slicer2.SetInput(&data_set);
  slicer2.Update();
  vtkh::DataSet *slice2  = slicer2.GetOutput();
  domain_ids = slice2->GetDomainIds();
  bounds.Include(slice2->GetGlobalBounds());

  vtkh::RayTracer tracer2;
  tracer2.SetInput(slice2);
  tracer2.SetField("cell_data"); 

  scene.AddRenderer(&tracer2);

  float bg_color[4] = { 0.f, 0.f, 0.f, 1.f};

  vtkh::Render render = vtkh::MakeRender<vtkh::RayTracer>(512, 
                                                          512, 
                                                          bounds, 
                                                          domain_ids, 
                                                          "2slice",
                                                           bg_color);  
  
  scene.AddRender(render);
  scene.Render();
  
  delete slice1; 
  delete slice2; 
}
