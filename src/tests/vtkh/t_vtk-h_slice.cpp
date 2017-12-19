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
