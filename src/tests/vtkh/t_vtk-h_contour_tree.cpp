//-----------------------------------------------------------------------------
///
/// file: t_vtk-h_dataset.cpp
///
//-----------------------------------------------------------------------------

#include "gtest/gtest.h"

#include <vtkh/vtkh.hpp>
#include <vtkh/DataSet.hpp>
#include <vtkh/filters/ContourTree.hpp>
//#include <vtkh/filters/MarchingCubes.hpp>
#include <vtkh/rendering/RayTracer.hpp>
#include <vtkh/rendering/Scene.hpp>
#include "t_test_utils.hpp"

#include <iostream>



//----------------------------------------------------------------------------
TEST(vtkh_contour_tree, vtkh_contour_tree)
{
  vtkh::DataSet data_set;
 
  const int base_size = 32;
  const int num_blocks = 1; 
  
  for(int i = 0; i < num_blocks; ++i)
  {
    data_set.AddDomain(CreateTestData(i, num_blocks, base_size), i);
  }

  vtkh::ContourTree contour_tree;
  contour_tree.SetInput(&data_set);
  contour_tree.SetField("point_data"); 

  //const int num_vals = 2;
  //double iso_vals [num_vals];
  //iso_vals[0] = -1; // ask for something that does not exist
  //iso_vals[1] = (float)base_size * (float)num_blocks * 0.5f;

  //contour_tree.SetIsoValues(iso_vals, num_vals);
  contour_tree.Update();

  //vtkh::DataSet *iso_output = contour_tree.GetOutput();

  //vtkm::Bounds bounds = iso_output->GetGlobalBounds();
  //float bg_color[4] = { 0.f, 0.f, 0.f, 1.f};
  //vtkm::rendering::Camera camera;
  //camera.ResetToBounds(bounds);
  //vtkh::Render render = vtkh::MakeRender(512, 
  //                                       512, 
  //                                       camera, 
  //                                       *iso_output, 
  //                                       "iso",
  //                                        bg_color);  
  //vtkh::RayTracer tracer;
  //tracer.SetInput(iso_output);
  //tracer.SetField("cell_data"); 

  //vtkh::Scene scene;
  //scene.AddRenderer(&tracer);
  //scene.AddRender(render);
  //scene.Render();

  //delete iso_output; 
}
