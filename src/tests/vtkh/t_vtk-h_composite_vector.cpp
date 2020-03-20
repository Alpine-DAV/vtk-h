//-----------------------------------------------------------------------------
///
/// file: t_vtk-h_dataset.cpp
///
//-----------------------------------------------------------------------------

#include "gtest/gtest.h"

#include <vtkh/vtkh.hpp>
#include <vtkh/DataSet.hpp>
#include <vtkh/filters/CompositeVector.hpp>
#include <vtkh/filters/VectorMagnitude.hpp>
#include <vtkh/rendering/RayTracer.hpp>
#include <vtkh/rendering/Scene.hpp>
#include "t_test_utils.hpp"

#include <iostream>



//----------------------------------------------------------------------------
TEST(vtkh_composite_vector, vtkh_vector3d)
{
  vtkh::DataSet data_set;

  const int base_size = 32;
  const int num_blocks = 2;

  for(int i = 0; i < num_blocks; ++i)
  {
    data_set.AddDomain(CreateTestData(i, num_blocks, base_size), i);
  }

  vtkh::CompositeVector vec_maker;
  vec_maker.SetInput(&data_set);
  vec_maker.SetFields("point_data_Float64", "point_data_Float64", "point_data_Float64");
  vec_maker.SetResultField("my_vec3");

  vec_maker.Update();
  vtkh::DataSet *output = vec_maker.GetOutput();

  vtkh::VectorMagnitude mag;
  mag.SetInput(output);
  mag.SetField("my_vec3");

  mag.SetResultName("mag");
  mag.Update();

  vtkh::DataSet *mag_output = mag.GetOutput();
  vtkm::Bounds bounds = mag_output->GetGlobalBounds();

  vtkm::rendering::Camera camera;
  camera.SetPosition(vtkm::Vec<vtkm::Float64,3>(-16, -16, -16));
  camera.ResetToBounds(bounds);
  vtkh::Render render = vtkh::MakeRender(512,
                                         512,
                                         camera,
                                         *output,
                                         "vector_composite3d");
  vtkh::RayTracer tracer;
  tracer.SetInput(mag_output);
  tracer.SetField("mag");

  vtkh::Scene scene;
  scene.AddRender(render);
  scene.AddRenderer(&tracer);
  scene.Render();

  delete output;
  delete mag_output;
}
