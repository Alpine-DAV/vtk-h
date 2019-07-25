//-----------------------------------------------------------------------------
///
/// file: t_vtk-h_dataset.cpp
///
//-----------------------------------------------------------------------------

#include "gtest/gtest.h"

#include <mpi.h>
#include <vtkh/vtkh.hpp>
#include <vtkh/DataSet.hpp>
#include <vtkh/rendering/Scene.hpp>
#include <vtkh/rendering/VolumeRenderer.hpp>
#include "t_test_utils.hpp"

#include <iostream>
#include <mpi.h>


//----------------------------------------------------------------------------
TEST(vtkh_volume_renderer, vtkh_parallel_render)
{
  MPI_Init(NULL, NULL);
  int comm_size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  vtkh::SetMPICommHandle(MPI_Comm_c2f(MPI_COMM_WORLD));
  vtkh::DataSet data_set;

  const int base_size = 32;
  const int blocks_per_rank = 2;
  const int num_blocks = comm_size * blocks_per_rank;

  for(int i = 0; i < blocks_per_rank; ++i)
  {
    int domain_id = rank * blocks_per_rank + i;
    data_set.AddDomain(CreateTestData(domain_id, num_blocks, base_size), domain_id);
  }

  vtkm::Bounds bounds = data_set.GetGlobalBounds();

  vtkm::rendering::Camera camera;
  camera.SetPosition(vtkm::Vec<vtkm::Float64,3>(-16, -16, -16));
  camera.ResetToBounds(bounds);
  vtkh::Render render = vtkh::MakeRender(512,
                                         512,
                                         camera,
                                         data_set,
                                         "volume_par");

  vtkm::cont::ColorTable color_map("cool to warm");
  color_map.AddPointAlpha(0.0, .05);
  color_map.AddPointAlpha(1.0, .5);

  vtkh::VolumeRenderer tracer;
  tracer.SetColorTable(color_map);
  tracer.SetInput(&data_set);
  tracer.SetField("point_data_Float64");

  tracer.Update();

  vtkh::Scene scene;
  scene.AddRender(render);
  scene.AddRenderer(&tracer);
  scene.Render();

  MPI_Finalize();
}
