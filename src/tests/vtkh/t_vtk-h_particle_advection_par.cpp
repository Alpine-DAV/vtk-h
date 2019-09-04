//-----------------------------------------------------------------------------
///
/// file: t_vtk-h_particle_advection_par.cpp
///
//-----------------------------------------------------------------------------

#include "gtest/gtest.h"

#include <vtkh/vtkh.hpp>
#include <vtkh/DataSet.hpp>
#include <vtkh/filters/ParticleAdvection.hpp>
#include <vtkm/io/writer/VTKDataSetWriter.h>
#include <vtkm/cont/DataSet.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include "t_test_utils.hpp"
#include <iostream>
#include <mpi.h>

void checkValidity(vtkh::DataSet *data, const int maxSteps)
{
  int numDomains = data->GetNumberOfDomains();

  //Check all domains
  for(int i = 0; i < numDomains; i++)
  {
    auto currentDomain = data->GetDomain(i);
    vtkm::cont::CellSetExplicit<> cellSet =
          currentDomain.GetCellSet().Cast<vtkm::cont::CellSetExplicit<>>();

    //Ensure that streamlines took <= to the max number of steps
    for(int j = 0; j < cellSet.GetNumberOfCells(); j++)
    {
      EXPECT_LE(cellSet.GetNumberOfPointsInCell(j), maxSteps);
    }
  }
}

void writeDataSet(vtkh::DataSet *data, std::string fName, int rank)
{
  int numDomains = data->GetNumberOfDomains();
  std::cerr << "num domains " << numDomains << std::endl;
  for(int i = 0; i < numDomains; i++)
  {
    char fileNm[128];
    sprintf(fileNm, "%s.rank%d.domain%d.vtk", fName.c_str(), rank, i);
    vtkm::io::writer::VTKDataSetWriter write(fileNm);
    write.WriteDataSet(data->GetDomain(i));
  }
}

//----------------------------------------------------------------------------
TEST(vtkh_particle_advection, vtkh_serial_particle_advection)
{
  MPI_Init(NULL, NULL);
  int comm_size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  vtkh::SetMPICommHandle(MPI_Comm_c2f(MPI_COMM_WORLD));

  std::cout << "Running parallel Particle Advection, vtkh - with " << comm_size << " ranks" << std::endl;

  vtkh::DataSet data_set;
  const int base_size = 32;
  const int blocks_per_rank = 1;
  const int maxAdvSteps = 100;
  const int num_blocks = comm_size * blocks_per_rank;

  for(int i = 0; i < blocks_per_rank; ++i)
  {
    int domain_id = rank * blocks_per_rank + i;
    data_set.AddDomain(CreateTestDataRectilinear(domain_id, num_blocks, base_size), domain_id);
  }

  vtkh::ParticleAdvection streamline;
  streamline.SetInput(&data_set);
  streamline.SetField("vector_data_Float64");
  streamline.SetMaxSteps(maxAdvSteps);
  streamline.SetStepSize(0.1);
  streamline.SetSeedsRandomWhole(500);
  streamline.Update();
  vtkh::DataSet *streamline_output = streamline.GetOutput();

  checkValidity(streamline_output, maxAdvSteps);
  writeDataSet(streamline_output, "advection_SeedsRandomWhole", rank);

  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
}
