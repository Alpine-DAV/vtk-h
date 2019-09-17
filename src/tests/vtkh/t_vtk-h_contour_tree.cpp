//-----------------------------------------------------------------------------
///
/// file: t_vtk-h_dataset.cpp
///
//-----------------------------------------------------------------------------

#include "gtest/gtest.h"

#include <vtkh/vtkh.hpp>
#include <vtkh/DataSet.hpp>
#include <vtkh/filters/MarchingCubes.hpp>
#include <vtkh/rendering/RayTracer.hpp>
#include <vtkh/rendering/Scene.hpp>
#include "t_test_utils.hpp"

#include <vtkm/cont/DataSetBuilderUniform.h>
#include <vtkm/cont/DataSetFieldAdd.h>
#include <vtkm/cont/MultiBlock.h>

#include <iostream>
#ifdef VTKH_PARALLEL
#include <mpi.h>
#endif


using ValueType = vtkm::Float64;

#ifndef VTKH_PARALLEL
bool ReadTestData(const char* filename, vtkm::cont::DataSet& inDataSet)
{
  const int rank = 0;
#else
bool ReadTestData(const char* filename, vtkm::cont::MultiBlock& inDataSet,
                  int rank, int size)
{
#endif

  std::ifstream inFile(filename);
  if (inFile.bad())
    return false;

  // Read the dimensions of the mesh, i.e,. number of elementes in x, y, and z
  std::vector<std::size_t> dims;
  std::string line;
  getline(inFile, line);
  std::istringstream linestream(line);
  std::size_t dimVertices;
  while (linestream >> dimVertices)
  {
    dims.push_back(dimVertices);
  }

  // Compute the number of vertices, i.e., xdim * ydim * zdim
  unsigned short nDims = static_cast<unsigned short>(dims.size());
  std::size_t nVertices = static_cast<std::size_t>(
    std::accumulate(dims.begin(), dims.end(), 1, std::multiplies<std::size_t>()));

  // Print the mesh metadata
  if (rank == 0)
  {
    std::cout << "Number of dimensions: " << nDims << std::endl;
    std::cout << "Number of mesh vertices: " << nVertices << std::endl;
  }
  // Check the the number of dimensiosn is either 2D or 3D
  bool invalidNumDimensions = (nDims < 2 || nDims > 3);
  if (rank == 0)
  {
    if (invalidNumDimensions)
    {
      std::cout << "The input mesh is " << nDims << "D. Input data must be either 2D or 3D."
                << std::endl;
    }
  }
  if (invalidNumDimensions)
  {
    return false;
  }

  // Read data
  std::vector<ValueType> values(nVertices);
  for (std::size_t vertex = 0; vertex < nVertices; ++vertex)
  {
    inFile >> values[vertex];
  }

  // finish reading the data
  inFile.close();

  vtkm::cont::DataSetBuilderUniform dsb;
#ifndef VTKH_PARALLEL
  {
    // build the input dataset
    // 2D data
    if (nDims == 2)
    {
      vtkm::Id2 vdims;
      vdims[0] = static_cast<vtkm::Id>(dims[0]);
      vdims[1] = static_cast<vtkm::Id>(dims[1]);
      inDataSet = dsb.Create(vdims);
    }
    // 3D data
    else
    {
      vtkm::Id3 vdims;
      vdims[0] = static_cast<vtkm::Id>(dims[0]);
      vdims[1] = static_cast<vtkm::Id>(dims[1]);
      vdims[2] = static_cast<vtkm::Id>(dims[2]);
      inDataSet = dsb.Create(vdims);
    }
    vtkm::cont::DataSetFieldAdd dsf;
    dsf.AddPointField(inDataSet, "values", values);
  }
#else // VTKH_PARALLEL
  int numBlocks = size;
  vtkm::Id3 blocksPerDim =
    nDims == 3 ? vtkm::Id3(1, 1, numBlocks) : vtkm::Id3(1, numBlocks, 1); // Decompose the data into
  vtkm::Id3 globalSize = nDims == 3 ? vtkm::Id3(static_cast<vtkm::Id>(dims[0]),
                                                static_cast<vtkm::Id>(dims[1]),
                                                static_cast<vtkm::Id>(dims[2]))
                                    : vtkm::Id3(static_cast<vtkm::Id>(dims[0]),
                                                static_cast<vtkm::Id>(dims[1]),
                                                static_cast<vtkm::Id>(0));
  int blocksPerRank = 1;
  vtkm::cont::ArrayHandle<vtkm::Id3> localBlockIndices;
  vtkm::cont::ArrayHandle<vtkm::Id3> localBlockOrigins;
  vtkm::cont::ArrayHandle<vtkm::Id3> localBlockSizes;
  localBlockIndices.Allocate(blocksPerRank);
  localBlockOrigins.Allocate(blocksPerRank);
  localBlockSizes.Allocate(blocksPerRank);
  auto localBlockIndicesPortal = localBlockIndices.GetPortalControl();
  auto localBlockOriginsPortal = localBlockOrigins.GetPortalControl();
  auto localBlockSizesPortal = localBlockSizes.GetPortalControl();
  {
    vtkm::Id lastDimSize =
      (nDims == 2) ? static_cast<vtkm::Id>(dims[1]) : static_cast<vtkm::Id>(dims[2]);
    if (size > (lastDimSize / 2.))
    {
      if (rank == 0)
      {
        std::cout << "Number of ranks to large for data. Use " << lastDimSize / 2
                  << "or fewer ranks" << std::endl;
      }
      return false;
    }
    vtkm::Id standardBlockSize = (vtkm::Id)(lastDimSize / numBlocks);
    vtkm::Id blockSize = standardBlockSize;
    vtkm::Id blockSliceSize =
      nDims == 2 ? static_cast<vtkm::Id>(dims[0]) : static_cast<vtkm::Id>((dims[0] * dims[1]));
    vtkm::Id blockNumValues = blockSize * blockSliceSize;

    vtkm::Id startBlock = blocksPerRank * rank;
    vtkm::Id endBlock = startBlock + blocksPerRank;
    for (vtkm::Id blockIndex = startBlock; blockIndex < endBlock; ++blockIndex)
    {
      vtkm::Id localBlockIndex = blockIndex - startBlock;
      vtkm::Id blockStart = blockIndex * blockNumValues;
      vtkm::Id blockEnd = blockStart + blockNumValues;
      if (blockIndex < (numBlocks - 1)) // add overlap between regions
      {
        blockEnd += blockSliceSize;
      }
      else
      {
        blockEnd = lastDimSize * blockSliceSize;
      }
      vtkm::Id currBlockSize = (vtkm::Id)((blockEnd - blockStart) / blockSliceSize);

      vtkm::cont::DataSet ds;

      // 2D data
      if (nDims == 2)
      {
        vtkm::Id2 vdims;
        vdims[0] = static_cast<vtkm::Id>(currBlockSize);
        vdims[1] = static_cast<vtkm::Id>(dims[0]);
        vtkm::Vec<ValueType, 2> origin(0, blockIndex * blockSize);
        vtkm::Vec<ValueType, 2> spacing(1, 1);
        ds = dsb.Create(vdims, origin, spacing);

        localBlockIndicesPortal.Set(localBlockIndex, vtkm::Id3(blockIndex, 0, 0));
        localBlockOriginsPortal.Set(localBlockIndex,
                                    vtkm::Id3((blockStart / blockSliceSize), 0, 0));
        localBlockSizesPortal.Set(localBlockIndex,
                                  vtkm::Id3(currBlockSize, static_cast<vtkm::Id>(dims[0]), 0));
      }
      // 3D data
      else
      {
        vtkm::Id3 vdims;
        vdims[0] = static_cast<vtkm::Id>(dims[0]);
        vdims[1] = static_cast<vtkm::Id>(dims[1]);
        vdims[2] = static_cast<vtkm::Id>(currBlockSize);
        vtkm::Vec<ValueType, 3> origin(0, 0, (blockIndex * blockSize));
        vtkm::Vec<ValueType, 3> spacing(1, 1, 1);
        ds = dsb.Create(vdims, origin, spacing);

        localBlockIndicesPortal.Set(localBlockIndex, vtkm::Id3(0, 0, blockIndex));
        localBlockOriginsPortal.Set(localBlockIndex,
                                    vtkm::Id3(0, 0, (blockStart / blockSliceSize)));
        localBlockSizesPortal.Set(
          localBlockIndex,
          vtkm::Id3(static_cast<vtkm::Id>(dims[0]), static_cast<vtkm::Id>(dims[1]), currBlockSize));
      }

      std::vector<vtkm::Float32> subValues((values.begin() + blockStart),
                                           (values.begin() + blockEnd));

      vtkm::cont::DataSetFieldAdd dsf;
      dsf.AddPointField(ds, "values", subValues);
      inDataSet.AddBlock(ds);
      std::cout << "blockPerDim: " << blocksPerDim << std::endl
                << "globalSize: " << globalSize << std::endl
                << "localblockIndices: " << localBlockIndicesPortal.Get(localBlockIndex) << std::endl
                << "localBlockOrigins: " << localBlockOriginsPortal.Get(localBlockIndex) << std::endl
                << "localBlockSizes: " << localBlockSizesPortal.Get(localBlockIndex) << std::endl
                << std::endl;

    }
  }
#endif // VTKH_PARALLEL
  return true;
}


//----------------------------------------------------------------------------
TEST(vtkh_contour_tree, vtkh_contour_tree)
{
  vtkh::DataSet data_set;
#ifndef VTKH_PARALLEL
  vtkm::cont::DataSet vtkmData;
  ReadTestData("fuel.txt", vtkmData);
  data_set.AddDomain(vtkmData, 0);
#else
  MPI_Init(NULL, NULL);
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  vtkh::SetMPICommHandle(MPI_COMM_WORLD);
  vtkm::cont::MultiBlock mb;
  ReadTestData("fuel.txt", mb, rank, size);
  for (vtkm::Id id = 0; id < mb.GetNumberOfBlocks(); ++id)
  {
    data_set.AddDomain(mb.GetBlock(id), id);
  }
#endif
  vtkh::MarchingCubes marcher;
  marcher.SetInput(&data_set);
  marcher.SetField("values");

  int num_levels = 5;
  marcher.SetLevels(num_levels);
  marcher.SetUseContourTree(true);
  marcher.AddMapField("values");
  marcher.Update();
  std::vector<double> isoValues = marcher.GetIsoValues();
  std::sort(isoValues.begin(), isoValues.end());

  EXPECT_FLOAT_EQ(isoValues[0], 1e-05);
  EXPECT_FLOAT_EQ(isoValues[1], 82);
  EXPECT_FLOAT_EQ(isoValues[2], 133);
  EXPECT_FLOAT_EQ(isoValues[3], 168);
  EXPECT_FLOAT_EQ(isoValues[4], 177);
#ifdef VTKH_PARALLEL
  MPI_Finalize();
#endif
}
