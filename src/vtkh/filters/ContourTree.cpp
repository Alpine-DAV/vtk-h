#include <vtkh/filters/ContourTree.hpp>
// vtkm includes
#include <vtkm/cont/DeviceAdapter.h>
#include <vtkm/cont/Storage.h>
#include <vtkm/internal/Configure.h>
#include <vtkm/filter/ContourTreeUniformAugmented.h>
#include <vtkm/worklet/contourtree_augmented/PrintVectors.h>
#include <vtkm/worklet/contourtree_augmented/ProcessContourTree.h>
#include <vtkm/worklet/contourtree_augmented/processcontourtree/Branch.h>
#include <vtkm/worklet/contourtree_augmented/processcontourtree/PiecewiseLinearFunction.h>

namespace caugmented_ns = vtkm::worklet::contourtree_augmented;
using DataValueType = vtkm::Float64;
using ValueArray = vtkm::cont::ArrayHandle<DataValueType>;
using BranchType = vtkm::worklet::contourtree_augmented::process_contourtree_inc::Branch<DataValueType>;
using PLFType = vtkm::worklet::contourtree_augmented::process_contourtree_inc::PiecewiseLinearFunction<DataValueType>;

#ifdef VTKH_PARALLEL
#include <mpi.h>
#endif


namespace vtkh
{

ContourTree::ContourTree()
  : m_levels(5)
{

}

ContourTree::~ContourTree()
{

}

void
ContourTree::SetField(const std::string &field_name)
{
  m_field_name = field_name;
}

void
ContourTree::SetNumLevels(int levels)
{
  m_levels = levels;
}

std::vector<double>
ContourTree::GetIsoValues()
{
  return m_iso_values;
}

void ContourTree::PreExecute()
{
  Filter::PreExecute();
}

void ContourTree::PostExecute()
{
  Filter::PostExecute();
}

void ContourTree::DoExecute()
{
  int rank = 0;
  this->m_output = new DataSet();
  const int num_domains = this->m_input->GetNumberOfDomains();
#ifndef VTKH_PARALLEL
  assert(num_domains == 1);
  vtkm::cont::DataSet inDataSet;
  vtkm::Id domain_id;
  this->m_input->GetDomain(0, inDataSet, domain_id);
  this->m_output->AddDomain(inDataSet, domain_id);
#else // VTKH_PARALLEL
  int size;
  MPI_Comm mpi_comm = MPI_Comm_f2c(vtkh::GetMPICommHandle());
  vtkm::cont::EnvironmentTracker::SetCommunicator(vtkmdiy::mpi::communicator(mpi_comm));
  MPI_Comm_size(mpi_comm, &size);
  MPI_Comm_rank(mpi_comm, &rank);
  vtkm::cont::PartitionedDataSet inDataSet;
  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);
    inDataSet.AppendPartition(dom);
    std::ostringstream ostr;
    ostr << "rank: " << rank
         << " coord system range: " << dom.GetCoordinateSystem(0).GetRange() << std::endl;
    std::cout << ostr.str();
  }

  // TODO: change hardcode to computation.
  vtkm::Id3 blocksPerDim = vtkm::Id3 (1, 1, 2);
  vtkm::Id3 globalSize = vtkm::Id3 (64, 64, 64);
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
  if (rank == 0)
  {
    localBlockIndicesPortal.Set(0, vtkm::Id3 (0, 0, 0));
    localBlockOriginsPortal.Set(0, vtkm::Id3 (0, 0, 0));
    localBlockSizesPortal.Set(0, vtkm::Id3(64, 64, 33));
  }
  else
  {
    localBlockIndicesPortal.Set(0, vtkm::Id3 (0, 0, 1));
    localBlockOriginsPortal.Set(0, vtkm::Id3 (0, 0, 32));
    localBlockSizesPortal.Set(0, vtkm::Id3(64, 64, 32));
  }

  // compute

#endif // VTKH_PARALLEL
  std::cout<<"RUNNING TREE\n";
  bool useMarchingCubes = false;
  // Compute the fully augmented contour tree.
  // This should always be true for now in order for the isovalue selection to work.
  bool computeRegularStructure = true;
  //Convert the mesh of values into contour tree, pairs of vertex ids
#ifdef VTKH_PARALLEL
    {
      std::ostringstream ostr;
      ostr << "global bounds: " << vtkm::cont::BoundsGlobalCompute(inDataSet) << std::endl
           << "useMarchingCubes: " << useMarchingCubes << std::endl
           << "computeRegularStructure: " << computeRegularStructure << std::endl;
      std::cout << ostr.str();

      ostr.str("");
      ostr << "blockPerDim: " << blocksPerDim << std::endl
           << "globalSize: " << globalSize << std::endl
           << "localblockIndices: " << localBlockIndicesPortal.Get(0) << std::endl
           << "localBlockOrigins: " << localBlockOriginsPortal.Get(0) << std::endl
           << "localBlockSizes: " << localBlockSizesPortal.Get(0) << std::endl
           << std::endl;
      std::cout << ostr.str();
    }
#endif
  vtkm::filter::ContourTreePPP2 filter(useMarchingCubes,
                                       computeRegularStructure);
  std::vector<DataValueType> iso_values;
#ifdef VTKH_PARALLEL
  filter.SetSpatialDecomposition(
    blocksPerDim, globalSize, localBlockIndices, localBlockOrigins, localBlockSizes);
#endif // VTKH_PARALLEL
  filter.SetActiveField(m_field_name);
  auto result = filter.Execute(inDataSet);

  m_iso_values.resize(m_levels);
  if (rank == 0)
  {
    DataValueType eps = 0.00001;        // Error away from critical point
    vtkm::Id numComp = m_levels + 1;    // Number of components the tree should be simplified to
    vtkm::Id contourType = 0;           // Approach to be used to select contours based on the tree
    vtkm::Id contourSelectMethod = 0;   // Method to be used to compute the relevant iso values
    bool usePersistenceSorter = true;
    ////////////////////////////////////////////
    // Compute the branch decomposition
    ////////////////////////////////////////////
    // compute the volume for each hyperarc and superarc
    caugmented_ns::IdArrayType superarcIntrinsicWeight;
    caugmented_ns::IdArrayType superarcDependentWeight;
    caugmented_ns::IdArrayType supernodeTransferWeight;
    caugmented_ns::IdArrayType hyperarcDependentWeight;

    caugmented_ns::ProcessContourTree::ComputeVolumeWeights(
      filter.GetContourTree(),
      filter.GetNumIterations(),
      superarcIntrinsicWeight,  // (output)
      superarcDependentWeight,  // (output)
      supernodeTransferWeight,  // (output)
      hyperarcDependentWeight); // (output)

    // compute the branch decomposition by volume
    caugmented_ns::IdArrayType whichBranch;
    caugmented_ns::IdArrayType branchMinimum;
    caugmented_ns::IdArrayType branchMaximum;
    caugmented_ns::IdArrayType branchSaddle;
    caugmented_ns::IdArrayType branchParent;

    caugmented_ns::ProcessContourTree::ComputeVolumeBranchDecomposition(
      filter.GetContourTree(),
      superarcDependentWeight,
      superarcIntrinsicWeight,
      whichBranch,               // (output)
      branchMinimum,             // (output)
      branchMaximum,             // (output)
      branchSaddle,              // (output)
      branchParent);             // (output)

    // create explicit representation of the branch decompostion from the array representation
#ifndef VTKH_PARALLEL
    ValueArray dataField;
    result.GetField(0).GetData().CopyTo(dataField);
    bool dataFieldIsSorted = false;
#else // VTKH_PARALLEL
    ValueArray dataField;
    result.GetPartition(0).GetField(0).GetData().CopyTo(dataField);
    bool dataFieldIsSorted = true;
#endif // VTKH_PARALLEL
    BranchType* branchDecompostionRoot = caugmented_ns::ProcessContourTree::ComputeBranchDecomposition<DataValueType>(
      filter.GetContourTree().superparents,
      filter.GetContourTree().supernodes,
      whichBranch,
      branchMinimum,
      branchMaximum,
      branchSaddle,
      branchParent,
      filter.GetSortOrder(),
      dataField,
      dataFieldIsSorted
    );

    // Simplify the contour tree of the branch decompostion
    branchDecompostionRoot->simplifyToSize(numComp, usePersistenceSorter);

    // Compute the relevant iso-values
    switch(contourSelectMethod)
    {
    default:
    case 0:
      {
        branchDecompostionRoot->getRelevantValues(contourType, eps, iso_values);
      }
      break;
    case 1:
      {
        PLFType plf;
        branchDecompostionRoot->accumulateIntervals(contourType, eps, plf);
        iso_values = plf.nLargest(m_levels);
      }
      break;
    }

    // Print the compute iso values
    std::sort(iso_values.begin(), iso_values.end());
    std::cout << "Isovalues: ";
    for (DataValueType val : iso_values) std::cout << val << " ";
    std::cout << std::endl;

    // Unique isovalues
    std::vector<DataValueType>::iterator it = std::unique (iso_values.begin(), iso_values.end());
    iso_values.resize( std::distance(iso_values.begin(),it) );

    std::cout << iso_values.size() << "  Unique Isovalues: ";
    for (DataValueType val : iso_values) std::cout << val << " ";
    std::cout << std::endl;

    std::cout<<"Acrs : " << filter.GetContourTree().arcs.GetNumberOfValues() <<std::endl;




    for(size_t x = 0; x < iso_values.size(); ++x)
    {
      m_iso_values[x] = iso_values[x];
    }
  }
#ifdef VTKH_PARALLEL
  MPI_Bcast(&m_iso_values[0], m_levels, MPI_DOUBLE, 0, mpi_comm);
#endif // VTKH_PARALLEL
}

std::string
ContourTree::GetName() const
{
  return "vtkh::ContourTree";
}

} //  namespace vtkh
