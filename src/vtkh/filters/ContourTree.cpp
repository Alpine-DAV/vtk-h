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
//using DataValueType = vtkm::Float32;
using DataValueType = vtkm::Float64;
using ValueArray = vtkm::cont::ArrayHandle<DataValueType, vtkm::cont::StorageTagBasic>;
using BranchType = vtkm::worklet::contourtree_augmented::process_contourtree_inc::Branch<DataValueType>;
using PLFType = vtkm::worklet::contourtree_augmented::process_contourtree_inc::PiecewiseLinearFunction<DataValueType>;

namespace vtkh
{

namespace detail
{

template<typename DeviceAdapter>
void ComputeContourValues(vtkm::cont::DataSet &inDataSet,
                          std::string fieldName,
                          vtkm::filter::ContourTreePPP2 &filter,
                          vtkm::Id numLevels,
                          std::vector<DataValueType> &iso_values)
{

  DataValueType eps = 0.00001;        // Error away from critical point
  vtkm::Id numComp = numLevels + 1;   // Number of components the tree should be simplified to
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
  ValueArray vtkmValues = inDataSet.GetField(fieldName).GetData().CastToTypeStorage<DataValueType, vtkm::cont::StorageTagBasic>();
  BranchType* branchDecompostionRoot = caugmented_ns::ProcessContourTree::ComputeBranchDecomposition(
          filter.GetContourTree().superparents,
          filter.GetContourTree().supernodes,
          whichBranch,
          branchMinimum,
          branchMaximum,
          branchSaddle,
          branchParent,
          filter.GetSortOrder(),
          vtkmValues
    );

  // Simplify the contour tree of the branch decompostion
  branchDecompostionRoot->simplifyToSize(numComp, usePersistenceSorter);
#ifdef DEBUG_PRINT
  branchDecompostionRoot->print(std::cout);
#endif

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
        iso_values = plf.nLargest(numLevels);
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
};

struct ComputeCaller
{

  template <typename Device>
  VTKM_CONT bool operator()(Device,
                            vtkm::cont::DataSet &in,
                            std::string fieldName,
                            vtkm::filter::ContourTreePPP2 &filter,
                            vtkm::Id numLevels,
                            std::vector<DataValueType> &iso_values) const
  {
    VTKM_IS_DEVICE_ADAPTER_TAG(Device);
    ComputeContourValues<Device>(in, fieldName, filter, numLevels, iso_values);
    return true;
  }
};

} // namespace detail

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
  this->m_output = new DataSet();
  const int num_domains = this->m_input->GetNumberOfDomains();

  for(int i = 0; i < num_domains; ++i)
  {
    std::cout<<"RUNNING TREE\n";
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);
    // insert interesting stuff
    m_output->AddDomain(dom, domain_id);
    bool useMarchingCubes = true;
    // Compute the fully augmented contour tree.
    // This should always be true for now in order for the isovalue selection to work.
    bool computeRegularStructure = true;
    vtkm::cont::DataSet result;
    //Convert the mesh of values into contour tree, pairs of vertex ids
    vtkm::filter::ContourTreePPP2 filter(useMarchingCubes,
                                         computeRegularStructure);
    std::vector<DataValueType> iso_values;
    filter.SetActiveField(m_field_name);
    result = filter.Execute(dom);
    vtkm::cont::TryExecute(detail::ComputeCaller(), result, m_field_name, filter, m_levels, iso_values);
    m_iso_values.resize(iso_values.size());
    for(size_t x = 0; x < iso_values.size(); ++x)
    {
      m_iso_values[x] = iso_values[x];
    }
  }
}

std::string
ContourTree::GetName() const
{
  return "vtkh::ContourTree";
}

} //  namespace vtkh
