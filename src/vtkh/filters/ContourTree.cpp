#include <vtkh/filters/ContourTree.hpp>
// vtkm includes
#include <vtkm/cont/DeviceAdapter.h>
#include <vtkm/cont/Storage.h>
#include <vtkm/internal/Configure.h>
#include <vtkm/filter/ContourTreeUniformPPP2.h>
#include <vtkm/worklet/contourtree_ppp2/PrintVectors.h>
#include <vtkm/worklet/contourtree_ppp2/ProcessContourTree.h>
#include <vtkm/worklet/contourtree_ppp2/ProcessContourTree_Inc/Branch.h>
#include <vtkm/worklet/contourtree_ppp2/ProcessContourTree_Inc/PiecewiseLinearFunction.h>

namespace cppp2_ns = vtkm::worklet::contourtree_ppp2;
//using DataValueType = vtkm::Float32;
using DataValueType = vtkm::Float64;
using ValueArray = vtkm::cont::ArrayHandle<DataValueType, vtkm::cont::StorageTagBasic>;
using BranchType = vtkm::worklet::contourtree_ppp2::process_contourtree_inc::Branch<DataValueType>;
using PLFType = vtkm::worklet::contourtree_ppp2::process_contourtree_inc::PiecewiseLinearFunction<DataValueType>;

namespace vtkh 
{

namespace detail
{

template<typename DeviceAdapter>
void ComputeContourValues(vtkm::cont::DataSet &inDataSet, 
                          std::string fieldName,                       
                          vtkm::filter::ContourTreePPP2 &filter)
{

  vtkm::Id numLevels = 5;             // Number of iso levels to be selected
  DataValueType eps = 0.00001;        // Error away from critical point
  vtkm::Id numComp = numLevels + 1;   // Number of components the tree should be simplified to
  vtkm::Id contourType = 0;           // Approach to be used to select contours based on the tree
  vtkm::Id contourSelectMethod = 0;          // Method to be used to compute the relevant iso values
  bool usePersistenceSorter = true;
  ////////////////////////////////////////////
  // Compute the branch decomposition
  ////////////////////////////////////////////
  // compute the volume for each hyperarc and superarc
  cppp2_ns::IdArrayType superarcIntrinsicWeight;
  cppp2_ns::IdArrayType superarcDependentWeight;
  cppp2_ns::IdArrayType supernodeTransferWeight;
  cppp2_ns::IdArrayType hyperarcDependentWeight;

  cppp2_ns::ProcessContourTree::ComputeVolumeWeights<DeviceAdapter>(
          filter.GetContourTree(),
          filter.GetNumIterations(),
          superarcIntrinsicWeight,  // (output)
          superarcDependentWeight,  // (output)
          supernodeTransferWeight,  // (output)
          hyperarcDependentWeight); // (output)

  // compute the branch decomposition by volume
  cppp2_ns::IdArrayType whichBranch;
  cppp2_ns::IdArrayType branchMinimum;
  cppp2_ns::IdArrayType branchMaximum;
  cppp2_ns::IdArrayType branchSaddle;
  cppp2_ns::IdArrayType branchParent;

  cppp2_ns::ProcessContourTree::ComputeVolumeBranchDecomposition<DeviceAdapter>(
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
  BranchType* branchDecompostionRoot = cppp2_ns::ProcessContourTree::ComputeBranchDecomposition<DataValueType>(
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
  std::vector<DataValueType> isoValues;
  switch(contourSelectMethod)
  {
    default:
    case 0:
      {
        branchDecompostionRoot->getRelevantValues(contourType, eps, isoValues);
      }
      break;
    case 1:
      {
        PLFType plf;
        branchDecompostionRoot->accumulateIntervals(contourType, eps, plf);
        isoValues = plf.nLargest(numLevels);
      }
      break;
  }

  // Print the compute iso values
  std::sort(isoValues.begin(), isoValues.end());
  std::cout << "Isovalues: ";
  for (DataValueType val : isoValues) std::cout << val << " ";
  std::cout << std::endl;

  // Unique isovalues
  std::vector<DataValueType>::iterator it = std::unique (isoValues.begin(), isoValues.end());
  isoValues.resize( std::distance(isoValues.begin(),it) );

  std::cout << isoValues.size() << "  Unique Isovalues: ";
  for (DataValueType val : isoValues) std::cout << val << " ";
  std::cout << std::endl;

  std::cout<<"Acrs : " << filter.GetContourTree().arcs.GetNumberOfValues() <<std::endl;
};

struct ComputeCaller
{
   
  template <typename Device>
  VTKM_CONT bool operator()(Device, 
                            vtkm::cont::DataSet &in,
                            std::string fieldName,
                            vtkm::filter::ContourTreePPP2 &filter) const
  {
    VTKM_IS_DEVICE_ADAPTER_TAG(Device);
    ComputeContourValues<Device>(in, fieldName, filter);
    return true;
  }
};

} // namespace detail

ContourTree::ContourTree()
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

    filter.SetActiveField(m_field_name);
    result = filter.Execute(dom);
    vtkm::cont::TryExecute(detail::ComputeCaller(), result, m_field_name, filter);
  }
}

std::string
ContourTree::GetName() const
{
  return "vtkh::ContourTree";
}

} //  namespace vtkh
