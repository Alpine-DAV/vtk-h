#include "Threshold.hpp"
#include <vtkh/Error.hpp>


#include <vtkm/cont/CellSetPermutation.h>
#include <vtkm/cont/TryExecute.h>
#include <vtkm/filter/Threshold.h>
#include <vtkm/worklet/CellDeepCopy.h>

namespace vtkh 
{

namespace detail
{
typedef vtkm::cont::CellSetPermutation<vtkm::cont::CellSetStructured<2>>
  PermStructured2d; 

typedef vtkm::cont::CellSetPermutation<vtkm::cont::CellSetStructured<3>>
  PermStructured3d;

typedef vtkm::cont::CellSetPermutation<vtkm::cont::CellSetExplicit<>>
  PermExplicit;

typedef  vtkm::cont::CellSetPermutation<vtkm::cont::CellSetSingleType<>>
  PermExplicitSingle;
//
// Theshold outputs CellSetPermutations which cannot
// be consumed by anything else in vtkm, so we need 
// to explicitly do a deep copy and make the cell set 
// explicit
//
template <typename CellSetType>
struct DeepCopy
{
  const CellSetType &m_input;
  vtkm::cont::CellSetExplicit<> &m_output;

  DeepCopy(CellSetType &input,
           vtkm::cont::CellSetExplicit<> &output)
    : m_input(input)
    , m_output(output)
  {
  }

  template <typename Device>
  bool operator()(Device device)
  { 
    m_output = vtkm::worklet::CellDeepCopy::Run(m_input, Device());
    return true;
  }
};

void StripPermutation(vtkm::cont::DataSet &data_set)
{
  vtkm::cont::DynamicCellSet cell_set = data_set.GetCellSet(); 
  vtkm::cont::DataSet result;
  vtkm::cont::CellSetExplicit<> explicit_cells;

  if(cell_set.IsSameType(PermStructured2d()))
  {
    PermStructured2d perm = cell_set.Cast<PermStructured2d>();
    DeepCopy<PermStructured2d> functor(perm, explicit_cells);
    vtkm::cont::TryExecute(functor);
  }
  else if(cell_set.IsSameType(PermStructured3d()))
  {
    PermStructured3d perm = cell_set.Cast<PermStructured3d>();
    DeepCopy<PermStructured3d> functor(perm, explicit_cells);
    vtkm::cont::TryExecute(functor);
  }
  else if(cell_set.IsSameType(PermExplicit()))
  {
    PermExplicit perm = cell_set.Cast<PermExplicit>();
    DeepCopy<PermExplicit> functor(perm, explicit_cells);
    vtkm::cont::TryExecute(functor);
  }
  else if(cell_set.IsSameType(PermExplicitSingle()))
  {
    PermExplicitSingle perm = cell_set.Cast<PermExplicitSingle>();
    DeepCopy<PermExplicitSingle> functor(perm, explicit_cells);
    vtkm::cont::TryExecute(functor);
  }
  
  result.AddCellSet(explicit_cells);

  vtkm::Id num_coords = data_set.GetNumberOfCoordinateSystems();
  for(vtkm::Id i = 0; i < num_coords; ++i)
  {
    result.AddCoordinateSystem(data_set.GetCoordinateSystem(i));
  }

  vtkm::Id num_fields = data_set.GetNumberOfFields();
  for(vtkm::Id i = 0; i < num_fields; ++i)
  {
    result.AddField(data_set.GetField(i));
  }
   
  data_set = result;
}

} // namespace detail

Threshold::Threshold()
{
}

Threshold::~Threshold()
{

}

void 
Threshold::SetUpperThreshold(const double &value)
{
  m_range.Max = value;
}

void 
Threshold::SetLowerThreshold(const double &value)
{
  m_range.Min = value;
}

void 
Threshold::SetField(const std::string &field_name)
{
  m_field_name = field_name;
}

double 
Threshold::GetUpperThreshold() const
{
  return m_range.Max;
}

double 
Threshold::GetLowerThreshold() const
{
  return m_range.Min;
}

std::string
Threshold::GetField() const
{
  return m_field_name;
}

void Threshold::PreExecute() 
{
  Filter::PreExecute();
}

void Threshold::PostExecute()
{
  Filter::PostExecute();
}

void Threshold::DoExecute()
{
  
  this->m_output = new DataSet();

  const int num_domains = this->m_input->GetNumberOfDomains(); 

  vtkm::filter::Threshold thresholder;
  thresholder.SetUpperThreshold(m_range.Max);
  thresholder.SetLowerThreshold(m_range.Min);

  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);

    vtkm::filter::Result res = thresholder.Execute(dom, m_field_name);

    for(size_t f = 0; f < m_map_fields.size(); ++f)
    {
      thresholder.MapFieldOntoOutput(res, dom.GetField(m_map_fields[f]));
    }

    vtkm::cont::DataSet data_set = res.GetDataSet();
    detail::StripPermutation(data_set);
    this->m_output->AddDomain(data_set, domain_id);
    
  }
}

} //  namespace vtkh
