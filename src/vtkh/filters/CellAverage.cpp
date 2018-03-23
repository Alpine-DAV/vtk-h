#include <vtkh/filters/CellAverage.hpp>
#include <vtkm/filter/CellAverage.h>

namespace vtkh 
{

CellAverage::CellAverage()
{

}

CellAverage::~CellAverage()
{

}

void 
CellAverage::SetField(const std::string &field_name)
{
  m_field_name = field_name;
}

void 
CellAverage::SetOutputField(const std::string &field_name)
{
  m_output_field_name = field_name;
}    

void CellAverage::PreExecute() 
{
  Filter::PreExecute();
  assert(m_field_name != "");
  assert(m_output_field_name != "");
}

void CellAverage::PostExecute()
{
  Filter::PostExecute();
}

void CellAverage::DoExecute()
{
  this->m_output = new DataSet();
  const int num_domains = this->m_input->GetNumberOfDomains();
  
  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);

    if(!dom.HasField(m_field_name))
    {
      continue;
    }

    vtkm::filter::CellAverage avg;
    avg.SetOutputFieldName(m_output_field_name);
    avg.SetFieldsToPass(this->GetFieldSelection());
    avg.SetActiveField(m_field_name);

    auto dataset = avg.Execute(dom);
    m_output->AddDomain(dataset, domain_id);
  }
}

std::string
CellAverage::GetName() const
{
  return "vtkh::CellAverage";
}

} //  namespace vtkh
