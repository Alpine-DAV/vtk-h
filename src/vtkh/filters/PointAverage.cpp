#include <vtkh/filters/PointAverage.hpp>
#include <vtkm/filter/PointAverage.h>

namespace vtkh 
{

PointAverage::PointAverage()
{

}

PointAverage::~PointAverage()
{

}

void 
PointAverage::SetField(const std::string &field_name)
{
  m_field_name = field_name;
}

void 
PointAverage::SetOutputField(const std::string &field_name)
{
  m_output_field_name = field_name;
}    

void PointAverage::PreExecute() 
{
  Filter::PreExecute();
  assert(m_field_name != "");
  assert(m_output_field_name != "");
}

void PointAverage::PostExecute()
{
  Filter::PostExecute();
}

void PointAverage::DoExecute()
{
  this->m_output = new DataSet();
  const int num_domains = this->m_input->GetNumberOfDomains();
  
  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);

    vtkm::filter::PointAverage avg;
    avg.SetOutputFieldName(m_output_field_name);
    vtkm::filter::Result res = avg.Execute(dom, m_field_name);
    
    m_output->AddDomain(res.GetDataSet(), domain_id);
  }
}

std::string
PointAverage::GetName() const
{
  return "vtkh::PointAverage";
}

} //  namespace vtkh
