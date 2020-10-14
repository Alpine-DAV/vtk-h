#include <vtkh/filters/ExtractStructured.hpp>
#include <vtkh/vtkm_filters/vtkmExtractStructured.hpp>

namespace vtkh
{

ExtractStructured::ExtractStructured()
{

}

ExtractStructured::~ExtractStructured()
{

}

void
ExtractStructured::SetField(const std::string &field_name)
{
  m_field_name = field_name;
}

void ExtractStructured::PreExecute()
{
  Filter::PreExecute();
  Filter::CheckForRequiredField(m_field_name);
}

void ExtractStructured::PostExecute()
{
  Filter::PostExecute();
}

void ExtractStructured::DoExecute()
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
  }
}

std::string
ExtractStructured::GetName() const
{
  return "vtkh::ExtractStructured";
}

} //  namespace vtkh
