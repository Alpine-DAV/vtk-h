#include "ClipField.hpp"

#include <vtkm/filter/ClipWithField.h>

namespace vtkh 
{

ClipField::ClipField()
  : m_clip_value(0.0)
{

}

ClipField::~ClipField()
{

}

void 
ClipField::SetClipValue(const vtkm::Float64 clip_value)
{
  m_clip_value = clip_value;
}

void 
ClipField::SetField(const std::string field_name)
{
  m_field_name = field_name;
}

void 
ClipField::PreExecute() 
{
  Filter::PreExecute();
}

void 
ClipField::PostExecute()
{
  Filter::PostExecute();
}

void ClipField::DoExecute()
{
  
  this->m_output = new DataSet();

  const int num_domains = this->m_input->GetNumberOfDomains(); 

  vtkm::filter::ClipWithField clipper;
  clipper.SetClipValue(m_clip_value);

  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);

    vtkm::filter::Result res = clipper.Execute(dom, m_field_name);

    for(size_t f = 0; f < m_map_fields.size(); ++f)
    {
      clipper.MapFieldOntoOutput(res, dom.GetField(m_map_fields[f]));
    }
    this->m_output->AddDomain(res.GetDataSet(), domain_id);
  }
}

std::string
ClipField::GetName() const
{
  return "vtkh::ClipField";
}

} //  namespace vtkh
