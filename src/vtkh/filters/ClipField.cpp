#include "ClipField.hpp"

#include <vtkm/filter/ClipWithField.h>
#include <vtkm/filter/PointAverage.h>

namespace vtkh 
{

ClipField::ClipField()
  : m_clip_value(0.0),
    m_invert(false)
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
ClipField::SetInvertClip(const bool invert)
{
  m_invert = invert;
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
  clipper.SetInvertClip(m_invert);

  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);

    if(!dom.HasField(m_field_name))
    {
      continue;
    }
    
    //bool isCellAssoc = dom.GetField(m_field_name).GetAssociation() == 
    //                   vtkm::cont::Field::ASSOC_CELL_SET; 
    //if(isCellAssoc)
    //{
    //  vtkm::cont::DataSet tmp;
    //  tmp.AddCellSet(dom.GetCellSet());
    //  //std::string field_name = m_field_name + "_clip_point_centered";
    //  vtkm::filter::PointAverage avg;
    //  avg.SetOutputFieldName(m_field_name);
    //  avg.SetActiveField(m_field_name);
    //  auto temp = avg.Execute(dom);
    //}
      

    clipper.SetActiveField(m_field_name);
    auto dataset = clipper.Execute(dom, this->GetFieldSelection());
    this->m_output->AddDomain(dataset, domain_id);
  }
}

std::string
ClipField::GetName() const
{
  return "vtkh::ClipField";
}

} //  namespace vtkh
