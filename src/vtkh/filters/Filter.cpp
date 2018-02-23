#include <vtkh/filters/Filter.hpp>
#include <vtkh/Error.hpp>

namespace vtkh
{

Filter::Filter() 
{ 
  m_input = nullptr; 
  m_output = nullptr; 
}

Filter::~Filter() 
{ 
};

void 
Filter::SetInput(DataSet *input)
{ 
  m_input = input; 
}

DataSet*
Filter::GetOutput()
{ 
  return m_output; 
}

DataSet* 
Filter::Update()
{
  PreExecute();
  DoExecute();
  PostExecute();
  return m_output;
}

void 
Filter::AddMapField(const std::string &field_name)
{
  m_map_fields.push_back(field_name);
}

void 
Filter::ClearMapFields()
{
  m_map_fields.clear();  
}

void 
Filter::PreExecute()
{
  if(m_input == nullptr)
  {
    std::stringstream msg;
    msg<<"Input for filter '"<<this->GetName()<<"' is null.";
    throw Error(msg.str());
  }

  if(m_map_fields.size() == 0)
  {
    this->MapAllFields(); 
  }
  
};

void 
Filter::PostExecute()
{
  this->PropagateMetadata();
};

void 
Filter::MapAllFields()
{
  if(m_input->GetNumberOfDomains() > 0)
  {
    vtkm::cont::DataSet dom = m_input->GetDomain(0);
    vtkm::IdComponent num_fields = dom.GetNumberOfFields();  
    for(vtkm::IdComponent i = 0; i < num_fields; ++i)
    {
      std::string field_name = dom.GetField(i).GetName();
      m_map_fields.push_back(field_name);
    }
  }
}

void 
Filter::PropagateMetadata()
{
  m_output->SetCycle(m_input->GetCycle());
}

} //namespace vtkh
