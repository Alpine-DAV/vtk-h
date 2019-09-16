#include <vtkh/filters/Filter.hpp>
#include <vtkh/Error.hpp>
#include <vtkh/Logger.hpp>

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
  VTKH_DATA_OPEN(this->GetName());
#ifdef VTKH_ENABLE_LOGGING
#endif
  PreExecute();
  DoExecute();
  PostExecute();
  VTKH_DATA_CLOSE();
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
    msg<<"Input for vtkh filter '"<<this->GetName()<<"' is null.";
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
Filter::CheckForRequiredField(const std::string &field_name)
{
  if(m_input == nullptr)
  {
    std::stringstream msg;
    msg<<"Cannot verify required field '"<<field_name;
    msg<<"' for vkth filter '"<<this->GetName()<<"' because input is null.";
    throw Error(msg.str());
  }

  if(!m_input->GlobalFieldExists(field_name))
  {
    std::stringstream msg;
    msg<<"Required field '"<<field_name;
    msg<<"' for vkth filter '"<<this->GetName()<<"' does not exist";
    throw Error(msg.str());
  }
}

void
Filter::PropagateMetadata()
{
  m_output->SetCycle(m_input->GetCycle());
}


vtkm::filter::FieldSelection
Filter::GetFieldSelection() const
{
  vtkm::filter::FieldSelection sel;
  for (const auto& str : this->m_map_fields)
  {
    sel.AddField(str);
  }
  return sel;
}


} //namespace vtkh
