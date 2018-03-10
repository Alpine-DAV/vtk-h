
#include <vtkh/filters/CleanGrid.hpp>
#include <vtkh/Error.hpp>

#include <vtkm/filter/CleanGrid.h>

namespace vtkh
{


CleanGrid::CleanGrid()
{

}

CleanGrid::~CleanGrid()
{

}

void
CleanGrid::PreExecute()
{
  Filter::PreExecute(); 
}

void
CleanGrid::DoExecute()
{
  this->m_output = new DataSet();
  vtkm::filter::CleanGrid cleaner;

  const int num_domains = this->m_input->GetNumberOfDomains(); 
  int valid = 0;
  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);
    auto dataset = cleaner.Execute(dom, this->GetFieldSelection());
    this->m_output->AddDomain(dataset, domain_id);
  }

}

void
CleanGrid::PostExecute()
{
  Filter::PostExecute();
}

std::string 
CleanGrid::GetName() const
{
  return "vtkh::CleanGrid";
}

} // namespace vtkh
