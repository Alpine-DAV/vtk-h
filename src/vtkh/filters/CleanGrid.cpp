
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
  this->MapAllFields();
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

    vtkm::filter::Result res = cleaner.Execute(dom);
    for(size_t f = 0; f < m_map_fields.size(); ++f)
    {
      cleaner.MapFieldOntoOutput(res, dom.GetField(m_map_fields[f]));
    }
    this->m_output->AddDomain(res.GetDataSet(), domain_id);
  }

  this->PropagateMetadata();
}

void
CleanGrid::PostExecute()
{

}

} // namespace vtkh
