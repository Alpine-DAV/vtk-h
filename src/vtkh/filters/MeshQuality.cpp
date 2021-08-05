#include <vtkh/filters/MeshQuality.hpp>
#include <vtkh/vtkm_filters/vtkmMeshQuality.hpp>
#include <vtkh/Error.hpp>

namespace vtkh
{

MeshQuality::MeshQuality()
  : m_metric(vtkm::filter::CellMetric::VOLUME)
{

}

MeshQuality::~MeshQuality()
{

}

void MeshQuality::cell_metric(vtkm::filter::CellMetric metric)
{
  m_metric = metric;
}

void MeshQuality::PreExecute()
{
  Filter::PreExecute();
  if(!m_input->IsUnstructured())
  {
    throw Error("Mesh quality requires that meshes be completely unstructured");
  }
}

void MeshQuality::PostExecute()
{
  Filter::PostExecute();
}

void MeshQuality::DoExecute()
{
  this->m_output = new DataSet();
  const int num_domains = this->m_input->GetNumberOfDomains();

  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);
    vtkmMeshQuality quali;
    vtkm::cont::DataSet res = quali.Run(dom, m_metric, this->GetFieldSelection());
    m_output->AddDomain(res, domain_id);
  }
}

std::string
MeshQuality::GetName() const
{
  return "vtkh::MeshQuality";
}

} //  namespace vtkh
