#include <vtkh/filters/CommTest.hpp>

#include "communication/BoundsMap.hpp"
#ifdef VTKH_PARALLEL
#include "communication/avtParICAlgorithm.hpp"
#endif

namespace vtkh
{

CommTest::CommTest()
{

}

CommTest::~CommTest()
{

}

void
CommTest::SetField(const std::string &field_name)
{
  m_field_name = field_name;
}

void CommTest::PreExecute()
{
  Filter::PreExecute();
}

void CommTest::PostExecute()
{
  Filter::PostExecute();
}

void CommTest::DoExecute()
{
  this->m_output = new DataSet();
  const int num_domains = this->m_input->GetNumberOfDomains();
  BoundsMap bmap;
#ifdef VTKH_PARALLEL
#endif
  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);
    bmap.AddBlock(domain_id, dom.GetCoordinateSystem().GetBounds());
  }
  bmap.Build();


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
CommTest::GetName() const
{
  return "vtkh::CommTest";
}

} //  namespace vtkh
