//#include <vtkm/filter/your_vtkm_filter.h>
#include <vtkh/filters/PathTrace.hpp>

#ifdef VTKH_PARALLEL
#include <mpi.h>
#endif

namespace vtkh
{

PathTrace::PathTrace()
{

}

PathTrace::~PathTrace()
{

}

void
PathTrace::SetField(const std::string &field_name)
{
  m_field_name = field_name;
}

void
PathTrace::SetCamera(const vtkmCamera &camera)
{
  m_camera= camera;
}

void PathTrace::PreExecute()
{
  Filter::PreExecute();

  // Build the bounds map so we can ask questions
  m_bounds_map.Clear();
  const int num_domains = this->m_input->GetNumberOfDomains();
  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);
    vtkm::Bounds bounds = dom.GetCoordinateSystem().GetBounds();
    m_bounds_map.AddBlock(domain_id, bounds);

  }

}

void PathTrace::PostExecute()
{
  Filter::PostExecute();
}

void PathTrace::DoExecute()
{
  int rank = 0;
  int procs = 1;

#ifdef VTKH_PARALLEL
  MPI_Comm comm = MPI_Comm_f2c(vtkh::GetMPICommHandle());
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &procs);
#endif

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

void PathTrace::CreateRays(vtkmRay &rays)
{

}

std::string
PathTrace::GetName() const
{
  return "vtkh::PathTrace";
}

} //  namespace vtkh
