//#include <vtkm/filter/your_vtkm_filter.h>
#include <vtkh/filters/PathTrace.hpp>

#include <vtkm/rendering/raytracing/Camera.h>

#ifdef VTKH_PARALLEL
#include <mpi.h>
#endif

namespace vtkh
{

PathTrace::PathTrace()
  : m_rank(0)
   ,m_procs(1)
#ifdef VTKH_PARALLEL
   ,m_messenger(MPI_Comm_f2c(vtkh::GetMPICommHandle()))
#endif
{
#ifdef VTKH_PARALLEL
  MPI_Comm comm = MPI_Comm_f2c(vtkh::GetMPICommHandle());
  MPI_Comm_rank(comm, &m_rank);
  MPI_Comm_size(comm, &m_procs);
  m_messenger.RegisterMessages(2, m_procs, m_procs);
#endif
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

  m_bounds_map.Build();
  m_spatial_query.SetBoundsMap(m_bounds_map);
  m_bounds = this->m_input->GetGlobalBounds();

  m_camera.ResetToBounds(m_bounds);
}

void PathTrace::PostExecute()
{
  Filter::PostExecute();
}

void PathTrace::PackOutgoing(vtkmRay &rays)
{
  SpatialQuery::Result dests = m_spatial_query.IntersectRays(rays);

  std::cout<<"Packing\n";
  std::vector<int> ray_dests(10); // destination for a single ray

  for(int i = 0; i < rays.NumRays; ++i)
  {
    dests.GetDomains(i, ray_dests);
    if(ray_dests.size() == 0) continue;
    Ray ray;

    ray.m_origin = rays.Origin.GetPortalConstControl().Get(i);
    ray.m_dir = rays.Dir.GetPortalConstControl().Get(i);
    const vtkm::Id c_offset = i * 4;
    ray.m_color[0] = rays.Buffers.at(0).Buffer.GetPortalConstControl().Get(c_offset + 0);
    ray.m_color[1] = rays.Buffers.at(0).Buffer.GetPortalConstControl().Get(c_offset + 1);
    ray.m_color[2] = rays.Buffers.at(0).Buffer.GetPortalConstControl().Get(c_offset + 2);

    ray.m_status = rays.Status.GetPortalConstControl().Get(i);
    //ray.m_depth = rays.;
    ray.m_distance = rays.Distance.GetPortalConstControl().Get(i);
    ray.m_max_distance = rays.MaxDistance.GetPortalConstControl().Get(i);
    ray.m_pixel_id = rays.PixelIdx.GetPortalConstControl().Get(i);

    //const vtkm::Id t_offset = i * 3;
    //ray.m_throughput[0] = rays.GetBuffer("throughput").Buffer.GetPortalControl().Get(t_offset + 0);
    //ray.m_throughput[1] = rays.GetBuffer("throughput").Buffer.GetPortalControl().Get(t_offset + 1);
    //ray.m_throughput[2] = rays.GetBuffer("throughput").Buffer.GetPortalControl().Get(t_offset + 2);
    for(int d = 0; d < ray_dests.size(); ++d)
    {
      const int domain = ray_dests[d];
      ray.m_dest_dom = domain;
      int dest_rank = m_bounds_map.m_rank_map[domain];
      //std::cout<<"Dest rank "<<dest_rank<<" domain "<<domain<<"\n";
      m_out_q[dest_rank].push_back(ray);
    }
  }
}

void PathTrace::RouteIncoming(std::vector<Ray> &in_rays)
{
  const int size = in_rays.size();
  for(int i = 0; i < size; ++i)
  {
    const int dom = in_rays[i].m_dest_dom;
    m_in_q[dom].push_back(in_rays[i]);
  }
}

void PathTrace::Recv()
{

#ifdef VTKH_PARALLEL
  std::vector<Ray> ray_data;
  for(int i = 0; i <1000000; ++i)
  {
    m_messenger.RecvRays(ray_data);
    if(ray_data.size() > 0 )
    {
      RouteIncoming(ray_data);
    }
  }
#endif
}

void PathTrace::RouteOutgoing()
{
  for(auto it = m_out_q.begin(); it != m_out_q.end(); ++it)
  {
    const int dest = it->first;

    if(m_rank == dest)
    {
      RouteIncoming(it->second);
      continue;
    }

#ifdef VTKH_PARALLEL
    std::cout<<"["<<m_rank<<"] --> ["<<dest<<"] "<<it->second.size()<<"\n";
    m_messenger.SendRays(dest, it->second);
#else
  // non-mpi should never get here
  std::cout<<"Non-mpi: can't route to rank "<<dest<<"\n";
#endif
  }
}

void PathTrace::DoExecute()
{
  int rank = 0;
  int procs = 1;

  vtkmRay rays;
  if(rank == 0)
  {
    CreateRays(rays);
  }

  SpatialQuery::Result dests = m_spatial_query.IntersectRays(rays);

  std::cout<<"Packing\n";
  PackOutgoing(rays);
  std::cout<<"Routing\n";
  RouteOutgoing();
  Recv();

  this->m_output = new DataSet();
  //const int num_domains = this->m_input->GetNumberOfDomains();

  //for(int i = 0; i < num_domains; ++i)
  //{
  //  vtkm::Id domain_id;
  //  vtkm::cont::DataSet dom;
  //  this->m_input->GetDomain(i, dom, domain_id);
  //  // insert interesting stuff
  //  m_output->AddDomain(dom, domain_id);
  //}
}

void PathTrace::CreateRays(vtkmRay &rays)
{
  vtkm::rendering::raytracing::Camera camera;
  int height = 100;
  int width = 100;
  vtkm::rendering::CanvasRayTracer canvas(width, height);
  camera.SetParameters(m_camera, canvas);
  camera.CreateRays(rays, m_bounds);
  std::cout<<"created rays\n";
}

std::string
PathTrace::GetName() const
{
  return "vtkh::PathTrace";
}

} //  namespace vtkh
