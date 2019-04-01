//#include <vtkm/filter/your_vtkm_filter.h>
#include <vtkh/filters/PathTrace.hpp>

#include <vtkm/rendering/raytracing/Camera.h>

#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>

#ifdef VTKH_PARALLEL
#include <mpi.h>
#endif

namespace vtkh
{

class CalcRayStart : public vtkm::worklet::WorkletMapField
{
  vtkm::Float32 Xmin;
  vtkm::Float32 Ymin;
  vtkm::Float32 Zmin;
  vtkm::Float32 Xmax;
  vtkm::Float32 Ymax;
  vtkm::Float32 Zmax;

public:
  VTKM_CONT
  CalcRayStart(const vtkm::Bounds boundingBox)
  {
    Xmin = static_cast<vtkm::Float32>(boundingBox.X.Min);
    Xmax = static_cast<vtkm::Float32>(boundingBox.X.Max);
    Ymin = static_cast<vtkm::Float32>(boundingBox.Y.Min);
    Ymax = static_cast<vtkm::Float32>(boundingBox.Y.Max);
    Zmin = static_cast<vtkm::Float32>(boundingBox.Z.Min);
    Zmax = static_cast<vtkm::Float32>(boundingBox.Z.Max);
  }

  VTKM_EXEC
  vtkm::Float32 rcp(vtkm::Float32 f) const { return 1.0f / f; }

  VTKM_EXEC
  vtkm::Float32 rcp_safe(vtkm::Float32 f) const { return rcp((fabs(f) < 1e-8f) ? 1e-8f : f); }

  using ControlSignature = void(FieldIn, FieldOut, FieldInOut, FieldInOut, FieldIn);
  using ExecutionSignature = void(_1, _2, _3, _4, _5);
  template <typename Precision>
  VTKM_EXEC void operator()(const vtkm::Vec<Precision, 3>& rayDir,
                            vtkm::Float32& minDistance,
                            vtkm::Float32& distance,
                            vtkm::Float32& maxDistance,
                            const vtkm::Vec<Precision, 3>& rayOrigin) const
  {
    vtkm::Float32 dirx = static_cast<vtkm::Float32>(rayDir[0]);
    vtkm::Float32 diry = static_cast<vtkm::Float32>(rayDir[1]);
    vtkm::Float32 dirz = static_cast<vtkm::Float32>(rayDir[2]);
    vtkm::Float32 origx = static_cast<vtkm::Float32>(rayOrigin[0]);
    vtkm::Float32 origy = static_cast<vtkm::Float32>(rayOrigin[1]);
    vtkm::Float32 origz = static_cast<vtkm::Float32>(rayOrigin[2]);

    vtkm::Float32 invDirx = rcp_safe(dirx);
    vtkm::Float32 invDiry = rcp_safe(diry);
    vtkm::Float32 invDirz = rcp_safe(dirz);

    vtkm::Float32 odirx = origx * invDirx;
    vtkm::Float32 odiry = origy * invDiry;
    vtkm::Float32 odirz = origz * invDirz;

    vtkm::Float32 xmin = Xmin * invDirx - odirx;
    vtkm::Float32 ymin = Ymin * invDiry - odiry;
    vtkm::Float32 zmin = Zmin * invDirz - odirz;
    vtkm::Float32 xmax = Xmax * invDirx - odirx;
    vtkm::Float32 ymax = Ymax * invDiry - odiry;
    vtkm::Float32 zmax = Zmax * invDirz - odirz;


    minDistance = vtkm::Max(
      vtkm::Max(vtkm::Max(vtkm::Min(ymin, ymax), vtkm::Min(xmin, xmax)), vtkm::Min(zmin, zmax)),
      minDistance);
    vtkm::Float32 exitDistance =
      vtkm::Min(vtkm::Min(vtkm::Max(ymin, ymax), vtkm::Max(xmin, xmax)), vtkm::Max(zmin, zmax));
    //maxDistance = vtkm::Min(maxDistance, exitDistance);
    if (exitDistance < minDistance)
    {
      minDistance = -1.f; //flag for miss
    }
    else
    {
      // set the ray to leave the domain
      minDistance = exitDistance;
      distance = exitDistance;
    }
  }
}; //class CalcRayStart
//
// stand in code for tracing
//
void PathTrace::ForwardRays()
{

  const int num_domains = this->m_input->GetNumberOfDomains();

  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);
    vtkm::Bounds bounds = dom.GetCoordinateSystem().GetBounds();

    vtkmRay &rays = m_dom_rays[domain_id];

    CalcRayStart ttt(bounds);
    vtkm::worklet::DispatcherMapField<CalcRayStart> calcDispatcher(ttt);

    calcDispatcher.Invoke(rays.Dir, rays.MinDistance, rays.Distance, rays.MaxDistance, rays.Origin);

    //vtkm::cont::Algorithm::Copy(rays.Distance, rays.MinDistance);

  }
}

PathTrace::PathTrace()
  : m_rank(0)
   ,m_procs(1)
   ,m_killed(false)
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

static int out = 0;
static int total_death = 0;
static int in = 0;

void PathTrace::PackRays(vtkmRay &rays)
{
  SpatialQuery::Result dests = m_spatial_query.IntersectRays(rays);

  std::vector<int> ray_dests(10); // destination for a single ray
  for(int i = 0; i < rays.NumRays; ++i)
  {
    dests.GetDomains(i, ray_dests);

    if(ray_dests.size() > 1) std::cout<<"spec ";

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
    if(ray_dests.size() == 0)
    {
      m_out_q[-1].push_back(ray);
    }

    for(int d = 0; d < ray_dests.size(); ++d)
    {
      const int domain = ray_dests[d];
      ray.m_dest_dom = domain;
      int dest_rank = m_bounds_map.m_rank_map[domain];
      //std::cout<<"Dest rank "<<dest_rank<<" domain "<<domain<<"\n";
      if(dest_rank < 0 || dest_rank >= m_procs) std::cout<<"Error: bad rank "<<dest_rank<<"\n";
      m_out_q[dest_rank].push_back(ray);
    }
  }


}

void PathTrace::PackOutgoing()
{
  for(auto it = m_dom_rays.begin(); it != m_dom_rays.end(); ++it)
  {
    vtkmRay &rays = it->second;
    PackRays(rays);
    rays.Resize(0);
  }
}

void PathTrace::RouteToDomains(std::vector<Ray> &in_rays)
{
  const int size = in_rays.size();
  for(int i = 0; i < size; ++i)
  {
    const int dom = in_rays[i].m_dest_dom;
    m_in_q[dom].push_back(in_rays[i]);
  }
}

int PathTrace::ActiveRays()
{

  int sum = 0;
  for(auto it = m_dom_rays.begin(); it != m_dom_rays.end(); ++it)
  {
    vtkmRay &rays = it->second;
    sum += rays.NumRays;
  }
  return sum;
}

void PathTrace::Recv()
{

#ifdef VTKH_PARALLEL
  std::vector<Ray> ray_data;
  m_messenger.RecvRays(ray_data);
  if(ray_data.size() > 0 )
  {
    RouteToDomains(ray_data);
  }

  if(m_rank == 0)
  {
    // check death signals
    std::vector<MsgCommData> msgs;
    m_messenger.RecvMsg(msgs);
    if(msgs.size() != 0)
    {
      for(int i = 0; i < msgs.size(); ++i)
      {
        if(msgs[i].message[0] == MessageType::DEATH)
        {
          const int count = msgs[i].message[1];
          m_total_rays -= count;
          std::cout<<"[0] Death message "<<count<<" remaining "<<m_total_rays<<"\n";
        }
      }
    }
  }
  else
  {
    // all other ranks check for control signals
    std::vector<MsgCommData> msgs;
    m_messenger.RecvMsg(msgs);
    if(msgs.size() != 0)
    {
      for(int i = 0; i < msgs.size(); ++i)
      {
        if(msgs[i].message[0] == MessageType::KILL)
        {
          m_killed = true;
          std::cout<<"["<<m_rank<<"] killed\n";
        }
      }
    }
  }
#endif
}

void PathTrace::Send()
{
  for(auto it = m_out_q.begin(); it != m_out_q.end(); ++it)
  {
    const int dest = it->first;

    out += it->second.size();
    // send death notices
    if(dest == -1)
    {
      const int death_count = it->second.size();
      if(death_count != 0)
      {
        if(m_rank != 0)
        {
#ifdef VTKH_PARALLEL
          std::vector<int> msg(2);
          msg[0] = MessageType::DEATH;
          msg[1] = death_count;
          m_messenger.SendMsg(0, msg);
#endif
        }
        else
        {
          m_total_rays -= death_count;
          std::cout<<"["<<m_rank<<"]  remaining "<< m_total_rays<<"\n";
        }

        std::cout<<"["<<m_rank<<"] death "<<death_count<<"\n";

        total_death += death_count;
      }

    }
    // route to self
    if(m_rank == dest)
    {
      RouteToDomains(it->second);
    }
#ifdef VTKH_PARALLEL
    else if(it->second.size() > 0)
    {
      std::cout<<"["<<m_rank<<"] --> ["<<dest<<"] "<<it->second.size()<<"\n";
      m_messenger.SendRays(dest, it->second);
    }
#else
  // non-mpi should never get here
    std::cout<<"Non-mpi: can't route to rank "<<dest<<"\n";
#endif
    it->second.clear();
  }
}

void PathTrace::UnpackIncoming()
{
  for(auto it = m_in_q.begin(); it != m_in_q.end(); ++it)
  {
    const vtkm::Id domain_id = it->first;
    if(!this->m_input->HasDomainId(domain_id))
    {
      std::cout<<"["<<m_rank<<"] ERROR: does not have domain id "<<domain_id<<"\n";
    }

    const int num_rays = it->second.size();
    vtkmRay &rays = m_dom_rays[domain_id];
    rays.Resize(num_rays);
    // TODO: add throughput buffer
    // or setup buffers at the begging and let
    //i resize take care of it
    if(num_rays > 0)
    {
      std::cout<<"["<<m_rank<<"] unpacking domain "<<domain_id<<" size "<<num_rays<<"\n";
    }
    in += num_rays;
    for(int i = 0; i < num_rays; ++i)
    {
      Ray &ray = it->second[i];
      rays.Origin.GetPortalControl().Set(i, ray.m_origin);
      rays.Dir.GetPortalControl().Set(i, ray.m_dir);

      const vtkm::Id c_offset = i * 4;
      rays.Buffers.at(0).Buffer.GetPortalControl().Set(c_offset + 0, ray.m_color[0]);
      rays.Buffers.at(0).Buffer.GetPortalControl().Set(c_offset + 1, ray.m_color[1]);
      rays.Buffers.at(0).Buffer.GetPortalControl().Set(c_offset + 2, ray.m_color[2]);

      rays.Status.GetPortalControl().Set(i, ray.m_status);
      rays.Distance.GetPortalControl().Set(i, ray.m_distance);
      rays.MaxDistance.GetPortalControl().Set(i, ray.m_max_distance);
      rays.PixelIdx.GetPortalControl().Set(i, ray.m_pixel_id);
      ////ray.m_depth = rays.;
    }
    // clear the incoming q
    it->second.clear();
  }
}

void PathTrace::Kill()
{
#ifdef VTKH_PARALLEL
  std::vector<int> msg(1);
  msg[0] = MessageType::KILL;
  m_messenger.SendAllMsg(msg);
  m_killed = true;
#endif
}

void PathTrace::DoExecute()
{
  vtkmRay rays;
  if(m_rank == 0)
  {
    CreateRays(rays);
  }

  m_total_rays = rays.NumRays;
  int t_ray = m_total_rays;

  PackRays(rays);
  std::cout<<"init rays "<<rays.NumRays<<"\n";
  Send();

  in = 0;
  out = 0;

  // TODO: filter spatial query for termination status
  // if a ray count is 0 it has exited the distrubuted mesh.
  // Add the light / background contribution and send the dead
  // ray home


  while(!m_killed)
  {

    Recv();
    UnpackIncoming();

    if(ActiveRays() > 0) std::cout<<"["<<m_rank<<"] tracing "<<ActiveRays()<<"\n";
    ForwardRays(); // stand in for trace

    PackOutgoing();
    Send();

    if(m_rank == 0 && m_total_rays < 1)
    {
      Kill();
    }

    if(m_rank == 0)
    {
      static int stagnation = 0;
      static int count = 0;
      if(m_total_rays == count)
      {
        stagnation++;
      }
      count = m_total_rays;
      if(stagnation > 10000) Kill();
      //std::cout<<" ***** outstanding "<<m_total_rays<<"\n";
    }
  }

  if(in != out)
  {
    std::cout<<"["<<m_rank<<"] missing "<<in<<" != "<<out<<"\n";
  }
  if(ActiveRays() > 0)
  {
    std::cout<<"["<<m_rank<<"] still active "<<ActiveRays()<<"\n";
  }

  if(m_rank == 0 && m_total_rays > 0 )
  {
    std::cout<<" ***** outstanding "<<m_total_rays<<"\n";
  }
#ifdef VTKH_PARALLEL
  MPI_Comm comm = MPI_Comm_f2c(vtkh::GetMPICommHandle());
  int res = 0;
  MPI_Reduce(&total_death, &res, 1, MPI_INT, MPI_SUM, 0, comm);
  if(m_rank == 0)
  {
    std::cout<<"TOTAL reported death "<<res<<" rays "<<t_ray<<"\n";
  }
#endif


  // forward rays


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
