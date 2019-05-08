//#include <vtkm/filter/your_vtkm_filter.h>
#include <vtkh/filters/PathTrace.hpp>

#include <vtkm/rendering/raytracing/Camera.h>

#include <vtkm/cont/ArrayCopy.h>
#include <vtkm/cont/ColorTable.hxx>
#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/worklet/WorkletMapField.h>


#ifdef VTKH_PARALLEL
#include <mpi.h>
#endif



#include <vtkh/filters/communication/DebugMeowMeow.hpp>

#ifdef  DEBUG_LOG
#define DLOG(x) m_log<<x
#else
#define DLOG(x)
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
  std::stringstream logname;
  logname<<"log_"<<m_rank<<".log";
  m_log.open(logname.str().c_str(), std::ofstream::out);
}

PathTrace::~PathTrace()
{
  m_log.close();
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

  vtkm::cont::ArrayHandle<vtkm::Range> ranges = m_input->GetGlobalRange(m_field_name);
  int num_components = ranges.GetPortalControl().GetNumberOfValues();

  //
  // current vtkm renderers only supports single component scalar fields
  //
  assert(num_components == 1);
  vtkm::Range global_range = ranges.GetPortalControl().Get(0);

  // Build the bounds map so we can ask questions
  m_bounds_map.Clear();
  const int num_domains = this->m_input->GetNumberOfDomains();
  m_trackers.resize(num_domains);

  // color table
  vtkm::cont::ColorTable ct("cool to warm");
  constexpr vtkm::Float32 conversionToFloatSpace = (1.0f / 255.0f);
  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::UInt8, 4>> temp;
  ct.Sample(1024, temp);

  vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32,4>> color_map;
  color_map.Allocate(1024);
  auto portal = color_map.GetPortalControl();
  auto colorPortal = temp.GetPortalConstControl();
  for (vtkm::Id i = 0; i < 1024; ++i)
  {
    auto color = colorPortal.Get(i);
    vtkm::Vec<vtkm::Float32, 4> t(color[0] * conversionToFloatSpace,
                                  color[1] * conversionToFloatSpace,
                                  color[2] * conversionToFloatSpace,
                                  color[3] * conversionToFloatSpace);
    portal.Set(i, t);
  }


  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);
    vtkm::Bounds bounds = dom.GetCoordinateSystem().GetBounds();
    m_bounds_map.AddBlock(domain_id, bounds);

    auto cellset = dom.GetCellSet();
    auto coords = dom.GetCoordinateSystem();
    vtkm::cont::Field field = dom.GetField(m_field_name);

    m_trackers[i].SetData(coords, cellset.Cast<vtkm::cont::CellSetStructured<3>>());
    m_trackers[i].SetScalarField(field, global_range);
    m_trackers[i].SetColorMap(color_map);
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

struct CompletedFunctor
{
  VTKM_EXEC_CONT
  template< typename T>
  bool operator()(const T& dests)
  {
    bool is_complete = false;
    if(dests == 0)
    {
      is_complete = true;
    }
    return is_complete;
  }

};

void PathTrace::PackRays(vtkmRay &rays)
{
  SpatialQuery::Result dests = m_spatial_query.IntersectRays(rays);

  vtkmRay completed_rays = rays.CopyIf(dests.m_counts,CompletedFunctor());

  std::vector<int> ray_dests(10); // destination for a single ray
  int scount = 0;
  int tcount = 0;
  for(int i = 0; i < rays.NumRays; ++i)
  {
    dests.GetDomains(i, ray_dests);

    if(ray_dests.size() > 1) std::cout<<"spec ";
    if(ray_dests.size() == 0)
    {
      tcount++;
      continue;
    }
    scount++;

    Ray ray;

    ray.m_origin = rays.Origin.GetPortalConstControl().Get(i);
    ray.m_dir = rays.Dir.GetPortalConstControl().Get(i);
    const vtkm::Id c_offset = i * 4;
    ray.m_color[0] = rays.Buffers.at(0).Buffer.GetPortalConstControl().Get(c_offset + 0);
    ray.m_color[1] = rays.Buffers.at(0).Buffer.GetPortalConstControl().Get(c_offset + 1);
    ray.m_color[2] = rays.Buffers.at(0).Buffer.GetPortalConstControl().Get(c_offset + 2);
    ray.m_color[3] = 1.f;

    ray.m_status = rays.Status.GetPortalConstControl().Get(i);
    //ray.m_depth = rays.;
    ray.m_distance = rays.Distance.GetPortalConstControl().Get(i);
    ray.m_max_distance = rays.MaxDistance.GetPortalConstControl().Get(i);
    ray.m_pixel_id = rays.PixelIdx.GetPortalConstControl().Get(i);

    const vtkm::Id t_offset = i * 3;
    ray.m_throughput[0] = rays.GetBuffer("throughput").Buffer.GetPortalControl().Get(t_offset + 0);
    ray.m_throughput[1] = rays.GetBuffer("throughput").Buffer.GetPortalControl().Get(t_offset + 1);
    ray.m_throughput[2] = rays.GetBuffer("throughput").Buffer.GetPortalControl().Get(t_offset + 2);

    for(int d = 0; d < ray_dests.size(); ++d)
    {
      const int domain = ray_dests[d];
      ray.m_dest_dom = domain;
      int dest_rank = m_bounds_map.m_rank_map[domain];
      if(dest_rank < 0 || dest_rank >= m_procs) std::cout<<"Error: bad rank "<<dest_rank<<"\n";
      m_out_q[dest_rank].push_back(ray);
      DLOG("sending "<<ray.m_pixel_id<<" "<<dest_rank<<"\n");
    }
  }

  std::cout<<"Total "<<rays.NumRays<<" Completed rays "<<completed_rays.NumRays<<" active "<<scount<<"\n";
  std::cout<<"tcout "<<tcount<<"\n";
  if(completed_rays.NumRays > 0)
  {
    vtkm::Vec<vtkm::Float32, 3> bg_color(1.f, 1.f, 1.f);
    assert(m_trackers.size() > 0);
    // All trackers have the same lights
    m_trackers[0].AddLightContribution(completed_rays, bg_color);
  }

  for(int i = 0; i < completed_rays.NumRays; ++i)
  {
    RayResult res;
    //std::cout<<" "<<ray_dests.size();
    //std::cout<<" "<<(int)completed_rays.Status.GetPortalConstControl().Get(i);

    const vtkm::Id c_offset = i * 4;
    res.m_color[0] = completed_rays.Buffers.at(0).Buffer.GetPortalConstControl().Get(c_offset + 0);
    res.m_color[1] = completed_rays.Buffers.at(0).Buffer.GetPortalConstControl().Get(c_offset + 1);
    res.m_color[2] = completed_rays.Buffers.at(0).Buffer.GetPortalConstControl().Get(c_offset + 2);

    res.m_pixel_id = completed_rays.PixelIdx.GetPortalConstControl().Get(i);

    m_result_q.push_back(res);
    DLOG("sending "<<res.m_pixel_id<<" "<<"-1"<<"\n");
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

void PathTrace::TraceRays()
{
  const int num_domains = this->m_input->GetNumberOfDomains();

  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);

    vtkmRay &rays = m_dom_rays[domain_id];
    if(rays.NumRays > 0)
    {
      vtkm::rendering::raytracing::ResidualTracker::RandomState rngs;
      rngs.Allocate(rays.NumRays);
      vtkm::rendering::raytracing::seedRng(rngs, true);

      m_trackers[i].Trace(rays, rngs);
    }

  }
  //ForwardRays();
}

void PathTrace::RouteToDomains(std::vector<Ray> &in_rays)
{
  const int size = in_rays.size();
  for(int i = 0; i < size; ++i)
  {
    const int dom = in_rays[i].m_dest_dom;

    DLOG("recieved "<<in_rays[i].m_pixel_id<<"\n");

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
    //// check signals
    //std::vector<MsgCommData> msgs;
    //m_messenger.RecvMsg(msgs);
    //if(msgs.size() != 0)
    //{
    //  for(int i = 0; i < msgs.size(); ++i)
    //  {
    //    if(msgs[i].message[0] == MessageType::DEATH)
    //    {
    //      const int count = msgs[i].message[1];
    //      m_total_rays -= count;
    //      DPRINT("[0] Death message "<<count<<" remaining "<<m_total_rays<<"\n");
    //    }
    //  }
    //}
    // check signals
    std::vector<RayResult> results;
    m_messenger.RecvResults(results);
    if(results.size() != 0)
    {
      m_total_rays -= results.size();
      DPRINT("[0] Death message "<<results.size()<<" remaining "<<m_total_rays<<"\n");
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
    // route to self
    if(m_rank == dest)
    {
      RouteToDomains(it->second);
    }
#ifdef VTKH_PARALLEL
    else if(it->second.size() > 0)
    {
      DPRINT("["<<m_rank<<"] --> ["<<dest<<"] "<<it->second.size()<<"\n");
      m_messenger.SendRays(dest, it->second);
    }
#else
  // non-mpi should never get here
    std::cout<<"Non-mpi: can't route to rank "<<dest<<"\n";
#endif
    it->second.clear();
  }

  // send death notices
  const int death_count = m_result_q.size();
  if(death_count != 0)
  {
    out += death_count;
    for(int i =0; i < death_count; ++i)
    {
      DLOG("death "<<m_result_q[i].m_pixel_id<<"\n");
    }
    if(m_rank != 0)
    {
#ifdef VTKH_PARALLEL
      m_messenger.SendResults(0, m_result_q);
      std::cout<<"SEND death\n";
#endif
   }
    else
    {
      m_total_rays -= death_count;
      //std::cout<<"["<<m_rank<<"]  remaining "<< m_total_rays<<"\n";
    }

    DPRINT("["<<m_rank<<"] death "<<death_count<<"\n");

    total_death += death_count;
    m_result_q.clear();
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
    rays.AddBuffer(3, "throughput");

    if(num_rays > 0)
    {
      DPRINT("["<<m_rank<<"] unpacking domain "<<domain_id<<" size "<<num_rays<<"\n");
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
      rays.Buffers.at(0).Buffer.GetPortalControl().Set(c_offset + 3, 1.f);

      const vtkm::Id t_offset = i * 3;
      rays.GetBuffer("throughput").Buffer.GetPortalControl().Set(t_offset + 0, ray.m_throughput[0]);
      rays.GetBuffer("throughput").Buffer.GetPortalControl().Set(t_offset + 1, ray.m_throughput[1]);
      rays.GetBuffer("throughput").Buffer.GetPortalControl().Set(t_offset + 2, ray.m_throughput[2]);

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

  rays.AddBuffer(3, "throughput");
  rays.GetBuffer("throughput").InitConst(1.f);
  rays.Buffers.at(0).InitConst(0.f);

  m_total_rays = rays.NumRays;
  int t_ray = m_total_rays;

  PackRays(rays);
  DPRINT("init rays "<<rays.NumRays<<"\n");
  if(m_rank == 0)
  {
    std::ofstream outrays("init_rays", std::ofstream::out);
    for(auto it = m_out_q.begin(); it != m_out_q.end(); ++it)
    {
      const int dest = it->first;
      for(int i = 0; i < it->second.size(); ++i)
      {
        outrays<<it->second[i].m_pixel_id<<"\n";
      }
    }
    outrays.close();
  }
  Send();

  DPRINT("BEGIN\n");
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

    if(ActiveRays() > 0)
    {
      DPRINT("["<<m_rank<<"] tracing "<<ActiveRays()<<"\n");
    }

    //ForwardRays(); // stand in for trace
    TraceRays();

    PackOutgoing();
    Send();

    if(m_rank == 0 && m_total_rays < 1)
    {
      Kill();
    }

    //if(m_rank == 0)
    //{
    //  static int stagnation = 0;
    //  static int count = 0;
    //  if(m_total_rays == count)
    //  {
    //    stagnation++;
    //  }
    //  count = m_total_rays;
    //  if(stagnation > 10000)
    //  {
    //    std::cout<<"Stagnation!!!\n";
    //    Kill();
    //  }
    //  //std::cout<<" ***** outstanding "<<m_total_rays<<"\n";
    //}
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
  int height = 1024;
  int width = 1024;
  //int height = 512;
  //int width = 512;
  vtkm::rendering::CanvasRayTracer canvas(width, height);
  camera.SetParameters(m_camera, canvas);
  camera.CreateRays(rays, m_bounds);

}

std::string
PathTrace::GetName() const
{
  return "vtkh::PathTrace";
}

} //  namespace vtkh
