#include <vtkh/filters/CommTest.hpp>

#include "communication/BoundsMap.hpp"
#include "communication/Ray.hpp"
#ifdef VTKH_PARALLEL
#include "communication/avtParICAlgorithm.hpp"
#include "communication/RayMessenger.hpp"
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
  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);
    bmap.AddBlock(domain_id, dom.GetCoordinateSystem().GetBounds());
  }

  bmap.Build(); // does parallel comm and builds lookups

#ifdef VTKH_PARALLEL
  MPI_Comm comm = MPI_Comm_f2c(vtkh::GetMPICommHandle());
  //MPICommunicator communicator(comm);
  int rank;
  int procs;
  MPI_Comm_rank(comm, &rank);
  MPI_Comm_size(comm, &procs);

  int dest = (rank + 1) % procs;



  Particle p;
  p.coords[0] = rank;
  p.coords[1] = rank;
  p.coords[2] = rank;
  p.id = rank;
  p.nSteps = rank;
  p.blockId = rank;
  p.status = Particle::ACTIVE;

  std::list<Particle> plist;
  plist.push_back(p);

  avtParICAlgorithm pcomm(comm);
  pcomm.RegisterMessages(2, procs, procs);
  pcomm.SendICs(dest, plist);

  bool expect_message = true;
  list<ParticleCommData<Particle>> particleData;
  vector<MsgCommData> msgData;
  bool blockAndWait = false;
  DBG(" ---RecvAny..."<<endl);
  //pcomm.RecvAny(&msgData, &particleData, NULL, blockAndWait);
  while(expect_message)
  {
    pcomm.RecvICs(particleData);
    if(particleData.size() != 0) expect_message = false;
  }

  std::cout<<"done\n";
  Ray ray;
  ray.m_origin[0] = rank;
  ray.m_origin[1] = rank;
  ray.m_origin[2] = rank;
  std::vector<Ray> rays;
  rays.push_back(ray);

  RayMessenger rmessenger(comm);
  rmessenger.RegisterMessages(2, procs, procs);
  rmessenger.SendRays(dest, rays);
  std::vector<ParticleCommData<Ray>> ray_data;;
  expect_message = true;

  while(expect_message)
  {
    rmessenger.RecvRays(ray_data);
    if(ray_data.size() != 0) expect_message = false;
  }

#endif

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
