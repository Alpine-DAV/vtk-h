#include <iostream>
#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkh/filters/ParticleAdvection.hpp>
#include <vtkh/filters/Communicator.hpp>
#include <vtkh/filters/util.hpp>
#include <vtkh/vtkh.hpp>
#include <vtkh/Error.hpp>

#include <vtkm/filter/Streamline.h>

#ifdef VTKH_PARALLEL
#include <mpi.h>
#endif


#include <vtkh/filters/DebugMeowMeow.hpp>
ofstream dbg;

using namespace std;

namespace vtkh
{
static inline float
rand01()
{
  return (float)rand() / (RAND_MAX+1.0f);
}

static inline float
randRange(const float &a, const float &b)
{
    return a + (b-a)*rand01();
}

ParticleAdvection::ParticleAdvection()
    : rank(0), numRanks(1), seedMethod(RANDOM),
      numSeeds(1000), totalNumSeeds(-1), randSeed(314),
      stepSize(.01),
      maxSteps(1000)
{
#ifdef VTKH_PARALLEL
  mpiComm = MPI_Comm_f2c(vtkh::GetMPICommHandle());
  MPI_Comm_rank(mpiComm, &rank);
  MPI_Comm_size(mpiComm, &numRanks);
#endif

#ifdef TRACE_DEBUG
  dbg = ofstream();
  char nm[32];
  sprintf(nm, "%d.out", rank);
  dbg.open(nm, ofstream::out);
#endif
}

ParticleAdvection::~ParticleAdvection()
{
}

void ParticleAdvection::PreExecute()
{
  Filter::PreExecute();
}

void ParticleAdvection::PostExecute()
{
  Filter::PostExecute();
}

void ParticleAdvection::DoExecute()
{
  cout<<__FILE__<<"  ParticleAdvection::DoExecute()"<<endl;

  //steps:
  //- compute bounds
  //- load data
  //- create seeds, and distribute
  //- iterate
  //- finish

  this->Init();
  this->CreateSeeds();
  cout<<rank<<" dom2rank: "<<domToRank<<endl;
  cout<<rank<<": "<<totalNumSeeds<<" "<<active.size()<<endl;

  bool done = false;
  int N = totalNumSeeds;
#ifdef VTKH_PARALLEL
  MPICommunicator communicator(mpiComm, domToRank);
  int cnt = 0;

  while (N > 0)
  {
      DBG("**** Loop: N= "<<N<<" AIT: "<<active.size()<<" "<<inactive.size()<<" "<<terminated.size()<<endl);

      vector<Particle> v;
      list<Particle> I, T, A;
      int numTerm = 0;
      if (GetActiveParticles(v))
      {
          DataBlock *blk = GetBlock(v[0].blockId);

          DBG("Integrate: "<<v.size()<<endl);
          blk->integrator.Go(false, v, maxSteps, I, T, A);
          DBG("  -----Integrate:  ITA: "<<I.size()<<" "<<T.size()<<" "<<A.size()<<endl);
          DBG("                   I= "<<I<<endl);
          numTerm = T.size();
          if (numTerm > 0)
              terminated.insert(terminated.end(), T.begin(), T.end());
      }

      bool haveActive = !active.empty();
      list<Particle> in, term;
      int totalTerm = communicator.Exchange(haveActive, I, in, term, boundsMap, numTerm);
      N -= (numTerm+totalTerm);

      DBG(" Exchange done: in="<<in<<" numTerm: "<<numTerm<<" totalTerm: "<<totalTerm<<endl);
      if (!term.empty())
          terminated.insert(terminated.end(), term.begin(), term.end());
      if (!in.empty())
          active.insert(active.end(), in.begin(), in.end());
      cnt++;
//      if (cnt > 10)
//          break;

      // no work available, take a snooze.
      if (active.empty())
      {
          DBG("    sleep"<<endl);
          usleep(1000);
      }
      DBG("   N --> "<<N<<endl);
      DBG(endl<<endl<<endl);
  }
  DBG("All done"<<endl);
  MPI_Barrier(mpiComm);
  exit(-1);

#endif



  this->m_output = new DataSet();
  const int nDoms = this->m_input->GetNumberOfDomains();

  std::vector<vtkm::cont::DataSet> datasets;
  std::vector<vtkm::Id> blockIds;
  for (int i = 0; i < nDoms; i++)
  {
    vtkm::Id domId;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domId);
    if(dom.HasField(m_field_name))
    {
      using vectorField_d = vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float64, 3>>;
      using vectorField_f = vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32, 3>>;
      auto field = dom.GetField(m_field_name).GetData();
      if(!field.IsSameType(vectorField_d()) && !field.IsSameType(vectorField_f()))
      {
        throw Error("Vector field type does not match <vtkm::Vec<vtkm::Float32,3>> or <vtkm::Vec<vtkm::Float64,3>>");
      }
    }
    else
    {
      throw Error("Domain does not contain specified vector field for ParticleAdvection analysis.");
    }
    datasets.push_back(dom);
    blockIds.push_back(domId);

    vector<vtkm::Vec<vtkm::FloatDefault, 3>> pts;
    vtkm::Bounds bounds = dom.GetCoordinateSystem().GetBounds();
    int numPts = 100;
    float dx = bounds.X.Max - bounds.X.Min;
    float dy = bounds.Y.Max - bounds.Y.Min;
    float dz = bounds.Z.Max - bounds.Z.Min;
    pts.resize(numPts);
    for (int i = 0; i < numPts; i++)
    {
        float x = bounds.X.Min + rand01()*dx;
        float y = bounds.Y.Min + rand01()*dy;
        float z = bounds.Z.Min + rand01()*dz;
        pts[i] = vtkm::make_Vec(x,y,z);
    }

    vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault, 3>> seeds;
    seeds = vtkm::cont::make_ArrayHandle(pts);

    vtkm::filter::Streamline streamline;
    streamline.SetStepSize(0.01);
    streamline.SetNumberOfSteps(100);
    streamline.SetSeeds(seeds);

    streamline.SetActiveField(m_field_name);
    vtkm::cont::DataSet res = streamline.Execute(dom);
    res.PrintSummary(cout);

    m_output->AddDomain(res, domId);
  }

  cout<<__FILE__<<"  ParticleAdvection::DoExecute() DONE"<<endl;
}


bool
ParticleAdvection::GetActiveParticles(std::vector<Particle> &v)
{
    v.clear();
    if (active.empty())
        return false;

    int blockId = active.front().blockId;
    while (!active.empty())
    {
        Particle p = active.front();
        v.push_back(p);
        active.pop_front();
    }

    return !v.empty();
}

std::string
ParticleAdvection::GetName() const
{
  return "vtkh::ParticleAdvection";
}


void
ParticleAdvection::Init()
{
  int totalNumDoms = this->m_input->GetGlobalNumberOfDomains();
  int nDoms = this->m_input->GetNumberOfDomains();

  vector<double> domBounds(6*totalNumDoms, 0.0), locBounds(6*totalNumDoms, 0.0);
  cout<<"TOTAL NUM DOMAINS= "<<totalNumDoms<<" nDoms= "<<nDoms<<endl;

  for (int i = 0; i < nDoms; i++)
  {
    vtkm::Id domId;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domId);
    dom.PrintSummary(cout);
    dataBlocks.push_back(new DataBlock(domId, dom, m_field_name, stepSize));

    vtkm::Bounds b = dom.GetCoordinateSystem().GetBounds();
    cout<<i<<" domId= "<<domId<<endl;
    cout<<i<<" bounds: ("<<b.X.Min<<" "<<b.Y.Min<<" "<<b.Z.Min<<") ("<<b.X.Max<<" "<<b.Y.Max<<" "<<b.Z.Max<<endl;

    locBounds[domId*6 + 0] = b.X.Min;
    locBounds[domId*6 + 1] = b.X.Max;
    locBounds[domId*6 + 2] = b.Y.Min;
    locBounds[domId*6 + 3] = b.Y.Max;
    locBounds[domId*6 + 4] = b.Z.Min;
    locBounds[domId*6 + 5] = b.Z.Max;
  }

#ifdef VTKH_PARALLEL
  MPI_Allreduce(&locBounds[0], &domBounds[0], 6*totalNumDoms, MPI_DOUBLE, MPI_SUM, mpiComm);
#endif

  for (int i = 0; i < totalNumDoms; i++)
  {
    vtkm::Bounds b(domBounds[i*6 + 0],
                   domBounds[i*6 + 1],
                   domBounds[i*6 + 2],
                   domBounds[i*6 + 3],
                   domBounds[i*6 + 4],
                   domBounds[i*6 + 5]);
    boundsMap.AddBlock(i, b);
    if (i == 0)
        globalBounds = b;
    else
        globalBounds.Include(b);
  }

  domToRank.resize(totalNumDoms, 0);
#ifdef VTKH_PARALLEL
  vector<int> tmp(totalNumDoms, 0);
  for (int i = 0; i < nDoms; i++)
      tmp[this->m_input->GetDomainIds()[i]] = rank;
  MPI_Allreduce(&tmp[0], &domToRank[0], totalNumDoms, MPI_INT, MPI_SUM, mpiComm);
#endif
}
DataBlock *
ParticleAdvection::GetBlock(int blockId)
{
    for (auto &d : dataBlocks)
        if (d->id == blockId)
            return d;
    return NULL;
}

void
ParticleAdvection::BoxOfSeeds(const vtkm::Bounds &box,
                              std::vector<Particle> &seeds,
                              vtkm::Id domId,
                              bool shrink)
{
  float boxRange[6] = {(float)box.X.Min, (float)box.X.Max,
                       (float)box.Y.Min, (float)box.Y.Max,
                       (float)box.Z.Min, (float)box.Z.Max};
  //shrink by 5%
  if (shrink)
  {
    float d[3] = {(boxRange[1]-boxRange[0]) * 0.025f,
                  (boxRange[3]-boxRange[2]) * 0.025f,
                  (boxRange[5]-boxRange[4]) * 0.025f};
    boxRange[0] += d[0];
    boxRange[1] -= d[0];
    boxRange[2] += d[1];
    boxRange[3] -= d[1];
    boxRange[4] += d[2];
    boxRange[5] -= d[2];
  }

  srand(randSeed);
  //TODO: Take care of idselect
  int N = numSeeds;
  for (int i = 0; i < N; i++)
  {
      int id = i;
      if (seedMethod == RANDOM_BLOCK)
          id = rank*N+i;

      Particle p;
      p.blockId = domId;
      p.id = id;
      p.coords[0] = randRange(boxRange[0], boxRange[1]);
      p.coords[1] = randRange(boxRange[2], boxRange[3]);
      p.coords[2] = randRange(boxRange[4], boxRange[5]);
      seeds.push_back(p);
  }
  cout<<rank<<" boxof seeds: "<<seeds.size()<<" "<<box<<" "<<seedMethod<<endl;
}

void
ParticleAdvection::CreateSeeds()
{
    active.clear();
    inactive.clear();
    terminated.clear();

    std::vector<Particle> seeds;
    if (seedMethod == RANDOM_BLOCK)
    {
        const int nDoms = this->m_input->GetNumberOfDomains();
        for (int i = 0; i < nDoms; i++)
        {
            vtkm::Id domId;
            vtkm::cont::DataSet dom;
            this->m_input->GetDomain(i, dom, domId);
            vtkm::Bounds b = dom.GetCoordinateSystem().GetBounds();
            BoxOfSeeds(b, seeds, domId);
        }
    }
    else if (seedMethod == RANDOM)
        BoxOfSeeds(globalBounds, seeds);
    else if (seedMethod == RANDOM_BOX)
        BoxOfSeeds(globalBounds, seeds);

    if (seedMethod == RANDOM_BLOCK)
        active.insert(active.end(), seeds.begin(), seeds.end());
    else
    {
        vector<int> domainIds;
        boundsMap.FindBlockIDs(seeds, domainIds);
        for (int i = 0; i < seeds.size(); i++)
            if (DomainToRank(domainIds[i]) == rank)
            {
                seeds[i].blockId = domainIds[i];
                active.push_back(seeds[i]);
            }
    }

    totalNumSeeds = active.size() + inactive.size();

#ifdef VTKH_PARALLEL
    MPI_Allreduce(MPI_IN_PLACE, &totalNumSeeds, 1, MPI_INT, MPI_SUM, mpiComm);
#endif
}

} //  namespace vtkh
