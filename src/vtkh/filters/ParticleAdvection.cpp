#include <iostream>
#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkh/filters/ParticleAdvection.hpp>
#include <vtkh/vtkh.hpp>
#include <vtkh/Error.hpp>

#include <vtkm/filter/Streamline.h>
#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>
#include <vtkm/cont/Algorithm.h>

#include <vtkh/filters/util.hpp>
#include <vector>
#include <omp.h>
#include <vtkh/filters/communication/ThreadSafeContainer.hpp>

#ifdef VTKH_PARALLEL
#include <mpi.h>
#include <vtkh/filters/communication/ParticleMessenger.hpp>
#endif

#include <vtkh/filters/communication/DebugMeowMeow.hpp>

//#define VTKH_MULTI_THREAD_PARTICLE_ADVECTION

std::ofstream dbg;

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

class OpenMPLock
{
public:
    OpenMPLock()
    {
        omp_init_lock(&lock);
    }
    ~OpenMPLock()
    {
        omp_destroy_lock(&lock);
    }

    void set() {omp_set_lock(&lock);}
    void unset() {omp_unset_lock(&lock);}

private:
    omp_lock_t lock;
};

#ifdef VTKH_PARALLEL
class Task
{
public:
    Task(MPI_Comm comm, const vtkh::BoundsMap &bmap, ParticleAdvection *pa) :
        managerThread(false), numWorkerThreads(-1), done(false), begin(false),
        communicator(comm, bmap), boundsMap(bmap), filter(pa)
    {
        m_Rank = vtkh::GetMPIRank();
        m_NumRanks = vtkh::GetMPISize();
        communicator.RegisterMessages(2, m_NumRanks, m_NumRanks);
    }

    void Init(const std::list<Particle> &particles, int N)
    {
        managerThread = true;
        numWorkerThreads = 1;
        TotalNumParticles = N;
        active.Assign(particles);
        inactive.Clear();
        terminated.Clear();
    }

    bool CheckDone()
    {
        bool val;
        stateLock.set();
        val = done;
        stateLock.unset();
        return val;
    }
    void SetDone()
    {
        stateLock.set();
        done = true;
        stateLock.unset();
    }

    bool GetBegin()
    {
        bool val;
        stateLock.set();
        val = begin;
        stateLock.unset();
        return val;
    }

    void SetBegin()
    {
        stateLock.set();
        begin = true;
        stateLock.unset();
    }

    void Go()
    {
        DBG("Go_bm: "<<boundsMap<<std::endl);
        DBG("actives= "<<active<<std::endl);

        #pragma omp parallel sections num_threads(2)
        {
            #pragma omp section
            #pragma omp parallel num_threads(1)
            {
                this->Manage();
            }

            #pragma omp section
            #pragma omp parallel num_threads(numWorkerThreads)
            #pragma omp master
            {
                this->Work();
            }
        }
    }

    void Work()
    {
        vector<vtkm::worklet::StreamlineResult> traces;
        DBG("work_bm: "<<boundsMap<<std::endl);

        while (!CheckDone())
        {
            std::vector<Particle> particles;
            if (active.Get(particles))
            {
                std::list<Particle> I, T, A;

                DataBlock *blk = filter->GetBlock(particles[0].blockId);

                blk->integrator.Trace(particles, filter->GetMaxSteps(), I, T, A, &traces);

                terminated.Insert(T);
                worker_active.Insert(A);
                worker_inactive.Insert(I);
            }
            else
            {
                usleep(1000);
            }
        }
        DBG("WORKER is DONE"<<std::endl);
        results.Insert(traces);
    }

    void Manage_works()
    {
        DBG("manage_bm: "<<boundsMap<<std::endl);
        int N = 0;

        int prevTermCount = 0;

        while (true)
        {
            DBG("MANAGE: termCount= "<<terminated.Size()<<std::endl<<std::endl);
            std::list<Particle> out, in, term;
            worker_inactive.Get(out);

            if (!out.empty())
                DBG("M-Send: "<<out<<std::endl);
            communicator.ExchangeParticles(out, in, term);
            if (!in.empty())
                DBG("M-Recv: "<<in<<std::endl);

            if (!in.empty())
                active.Insert(in);
            if (!term.empty())
                terminated.Insert(term);

            int numTerm = terminated.Size();
            int dT = numTerm - prevTermCount;
            int x = communicator.ExchangeCounter(dT);
            dT += x;
            N += dT;

            prevTermCount = numTerm;

            DBG("Manage: N= "<<N<<std::endl);
            if (N >= TotalNumParticles)
                break;

            if (active.Empty())
                usleep(1e6);
        }
        DBG("TIA: "<<terminated.Size()<<" "<<inactive.Size()<<" "<<active.Size()<<std::endl);
        DBG("RESULTS= "<<results.Size()<<std::endl);

        SetDone();
    }


    void Manage()
    {
        DBG("manage_bm: "<<boundsMap<<std::endl);

        int N = 0;
        int prevTermCount = 0;

        while (true)
        {
            DBG("MANAGE: termCount= "<<terminated.Size()<<std::endl<<std::endl);
            std::list<Particle> out, in, term;
            worker_inactive.Get(out);

            int numTerm = terminated.Size();
            int dT = numTerm - prevTermCount;
            int val = communicator.Exchange2(out, in, term, dT);

            if (!in.empty())
                active.Insert(in);
            if (!term.empty())
                terminated.Insert(term);

            N += (dT+val);
            prevTermCount = numTerm;

            DBG("Manage: N= "<<N<<std::endl);
            if (N >= TotalNumParticles)
                break;

            if (active.Empty())
                usleep(1e6);
        }
        DBG("TIA: "<<terminated.Size()<<" "<<inactive.Size()<<" "<<active.Size()<<std::endl);
        DBG("RESULTS= "<<results.Size()<<std::endl);

        SetDone();
    }

    int m_Rank, m_NumRanks;
    int TotalNumParticles;

    using ParticleList = vtkh::ThreadSafeContainer<Particle, std::list, vtkh::OpenMPLock>;
    using ResultsVec = vtkh::ThreadSafeContainer<vtkm::worklet::StreamlineResult, std::vector, vtkh::OpenMPLock>;

    ParticleMessenger communicator;
    ParticleList active, inactive, terminated;
    ParticleList worker_active, worker_inactive;
    ResultsVec results;

    bool managerThread;
    int numWorkerThreads;

    bool done, begin;
    vtkh::OpenMPLock stateLock;
    BoundsMap boundsMap;
    ParticleAdvection *filter;
};
#endif

ParticleAdvection::ParticleAdvection()
    : rank(0), numRanks(1), seedMethod(RANDOM),
      numSeeds(1000), totalNumSeeds(-1), randSeed(314),
      stepSize(.01),
      maxSteps(1000)
{
#ifdef VTKH_PARALLEL
  rank = vtkh::GetMPIRank();
  numRanks = vtkh::GetMPISize();
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

  //Create the bounds map and dataBlocks list.
  boundsMap.Clear();
  const int nDoms = this->m_input->GetNumberOfDomains();
  for(int i = 0; i < nDoms; i++)
  {
    vtkm::Id id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, id);

    dataBlocks.push_back(new DataBlock(id, dom, m_field_name, stepSize));

    boundsMap.AddBlock(id, dom.GetCoordinateSystem().GetBounds());
  }

  boundsMap.Build();
  DBG("PreExecute: "<<boundsMap<<std::endl);
}

void ParticleAdvection::PostExecute()
{
  Filter::PostExecute();
}

void ParticleAdvection::TraceSeeds(vector<vtkm::worklet::StreamlineResult> &traces)
{
#ifdef VTKH_PARALLEL
  MPI_Comm mpiComm = MPI_Comm_f2c(vtkh::GetMPICommHandle());
  std::cout<<rank<<": "<<totalNumSeeds<<" "<<active.size()<<std::endl;

#ifdef VTKH_MULTI_THREAD_PARTICLE_ADVECTION

  Task *task = new Task(mpiComm, boundsMap, this);

  std::cout<<boundsMap<<std::endl;

  task->Init(active, totalNumSeeds);
  task->Go();
  task->results.Get(traces);
  DBG("We have some results: "<<traces.size()<<std::endl);

#else

  ParticleMessenger communicator(mpiComm, boundsMap);
  communicator.RegisterMessages(2, numRanks, numRanks);

  int N = 0;
  int prevTermCount = 0;

  while (true)
  {
      DBG("MANAGE: termCount= "<<terminated.size()<<std::endl<<std::endl);
      std::vector<Particle> v;
      std::list<Particle> I;

      if (GetActiveParticles(v))
      {
          DBG("GetActiveParticles: "<<v<<std::endl);
          DataBlock *blk = GetBlock(v[0].blockId);
          DBG("Integrate: "<<v<<std::endl);

          std::list<Particle> T, A;
          blk->integrator.Trace(v, maxSteps, I, T, A, &traces);
          DBG("--Integrate:  ITA: "<<I<<" "<<T<<" "<<A<<std::endl);
          DBG("                   I= "<<I<<std::endl);
          if (!T.empty())
              terminated.insert(terminated.end(), T.begin(), T.end());
          if (!A.empty())
              active.insert(active.end(), A.begin(), A.end());
      }

      int numTerm = terminated.size();
      int dT = numTerm - prevTermCount;
      std::list<Particle> in, term;
      int val = communicator.Exchange2(I, in, term, dT);

      if (!in.empty())
          active.insert(active.end(), in.begin(), in.end());
      if (!term.empty())
          terminated.insert(terminated.end(), term.begin(), term.end());
      numTerm = terminated.size();
      dT = numTerm - prevTermCount;

      N += (dT+val);
      prevTermCount = numTerm;

      DBG("Manage: N= "<<N<<std::endl);
      if (N >= totalNumSeeds)
          break;

      if (active.empty())
          usleep(1e6);
  }
  DBG("TIA: "<<terminated.size()<<" "<<inactive.size()<<" "<<active.size()<<std::endl);
  DBG("RESULTS= "<<traces.size()<<std::endl);

//  DumpTraces(rank, traces);
  std::cout<<rank<<": done tracing. # of traces= "<<traces.size()<<std::endl;
  DBG("All done"<<std::endl);

#endif

#endif
}

void ParticleAdvection::DoExecute()
{
  this->Init();
  this->CreateSeeds();
  this->DumpDS();

  vector<vtkm::worklet::StreamlineResult> particleTraces;
  this->TraceSeeds(particleTraces);

  this->m_output = new DataSet();

  //Compact all the traces into a single dataset.
  if (!particleTraces.empty())
  {
      //Collapse particle trace data into one set of polylines.
      vtkm::Id totalNumPts = 0, totalNumCells = 0;
      std::vector<vtkm::Id> offsets(particleTraces.size(), 0), numLines(particleTraces.size(),0);
      for (int i = 0; i < particleTraces.size(); i++)
      {
          vtkm::Id n = particleTraces[i].positions.GetNumberOfValues();
          offsets[i] = n;
          totalNumPts += n;
          n = particleTraces[i].polyLines.GetNumberOfCells();
          totalNumCells += n;
          numLines[i] = n;

      }

      //Append all the positions into one array.
      vtkm::cont::ArrayHandle<vtkm::Vec<double, 3>> positions;
      positions.Allocate(totalNumPts);
      auto posPortal = positions.GetPortalControl();
      vtkm::Id idx = 0;
      for (int i = 0; i < particleTraces.size(); i++)
      {
          auto inP = particleTraces[i].positions.GetPortalConstControl();
          for (int j = 0; j < offsets[i]; j++, idx++)
              posPortal.Set(idx, inP.Get(j));
      }

      //Cell types are all lines...
      vtkm::cont::ArrayHandle<vtkm::UInt8> cellTypes;
      cellTypes.Allocate(totalNumCells);
      vtkm::cont::ArrayHandleConstant<vtkm::UInt8> polyLineShape(vtkm::CELL_SHAPE_LINE, totalNumCells);
      vtkm::cont::Algorithm::Copy(polyLineShape, cellTypes);

      //Append all the conn and cellCounts.
      vtkm::cont::ArrayHandle<vtkm::Id> connectivity;
      vtkm::cont::ArrayHandle<vtkm::IdComponent> cellCounts;
      connectivity.Allocate(totalNumPts);
      cellCounts.Allocate(totalNumCells);
      auto connPortal = connectivity.GetPortalControl();
      auto cntPortal = cellCounts.GetPortalControl();

      vtkm::Id offset = 0, connIdx = 0, cntIdx = 0;
      for (int i = 0; i < particleTraces.size(); i++)
      {
          if (i > 0)
              offset = offsets[i-1];

          vtkm::Id n = particleTraces[i].polyLines.GetNumberOfCells();
          vtkm::cont::ArrayHandle<vtkm::Id> ids;

          for (vtkm::Id j = 0; j < n; j++)
          {
              particleTraces[i].polyLines.GetIndices(j, ids);
              vtkm::Id nids = ids.GetNumberOfValues();
              auto idsPortal = ids.GetPortalControl();
              for (vtkm::Id k = 0; k < nids; k++, connIdx++)
                  connPortal.Set(connIdx, idsPortal.Get(k)+offset);
              cntPortal.Set(cntIdx, nids);
              cntIdx++;
          }
      }

      //Create a single polyLines cell set.
      vtkm::cont::CellSetExplicit<> polyLines;
      polyLines.Fill(positions.GetNumberOfValues(), cellTypes, cellCounts, connectivity);

      vtkm::cont::DataSet ds;
      vtkm::cont::CoordinateSystem outputCoords("coordinates", positions);
      ds.AddCoordinateSystem(outputCoords);
      ds.AddCellSet(polyLines);
      this->m_output->AddDomain(ds, rank);

      this->DumpSLOutput(ds, rank);
  }
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
ParticleAdvection::DumpSLOutput(const vtkm::cont::DataSet &ds, int domId)
{
  char nm[128];
  sprintf(nm, "ds.%03d.vtk", domId);
  vtkm::io::writer::VTKDataSetWriter writer(nm);
  writer.WriteDataSet(ds);
  if (rank == 0)
  {
    ofstream output;
    output.open("ds.visit", ofstream::out);
    output<<"!NBLOCKS "<<numRanks<<std::endl;
    for (int i = 0; i < numRanks; i++)
    {
        char nm[128];
        sprintf(nm, "ds.%03d.vtk", i);
        output<<nm<<std::endl;
    }
  }
}

void
ParticleAdvection::DumpDS()
{
  int totalNumDoms = this->m_input->GetGlobalNumberOfDomains();
  int nDoms = this->m_input->GetNumberOfDomains();

  for (int i = 0; i < nDoms; i++)
  {
    vtkm::cont::DataSet dom;
    vtkm::Id domId;
    this->m_input->GetDomain(i, dom, domId);

    char nm[128];
    sprintf(nm, "dom.%03lld.vtk", domId);

    vtkm::io::writer::VTKDataSetWriter writer(nm);
    writer.WriteDataSet(dom);
  }

  if (rank == 0)
  {
    ofstream output;
    output.open("dom.visit", ofstream::out);
    output<<"!NBLOCKS "<<numRanks<<std::endl;
    for (int i = 0; i < totalNumDoms; i++)
    {
      char nm[128];
      sprintf(nm, "dom.%03d.vtk", i);
      output<<nm<<std::endl;
    }
  }
}

void
ParticleAdvection::Init()
{
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
  DBG("Box of Seeds: N= "<<numSeeds<<" "<<box<<std::endl);
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
}

void
ParticleAdvection::CreateSeeds()
{
    std::cout<<"CreateSeeds: "<<seedMethod<<std::endl;
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
        BoxOfSeeds(boundsMap.globalBounds, seeds);
    else if (seedMethod == RANDOM_BOX)
        BoxOfSeeds(seedBox, seeds);
    else if (seedMethod == POINT)
    {
        Particle p(seedPoint, 0);
        seeds.push_back(p);
    }

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
    MPI_Comm mpiComm = MPI_Comm_f2c(vtkh::GetMPICommHandle());
    MPI_Allreduce(MPI_IN_PLACE, &totalNumSeeds, 1, MPI_INT, MPI_SUM, mpiComm);
#endif
}

void
ParticleAdvection::DumpTraces(int idx, const vector<vtkm::Vec<double,4>> &particleTraces)
{
    ofstream output;
    char nm[128];
    sprintf(nm, "output.%03d.txt", idx);
    output.open(nm, ofstream::out);

    output<<"X,Y,Z,ID"<<std::endl;
    for (auto &p : particleTraces)
        output<<p[0]<<", "<<p[1]<<", "<<p[2]<<", "<<(int)p[3]<<std::endl;

    if (idx == 0)
    {
        ofstream output;
        output.open("output.visit", ofstream::out);
        output<<"!NBLOCKS "<<numRanks<<std::endl;
        for (int i = 0; i < numRanks; i++)
        {
            char nm[128];
            sprintf(nm, "output.%03d.txt", i);
            output<<nm<<std::endl;
        }
    }
}

} //  namespace vtkh
