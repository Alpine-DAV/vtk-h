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
        numWorkerThreads(-1),
        done(false),
        begin(false),
        communicator(comm, bmap, pa->GetStats()),
        boundsMap(bmap),
        filter(pa),
        stats(pa->GetStats()),
        sleepUS(100)
    {
        m_Rank = vtkh::GetMPIRank();
        m_NumRanks = vtkh::GetMPISize();
        communicator.RegisterMessages(2, m_NumRanks, m_NumRanks);
        stats->addTimer("worker_sleep");
        stats->addCounter("worker_naps");
    }

    void Init(const std::list<Particle> &particles, int N, int _sleepUS)
    {
        numWorkerThreads = 1;
        TotalNumParticles = N;
        sleepUS = _sleepUS;
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

                DataBlock *blk = filter->GetBlock(particles[0].blockIds[0]);

                stats->start("advect");
                int n = blk->integrator.Trace(particles, filter->GetMaxSteps(), I, T, A, stats, &traces);
                stats->stop("advect");
                stats->increment("advectSteps", n);

                terminated.Insert(T);
                worker_active.Insert(A);
                worker_inactive.Insert(I);
            }
            else
            {
                stats->start("worker_sleep");
                usleep(sleepUS);
                stats->stop("worker_sleep");
                stats->increment("worker_naps");
            }
        }
        DBG("WORKER is DONE"<<std::endl);
        results.Insert(traces);
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
            int val = communicator.Exchange(out, in, term, dT);

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
            {
                stats->start("sleep");
                usleep(sleepUS);
                stats->stop("sleep");
                stats->increment("naps");
            }
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

    int numWorkerThreads;
    int sleepUS;

    bool done, begin;
    vtkh::OpenMPLock stateLock;
    BoundsMap boundsMap;
    ParticleAdvection *filter;
    StatisticsDB *stats;
};
#endif

static int currentStep = 0;   
static std::ofstream *timingInfo = NULL;    
void RecordTime(const std::string &nm, double time)
{
    int rank = 0, numRanks = 0;
#ifdef VTKH_PARALLEL
    rank = vtkh::GetMPIRank();
    numRanks = vtkh::GetMPISize();
#endif
    
    if (timingInfo == NULL)
    {
        timingInfo = new ofstream;
        char nm[32];
        sprintf(nm, "timing.vtkh.%d.out", rank);
        timingInfo->open(nm, ofstream::out);
    }
    (*timingInfo)<<currentStep<<", VTKHTiming_"<<rank<<"_"<<numRanks<<", "<<nm<<", "<<time<<endl;
    //cout<<nm<<" rank "<<rank<<" time "<<time<<endl;
}


ParticleAdvection::ParticleAdvection()
    : rank(0), numRanks(1), seedMethod(RANDOM),
      numSeeds(1000), totalNumSeeds(-1), randSeed(314),
      stepSize(.01),
      maxSteps(1000),
      useThreadedVersion(false),
      gatherTraces(true),
      dumpOutputFiles(false),
      sleepUS(100)
{
#ifdef VTKH_PARALLEL
  rank = vtkh::GetMPIRank();
  numRanks = vtkh::GetMPISize();
#endif

#ifdef TRACE_DEBUG
  dbg = ofstream();
  char nm[32];
  sprintf(nm, "%d.ts%d.out", rank, currentStep);
  dbg.open(nm, ofstream::out);
#endif
}

void
ParticleAdvection::InitStats()
{
  stats.addTimer("total");
  stats.addTimer("sleep");
  stats.addTimer("advect");
  stats.addTimer("VTKmAdvect");
  stats.addCounter("advectSteps");
  stats.addCounter("myParticles");
  stats.addCounter("naps");
}
 
static std::ofstream *out = NULL;
void
ParticleAdvection::DumpStats(const std::string &fname)
{
    stats.calcStats();

    if (rank != 0)
        return;

    if (out == NULL)
    {
      out = new ofstream;
      out->open(fname, ofstream::out);
    }

    (*out)<<std::endl<<std::endl;
    (*out)<<"Step "<<currentStep<<std::endl;
    if (!stats.timers.empty())
    {
        (*out)<<"TIMERS:"<<std::endl;
        for (auto &ti : stats.timers)
            (*out)<<ti.first<<": "<<ti.second.GetTime()<<std::endl;
        (*out)<<std::endl;
        (*out)<<"TIMER_STATS"<<std::endl;
        for (auto &ti : stats.timers)
            (*out)<<ti.first<<" "<<stats.timerStat(ti.first)<<std::endl;
    }
    if (!stats.counters.empty())
    {
        (*out)<<std::endl;
        (*out)<<"COUNTERS:"<<std::endl;
        for (auto &ci : stats.counters)
            (*out)<<ci.first<<" "<<stats.totalVal(ci.first)<<std::endl;
        (*out)<<std::endl;
        (*out)<<"COUNTER_STATS"<<std::endl;
        for (auto &ci : stats.counters)
            (*out)<<ci.first<<" "<<stats.counterStat(ci.first)<<std::endl;
    }
}

ParticleAdvection::~ParticleAdvection()
{
  for (auto p : dataBlocks)
  {
    delete p;
  } 
  dataBlocks.clear();
  currentStep++;
}

void ParticleAdvection::PreExecute()
{
  auto startT = std::chrono::steady_clock::now();

  Filter::PreExecute();

  //RecordTime("ParticleAdevection-FilterPreExecute", std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now()-startT).count());

  //startT = std::chrono::steady_clock::now();
  //Create the bounds map and dataBlocks list.
  boundsMap.Clear();
  //RecordTime("ParticleAdevection-ClearBoundsMapPreExecute", std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now()-startT).count());
  
  //startT = std::chrono::steady_clock::now();
  const int nDoms = this->m_input->GetNumberOfDomains();
  //RecordTime("ParticleAdevection-GetNumDomainsPreExecute", std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now()-startT).count());
  
  auto s = std::chrono::steady_clock::now();
  //double getCount = 0;
  //double dataBlockCount = 0;
  //double boundsMapAddCount = 0;
  for(int i = 0; i < nDoms; i++)
  {
    s = std::chrono::steady_clock::now();
    vtkm::Id id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, id);
    //getCount += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now()-s).count();

    s = std::chrono::steady_clock::now();
    dataBlocks.push_back(new DataBlock(id, &dom, m_field_name, stepSize));
    //dataBlockCount += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now()-s).count();
    
    s = std::chrono::steady_clock::now();
    boundsMap.AddBlock(id, dom.GetCoordinateSystem().GetBounds());
    //boundsMapAddCount += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now()-s).count();
  }
  
  //RecordTime("ParticleAdevection-AddingBlocksFilterPreExecute", std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now()-startT).count());

//RecordTime("ParticleAdevection-GetCountFilterPreExecute", getCount);
//RecordTime("ParticleAdevection-DataBlockCountFilterPreExecute", dataBlockCount);
//RecordTime("ParticleAdevection-BoundsMapAddCountFilterPreExecute", boundsMapAddCount);

  //startT = std::chrono::steady_clock::now();
  boundsMap.Build();
  DBG("PreExecute: "<<boundsMap<<std::endl);

  RecordTime("PreExecute", std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now()-startT).count());
}

void ParticleAdvection::PostExecute()
{
  auto startT = std::chrono::steady_clock::now();

  Filter::PostExecute();

  RecordTime("PostExecute", std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now()-startT).count());
}


/*template <typename ResultT>
void ParticleAdvection::TraceMultiThread(vector<ResultT> &traces)
{
#ifdef VTKH_PARALLEL
  MPI_Comm mpiComm = MPI_Comm_f2c(vtkh::GetMPICommHandle());

  Task *task = new Task(mpiComm, boundsMap, this);

  task->Init(active, totalNumSeeds, sleepUS);
  task->Go();
  task->results.Get(traces);
#endif
}*/

template<>
int
ParticleAdvection::InternalIntegrate<vtkm::worklet::ParticleAdvectionResult>(DataBlock &blk,
                                     std::vector<Particle> &v, 
                                     std::list<Particle> &I,
                                     std::list<Particle> &T,
                                     std::list<Particle> &A,
                                     vector<vtkm::worklet::ParticleAdvectionResult> &traces
                                     )
{
  return blk.integrator.Advect(v, maxSteps, I, T, A, &traces);
}

template<>
int
ParticleAdvection::InternalIntegrate<vtkm::worklet::StreamlineResult>(DataBlock &blk,
                                     std::vector<Particle> &v, 
                                     std::list<Particle> &I,
                                     std::list<Particle> &T,
                                     std::list<Particle> &A,
                                     vector<vtkm::worklet::StreamlineResult> &traces
                                     )
{
  return blk.integrator.Trace(v, maxSteps, I, T, A, GetStats(), &traces);
}

template <typename ResultT>
void ParticleAdvection::TraceSingleThread(vector<ResultT> &traces)
{
#ifdef VTKH_PARALLEL
  MPI_Comm mpiComm = MPI_Comm_f2c(vtkh::GetMPICommHandle());

  ParticleMessenger communicator(mpiComm, boundsMap, GetStats());
  communicator.RegisterMessages(2, 64, 64);
  //Trying smaller number as to not overwhelm system at larger scale
  //communicator.RegisterMessages(2, numRanks, numRanks);
  const int nDoms = this->m_input->GetNumberOfDomains();
  for (int i = 0; i < nDoms; i++)
  {
    vtkm::Id id;
    vtkm::cont::DataSet ds;
    this->m_input->GetDomain(i, ds, id);
//auto startT = std::chrono::steady_clock::now();
    communicator.AddLocator(id, ds);
//RecordTime("TraceSingleThread-AddLocator", std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now()-startT).count());
  }

  int N = 0;
  int prevTermCount = 0;

  while (true)
  {
      DBG("MANAGE: termCount= "<<terminated.size()<<std::endl<<std::endl);
      std::vector<Particle> v;
      std::list<Particle> I;

      if (GetActiveParticles(v))
      {
          stats.increment("myParticles", v.size());
          DBG("GetActiveParticles: "<<v<<std::endl);
          DataBlock *blk = GetBlock(v[0].blockIds[0]);
          DBG("Integrate: "<<v<<std::endl);
          DBG("Loading Block: "<<v[0].blockIds[0]<<std::endl);

          std::list<Particle> T, A;
          stats.start("advect");
          int n;
          n = InternalIntegrate<ResultT>(*blk, v, I, T, A, traces);
          stats.stop("advect");
          stats.increment("advectSteps", n);
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
      int val = communicator.Exchange(I, in, term, dT);

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
      {
          stats.start("sleep");
          usleep(sleepUS);
          stats.stop("sleep");
          stats.increment("naps");
      }
  }
  DBG("TIA: "<<terminated.size()<<" "<<inactive.size()<<" "<<active.size()<<std::endl);
  DBG("RESULTS= "<<traces.size()<<std::endl);

//  DumpTraces(rank, traces);
//  std::cout<<rank<<": done tracing. # of traces= "<<traces.size()<<std::endl;
  MPI_Barrier(mpiComm);
  DBG("All done"<<std::endl);

if(0)
{
//save number of steps for each individual particle
  vector<int> outStuff(totalNumSeeds, 0);
  vector<int> inStuff(totalNumSeeds, 0);
  std::list<Particle>::iterator listIt = terminated.begin();
  while (listIt != terminated.end())
  {
    outStuff[(*listIt).id] = (*listIt).nSteps;
    listIt++;
  }
  MPI_Reduce(outStuff.data(), inStuff.data(), totalNumSeeds, MPI_INT, MPI_SUM, 0, mpiComm);

  if(rank == 0)
  {
    ofstream output;
    char nm[128];
    sprintf(nm, "stepCounts.ts%03i.txt", currentStep);
    output.open(nm, ofstream::out);
    output << "seed id's and counts" << endl;
    for(int i = 0; i < totalNumSeeds; i++)
    {
      output<<"["<<i<<"] "<<inStuff[i] << endl;
    }
  }
}

#endif
}

template <typename ResultT>
void ParticleAdvection::TraceSeeds(vector<ResultT> &traces)
{
  stats.start("total");

//  if (useThreadedVersion)
//      //TraceMultiThread<ResultT>(traces);
//  else
      TraceSingleThread<ResultT>(traces);

  stats.stop("total");
  DumpStats("particleAdvection.stats.txt");
}

void ParticleAdvection::DoExecute()
{
  auto startT = std::chrono::steady_clock::now();

  this->Init();
  this->CreateSeeds();
  if (this->dumpOutputFiles)
    this->DumpDS();

  if(!gatherTraces)
  {
    vector<vtkm::worklet::ParticleAdvectionResult> particleTraces;
    this->TraceSeeds<vtkm::worklet::ParticleAdvectionResult>(particleTraces);
    this->m_output = new DataSet();
  }
  else
  {
    vector<vtkm::worklet::StreamlineResult> particleTraces;
    this->TraceSeeds<vtkm::worklet::StreamlineResult>(particleTraces);

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
        if (this->dumpOutputFiles)
            this->DumpSLOutput(&ds, rank);
    }
    else
    {
        if (this->dumpOutputFiles)
            this->DumpSLOutput(NULL, rank);
    }
  }

  RecordTime("DoExecute", std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now()-startT).count());
}


bool
ParticleAdvection::GetActiveParticles(std::vector<Particle> &v)
{
    v.clear();
    if (active.empty())
        return false;

    int workingBlockID = active.front().blockIds[0];

    std::list<Particle>::iterator listIt = active.begin();
    while (listIt != active.end())
    {
      Particle p = *listIt;
      if(workingBlockID == p.blockIds[0])
      {
        v.push_back(p);
        listIt = active.erase(listIt);
      }
      else
      {
        listIt++;
      }
    }
    return !v.empty();
}

std::string
ParticleAdvection::GetName() const
{
  return "vtkh::ParticleAdvection";
}

void
ParticleAdvection::DumpSLOutput(vtkm::cont::DataSet *ds, int domId)
{
  char nm[128];
  if (ds)
  {
      sprintf(nm, "ds.ts%03i.block%03d.vtk", currentStep, domId);
      vtkm::io::writer::VTKDataSetWriter writer(nm);
      writer.WriteDataSet(*ds);
  }
  if (rank == 0)
  {
    ofstream output;
    output.open("ds.visit", ofstream::out);
    output<<"!NBLOCKS "<<numRanks<<std::endl;
    for (int i = 0; i < numRanks; i++)
    {
        char nm[128];
        sprintf(nm, "ds.ts%03i.block%03d.vtk", currentStep, i);
        output<<nm<<std::endl;
    }
  }

  sprintf(nm, "pts.ts%03i.block%03d.txt", currentStep, domId);
  ofstream pout;
  pout.open(nm, ofstream::out);

  if (ds)
  {
      int nPts = ds->GetCoordinateSystem(0).GetNumberOfPoints();
      auto portal = ds->GetCoordinateSystem(0).GetData().GetPortalConstControl();
      for (int i = 0; i < nPts; i++)
      {
          vtkm::Vec<float,3> pt;
          pt = portal.Get(i);
          pout<<pt[0]<<","<<pt[1]<<","<<pt[2]<<","<<i<<std::endl;
      }
  }
  //else
  //    pout<<"-12,-12,-12,-1"<<std::endl;

  if (rank == 0)
  {
    ofstream output;
    output.open("pts.visit", ofstream::out);
    output<<"!NBLOCKS "<<numRanks<<std::endl;
    for (int i = 0; i < numRanks; i++)
    {
        char nm[128];
        sprintf(nm, "pts.ts%03i.block%03d.txt", currentStep, i);
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
    sprintf(nm, "dom.ts%03i.block%03d.vtk", currentStep, domId);

    vtkm::io::writer::VTKDataSetWriter writer(nm);
    writer.WriteDataSet(dom);

  }

  if (rank == 0)
  {
    ofstream output;
    output.open("dom.visit", ofstream::out);
    output<<"!NBLOCKS "<<totalNumDoms<<std::endl;
    for (int i = 0; i < totalNumDoms; i++)
    {
      char nm[128];
      sprintf(nm, "dom.ts%03i.block%03d.vtk", currentStep, i);
      output<<nm<<std::endl;
    }
  }
}

void
ParticleAdvection::Init()
{
    InitStats();
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
    const float factor = 0.025;
    float d[3] = {(boxRange[1]-boxRange[0]) * factor,
                  (boxRange[3]-boxRange[2]) * factor,
                  (boxRange[5]-boxRange[4]) * factor};
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
      p.id = id;
      p.coords[0] = randRange(boxRange[0], boxRange[1]);
      p.coords[1] = randRange(boxRange[2], boxRange[3]);
      p.coords[2] = randRange(boxRange[4], boxRange[5]);

      /*
      if (N == 1)
      {
          p.coords[0] = -1.0;
          p.coords[1] = -3.0;
          p.coords[2] = -5.0;

          p.coords[0] = -1.0;
          p.coords[1] = -0.6;
          p.coords[2] =  0.1;
      }
      */

      seeds.push_back(p);
  }
  DBG("Seeds: "<<seeds<<std::endl);
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
        BoxOfSeeds(boundsMap.globalBounds, seeds);
    else if (seedMethod == RANDOM_BOX)
        BoxOfSeeds(seedBox, seeds);
    else if (seedMethod == POINT)
    {
        Particle p(seedPoint, 0);
        seeds.push_back(p);
    }

    //Set the blockIds for each seed.
    std::vector<std::vector<int>> domainIds;
    boundsMap.FindBlockIDs(seeds, domainIds);
    for (int i = 0; i < seeds.size(); i++)
    {
        if (!domainIds[i].empty() && DomainToRank(domainIds[i][0]) == rank)
        {
            seeds[i].blockIds = domainIds[i];
            active.push_back(seeds[i]);
            if (domainIds[i].size() > 1) DBG("WE have a DUP: "<<seeds[i]<<std::endl);

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
    sprintf(nm, "output.ts%03i.block%03d.txt", currentStep, idx);
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
            sprintf(nm, "output.ts%03i.block%03d.txt", currentStep, i);
            output<<nm<<std::endl;
        }
    }
}

} //  namespace vtkh
