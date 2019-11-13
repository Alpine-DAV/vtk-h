#include <iostream>
#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkh/filters/ParticleAdvection.hpp>
#include <vtkm/filter/Streamline.h>
#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>
#include <vtkm/cont/Algorithm.h>

#include <vtkh/vtkh.hpp>
#include <vtkh/Error.hpp>
#include <vtkh/Logger.hpp>
#include <vtkh/utils/StreamUtil.hpp>
#include <vtkh/utils/ThreadSafeContainer.hpp>

#include <vector>
#include <thread>

#ifdef VTKH_PARALLEL
#include <mpi.h>
#include <vtkh/filters/communication/ParticleMessenger.hpp>
#include <vtkh/filters/ParticleAdvectionTask.hpp>
#endif

#ifdef VTKH_ENABLE_LOGGING
#define DBG(msg) vtkh::Logger::GetInstance("out")->GetStream()<<msg
#define WDBG(msg) vtkh::Logger::GetInstance("wout")->GetStream()<<msg
#else
#define DBG(msg)
#define WDBG(msg)
#endif

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
}

ParticleAdvection::~ParticleAdvection()
{
  for (auto p : dataBlocks)
    delete p;
  dataBlocks.clear();
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

    dataBlocks.push_back(new DataBlockIntegrator(id, &dom, m_field_name, stepSize));
    boundsMap.AddBlock(id, dom.GetCoordinateSystem().GetBounds());
  }

  boundsMap.Build();
}

void ParticleAdvection::PostExecute()
{
  Filter::PostExecute();
}

template <typename ResultT>
void ParticleAdvection::TraceMultiThread(std::vector<ResultT> &traces)
{
#ifdef VTKH_PARALLEL
  MPI_Comm mpiComm = MPI_Comm_f2c(vtkh::GetMPICommHandle());

  vtkh::ParticleAdvectionTask<ResultT> *task = new vtkh::ParticleAdvectionTask<ResultT>(mpiComm, boundsMap, this);

  task->Init(active, totalNumSeeds, sleepUS);
  task->Go();
  task->results.Get(traces);
#endif
}

template<>
int
ParticleAdvection::InternalIntegrate<vtkm::worklet::ParticleAdvectionResult>(DataBlockIntegrator &blk,
                                     std::vector<Particle> &v,
                                     std::vector<Particle> &I,
                                     std::vector<Particle> &T,
                                     std::vector<Particle> &A,
                                     std::vector<vtkm::worklet::ParticleAdvectionResult> &traces
                                     )
{
  return blk.integrator.Advect(v, maxSteps, I, T, A, &traces);
}

template<>
int
ParticleAdvection::InternalIntegrate<vtkm::worklet::StreamlineResult>(DataBlockIntegrator &blk,
                                     std::vector<Particle> &v,
                                     std::vector<Particle> &I,
                                     std::vector<Particle> &T,
                                     std::vector<Particle> &A,
                                     std::vector<vtkm::worklet::StreamlineResult> &traces
                                     )
{
  return blk.integrator.Trace(v, maxSteps, I, T, A, &traces);
}

template <typename ResultT>
void ParticleAdvection::TraceSingleThread(std::vector<ResultT> &traces)
{
#ifdef VTKH_PARALLEL
  MPI_Comm mpiComm = MPI_Comm_f2c(vtkh::GetMPICommHandle());

  ParticleMessenger communicator(mpiComm, boundsMap);
  communicator.RegisterMessages(2, std::min(64, numRanks-1), 128, std::min(64, numRanks-1));

  int N = 0;
  while (true)
  {
      DBG("MANAGE: termCount= "<<terminated.size()<<std::endl<<std::endl);
      std::vector<Particle> v, I, T, A;

      if (GetActiveParticles(v))
      {
          COUNTER_INC("myParticles", v.size());
          DBG("GetActiveParticles: "<<v<<std::endl);
          DataBlockIntegrator *blk = GetBlock(v[0].blockIds[0]);
          DBG("Integrate: "<<v<<std::endl);
          DBG("Loading Block: "<<v[0].blockIds[0]<<std::endl);

          TIMER_START("advect");
          int n = InternalIntegrate<ResultT>(*blk, v, I, T, A, traces);
          TIMER_STOP("advect");
          COUNTER_INC("advectSteps", n);
          DBG("--Integrate:  ITA: "<<I<<" "<<T<<" "<<A<<std::endl);
          DBG("                   I= "<<I<<std::endl);
          if (!A.empty())
              active.insert(active.end(), A.begin(), A.end());
      }

      std::vector<Particle> in;
      int numTermMessages;
      communicator.Exchange(I, in, T, numTermMessages);
      int numTerm = T.size() + numTermMessages;

      if (!in.empty())
          active.insert(active.end(), in.begin(), in.end());
      if (!T.empty())
          terminated.insert(terminated.end(), T.begin(), T.end());

      N += numTerm;

      DBG("Manage: N= "<<N<<std::endl);
      if (N > totalNumSeeds)
          throw "Particle count error";

      if (N == totalNumSeeds)
          break;

      if (active.empty())
      {
          TIMER_START("sleep");
          usleep(sleepUS);
          TIMER_STOP("sleep");
          COUNTER_INC("naps", 1);
      }
  }
  DBG("TIA: "<<terminated.size()<<" "<<inactive.size()<<" "<<active.size()<<std::endl);
  DBG("RESULTS= "<<traces.size()<<std::endl);

  DBG("All done"<<std::endl);
#endif
}

template <typename ResultT>
void ParticleAdvection::TraceSeeds(std::vector<ResultT> &traces)
{
  TIMER_START("total");

  if (useThreadedVersion)
      TraceMultiThread<ResultT>(traces);
  else
      TraceSingleThread<ResultT>(traces);

  TIMER_STOP("total");
  DUMP_STATS("particleAdvection.stats.txt");
}

void ParticleAdvection::DoExecute()
{
  this->Init();
  this->CreateSeeds();
  if (this->dumpOutputFiles)
    this->DumpDS(0);

  if(!gatherTraces)
  {
    std::vector<vtkm::worklet::ParticleAdvectionResult> particleTraces;
    this->TraceSeeds<vtkm::worklet::ParticleAdvectionResult>(particleTraces);
    this->m_output = new DataSet();
  }
  else
  {
    std::vector<vtkm::worklet::StreamlineResult> particleTraces;
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
        ds.SetCellSet(polyLines);
        this->m_output->AddDomain(ds, rank);
        if (this->dumpOutputFiles)
            this->DumpSLOutput(&ds, rank, 0);
    }
    else
    {
        if (this->dumpOutputFiles)
            this->DumpSLOutput(NULL, rank, 0);
    }
  }
}


bool
ParticleAdvection::GetActiveParticles(std::vector<Particle> &v)
{
    v.clear();
    if (active.empty())
        return false;

    int workingBlockID = active.front().blockIds[0];

    std::vector<Particle>::iterator listIt = active.begin();
    while (listIt != active.end())
    {
      Particle p = *listIt;
      if(workingBlockID == p.blockIds[0])
      {
        v.push_back(p);
        listIt = active.erase(listIt);
      }
      else
        listIt++;
    }
    return !v.empty();
}

std::string
ParticleAdvection::GetName() const
{
  return "vtkh::ParticleAdvection";
}

void
ParticleAdvection::DumpSLOutput(vtkm::cont::DataSet *ds, int domId, int ts)
{
  char nm[128];
  if (ds)
  {
      sprintf(nm, "ds.ts%03i.block%03d.vtk", ts, domId);
      vtkm::io::writer::VTKDataSetWriter writer(nm);
      writer.WriteDataSet(*ds);
  }
  if (rank == 0)
  {
    std::ofstream output;
    output.open("ds.visit", std::ofstream::out);
    output<<"!NBLOCKS "<<numRanks<<std::endl;
    for (int i = 0; i < numRanks; i++)
    {
        char nm[128];
        sprintf(nm, "ds.ts%03i.block%03d.vtk", ts, i);
        output<<nm<<std::endl;
    }
  }

  sprintf(nm, "pts.ts%03i.block%03d.txt", ts, domId);
  std::ofstream pout;
  pout.open(nm, std::ofstream::out);

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

  if (rank == 0)
  {
    std::ofstream output;
    output.open("pts.visit", std::ofstream::out);
    output<<"!NBLOCKS "<<numRanks<<std::endl;
    for (int i = 0; i < numRanks; i++)
    {
        char nm[128];
        sprintf(nm, "pts.ts%03d.block%03d.txt", ts, i);
        output<<nm<<std::endl;
    }
  }
}

void
ParticleAdvection::DumpDS(int ts)
{
  int totalNumDoms = this->m_input->GetGlobalNumberOfDomains();
  int nDoms = this->m_input->GetNumberOfDomains();

  for (int i = 0; i < nDoms; i++)
  {
    vtkm::cont::DataSet dom;
    vtkm::Id domId;
    this->m_input->GetDomain(i, dom, domId);

    char nm[128];
    sprintf(nm, "dom.ts%03d.block%03d.vtk", ts, (int)domId);

    vtkm::io::writer::VTKDataSetWriter writer(nm);
    writer.WriteDataSet(dom);

  }

  if (rank == 0)
  {
    std::ofstream output;
    output.open("dom.visit", std::ofstream::out);
    output<<"!NBLOCKS "<<totalNumDoms<<std::endl;
    for (int i = 0; i < totalNumDoms; i++)
    {
      char nm[128];
      sprintf(nm, "dom.ts%03d.block%03d.vtk", ts, i);
      output<<nm<<std::endl;
    }
  }
}

void
ParticleAdvection::Init()
{
  //Initialize timers/counters.
  ADD_TIMER("total");
  ADD_TIMER("sleep");
  ADD_TIMER("advect");
  ADD_COUNTER("advectSteps");
  ADD_COUNTER("myParticles");
  ADD_COUNTER("naps");

}

DataBlockIntegrator *
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
  //DBG("Box of Seeds: N= "<<numSeeds<<" "<<box<<std::endl);
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
      seeds.push_back(p);
  }
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
ParticleAdvection::DumpTraces(int ts, const std::vector<vtkm::Vec<double,4>> &particleTraces)
{
    std::ofstream output;
    char nm[128];
    sprintf(nm, "output.ts%03i.block%03d.txt", ts, rank);
    output.open(nm, std::ofstream::out);

    output<<"X,Y,Z,ID"<<std::endl;
    for (auto &p : particleTraces)
        output<<p[0]<<", "<<p[1]<<", "<<p[2]<<", "<<(int)p[3]<<std::endl;

    if (rank == 0)
    {
        std::ofstream output;
        output.open("output.visit", std::ofstream::out);
        output<<"!NBLOCKS "<<numRanks<<std::endl;
        for (int i = 0; i < numRanks; i++)
        {
            char nm[128];
            sprintf(nm, "output.ts%03i.block%03d.txt", ts, i);
            output<<nm<<std::endl;
        }
    }
}

} //  namespace vtkh
