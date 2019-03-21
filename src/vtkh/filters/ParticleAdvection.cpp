#include <iostream>
#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkh/filters/ParticleAdvection.hpp>
#include <vtkh/filters/Communicator.hpp>
#include <vtkh/filters/util.hpp>
#include <vtkh/vtkh.hpp>
#include <vtkh/Error.hpp>

#include <vtkm/filter/Streamline.h>
#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>
#include <vtkm/cont/Algorithm.h>

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
    cout<<"Streamline dtor"<<endl;
}

void ParticleAdvection::PreExecute()
{
  Filter::PreExecute();
}

void ParticleAdvection::PostExecute()
{
  Filter::PostExecute();
}

void ParticleAdvection::TraceSeeds(vector<vtkm::worklet::StreamlineResult> &traces)
{
#ifdef VTKH_PARALLEL
  MPI_Comm mpiComm = MPI_Comm_f2c(vtkh::GetMPICommHandle());
  cout<<rank<<" dom2rank: "<<domToRank<<endl;
  cout<<rank<<": "<<totalNumSeeds<<" "<<active.size()<<endl;

  int N = totalNumSeeds;
  int cnt = 0;
  MPICommunicator communicator(mpiComm, domToRank);

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
          blk->integrator.Trace(v, maxSteps, I, T, A, &traces);
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

      // no work available, take a snooze.
      if (active.empty())
      {
          DBG("    sleep"<<endl);
          //usleep(1000);
      }
      DBG("   N --> "<<N<<endl);
      DBG(endl<<endl<<endl);
  }
//  DumpTraces(rank, traces);
  cout<<rank<<": done tracing. # of traces= "<<traces.size()<<endl;
  DBG("All done"<<endl);
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
    output<<"!NBLOCKS "<<numRanks<<endl;
    for (int i = 0; i < numRanks; i++)
    {
        char nm[128];
        sprintf(nm, "ds.%03d.vtk", i);
        output<<nm<<endl;
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
    output<<"!NBLOCKS "<<numRanks<<endl;
    for (int i = 0; i < totalNumDoms; i++)
    {
      char nm[128];
      sprintf(nm, "dom.%03d.vtk", i);
      output<<nm<<endl;
    }
  }
}

void
ParticleAdvection::Init()
{
#ifdef VTKH_PARALLEL
    MPI_Comm mpiComm = MPI_Comm_f2c(vtkh::GetMPICommHandle());
#endif
  int totalNumDoms = this->m_input->GetGlobalNumberOfDomains();
  int nDoms = this->m_input->GetNumberOfDomains();

  vector<double> domBounds(6*totalNumDoms, 0.0), locBounds(6*totalNumDoms, 0.0);
  cout<<"TOTAL NUM DOMAINS= "<<totalNumDoms<<" nDoms= "<<nDoms<<endl;

  for (int i = 0; i < nDoms; i++)
  {
    vtkm::Id domId;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domId);
    //dom.PrintSummary(cout);
    dataBlocks.push_back(new DataBlock(domId, dom, m_field_name, stepSize));

    vtkm::Bounds b = dom.GetCoordinateSystem().GetBounds();
    cout<<i<<" domId= "<<domId<<endl;
    cout<<i<<" bounds: ("<<b.X.Min<<" "<<b.Y.Min<<" "<<b.Z.Min<<") ("<<b.X.Max<<" "<<b.Y.Max<<" "<<b.Z.Max<<")"<<endl;

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

    output<<"X,Y,Z,ID"<<endl;
    for (auto &p : particleTraces)
        output<<p[0]<<", "<<p[1]<<", "<<p[2]<<", "<<(int)p[3]<<endl;

    if (idx == 0)
    {
        ofstream output;
        output.open("output.visit", ofstream::out);
        output<<"!NBLOCKS "<<numRanks<<endl;
        for (int i = 0; i < numRanks; i++)
        {
            char nm[128];
            sprintf(nm, "output.%03d.txt", i);
            output<<nm<<endl;
        }
    }
}

} //  namespace vtkh
