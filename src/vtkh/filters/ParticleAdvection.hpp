#ifndef VTK_H_PARTICLE_ADVECTION_HPP
#define VTK_H_PARTICLE_ADVECTION_HPP

#include <list>
#include <vector>
#include <deque>
#include <map>
#include <algorithm>

#include <vtkh/vtkh_exports.h>
#include <vtkh/vtkh.hpp>
#include <vtkh/StatisticsDB.hpp>
#include <vtkh/filters/Filter.hpp>
#include <vtkh/filters/Particle.hpp>
#include <vtkh/filters/communication/BoundsMap.hpp>
#include <vtkh/filters/Integrator.hpp>
#include <vtkh/DataSet.hpp>

#ifdef VTKH_PARALLEL
#include <mpi.h>
#endif

namespace vtkh
{
class DataBlockIntegrator;

class VTKH_API ParticleAdvection : public Filter
{
public:
  enum SeedMethod {RANDOM=0, RANDOM_BLOCK, RANDOM_BOX, POINT};

  ParticleAdvection();
  virtual ~ParticleAdvection();
  std::string GetName() const override;

  DataSet * GetInput() const {return m_input;}

  void SetSeedPoint(const vtkm::Vec<double,3> &pt)
  {
    seedMethod = POINT;
    seedPoint = pt;
  }
  void SetSeedsRandomWhole(const int &n)
  {
    seedMethod = RANDOM;
    numSeeds = n;
  }
  void SetSeedsRandomBlock(const int &n)
  {
    seedMethod = RANDOM_BLOCK;
    numSeeds = n;
  }
  void SetSeedsRandomBox(const int &n, const vtkm::Bounds &box)
  {
    seedMethod = RANDOM_BOX;
    numSeeds = n;
    seedBox = box;
  }

  void SetUseThreadedVersion(bool useThreaded)
  {
    useThreadedVersion = useThreaded;
  }

  void SetGatherTraces(bool gTraces)
  {
    gatherTraces = gTraces;
  }

  void SetDumpOutputFiles(bool dumpOutput)
  {
    dumpOutputFiles = dumpOutput;
  }

  void SetField(const std::string &field_name) {m_field_name = field_name;}
  void SetStepSize(const double &v) { stepSize = v;}
  void SetMaxSteps(const int &n) { maxSteps = n;}
  int  GetMaxSteps() const { return maxSteps; }

  DataBlockIntegrator * GetBlock(int blockId);

  template <typename ResultT>
  int InternalIntegrate(DataBlockIntegrator &blk,
                        std::vector<Particle> &v,
                        std::vector<Particle> &I,
                        std::vector<Particle> &T,
                        std::vector<Particle> &A,
                        std::vector<ResultT> &traces);

protected:
  void PreExecute() override;
  void PostExecute() override;
  void DoExecute() override;

  void Init();
  void CreateSeeds();

  template <typename ResultT>
  void TraceSeeds(std::vector<ResultT> &traces);
  template <typename ResultT>
  void TraceMultiThread(std::vector<ResultT> &traces);
  template <typename ResultT>
  void TraceSingleThread(std::vector<ResultT> &traces);

  int DomainToRank(int blockId) {return boundsMap.GetRank(blockId);}
  void BoxOfSeeds(const vtkm::Bounds &box,
                  std::vector<Particle> &seeds,
                  vtkm::Id domId=-1,
                  bool shrink=true);

  bool useThreadedVersion;
  bool gatherTraces;
  bool dumpOutputFiles;
  int sleepUS;
  int rank, numRanks;
  std::string m_field_name;
  int numSeeds, totalNumSeeds;
  int randSeed;
  int seedMethod;
  int maxSteps;
  vtkm::Bounds seedBox;
  vtkm::Vec<double,3> seedPoint;

  float stepSize;

  BoundsMap boundsMap;
  std::vector<DataBlockIntegrator*> dataBlocks;

  //seed data
  std::vector<Particle> active, inactive, terminated;
  bool GetActiveParticles(std::vector<Particle> &v);

  void DumpTraces(int ts, const std::vector<vtkm::Vec<double,4>> &particleTraces);
  void DumpDS(int ts);
  void DumpSLOutput(vtkm::cont::DataSet *ds, int domId, int ts);
};


class DataBlockIntegrator
{
public:
    DataBlockIntegrator(int _id, vtkm::cont::DataSet *_ds, const std::string &fieldName, float advectStep)
        : id(_id), ds(_ds),
          integrator(_ds, fieldName, advectStep)
    {
    }
    ~DataBlockIntegrator() {}

    int id;
    vtkm::cont::DataSet *ds;
    Integrator integrator;

    friend std::ostream &operator<<(std::ostream &os, const DataBlockIntegrator &d)
    {
        os<<"DataBlockIntegrator {"<<std::endl;
        os<<"  id="<<d.id<<std::endl;
        d.ds->PrintSummary(os);
        os<<"} DataBlockIntegrator"<<std::endl;
        return os;
    }
};


} //namespace vtkh
#endif
