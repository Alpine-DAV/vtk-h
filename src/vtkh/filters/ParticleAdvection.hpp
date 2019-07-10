#ifndef VTK_H_PARTICLE_ADVECTION_HPP
#define VTK_H_PARTICLE_ADVECTION_HPP

#include <list>
#include <vector>
#include <deque>
#include <map>
#include <algorithm>

#include <vtkh/vtkh.hpp>
#include <vtkh/utils/StatisticsDB.hpp>
#include <vtkh/filters/Filter.hpp>
#include <vtkh/filters/Particle.hpp>
#include <vtkh/filters/communication/BoundsMap.hpp>
#include <vtkh/filters/Integrator.hpp>
#include <vtkh/utils/StatisticsDB.hpp>
#include <vtkh/DataSet.hpp>

#ifdef VTKH_PARALLEL
  #include <mpi.h>
#endif

namespace vtkh
{
class DataBlock;

class ParticleAdvection : public Filter
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

  DataBlock * GetBlock(int blockId);

  template <typename ResultT>
  int InternalIntegrate(DataBlock &blk,
                        std::vector<Particle> &v,
                        std::list<Particle> &I,
                        std::list<Particle> &T,
                        std::list<Particle> &A,
                        vector<ResultT> &traces);

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
  std::vector<DataBlock*> dataBlocks;

  //seed data
  std::list<Particle> active, inactive, terminated;
  bool GetActiveParticles(std::vector<Particle> &v);

  void DumpTraces(int ts, const vector<vtkm::Vec<double,4>> &particleTraces);
  void DumpDS(int ts);
  void DumpSLOutput(vtkm::cont::DataSet *ds, int domId, int ts);
};


class DataBlock
{
public:
    DataBlock(int _id, vtkm::cont::DataSet *_ds, const std::string &fieldName, float advectStep)
        : id(_id), ds(_ds),
          integrator(_ds, fieldName, advectStep)
          //refCount(0), used(false)

    {
        //dbg<<"DB ctor: "<<ds.use_count()<<endl;
    }
    ~DataBlock() {}//{ds=NULL; delete ds;}//cout<<"Delete datablock id= "<<id<<" cnt= "<<ds.use_count()<<endl;}

    int id;
    vtkm::cont::DataSet *ds;
    Integrator integrator;

    /*
    void getRefUsed(int &r, bool &u) {m.lock(); r=refCount; u=used; m.unlock();}
    bool getUsed() {bool v; m.lock(); v=used; m.unlock(); return v;}
    int getRefCount() {int v; m.lock(); v=refCount; m.unlock(); return v;}

    void addUsed() {m.lock(); used=true; m.unlock();}

    void addRef() {m.lock(); refCount++; m.unlock();}
    void release() {m.lock(); refCount--;  m.unlock(); if(refCount<0) throw std::runtime_error("Negative ref cnt");}

private:
    int refCount;
    bool used;
    mutex m;
    */

    friend std::ostream &operator<<(std::ostream &os, const DataBlock &d)
    {
        os<<"DataBlock {"<<std::endl;
        os<<"  id="<<d.id<<std::endl;
        d.ds->PrintSummary(os);
        os<<"} DataBlock"<<std::endl;
        return os;
    }
};


} //namespace vtkh
#endif
