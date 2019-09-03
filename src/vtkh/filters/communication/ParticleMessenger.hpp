#ifndef VTKH_PARTICLE_MESSENGER_H
#define VTKH_PARTICLE_MESSENGER_H

#include <mpi.h>
#include <list>
#include <vector>
#include <set>
#include <map>
#include "CommData.hpp"
#include <vtkh/filters/Particle.hpp>
#include <vtkh/filters/communication/Messenger.hpp>
#include <vtkh/filters/communication/BoundsMap.hpp>
#include <vtkm/cont/CellLocatorUniformBins.h>

namespace vtkh
{

class MemStream;

class ParticleMessenger : public Messenger
{
    const int MSG_TERMINATE = 1;
    const int MSG_DONE = 1;

  public:
    ParticleMessenger(MPI_Comm comm, const vtkh::BoundsMap &bm);
    ~ParticleMessenger() {}

    void RegisterMessages(int msgSz,
                          int nMsgRecvs,
                          int nParticles,
                          int nParticlesRecvs);

    void Exchange(std::list<vtkh::Particle> &outData,
                  std::list<vtkh::Particle> &inData,
                  std::list<vtkh::Particle> &term,
                  int &numTerminateMessages);

    // Send/Recv Integral curves.
    template <typename P, template <typename, typename> class Container,
              typename Allocator=std::allocator<P>>
    void SendParticles(int dst, const Container<P, Allocator> &c);

    template <typename P, template <typename, typename> class Container,
              typename Allocator=std::allocator<P>>
    void SendParticles(const std::map<int, Container<P, Allocator>> &m);

    template <typename P, template <typename, typename> class Container,
              typename Allocator=std::allocator<P>>
    bool RecvParticles(Container<P, Allocator> &recvICs);

    // Send/Recv messages.
    void SendMsg(int dst, const std::vector<int> &msg);
    void SendAllMsg(const std::vector<int> &msg);
    bool RecvMsg(std::vector<MsgCommData> &msgs);

    // Send/Recv datasets.
    bool RecvAny(std::vector<MsgCommData> *msgs,
                 std::list<ParticleCommData<Particle>> *recvParticles,
                 std::vector<DSCommData> *ds,
                 bool blockAndWait);

  void AddLocator(int domain, vtkm::cont::DataSet &ds)
  {
      vtkm::cont::CellLocatorUniformBins locator;
      locator.SetCoordinates(ds.GetCoordinateSystem());
      locator.SetCellSet(ds.GetCellSet());
      gridLocators.insert(std::pair<int,vtkm::cont::CellLocatorUniformBins>(domain, locator));
  }

  private:
    bool done;
    vtkh::BoundsMap boundsMap;

    template <typename P, template <typename, typename> class Container,
              typename Allocator=std::allocator<P>>
    bool RecvParticles(Container<ParticleCommData<P>, Allocator> &recvICs);

    void
    ParticleSorter(std::list<vtkh::Particle> &outData,
                   std::list<vtkh::Particle> &inData,
                   std::list<vtkh::Particle> &term,
                   std::map<int, std::list<Particle>> &sendData);

    enum
    {
        MESSAGE_TAG = 0x42000,
        PARTICLE_TAG = 0x42001
    };

    static int CalcParticleBufferSize(int nParticles, int numBlockIds=2);

    std::map<int, vtkm::cont::CellLocatorUniformBins> gridLocators;
};
} //namespace vtkh
#endif
