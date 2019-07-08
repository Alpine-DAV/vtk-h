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
#include <vtkh/utils/StatisticsDB.hpp>
#include <vtkm/worklet/particleadvection/GridEvaluators.h>
#include <vtkm/cont/CellLocatorTwoLevelUniformGrid.h>

namespace vtkh
{

class MemStream;

class ParticleMessenger : public Messenger
{
    const int MSG_TERMINATE = 1;
    const int MSG_DONE = 1;

  public:
    ParticleMessenger(MPI_Comm comm, const vtkh::BoundsMap &bm, vtkh::StatisticsDB *pSDB=NULL);
    ~ParticleMessenger() {}

    void RegisterMessages(int msgSize,
                          int numMsgRecvs,
                          int numICRecvs,
                          int numDSRecvs=0);

    int Exchange(std::list<vtkh::Particle> &outData,
                 std::list<vtkh::Particle> &inData,
                 std::list<vtkh::Particle> &term);

    // Send/Recv Integral curves.
    template <typename P, template <typename, typename> class Container,
              typename Allocator=std::allocator<P>>
    void SendICs(int dst, Container<P, Allocator> &c);

    template <typename P, template <typename, typename> class Container,
              typename Allocator=std::allocator<P>>
    void SendICs(std::map<int, Container<P, Allocator>> &m);

    template <typename P, template <typename, typename> class Container,
              typename Allocator=std::allocator<P>>
    bool RecvICs(Container<P, Allocator> &recvICs);

    template <typename P, template <typename, typename> class Container,
              typename Allocator=std::allocator<P>>
    bool RecvICs(Container<ParticleCommData<P>, Allocator> &recvICs);

    // Send/Recv messages.
    void SendMsg(int dst, std::vector<int> &msg);
    void SendAllMsg(std::vector<int> &msg);
    bool RecvMsg(std::vector<MsgCommData> &msgs);

    // Send/Recv datasets.
    bool RecvAny(std::vector<MsgCommData> *msgs,
                 std::list<ParticleCommData<Particle>> *recvICs,
                 std::vector<DSCommData> *ds,
                 bool blockAndWait);

  void AddLocator(int domain, vtkm::cont::DataSet &ds)
  {
      vtkm::cont::CellLocatorTwoLevelUniformGrid locator;
      locator.SetCoordinates(ds.GetCoordinateSystem());
      locator.SetCellSet(ds.GetCellSet());
      locator.Build();
      gridLocators.insert(std::pair<int,vtkm::cont::CellLocatorTwoLevelUniformGrid>(domain, locator));
  }

  private:
    template <typename P>
    bool DoSendICs(int dst, std::vector<P> &ics);
    bool done;

    vtkh::BoundsMap boundsMap;
    vtkh::StatisticsDB *stats;

    void
    CheckAllBlocks(Particle &p,
                  std::list<vtkh::Particle> &outData,
                  std::list<vtkh::Particle> &inData,
                  std::list<vtkh::Particle> &term,
                  int *earlyTerm,
                  map<int, list<Particle>> &sendData);

    enum
    {
        MESSAGE_TAG = 0x42000,
        PARTICLE_TAG = 0x42001
    };

    //Message headers.
    typedef struct
    {
        int rank, id, tag, numPackets, packet, packetSz, dataSz;
    } Header;

    std::map<int, vtkm::cont::CellLocatorTwoLevelUniformGrid> gridLocators;
};
} //namespace vtkh
#endif
