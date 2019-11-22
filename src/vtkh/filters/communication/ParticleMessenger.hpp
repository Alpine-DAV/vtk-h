#ifndef VTKH_PARTICLE_MESSENGER_H
#define VTKH_PARTICLE_MESSENGER_H

#include <mpi.h>
#include <list>
#include <vector>
#include <set>
#include <map>
#include <vtkh/vtkh_exports.h>
#include <vtkh/filters/Particle.hpp>
#include <vtkh/filters/communication/Messenger.hpp>
#include <vtkh/filters/communication/BoundsMap.hpp>
#include <vtkm/cont/CellLocatorUniformBins.h>

namespace vtkh
{

class MemStream;

class VTKH_API ParticleMessenger : public Messenger
{
    const int MSG_TERMINATE = 1;
    const int MSG_DONE = 1;

    using MsgCommType = std::pair<int, std::vector<int>>;
    using ParticleCommType = std::pair<int, std::vector<vtkh::Particle>>;

  public:
    ParticleMessenger(MPI_Comm comm, const vtkh::BoundsMap &bm);
    ~ParticleMessenger() {}

    void RegisterMessages(int msgSz,
                          int nMsgRecvs,
                          int nParticles,
                          int nParticlesRecvs);

    void Exchange(std::vector<vtkh::Particle> &outData,
                  std::vector<vtkh::Particle> &inData,
                  std::vector<vtkh::Particle> &term,
                  int &numTerminateMessages);

    // Send/Recv Integral curves.
    template <typename P, template <typename, typename> class Container,
              typename Allocator=std::allocator<P>>
    void SendParticles(int dst, const Container<P, Allocator> &c);

    template <typename P, template <typename, typename> class Container,
              typename Allocator=std::allocator<P>>
    void SendParticles(const std::map<int, Container<P, Allocator>> &m);

    // Send/Recv messages.
    void SendMsg(int dst, const std::vector<int> &msg);
    void SendAllMsg(const std::vector<int> &msg);
    bool RecvMsg(std::vector<MsgCommType> &msgs)
    {
      return RecvAny(&msgs, NULL, false);
    }

    // Send/Recv datasets.
    bool RecvAny(std::vector<MsgCommType> *msgs,
                 std::vector<ParticleCommType> *recvParticles,
                 bool blockAndWait);

  private:
    bool done;
    vtkh::BoundsMap boundsMap;

    void
    ParticleSorter(std::vector<vtkh::Particle> &outData,
                   std::vector<vtkh::Particle> &inData,
                   std::vector<vtkh::Particle> &term,
                   std::map<int, std::vector<Particle>> &sendData);

    enum
    {
        MESSAGE_TAG = 0x42000,
        PARTICLE_TAG = 0x42001
    };

    static int CalcParticleBufferSize(int nParticles, int numBlockIds=2);
};
} //namespace vtkh
#endif
