#ifndef VTKH_PAR_IC_ALGORITHM_H
#define VTKH_PAR_IC_ALGORITHM_H

#include <mpi.h>
#include <list>
#include <vector>
#include <set>
#include <map>
#include "CommData.hpp"
#include "../Particle.hpp"
#include "Messenger.hpp"
#include "BoundsMap.hpp"

namespace vtkh
{

class MemStream;

class avtParICAlgorithm : public Messenger
{
    const int MSG_TERMINATE = 1;
    const int MSG_DONE = 1;

  public:
    avtParICAlgorithm(MPI_Comm comm);
    ~avtParICAlgorithm() {}

    void RegisterMessages(int msgSize,
                          int numMsgRecvs,
                          int numICRecvs,
                          int numDSRecvs=0);

    size_t Exchange(bool haveWork,
                    list<Particle> &outData,
                    list<Particle> &inData,
                    list<Particle> &term,
                    vtkh::BoundsMap &boundsMap,
                    int numTerm);

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

  private:
    template <typename P>
    bool DoSendICs(int dst, std::vector<P> &ics);

    bool done;

    enum
    {
        MESSAGE_TAG = 0xbadbeef,
        PARTICLE_TAG = 0xfeebdab
    };

    //Message headers.
    typedef struct
    {
        int rank, id, tag, numPackets, packet, packetSz, dataSz;
    } Header;
};
} //namespace vtkh
#endif
