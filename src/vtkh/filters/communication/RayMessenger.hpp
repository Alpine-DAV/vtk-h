#ifndef VTKH_RAY_MESSENGER_H
#define VTKH_RAY_MESSENGER_H

#include <mpi.h>
#include <list>
#include <vector>
#include <set>
#include <map>
#include "CommData.hpp"
#include "Ray.hpp"
#include "Messenger.hpp"

namespace vtkh
{

class MemStream;

class RayMessenger : public Messenger
{
  public:
    RayMessenger(MPI_Comm comm);
    ~RayMessenger() {}

    void RegisterMessages(int msgSize,
                          int numMsgRecvs,
                          int numICRecvs,
                          int numDSRecvs=0);

    void SendRays(int dst, std::vector<Ray> &rays);

    void SendRays(std::map<int, std::vector<Ray>> &ray_map);

    bool RecvRays(std::vector<Ray> &rays);

    // Send/Recv messages.
    void SendMsg(int dst, std::vector<int> &msg);
    void SendAllMsg(std::vector<int> &msg);
    bool RecvMsg(std::vector<MsgCommData> &msgs);

    // Send/Recv datasets.
    bool RecvAny(std::vector<MsgCommData> *msgs,
                 std::vector<Ray> *rays,
                 bool blockAndWait);

  private:
    enum
    {
      MESSAGE_TAG = 0xbadbeef,
      RAY_TAG = 0xfeebdab
    };

    //Message headers.
    typedef struct
    {
        int rank, id, tag, numPackets, packet, packetSz, dataSz;
    } Header;
};

} //namespace vtkh
#endif
