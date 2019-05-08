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

// debuggin
#include <iostream>
#include <fstream>
#include <sstream>

#define LOG_MESSAGES 1

namespace vtkh
{

class MemStream;

class RayMessenger : public Messenger
{
  public:
    RayMessenger(MPI_Comm comm);
    ~RayMessenger();

    void RegisterMessages(int msgSize,
                          int numMsgRecvs,
                          int numICRecvs);

    void SendRays(std::map<int, std::vector<Ray>> &ray_map);
    void SendRays(int dst, std::vector<Ray> &rays);
    bool RecvRays(std::vector<Ray> &rays);

    void SendResults(int dst, std::vector<RayResult> &results);
    bool RecvResults(std::vector<RayResult> &results);

    // Send/Recv messages.
    void SendMsg(int dst, std::vector<int> &msg);
    void SendAllMsg(std::vector<int> &msg);
    bool RecvMsg(std::vector<MsgCommData> &msgs);

    bool RecvAny(std::vector<MsgCommData> *msgs,
                 std::vector<Ray> *rays,
                 std::vector<RayResult> *results,
                 bool blockAndWait);

  private:

#ifdef LOG_MESSAGES
    static int m_message_id;
    std::ofstream m_log;
#endif

    enum
    {
      MESSAGE_TAG = 0xbadbeef,
      RAY_TAG = 0xfeebdab,
      RESULT_TAG = 0xbadbad
    };

    //Message headers.
    typedef struct
    {
        int rank, id, tag, numPackets, packet, packetSz, dataSz;
    } Header;
};

} //namespace vtkh
#endif
