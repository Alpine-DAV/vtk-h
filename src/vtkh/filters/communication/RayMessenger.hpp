#ifndef VTKH_RAY_MESSENGER_H
#define VTKH_RAY_MESSENGER_H

#include <mpi.h>
#include <list>
#include <vector>
#include <set>
#include <map>
#include "Ray.hpp"
#include "Messenger.hpp"

// debuggin
#include <iostream>
#include <fstream>
#include <sstream>

//#define LOG_MESSAGES 1

namespace vtkh
{

class MemStream;

class RayMessenger : public Messenger
{
  using MsgCommType = std::pair<int, std::vector<int>>;

  public:
    RayMessenger(MPI_Comm comm);
    ~RayMessenger();

    void RegisterMessages(int msgSize,
                          int numMsgRecvs,
                          int nRays,
                          int nRaysRecvs);

    void SendRays(int dst, std::vector<vtkh::Ray> &rays);

    void SendRays(std::map<int, std::vector<vtkh::Ray>> &ray_map);

    bool RecvRays(std::vector<vtkh::Ray> &rays);

    // Send/Recv messages.
    void SendMsg(int dst, std::vector<int> &msg);
    void SendAllMsg(std::vector<int> &msg);
    bool RecvMsg(std::vector<MsgCommType> &msgs);

    // Send/Recv datasets.
    bool RecvAny(std::vector<MsgCommType> *msgs,
                 std::vector<vtkh::Ray> *rays,
                 bool blockAndWait);

  private:

#ifdef LOG_MESSAGES
    static int m_message_id;
    std::ofstream m_log;
#endif

    enum
    {
      MESSAGE_TAG = 0xbadbeef,
      RAY_TAG = 0xfeebdab
    };

    static int CalcRayBufferSize(int nRays);
};

} //namespace vtkh
#endif
