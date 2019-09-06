#include <iostream>
#include <string.h>
#include <vtkh/Logger.hpp>
#include <vtkh/filters/communication/MemStream.h>
#include <vtkh/filters/communication/RayMessenger.hpp>

using namespace std;
namespace vtkh
{
#ifdef LOG_MESSAGES
int RayMessenger::m_message_id = 0;
#endif

RayMessenger::RayMessenger(MPI_Comm comm)
  : Messenger(comm)
{
#ifdef LOG_MESSAGES
  std::stringstream ss;
  ss<<"msg_log_"<<this->rank<<".txt";
  m_log.open(ss.str());
#endif
}

RayMessenger::~RayMessenger()
{
#ifdef LOG_MESSAGES
  m_log.close();
#endif
}

int
RayMessenger::CalcRayBufferSize(int nRays)
{
    MemStream buff;
    int rank = 0;
    std::vector<vtkh::Ray> v(nRays);
    vtkh::write(buff, rank);
    vtkh::write(buff, v);

    return buff.len();
}

void
RayMessenger::RegisterMessages(int msgSize,
                               int numMsgRecvs,
                               int nRays,
                               int nRaysRecvs)
{
    int messageBuffSz = CalcMessageBufferSize(msgSize);
    int rayBuffSz = CalcRayBufferSize(nRays);

    this->RegisterTag(RayMessenger::MESSAGE_TAG, numMsgRecvs, messageBuffSz);
    this->RegisterTag(RayMessenger::RAY_TAG, nRaysRecvs, rayBuffSz);

    this->InitializeBuffers();
}

void
RayMessenger::SendMsg(int dst, vector<int> &msg)
{
    MemStream *buff = new MemStream;

    //Write data.
    vtkh::write(*buff, rank);
    vtkh::write(*buff, msg);
#ifdef LOG_MESSAGES
    vtkh::write(*buff, m_message_id);
    m_log<<rank<<" "<<m_message_id<<"\n";
    m_message_id++;
#endif
    SendData(dst, RayMessenger::MESSAGE_TAG, buff);
}

void
RayMessenger::SendAllMsg(vector<int> &msg)
{
    for (int i = 0; i < nProcs; i++)
        if (i != rank)
            SendMsg(i, msg);
}

bool
RayMessenger::RecvAny(vector<MsgCommType> *msgs,
                      std::vector<vtkh::Ray> *rays,
                      bool blockAndWait)
{
    set<int> tags;
    if(msgs)
    {
      tags.insert(RayMessenger::MESSAGE_TAG);
      msgs->resize(0);
    }
    if(rays)
    {
      tags.insert(RayMessenger::RAY_TAG);
      rays->resize(0);
    }

    if (tags.empty())
        return false;

    vector<pair<int, MemStream *> > buffers;
    if (! RecvData(tags, buffers, blockAndWait))
        return false;

    for (size_t i = 0; i < buffers.size(); i++)
    {
        if (buffers[i].first == RayMessenger::MESSAGE_TAG)
        {
          int sendRank;
          vector<int> m;
          vtkh::read(*buffers[i].second, sendRank);
          vtkh::read(*buffers[i].second, m);
#ifdef LOG_MESSAGES
          int message_id;
          vtkh::read(*buffers[i].second, message_id);
          m_log<<sendRank<<" "<<message_id<<"\n";
#endif
          msgs->push_back(std::make_pair(sendRank, m));
        }
        else if (buffers[i].first == RayMessenger::RAY_TAG)
        {
          int num, sendRank;
          vtkh::read(*buffers[i].second, sendRank);
          vtkh::read(*buffers[i].second, num);
#ifdef LOG_MESSAGES
          int message_id;
          vtkh::read(*buffers[i].second, message_id);
          m_log<<sendRank<<" "<<message_id<<"\n";
#endif

          rays->resize(num);
          for (int j = 0; j < num; j++)
          {
            vtkh::read(*(buffers[i].second), (*rays)[j]);
          }
          std::cout<<"["<<rank<<"] <-- ["<<sendRank<<"] "<<rays->size()<<"\n";
        }

        delete buffers[i].second;
    }

    return true;
}

bool
RayMessenger::RecvMsg(vector<MsgCommType> &msgs)
{
    return RecvAny(&msgs, NULL, false);
}

void RayMessenger::SendRays(int dst, std::vector<vtkh::Ray> &rays)
{
    if (dst == rank)
    {
        cerr<<"Error. Sending IC to yourself"<<endl;
        return;
    }
    if (rays.empty())
        return;

    MemStream *buff = new MemStream;
    vtkh::write(*buff, rank);
    const int num = rays.size();
    vtkh::write(*buff, num);

#ifdef LOG_MESSAGES
    vtkh::write(*buff, m_message_id);
    m_log<<rank<<" "<<m_message_id<<"\n";
    m_message_id++;
#endif

    for (auto &ray : rays)
    {
        //std::cout<<"writing "<<ray<<"\n";
        vtkh::write(*buff, ray);
    }
    SendData(dst, RayMessenger::RAY_TAG, buff);
    rays.clear();
}

void RayMessenger::SendRays(std::map<int, std::vector<vtkh::Ray>> &ray_map)
{
    for (auto mit = ray_map.begin(); mit != ray_map.end(); mit++)
        if (! mit->second.empty())
            SendRays(mit->first, mit->second);
}

bool RayMessenger::RecvRays(std::vector<vtkh::Ray> &rays)
{
    return RecvAny(NULL, &rays, false);
}

} // namespace vtkh
