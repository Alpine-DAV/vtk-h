#include <iostream>
#include <string.h>
#include <vtkh/utils/Logger.hpp>
#include <vtkh/utils/StatisticsDB.hpp>
#include <vtkh/filters/communication/MemStream.h>
#include <vtkh/filters/communication/ParticleMessenger.hpp>

#ifdef ENABLE_LOGGING
#define DBG(msg) vtkh::Logger::GetInstance("out")->GetStream()<<msg
#define WDBG(msg) vtkh::Logger::GetInstance("wout")->GetStream()<<msg
#else
#define DBG(msg)
#define WDBG(msg)
#endif

namespace vtkh
{

    ParticleMessenger::ParticleMessenger(MPI_Comm comm, const vtkh::BoundsMap &bm)
    : Messenger(comm),
      boundsMap(bm),
      done(false)
{
    ADD_TIMER("communication");
    ADD_TIMER("gridLocator");
    ADD_COUNTER("particlesSent");
    ADD_COUNTER("earlyTerm");
    ADD_COUNTER("messagesSent");
}

void
ParticleMessenger::RegisterMessages(int mSz,
                                    int nMsgRecvs,
                                    int nICRecvs,
                                    int nDSRecvs)
{
    numMsgRecvs = nMsgRecvs;
    numSLRecvs = nICRecvs;
    numDSRecvs = nDSRecvs;

    // Msgs are handled as vector<int>.
    // Serialization of msg consists: size_t (num elements) +
    // sender rank + message size.
    int msgSize = sizeof(size_t);
    msgSize += sizeof(int); // sender rank.
    msgSize += (mSz * sizeof(int));

    //During particle advection, the IC state is only serialized.
    slSize = 256;
    slsPerRecv = 64;

    int dsSize = 2*sizeof(int);

    this->RegisterTag(ParticleMessenger::MESSAGE_TAG, numMsgRecvs, msgSize);
    this->RegisterTag(ParticleMessenger::PARTICLE_TAG, numSLRecvs, slSize * slsPerRecv);

    this->InitializeBuffers();
}

void
ParticleMessenger::SendMsg(int dst, std::vector<int> &msg)
{
    MemStream *buff = new MemStream;

    //Write data.
    vtkh::write(*buff, rank);
    vtkh::write(*buff, msg);

    SendData(dst, ParticleMessenger::MESSAGE_TAG, buff);

    COUNTER_INC("messagesSent", 1);
}

void
ParticleMessenger::SendAllMsg(std::vector<int> &msg)
{
    for (int i = 0; i < nProcs; i++)
        if (i != rank)
        {
            DBG("          ***************** SendMsg to "<<i<<" "<<msg<<endl);
            SendMsg(i, msg);
        }
}

bool
ParticleMessenger::RecvAny(std::vector<MsgCommData> *msgs,
                           list<ParticleCommData<Particle>> *recvICs,
                           std::vector<DSCommData> *ds,
                           bool blockAndWait)
{
    set<int> tags;
    if (msgs)
    {
        tags.insert(ParticleMessenger::MESSAGE_TAG);
        msgs->resize(0);
    }
    if (recvICs)
    {
        tags.insert(ParticleMessenger::PARTICLE_TAG);
        recvICs->resize(0);
    }

    if (tags.empty())
        return false;

    std::vector<pair<int, MemStream *> > buffers;
    if (! RecvData(tags, buffers, blockAndWait))
        return false;

    for (size_t i = 0; i < buffers.size(); i++)
    {
        if (buffers[i].first == ParticleMessenger::MESSAGE_TAG)
        {
            int sendRank;
            std::vector<int> m;
            vtkh::read(*buffers[i].second, sendRank);
            vtkh::read(*buffers[i].second, m);

            MsgCommData msg(sendRank, m);

            msgs->push_back(msg);
        }
        else if (buffers[i].first == ParticleMessenger::PARTICLE_TAG)
        {
            int num, sendRank;
            vtkh::read(*buffers[i].second, sendRank);
            vtkh::read(*buffers[i].second, num);

            for (int j = 0; j < num; j++)
            {
                Particle recvP;
                vtkh::read(*(buffers[i].second), recvP);
                ParticleCommData<Particle> d(sendRank, recvP);
                recvICs->push_back(d);
            }
        }

        delete buffers[i].second;
    }

    return true;
}

bool
ParticleMessenger::RecvMsg(std::vector<MsgCommData> &msgs)
{
    return RecvAny(&msgs, NULL, NULL, false);
}

template <typename P, template <typename, typename> class Container,
          typename Allocator>
void ParticleMessenger::SendICs(int dst, Container<P, Allocator> &c)
{
    if (dst == rank)
    {
        std::cerr<<"Error. Sending IC to yourself"<<std::endl;
        return;
    }
    if (c.empty())
        return;

    MemStream *buff = new MemStream;
    vtkh::write(*buff, rank);
    int num = c.size();
    vtkh::write(*buff, num);
    for (auto &p : c)
        vtkh::write(*buff, p);
    SendData(dst, ParticleMessenger::PARTICLE_TAG, buff);

    COUNTER_INC("particlesSent", num);
    c.clear();
}

template <typename P, template <typename, typename> class Container,
          typename Allocator>
void ParticleMessenger::SendICs(std::map<int, Container<P, Allocator>> &m)
{
    for (auto mit = m.begin(); mit != m.end(); mit++)
        if (! mit->second.empty())
            SendICs(mit->first, mit->second);
}

template <typename P, template <typename, typename> class Container,
          typename Allocator>
bool ParticleMessenger::RecvICs(Container<ParticleCommData<P>, Allocator> &recvICs)
{
    return RecvAny(NULL, &recvICs, NULL, false);
}

template <typename P, template <typename, typename> class Container,
          typename Allocator>
bool ParticleMessenger::RecvICs(Container<P, Allocator> &recvICs)
{
    list<ParticleCommData<P>> incoming;

    if (RecvICs(incoming))
    {
        for (auto &it : incoming)
            recvICs.push_back(it.p);
        return true;
    }
    return false;
}

template <typename P>
bool ParticleMessenger::DoSendICs(int dst, std::vector<P> &ics)
{
    if (dst == rank)
        std::cerr << "Error in ParticleMessenger::DoSendICs() Sending ICs to yourself" << std::endl;

    if (ics.empty())
        return false;

    MemStream *buff = new MemStream;
    vtkh::write(*buff,rank);
    int num = ics.size();
    vtkh::write(*buff,num);
    for (size_t i = 0; i < ics.size(); i++)
        vtkh::write(*buff,ics[i]);
    SendData(dst, ParticleMessenger::PARTICLE_TAG, buff);

    return true;
}

//Check allof the blockIds that I could be in. If this rank owns that Id, see if it is really here,
//if not, send it on.
void
ParticleMessenger::CheckAllBlocks(Particle &p,
                                  std::list<vtkh::Particle> &outData,
                                  std::list<vtkh::Particle> &inData,
                                  std::list<vtkh::Particle> &term,
                                  int &earlyTerm,
                                  std::map<int, list<Particle>> &sendData)
{
  while (!p.blockIds.empty())
  {
      int dst = boundsMap.GetRank(p.blockIds[0]);
      if (dst != rank)
      {
        sendData[dst].push_back(p);
        break;
      }

      TIMER_START("gridLocator");
      auto loc = gridLocators[p.blockIds[0]];
      vtkm::cont::ArrayHandle<vtkm::Id> cellId;
      vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault,3>> point, pcoord;
      std::vector<vtkm::Vec<vtkm::FloatDefault,3>> pt = {p.coords};
      point = vtkm::cont::make_ArrayHandle(pt);
      loc.FindCells(point, cellId, pcoord);
      auto cellPortal = cellId.GetPortalConstControl();
      TIMER_STOP("gridLocator");

      //point NOT in this domain.
      if (cellPortal.Get(0) == -1)
      {
        p.blockIds.erase(p.blockIds.begin());
        if (p.blockIds.empty())
        {
            term.push_back(p);
            p.status = vtkh::Particle::TERMINATE;
            earlyTerm++;
            DBG("-----earlyterm: "<<p<<std::endl);
            break;
        }

        int dst = boundsMap.GetRank(p.blockIds[0]);
        if(dst != rank)
        {
          sendData[dst].push_back(p);
          break;
        }
      }
      else
      {
        inData.push_back(p);
        break;
      }
  }

  /*
  while (!p.blockIds.empty())
  {
      int dst = boundsMap.GetRank(p.blockIds[0]);

      //See if this
      if (dst == rank)
  }
  */

}


void
ParticleMessenger::Exchange(std::list<vtkh::Particle> &outData,
                            std::list<vtkh::Particle> &inData,
                            std::list<vtkh::Particle> &term,
                            int &numTerminatedMessages)
{
    DBG("----ExchangeParticles: O="<<outData<<" I="<<inData<<std::endl);
    std::map<int, list<Particle>> sendData;

    TIMER_START("communication");
    int earlyTerm = 0;
    if (!outData.empty())
    {
        std::vector<std::vector<int>> blockIds;
        boundsMap.FindBlockIDs(outData, blockIds, true);
        DBG("-----O.blockIds: "<<outData<<" "<<blockIds<<endl);

        auto bit = blockIds.begin();
        for (auto lit = outData.begin(); lit != outData.end(); lit++, bit++)
        {
            //No domains. it terminated.
            if ((*lit).blockIds.empty())
            {
                term.push_back(*lit);
                (*lit).status = vtkh::Particle::TERMINATE;
                earlyTerm++;
                DBG("-----earlyterm: "<<*lit<<std::endl);
            }
            else
            {
                int dst = boundsMap.GetRank((*lit).blockIds[0]);
                if (dst == -1)
                    throw "Block ID not found";

                if (dst != rank)
                    sendData[dst].push_back(*lit);
                else
                    CheckAllBlocks((*lit), outData, inData, term, earlyTerm, sendData);
//                    inData.push_back(*lit);
            }
        }
    }

    //Check if we have anything coming in.
    std::list<ParticleCommData<Particle>> particleData;
    std::vector<MsgCommData> msgData;
    numTerminatedMessages = 0;

    DBG("-----RecvAny..."<<endl);
    if (RecvAny(&msgData, &particleData, NULL, false))
    {
        DBG("-----Recv: M: "<<msgData<<" P: "<<particleData<<std::endl);
        for (auto &p : particleData)
        {
            CheckAllBlocks(p.p, outData, inData, term, earlyTerm, sendData);
 /*           if (stats) stats->start("gridLocator");
            auto loc = gridLocators[p.p.blockIds[0]];
            vtkm::cont::ArrayHandle<vtkm::Id> cellId;
            vtkm::cont::ArrayHandle<vtkm::Vec<double,3>> point, pcoord;
            std::vector<vtkm::Vec<double,3>> pt = {p.p.coords};
            point = vtkm::cont::make_ArrayHandle(pt);
            loc.FindCells(point, cellId, pcoord);
            auto cellPortal = cellId.GetPortalConstControl();
            if (stats) stats->stop("gridLocator");

            //point NOT in this domain.
            if (cellPortal.Get(0) == -1)
            {
                if (p.p.blockIds.size() == 1)
                {
                    p.p.blockIds.erase(p.p.blockIds.begin());
                    p.p.status = vtkh::Particle::TERMINATE;
                    term.push_back(p.p);
                    earlyTerm++;
                }
                else
                {
                    p.p.blockIds.erase(p.p.blockIds.begin());
                    int dst = boundsMap.GetRank(p.p.blockIds[0]);
                    if(dst != rank)
                      sendData[dst].push_back(p.p);
                    else
                      inData.push_back(p.p);
                }
            }
            else
                inData.push_back(p.p);
*/        }

        for (auto &m : msgData)
        {
            if (m.message[0] == MSG_TERMINATE)
            {
                numTerminatedMessages += m.message[1];
                DBG("-----TERMinate: Recv: "<<m.message[1]<<endl);
            }
            else if (m.message[0] == MSG_DONE)
            {
                DBG("-----DONE RECEIVED: "<<m.message[1]<<endl);
                done = true;
            }
        }
    }
    else
    {
        DBG("-----RecvAny --Nothing in the can"<<endl);
    }

    COUNTER_INC("earlyTerm", earlyTerm);

    //Do all the sending...
    if (!term.empty())
    {
        std::vector<int> msg = {MSG_TERMINATE, (int)term.size()};
        DBG("-----SendAllMsg: msg="<<msg<<endl);
        SendAllMsg(msg);
    }
    DBG("-----SendP: "<<sendData<<std::endl);
    if (!sendData.empty())
    {
        for (auto &i : sendData)
            SendICs(i.first, i.second);
    }

    CheckPendingSendRequests();
    DBG("----ExchangeParticles Done: I= "<<inData<<" T= "<<term<<std::endl<<std::endl);
    TIMER_STOP("communication");
}

}
