#include <iostream>
#include <string.h>
#include <vtkh/Logger.hpp>
#include <vtkh/StatisticsDB.hpp>
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

int
ParticleMessenger::CalcParticleBufferSize(int nParticles, int numBlockIds)
{
    MemStream buff;
    int rank = 0;

    //Make a vector of particles where each particle has 'numBlockIds' in the blockId array.
    std::vector<vtkh::Particle> v(nParticles);
    vtkh::Particle p;
    p.blockIds.resize(numBlockIds);
    for (int i = 0; i < nParticles; i++)
        v[i] = p;

    vtkh::write(buff, rank);
    vtkh::write(buff, v);

    return buff.len();
}

void
ParticleMessenger::RegisterMessages(int msgSz,
                                    int nMsgRecvs,
                                    int nParticles,
                                    int nParticlesRecvs)
{
    //Determine buffer size for msg and particle tags.
    int messageBuffSz = CalcMessageBufferSize(msgSz);
    int particleBuffSz = CalcParticleBufferSize(nParticles);

    this->RegisterTag(ParticleMessenger::MESSAGE_TAG, nMsgRecvs, messageBuffSz);
    this->RegisterTag(ParticleMessenger::PARTICLE_TAG, nParticlesRecvs, particleBuffSz);

    this->InitializeBuffers();
}

void
ParticleMessenger::SendMsg(int dst, const std::vector<int> &msg)
{
    MemStream *buff = new MemStream();

    //Write data.
    vtkh::write(*buff, rank);
    vtkh::write(*buff, msg);
    SendData(dst, ParticleMessenger::MESSAGE_TAG, buff);

    COUNTER_INC("messagesSent", 1);
}

void
ParticleMessenger::SendAllMsg(const std::vector<int> &msg)
{
    for (int i = 0; i < nProcs; i++)
      if (i != rank)
        SendMsg(i, msg);
}

bool
ParticleMessenger::RecvAny(std::vector<MsgCommData> *msgs,
                           std::list<ParticleCommData<Particle>> *recvParticles,
                           std::vector<DSCommData> *ds,
                           bool blockAndWait)
{
    std::set<int> tags;
    if (msgs)
    {
        tags.insert(ParticleMessenger::MESSAGE_TAG);
        msgs->resize(0);
    }
    if (recvParticles)
    {
        tags.insert(ParticleMessenger::PARTICLE_TAG);
        recvParticles->resize(0);
    }

    if (tags.empty())
        return false;

    std::vector<std::pair<int, MemStream *> > buffers;
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
                recvParticles->push_back(d);
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
void ParticleMessenger::SendParticles(int dst, const Container<P, Allocator> &c)
{
    if (dst == rank)
    {
        std::cerr<<"Error. Sending IC to yourself"<<std::endl;
        return;
    }
    if (c.empty())
        return;

    MemStream *buff = new MemStream();
    vtkh::write(*buff, rank);
    vtkh::write(*buff, c);
    SendData(dst, ParticleMessenger::PARTICLE_TAG, buff);

    COUNTER_INC("particlesSent", c.size());
}

template <typename P, template <typename, typename> class Container,
          typename Allocator>
void ParticleMessenger::SendParticles(const std::map<int, Container<P, Allocator>> &m)
{
    for (auto mit = m.begin(); mit != m.end(); mit++)
        if (! mit->second.empty())
            SendParticles(mit->first, mit->second);
}

template <typename P, template <typename, typename> class Container,
          typename Allocator>
bool ParticleMessenger::RecvParticles(Container<ParticleCommData<P>, Allocator> &recvParticles)
{
  return RecvAny(NULL, &recvParticles, NULL, false);
}

template <typename P, template <typename, typename> class Container,
          typename Allocator>
bool ParticleMessenger::RecvParticles(Container<P, Allocator> &recvParticles)
{
    std::list<ParticleCommData<P>> incoming;

    if (RecvParticles(incoming))
    {
        for (auto &it : incoming)
            recvParticles.push_back(it.p);
        return true;
    }
    return false;
}

//Check allof the blockIds that I could be in. If this rank owns that Id, see if it is really here,
//if not, send it on.
void
ParticleMessenger::ParticleBlockSorter(Particle &p,
                                       std::list<vtkh::Particle> &inData,
                                       std::list<vtkh::Particle> &term,
                                       std::map<int, std::list<Particle>> &sendData)
{
  std::vector<int> bids;
  //Examine each blockID for inclusion in MY blocks.
  for (auto &bid : p.blockIds)
  {
    int dst = boundsMap.GetRank(bid);
    //One of my data blocks. See if it's really in there.
    if (dst == rank)
    {
        //Dave: Fix this.
#if 0
      TIMER_START("gridLocator");
      auto loc = gridLocators[p.blockIds[0]];
      vtkm::cont::ArrayHandle<vtkm::Id> cellId;
      vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::FloatDefault,3>> point, pcoord;
      std::vector<vtkm::Vec<vtkm::FloatDefault,3>> pt = {p.coords};
      point = vtkm::cont::make_ArrayHandle(pt);
      loc.FindCells(point, cellId, pcoord);
      auto cellPortal = cellId.GetPortalConstControl();
      TIMER_STOP("gridLocator");

      //point in this domain. Set the domain and return.
      if (cellPortal.Get(0) != -1)
      {
#endif
          p.blockIds = {bid};
          inData.push_back(p);
          return;
#if 0
      }
#endif
    }
    else
      bids.push_back(bid);
  }

  p.blockIds = bids;

  //No blocks remain. Particle terminated.
  if (p.blockIds.empty())
  {
      p.status = vtkh::Particle::TERMINATE;
      term.push_back(p);
  }
  //Otherwise, not mine. Pass it along.
  else
      sendData[boundsMap.GetRank(p.blockIds[0])].push_back(p);
}

void
ParticleMessenger::Exchange(std::list<vtkh::Particle> &outData,
                            std::list<vtkh::Particle> &inData,
                            std::list<vtkh::Particle> &term,
                            int &numTerminatedMessages)
{
  DBG("----ExchangeParticles: O="<<outData<<" I="<<inData<<std::endl);
  std::map<int, std::list<Particle>> sendData;

  TIMER_START("communication");
  if (!outData.empty())
  {
    std::vector<std::vector<int>> blockIds;
    boundsMap.FindBlockIDs(outData, blockIds, true);
    DBG("-----O.blockIds: "<<outData<<" "<<blockIds<<endl);

    auto bit = blockIds.begin();
    for (auto lit = outData.begin(); lit != outData.end(); lit++, bit++)
      ParticleBlockSorter((*lit), inData, term, sendData);
  }

  //Check if we have anything coming in.
  std::list<ParticleCommData<Particle>> particleData;
  std::vector<MsgCommData> msgData;
  numTerminatedMessages = 0;

  if (RecvAny(&msgData, &particleData, NULL, false))
  {
    DBG("-----Recv: M: "<<msgData<<" P: "<<particleData<<std::endl);
    for (auto &p : particleData)
      ParticleBlockSorter(p.p, inData, term, sendData);

    for (auto &m : msgData)
    {
      if (m.message[0] == MSG_TERMINATE)
      {
          numTerminatedMessages += m.message[1];
          DBG("-----TERMinate: Recv: "<<m.message[1]<<std::endl);
      }
      else if (m.message[0] == MSG_DONE)
      {
        DBG("-----DONE RECEIVED: "<<m.message[1]<<std::endl);
        done = true;
      }
    }
  }

  //Do all the sending...
  if (!term.empty())
  {
    std::vector<int> msg = {MSG_TERMINATE, (int)term.size()};
    DBG("-----SendAllMsg: msg="<<msg<<std::endl);
    SendAllMsg(msg);
  }
  if (!sendData.empty())
  {
    for (auto &i : sendData)
      SendParticles(i.first, i.second);
    sendData.clear();
  }

  CheckPendingSendRequests();
  DBG("----ExchangeParticles Done: I= "<<inData<<" T= "<<term<<std::endl<<std::endl);
  TIMER_STOP("communication");
}

}
