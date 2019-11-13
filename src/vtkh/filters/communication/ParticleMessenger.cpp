#include <iostream>
#include <string.h>
#include <vtkh/Logger.hpp>
#include <vtkh/StatisticsDB.hpp>
#include <vtkh/filters/communication/MemStream.h>
#include <vtkh/filters/communication/ParticleMessenger.hpp>

#ifdef VTKH_ENABLE_LOGGING
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
ParticleMessenger::RecvAny(std::vector<MsgCommType> *msgs,
                           std::vector<ParticleCommType> *recvParticles,
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

            msgs->push_back(std::make_pair(sendRank, m));
        }
        else if (buffers[i].first == ParticleMessenger::PARTICLE_TAG)
        {
            int sendRank;
            std::size_t num;
            vtkh::read(*buffers[i].second, sendRank);
            vtkh::read(*buffers[i].second, num);
            if (num > 0)
            {
                std::vector<vtkh::Particle> particles(num);
                for (int j = 0; j < num; j++)
                    vtkh::read(*(buffers[i].second), particles[j]);
                recvParticles->push_back(std::make_pair(sendRank, particles));
            }
        }

        delete buffers[i].second;
    }

    return true;
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

void
ParticleMessenger::ParticleSorter(std::vector<vtkh::Particle> &outData,
                                  std::vector<vtkh::Particle> &inData,
                                  std::vector<vtkh::Particle> &term,
                                  std::map<int, std::vector<Particle>> &sendData)
{
    for (auto &p : outData)
    {
        //If particle in wrong domain (took no steps), remove this ID, and use what is left below.
        //otherwise, compute a new set of block ids.
        if (p.status == vtkh::Particle::WRONG_DOMAIN)
            p.blockIds.erase(p.blockIds.begin());
        else
            p.blockIds = boundsMap.FindBlock(p, true);

        //No blocks, it terminated
        if (p.blockIds.empty())
        {
            p.status = vtkh::Particle::TERMINATE;
            term.push_back(p);
        }
        else
        {
            //If we have more than blockId, we want to minimize communication
            //and put any blocks owned by this rank first.
            if (p.blockIds.size() > 1)
            {
                auto iter = p.blockIds.begin();
                for (auto iter = p.blockIds.begin(); iter != p.blockIds.end(); iter++)
                {
                    int dst = boundsMap.GetRank(*iter);
                    if (dst == rank)
                    {
                        int bid = *iter;
                        p.blockIds.erase(iter);
                        p.blockIds.insert(p.blockIds.begin(), bid);
                        break;
                    }
                }
            }

            //Particle is mine, or put it in the sendData.
            if (p.blockIds[0] == rank)
                inData.push_back(p);
            else
                sendData[boundsMap.GetRank(p.blockIds[0])].push_back(p);
        }
    }
    outData.clear();
}

void
ParticleMessenger::Exchange(std::vector<vtkh::Particle> &outData,
                            std::vector<vtkh::Particle> &inData,
                            std::vector<vtkh::Particle> &term,
                            int &numTerminatedMessages)
{
  DBG("----ExchangeParticles: O="<<outData<<" I="<<inData<<std::endl);
  std::map<int, std::vector<Particle>> sendData;

  TIMER_START("communication");

  if (!outData.empty())
    ParticleSorter(outData, inData, term, sendData);

  //Check if we have anything coming in.
  std::vector<ParticleCommType> particleData;
  std::vector<MsgCommType> msgData;
  numTerminatedMessages = 0;

  if (RecvAny(&msgData, &particleData, false))
  {
    DBG("-----Recv: M: "<<msgData<<" P: "<<particleData<<std::endl);
    for (auto &p : particleData)
        inData.insert(inData.end(), p.second.begin(), p.second.end());

    for (auto &m : msgData)
    {
      if (m.second[0] == MSG_TERMINATE)
      {
          numTerminatedMessages += m.second[1];
          DBG("-----TERMinate: Recv: "<<m.second[1]<<std::endl);
      }
      else if (m.second[0] == MSG_DONE)
      {
        DBG("-----DONE RECEIVED: "<<m.second[1]<<std::endl);
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
