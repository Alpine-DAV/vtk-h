#include <iostream>
#include <string.h>
#include "MemStream.h"
#include "DebugMeowMeow.hpp"
#include "ParticleMessenger.hpp"

using namespace std;
namespace vtkh
{

ParticleMessenger::ParticleMessenger(MPI_Comm comm, const vtkh::BoundsMap &bm, vtkh::StatisticsDB *pSDB)
    : Messenger(comm),
      boundsMap(bm),
      done(false),
      stats(pSDB)
{
    if (stats)
    {
        stats->addTimer("communication");
        stats->addCounter("particlesSent");
        stats->addCounter("messagesSent");
    }
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
ParticleMessenger::SendMsg(int dst, vector<int> &msg)
{
    MemStream *buff = new MemStream;

    //Write data.
    vtkh::write(*buff, rank);
    vtkh::write(*buff, msg);

    SendData(dst, ParticleMessenger::MESSAGE_TAG, buff);

    if (stats) stats->increment("messagesSent");
}

void
ParticleMessenger::SendAllMsg(vector<int> &msg)
{
    for (int i = 0; i < nProcs; i++)
        if (i != rank)
        {
            DBG("          ***************** SendMsg to "<<i<<" "<<msg<<endl);
            SendMsg(i, msg);
        }
}

bool
ParticleMessenger::RecvAny(vector<MsgCommData> *msgs,
                           list<ParticleCommData<Particle>> *recvICs,
                           vector<DSCommData> *ds,
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

    vector<pair<int, MemStream *> > buffers;
    if (! RecvData(tags, buffers, blockAndWait))
        return false;

//    int timerHandle = visitTimer->StartTimer();

    for (size_t i = 0; i < buffers.size(); i++)
    {
        if (buffers[i].first == ParticleMessenger::MESSAGE_TAG)
        {
            int sendRank;
            vector<int> m;
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

//    CommTime.value += visitTimer->StopTimer(timerHandle, "RecvAny");
    return true;
}

bool
ParticleMessenger::RecvMsg(vector<MsgCommData> &msgs)
{
    return RecvAny(&msgs, NULL, NULL, false);
}

template <typename P, template <typename, typename> class Container,
          typename Allocator>
void ParticleMessenger::SendICs(int dst, Container<P, Allocator> &c)
{
    if (dst == rank)
    {
        cerr<<"Error. Sending IC to yourself"<<endl;
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

    if (stats) stats->increment("particlesSent", num);
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
bool ParticleMessenger::DoSendICs(int dst, vector<P> &ics)
{
    if (dst == rank)
    {
        cerr << "Error in ParticleMessenger::DoSendICs() Sending ICs to yourself" << endl;
        for (size_t i = 0; i < ics.size(); i++)
            cerr << "Proc " << rank << "  "<<ics[i]<<endl;
    }

    if (ics.empty())
        return false;

    MemStream *buff = new MemStream;
    vtkh::write(*buff,rank);
    int num = ics.size();
    vtkh::write(*buff,num);
    //cout<<" ********************* DOSENDICS: "<<ics[0]<<endl;
    for (size_t i = 0; i < ics.size(); i++)
        vtkh::write(*buff,ics[i]);
    SendData(dst, ParticleMessenger::PARTICLE_TAG, buff);

    return true;
}

int
ParticleMessenger::Exchange(list<Particle> &outData,
                             list<Particle> &inData,
                             list<Particle> &term,
                             int increment)
{
    DBG("----ExchangeParticles: O="<<outData<<" I="<<inData<<std::endl);
    map<int, list<Particle>> sendData;

    if (stats) stats->start("communication");
    int earlyTerm = 0;
    if (!outData.empty())
    {
        vector<int> blockIds;
        boundsMap.FindBlockIDs(outData, blockIds);
        DBG("-----O.blockIds: "<<outData<<" "<<blockIds<<endl);

        auto bit = blockIds.begin();
        for (auto lit = outData.begin(); lit != outData.end(); lit++, bit++)
        {
            int id = *bit;
            lit->blockId = id;
            if (id == -1)
            {
                term.push_back(*lit);
                earlyTerm++;
                DBG("-----earlyterm: "<<*lit<<" id= "<<id<<endl);
            }
            else
            {
                int dst = boundsMap.GetRank(id);
                if(dst == -1)
                    throw "Block ID not found";

                if (dst == rank)
                    inData.push_back(*lit);
                else
                    sendData[dst].push_back(*lit);
            }
        }

        DBG("-----SendP: "<<sendData<<endl);
        for (auto &i : sendData)
            SendICs(i.first, i.second);
    }
    increment += earlyTerm;

    if (increment > 0)
    {
        std::vector<int> msg = {MSG_TERMINATE, increment};
        DBG("-----SendAllMsg: msg="<<msg<<endl);
        SendAllMsg(msg);
    }

    //Check if we have anything coming in.
    std::list<ParticleCommData<Particle>> particleData;
    std::vector<MsgCommData> msgData;
    int val = 0;

    DBG("-----RecvAny..."<<endl);
    if (RecvAny(&msgData, &particleData, NULL, false))
    {
        DBG("-----Recv: M: "<<msgData.size()<<" P: "<<particleData.size()<<std::endl);
        for (auto &p : particleData)
            inData.push_back(p.p);

        for (auto &m : msgData)
        {
            if (m.message[0] == MSG_TERMINATE)
            {
                val += m.message[1];
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

    CheckPendingSendRequests();
    DBG("----ExchangeParticles Done: I= "<<inData.size()<<" T= "<<term.size()<<endl<<endl);
    if (stats) stats->stop("communication");

    return val;
}

}
