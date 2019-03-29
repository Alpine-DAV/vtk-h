#include <iostream>
#include <string.h>
#include "MemStream.h"
#include "DebugMeowMeow.hpp"
#include "avtParICAlgorithm.hpp"

using namespace std;
namespace vtkh
{

avtParICAlgorithm::avtParICAlgorithm(MPI_Comm comm)
    : Messenger(comm),
      done(false)
{
}

void
avtParICAlgorithm::RegisterMessages(int mSz,
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

    this->RegisterTag(avtParICAlgorithm::MESSAGE_TAG, numMsgRecvs, msgSize);
    this->RegisterTag(avtParICAlgorithm::PARTICLE_TAG, numSLRecvs, slSize * slsPerRecv);

    this->InitializeBuffers();
}

void
avtParICAlgorithm::SendMsg(int dst, vector<int> &msg)
{
    MemStream *buff = new MemStream;

    //Write data.
    vtkh::write(*buff, rank);
    vtkh::write(*buff, msg);

    SendData(dst, avtParICAlgorithm::MESSAGE_TAG, buff);
//    MsgCnt.value++;
//    CommTime.value += visitTimer->StopTimer(timerHandle, "SendMsg");
}

void
avtParICAlgorithm::SendAllMsg(vector<int> &msg)
{
    for (int i = 0; i < nProcs; i++)
        if (i != rank)
        {
            DBG("          ***************** SendMsg to "<<i<<" "<<msg<<endl);
            SendMsg(i, msg);
        }
}

bool
avtParICAlgorithm::RecvAny(vector<MsgCommData> *msgs,
                           list<ParticleCommData<Particle>> *recvICs,
                           vector<DSCommData> *ds,
                           bool blockAndWait)
{
    set<int> tags;
    if (msgs)
    {
        tags.insert(avtParICAlgorithm::MESSAGE_TAG);
        msgs->resize(0);
    }
    if (recvICs)
    {
        tags.insert(avtParICAlgorithm::PARTICLE_TAG);
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
        if (buffers[i].first == avtParICAlgorithm::MESSAGE_TAG)
        {
            int sendRank;
            vector<int> m;
            vtkh::read(*buffers[i].second, sendRank);
            vtkh::read(*buffers[i].second, m);

            MsgCommData msg(sendRank, m);

            msgs->push_back(msg);
        }
        else if (buffers[i].first == avtParICAlgorithm::PARTICLE_TAG)
        {
            std::cout<<"particle tag\n";
            int num, sendRank;
            vtkh::read(*buffers[i].second, sendRank);
            vtkh::read(*buffers[i].second, num);

            for (int j = 0; j < num; j++)
            {
                Particle recvP;
                //buffers[i].second->read(recvP);
                vtkh::read(*(buffers[i].second), recvP);
                //Serialization<Particle>::read(*(buffers[i].second), recvP);
                std::cout<<"reading "<<recvP<<"\n";
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
avtParICAlgorithm::RecvMsg(vector<MsgCommData> &msgs)
{
    return RecvAny(&msgs, NULL, NULL, false);
}

template <typename P, template <typename, typename> class Container,
          typename Allocator>
void avtParICAlgorithm::SendICs(int dst, Container<P, Allocator> &c)
{
    std::cout<<"Sending to rank "<<dst<<"\n";
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
    {
        std::cout<<"writing "<<p<<"\n";
        vtkh::write(*buff, p);
        //buff->write(p);
    }
    SendData(dst, avtParICAlgorithm::PARTICLE_TAG, buff);
    c.clear();
}

template <typename P, template <typename, typename> class Container,
          typename Allocator>
void avtParICAlgorithm::SendICs(std::map<int, Container<P, Allocator>> &m)
{
    for (auto mit = m.begin(); mit != m.end(); mit++)
        if (! mit->second.empty())
            SendICs(mit->first, mit->second);
}

template <typename P, template <typename, typename> class Container,
          typename Allocator>
bool avtParICAlgorithm::RecvICs(Container<ParticleCommData<P>, Allocator> &recvICs)
{
    return RecvAny(NULL, &recvICs, NULL, false);
}

template <typename P, template <typename, typename> class Container,
          typename Allocator>
bool avtParICAlgorithm::RecvICs(Container<P, Allocator> &recvICs)
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
bool avtParICAlgorithm::DoSendICs(int dst, vector<P> &ics)
{
    if (dst == rank)
    {
        cerr << "Error in avtParICAlgorithm::DoSendICs() Sending ICs to yourself" << endl;
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
    SendData(dst, avtParICAlgorithm::PARTICLE_TAG, buff);

    return true;
}


size_t
avtParICAlgorithm::Exchange(bool haveWork,
                            list<vtkh::Particle> &outData,
                            list<vtkh::Particle> &inData,
                            list<vtkh::Particle> &term,
                            vtkh::BoundsMap &boundsMap,
                            int numTerm)
{
    DBG("----Exchange: O="<<outData<<" I="<<inData<<" NT= "<<numTerm<<endl);

    map<int, list<Particle>> sendData;

    int earlyTerm = 0;
    if (!outData.empty())
    {
        vector<int> blockIds;
        boundsMap.FindBlockIDs(outData, blockIds);
        DBG("          -O.blockIds: "<<outData<<" "<<blockIds<<endl);

        auto bit = blockIds.begin();
        for (auto lit = outData.begin(); lit != outData.end(); lit++, bit++)
        {
            int id = *bit;
            lit->blockId = id;
            if (id == -1)
            {
                term.push_back(*lit);
                earlyTerm++;
                DBG("    earlyterm: "<<*lit<<" id= "<<id<<endl);
            }
            else
            {
                int dst = boundsMap.GetRank(id);
                DBG("---- boundsMap.GetRank("<<id<<")= "<<dst<<std::endl);
                if(dst == -1) std::cout<<"COUND NOT FIND RANK for block "<<id<<"\n";
                if (dst == rank)
                    inData.push_back(*lit);
                else
                    sendData[dst].push_back(*lit);
            }
        }

        DBG("   ---SendP: "<<sendData<<endl);
        for (auto &i : sendData)
            SendICs(i.first, i.second);
//            sdb->increment("numMsgSent", sendData.size());
    }
//        if (earlyTerm > 0)
//            sdb->increment("earlyTerm", earlyTerm);


    int terminations = earlyTerm + numTerm;

    if (terminations > 0)
    {
        vector<int> msg = {MSG_TERMINATE, terminations};
        DBG("   ---SendAllMsg: msg="<<msg<<endl);
        SendAllMsg(msg);
    }

    //Check if we have anything coming in.
    vector<MsgCommData> msgData;
    list<ParticleCommData<Particle>> particleData;
    int incomingTerm = 0;

    bool blockAndWait = false;
    DBG(" ---RecvAny..."<<endl);
    if (RecvAny(&msgData, &particleData, NULL, blockAndWait))
    {
        DBG(" ---Recv: M: "<<msgData<<" P: "<<particleData<<endl);
//            sdb->increment("numMsgRecv", msgData.size());

        for (auto &m : msgData)
        {
            if (m.message[0] == MSG_TERMINATE)
            {
                incomingTerm += m.message[1];
                //numTermReceived += m.message[1];
                DBG("     TERMinate: Recv: "<<m.message[1]<<endl);
            }
            else if (m.message[0] == MSG_DONE)
                done = true;
        }
        for (auto &p : particleData)
            inData.push_back(p.p);
/*
  if (numTermReceived > 0)
  {
  if (rank == 0)
  termCounter += numTermReceived;
  else
  {
  vector<int> msg(2);
  msg[0] = MSG_TERMINATE;
  msg[1] = numTermReceived;
  a.SendMsg(termSendMap[rank], msg);
//                    sdb->increment("numTermSent");
//                    sdb->increment("numMsgSent");
}
}
*/
    }
    else
        DBG("  ---RecvAny --Nothing in the can"<<endl);


    CheckPendingSendRequests();
    int returnTerm = incomingTerm + earlyTerm;
    DBG(" ---ExchangeDone: nt= "<<returnTerm<<" I= "<<inData.size()<<" T= "<<term.size()<<endl);
    return returnTerm;
}
}
