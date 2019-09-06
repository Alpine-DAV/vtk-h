#include <iostream>
#include <string.h>
#include "MemStream.h"
#include "Messenger.hpp"
#include <vtkh/utils/StreamUtil.hpp>

namespace vtkh
{

Messenger::Messenger(MPI_Comm comm)
  : m_mpi_comm(comm)
{
    MPI_Comm_size(comm, &nProcs);
    MPI_Comm_rank(comm, &rank);
    msgID = 0;
}

int
Messenger::CalcMessageBufferSize(int msgSz)
{
    MemStream buff;
    int rank = 0;

    std::vector<int> m(msgSz);
    vtkh::write(buff, rank);
    vtkh::write(buff, m);

    return buff.len();
}

void
Messenger::RegisterTag(int tag, int num_recvs, int size)
{
  if (messageTagInfo.find(tag) != messageTagInfo.end() || tag == TAG_ANY)
  {
    std::stringstream msg;
    msg<<"Invalid message tag: "<<tag<<std::endl;
    throw msg.str();
  }

  messageTagInfo[tag] = std::pair<int,int>(num_recvs, size);
}

void
Messenger::InitializeBuffers()
{
    //Setup receive buffers.
    std::map<int, std::pair<int, int> >::const_iterator it;
    for (it = messageTagInfo.begin(); it != messageTagInfo.end(); it++)
    {
        int tag = it->first, num = it->second.first;
        for (int i = 0; i < num; i++)
            PostRecv(tag);
    }
}

void
Messenger::CleanupRequests(int tag)
{
    std::vector<RequestTagPair> delKeys;
    for (bufferIterator i = recvBuffers.begin(); i != recvBuffers.end(); i++)
    {
        if (tag == TAG_ANY || tag == i->first.second)
            delKeys.push_back(i->first);
    }

    if (! delKeys.empty())
    {
        std::vector<RequestTagPair>::const_iterator it;
        for (it = delKeys.begin(); it != delKeys.end(); it++)
        {
            RequestTagPair v = *it;

            unsigned char *buff = recvBuffers[v];
            MPI_Cancel(&(v.first));
            delete [] buff;
            recvBuffers.erase(v);
        }
    }
}

void
Messenger::PostRecv(int tag)
{
    auto it = messageTagInfo.find(tag);
    if (it != messageTagInfo.end())
        PostRecv(tag, it->second.second);
}

void
Messenger::PostRecv(int tag, int sz, int src)
{
    sz += sizeof(Messenger::Header);
    unsigned char *buff = new unsigned char[sz];
    memset(buff, 0, sz);

    MPI_Request req;
    if (src == -1)
        MPI_Irecv(buff, sz, MPI_BYTE, MPI_ANY_SOURCE, tag, m_mpi_comm, &req);
    else
        MPI_Irecv(buff, sz, MPI_BYTE, src, tag, m_mpi_comm, &req);

    RequestTagPair entry(req, tag);
    recvBuffers[entry] = buff;
}

void
Messenger::CheckPendingSendRequests()
{
    bufferIterator it;
    std::vector<MPI_Request> req, copy;
    std::vector<int> tags;

    for (it = sendBuffers.begin(); it != sendBuffers.end(); it++)
    {
        req.push_back(it->first.first);
        copy.push_back(it->first.first);
        tags.push_back(it->first.second);
    }

    if (req.empty())
        return;

    //See if any sends are done.
    int num = 0, *indices = new int[req.size()];
    MPI_Status *status = new MPI_Status[req.size()];
    int err = MPI_Testsome(req.size(), &req[0], &num, indices, status);
    if (err != MPI_SUCCESS)
    {
        std::cerr << "Err with MPI_Testsome in PARIC algorithm" << std::endl;
    }
    for (int i = 0; i < num; i++)
    {
        MPI_Request r = copy[indices[i]];
        int tag = tags[indices[i]];

        RequestTagPair k(r,tag);
        bufferIterator entry = sendBuffers.find(k);
        if (entry != sendBuffers.end())
        {
            delete [] entry->second;
            sendBuffers.erase(entry);
        }
    }

    delete [] indices;
    delete [] status;
}

bool
Messenger::PacketCompare(const unsigned char *a, const unsigned char *b)
{
    Messenger::Header ha, hb;
    memcpy(&ha, a, sizeof(ha));
    memcpy(&hb, b, sizeof(hb));

    return ha.packet < hb.packet;
}

void
Messenger::PrepareForSend(int tag, MemStream *buff, std::vector<unsigned char *> &buffList)
{
    auto  it = messageTagInfo.find(tag);
    if (it == messageTagInfo.end())
    {
      std::stringstream msg;
      msg<<"Message tag not found: "<<tag<<std::endl;
      throw msg.str();
    }

    int bytesLeft = buff->len();
    int maxDataLen = it->second.second;
    Messenger::Header header;
    header.tag = tag;
    header.rank = rank;
    header.id = msgID;
    header.numPackets = 1;
    if (buff->len() > (unsigned int)maxDataLen)
        header.numPackets += buff->len() / maxDataLen;

    header.packet = 0;
    header.packetSz = 0;
    header.dataSz = 0;
    msgID++;

    buffList.resize(header.numPackets);
    size_t pos = 0;
    for (int i = 0; i < header.numPackets; i++)
    {
        header.packet = i;
        if (i == (header.numPackets-1))
            header.dataSz = bytesLeft;
        else
            header.dataSz = maxDataLen;

        header.packetSz = header.dataSz + sizeof(header);
        unsigned char *b = new unsigned char[header.packetSz];

        //Write the header.
        unsigned char *bPtr = b;
        memcpy(bPtr, &header, sizeof(header));
        bPtr += sizeof(header);

        //Write the data.
        memcpy(bPtr, &buff->data()[pos], header.dataSz);
        pos += header.dataSz;

        buffList[i] = b;
        bytesLeft -= maxDataLen;
    }

}

void
Messenger::SendData(int dst, int tag, MemStream *buff)
{
    std::vector<unsigned char *> bufferList;

    //Add headers, break into multiple buffers if needed.
    PrepareForSend(tag, buff, bufferList);

    Messenger::Header header;
    for (size_t i = 0; i < bufferList.size(); i++)
    {
        memcpy(&header, bufferList[i], sizeof(header));
        MPI_Request req;
        int err = MPI_Isend(bufferList[i], header.packetSz, MPI_BYTE, dst,
                            tag, m_mpi_comm, &req);
        if (err != MPI_SUCCESS)
        {
            std::cerr << "Err with MPI_Isend in SendData algorithm" << std::endl;
        }

        //Add it to sendBuffers
        RequestTagPair entry(req, tag);
        sendBuffers[entry] = bufferList[i];
    }

    delete buff;
}

bool
Messenger::RecvData(int tag, std::vector<MemStream *> &buffers,
                            bool blockAndWait)
{
    std::set<int> setTag;
    setTag.insert(tag);
    std::vector<std::pair<int, MemStream *> > b;
    buffers.resize(0);
    if (RecvData(setTag, b, blockAndWait))
    {
        buffers.resize(b.size());
        for (size_t i = 0; i < b.size(); i++)
            buffers[i] = b[i].second;
        return true;
    }
    return false;
}

bool
Messenger::RecvData(std::set<int> &tags,
                    std::vector<std::pair<int, MemStream *> > &buffers,
                    bool blockAndWait)
{
    buffers.resize(0);

    //Find all recv of type tag.
    std::vector<MPI_Request> req, copy;
    std::vector<int> reqTags;
    for (bufferIterator i = recvBuffers.begin(); i != recvBuffers.end(); i++)
    {
        if (tags.find(i->first.second) != tags.end())
        {
            req.push_back(i->first.first);
            copy.push_back(i->first.first);
            reqTags.push_back(i->first.second);
        }
    }

    if (req.empty())
        return false;

    MPI_Status *status = new MPI_Status[req.size()];
    int *indices = new int[req.size()], num = 0;
    if (blockAndWait)
        MPI_Waitsome(req.size(), &req[0], &num, indices, status);
    else
        MPI_Testsome(req.size(), &req[0], &num, indices, status);

    if (num == 0)
    {
        delete [] status;
        delete [] indices;
        return false;
    }

    std::vector<unsigned char *> incomingBuffers(num);
    for (int i = 0; i < num; i++)
    {
        RequestTagPair entry(copy[indices[i]], reqTags[indices[i]]);
        bufferIterator it = recvBuffers.find(entry);
        if ( it == recvBuffers.end())
        {
            delete [] status;
            delete [] indices;
            throw "receive buffer not found";
        }

        incomingBuffers[i] = it->second;
        recvBuffers.erase(it);
    }

    ProcessReceivedBuffers(incomingBuffers, buffers);

    for (int i = 0; i < num; i++)
        PostRecv(reqTags[indices[i]]);

    delete [] status;
    delete [] indices;

    return ! buffers.empty();
}

void
Messenger::ProcessReceivedBuffers(std::vector<unsigned char*> &incomingBuffers,
                                  std::vector<std::pair<int, MemStream *> > &buffers)
{
    for (size_t i = 0; i < incomingBuffers.size(); i++)
    {
        unsigned char *buff = incomingBuffers[i];

        //Grab the header.
        Messenger::Header header;
        memcpy(&header, buff, sizeof(header));

        //Only 1 packet, strip off header and add to list.
        if (header.numPackets == 1)
        {
            MemStream *b = new MemStream(header.dataSz, (buff + sizeof(header)));
            b->rewind();
            std::pair<int, MemStream*> entry(header.tag, b);
            buffers.push_back(entry);
            delete [] buff;
        }

        //Multi packet....
        else
        {
            RankIdPair k(header.rank, header.id);
            packetIterator i2 = recvPackets.find(k);

            //First packet. Create a new list and add it.
            if (i2 == recvPackets.end())
            {
                std::list<unsigned char *> l;
                l.push_back(buff);
                recvPackets[k] = l;
            }
            else
            {
                i2->second.push_back(buff);

                // The last packet came in, merge into one MemStream.
                if (i2->second.size() == (size_t)header.numPackets)
                {
                    //Sort the packets into proper order.
                    i2->second.sort(Messenger::PacketCompare);

                    MemStream *mergedBuff = new MemStream;
                    std::list<unsigned char *>::iterator listIt;

                    for (listIt = i2->second.begin(); listIt != i2->second.end(); listIt++)
                    {
                        unsigned char *bi = *listIt;

                        Messenger::Header header;
                        memcpy(&header, bi, sizeof(header));
                        mergedBuff->write_binary((bi+sizeof(header)), header.dataSz);
                        delete [] bi;
                    }

                    mergedBuff->rewind();
                    std::pair<int, MemStream*> entry(header.tag, mergedBuff);
                    buffers.push_back(entry);
                    recvPackets.erase(i2);
                }
            }
        }
    }
}

}// namespace vtkh
