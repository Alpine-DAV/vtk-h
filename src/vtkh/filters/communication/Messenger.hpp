#ifndef VTKH_MESSENGER_H
#define VTKH_MESSENGER_H

#include <mpi.h>
#include <list>
#include <vector>
#include <set>
#include <map>

#include <vtkh/vtkh_exports.h>

namespace vtkh
{

class MemStream;

class VTKH_API Messenger
{
  public:
    Messenger(MPI_Comm comm);
    virtual ~Messenger() {Cleanup();}

    //Message headers.
    typedef struct
    {
        int rank, id, tag, numPackets, packet, packetSz, dataSz;
    } Header;

    static constexpr int TAG_ANY = -1;

    // Register message tags for this messenger
    // Must be called before InitializeBuffers
    void RegisterTag(int tag,         // Unique message tag
                     int num_recvs,   // number of receives to check each time
                     int size);       // size in bytes for each message

    // Creates receives buffers for all tags registered to this messenger
    void InitializeBuffers();

    void Cleanup() { CleanupRequests(); }


    //Manage communication.
    void CleanupRequests(int tag=TAG_ANY);
    void CheckPendingSendRequests();

  protected:
    void PostRecv(int tag);
    void PostRecv(int tag, int sz, int src=-1);
    void SendData(int dst, int tag, MemStream *buff);
    bool RecvData(std::set<int> &tags,
                  std::vector<std::pair<int,MemStream *>> &buffers,
                  bool blockAndWait=false);
    bool RecvData(int tag, std::vector<MemStream *> &buffers,
                  bool blockAndWait=false);
    void AddHeader(MemStream *buff);
    void RemoveHeader(MemStream *input, MemStream *header, MemStream *buff);

    template <typename P>
    bool DoSendICs(int dst, std::vector<P> &ics);
    void PrepareForSend(int tag, MemStream *buff, std::vector<unsigned char *> &buffList);
    static bool PacketCompare(const unsigned char *a, const unsigned char *b);
    void ProcessReceivedBuffers(std::vector<unsigned char*> &incomingBuffers,
                                std::vector<std::pair<int, MemStream *>> &buffers);

    // Send/Recv buffer management structures.
    typedef std::pair<MPI_Request, int> RequestTagPair;
    typedef std::pair<int, int> RankIdPair;
    typedef std::map<RequestTagPair, unsigned char *>::iterator bufferIterator;
    typedef std::map<RankIdPair, std::list<unsigned char *>>::iterator packetIterator;

    int rank, nProcs;
    MPI_Comm m_mpi_comm;
    std::map<RequestTagPair, unsigned char *> sendBuffers, recvBuffers;
    std::map<RankIdPair, std::list<unsigned char *>> recvPackets;

    // Maps MPI_TAG to pair(num buffers, data size).
    std::map<int, std::pair<int, int>> messageTagInfo;
    long msgID;

    static int CalcMessageBufferSize(int msgSz);
};

} // namespace vtkh
#endif
