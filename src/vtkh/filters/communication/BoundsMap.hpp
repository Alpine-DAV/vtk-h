#ifndef VTK_H_BOUNDS_MAP_HPP
#define VTK_H_BOUNDS_MAP_HPP

#include <vtkm/Bounds.h>

#include <vtkh/vtkh_exports.h>
#include <vtkh/vtkh.hpp>
#include <vtkh/filters/Filter.hpp>
#include <vtkh/DataSet.hpp>
#include <vtkh/filters/Particle.hpp>
#include <vtkh/utils/StreamUtil.hpp>

#include <string>
#include <vector>
#include <fstream>
#include <deque>
#include <algorithm>
#include <map>

#if VTKH_PARALLEL
#include <mpi.h>
#endif

namespace vtkh
{

class VTKH_API BoundsMap
{
public:
  BoundsMap() {}
  BoundsMap(const BoundsMap &_bm)
      : bm(_bm.bm), m_rank_map(_bm.m_rank_map), globalBounds(_bm.globalBounds)
  {
  }

  void Clear()
  {
    bm.clear();
    m_rank_map.clear();
  }

  void AddBlock(int id, const vtkm::Bounds &bounds)
  {
      if (bm.find(id) == bm.end())
          bm[id] = bounds;
      else
          throw "Duplicate block";
      m_rank_map[id] = -1;
  }

  template <template <typename, typename> class Container,
            typename Allocator=std::allocator<Particle>>
  void FindBlockIDs(const Container<Particle, Allocator> &particles,
                    std::vector<std::vector<int>> &blockIDs,
                    bool ignoreCurrentBlock=true) const
  {
      size_t sz = particles.size();
      blockIDs.resize(sz);
      auto pit = particles.begin();
      auto oit = blockIDs.begin();
      for ( ; pit != particles.end(); pit++, oit++)
          *oit = FindBlock(*pit, ignoreCurrentBlock && !pit->blockIds.empty());
  }

    std::vector<int> FindBlock(const vtkh::Particle &p,
                               bool ignoreCurrentBlock) const
  {
      std::vector<int> res;
      for (auto it = bm.begin(); it != bm.end(); it++)
      {
          if (ignoreCurrentBlock && !p.blockIds.empty() && p.blockIds[0] == it->first)
              continue;
          if (p.coords[0] >= it->second.X.Min &&
              p.coords[0] < it->second.X.Max &&
              p.coords[1] >= it->second.Y.Min &&
              p.coords[1] < it->second.Y.Max &&
              p.coords[2] >= it->second.Z.Min &&
              p.coords[2] < it->second.Z.Max)
          {
              res.push_back(it->first);
          }
      }
      return res;
  }

  int GetRank(const int &block_id)
  {
    auto it = m_rank_map.find(block_id);
    int rank = -1;
    if(it != m_rank_map.end())
      rank = m_rank_map[block_id];

    return rank;
  }

  void Build()
  {
    int size = bm.size();
#if VTKH_PARALLEL
    int rank;
    int procs;
    MPI_Comm comm = MPI_Comm_f2c(vtkh::GetMPICommHandle());
    MPI_Comm_rank(comm, &rank);
    MPI_Comm_size(comm, &procs);
    int *dom_counts = new int[procs];
    int *box_counts = new int[procs];
    MPI_Allgather(&size, 1, MPI_INT, dom_counts, 1, MPI_INT, comm);

    // prefix sum to build incoming buffers offsets
    int *box_offsets = new int[procs];
    int *dom_offsets = new int[procs];
    box_offsets[0] = 0;
    dom_offsets[0] = 0;
    box_counts[0] = dom_counts[0] * 6;
    for(int i = 1; i < procs; ++i)
    {
      box_offsets[i] = box_offsets[i-1] + dom_counts[i-1] * 6;
      dom_offsets[i] = dom_offsets[i-1] + dom_counts[i-1];
      box_counts[i] = dom_counts[i] * 6;
      //if(rank == 0) std::cout<<" box count "<<i<<" "<<box_counts[i]<<"\n";
      //if(rank == 0) std::cout<<" dom_offset "<<i<<" "<<dom_offsets[i]<<"\n";
      //if(rank == 0) std::cout<<" box_offsett "<<i<<" "<<box_offsets[i]<<"\n";
    }

    int total_boxs = dom_offsets[procs - 1] + dom_counts[procs - 1];
    double *box_send_buff = new double[size * 6];
    int *dom_send_buff = new int[size];

    int counter = 0;
    for (auto it = bm.begin(); it != bm.end(); it++)
    {
       const int offset = counter * 6;
       box_send_buff[offset + 0] = it->second.X.Min;
       box_send_buff[offset + 1] = it->second.X.Max;
       box_send_buff[offset + 2] = it->second.Y.Min;
       box_send_buff[offset + 3] = it->second.Y.Max;
       box_send_buff[offset + 4] = it->second.Z.Min;
       box_send_buff[offset + 5] = it->second.Z.Max;
       dom_send_buff[counter] = it->first;
       counter++;
    }
    //std::cout<<"total boxs "<<total_boxs<<"\n";
    double *box_rec_buff = new double[total_boxs * 6];
    int *dom_rec_buff = new int[total_boxs];
    MPI_Allgatherv(box_send_buff, size*6, MPI_DOUBLE, box_rec_buff, box_counts, box_offsets, MPI_DOUBLE, comm);
    MPI_Allgatherv(dom_send_buff, size, MPI_INT, dom_rec_buff, dom_counts, dom_offsets, MPI_INT, comm);

    bm.clear();
    m_rank_map.clear();

    //build a map of rank that handles empty counts
    int *rank_map = new int[total_boxs];
    int idx = 0;
    for(int i = 0; i < procs; ++i)
    {
      for(int d = 0; d < dom_counts[i]; ++d)
      {
        rank_map[idx] = i;
        ++idx;
      }
    }

    for(int i = 0; i < total_boxs; ++i)
    {
      const int offset = i * 6;
      int dom_id = dom_rec_buff[i];
      vtkm::Bounds &bounds = bm[dom_id];
      bounds.X.Min = box_rec_buff[offset + 0];
      bounds.X.Max = box_rec_buff[offset + 1];
      bounds.Y.Min = box_rec_buff[offset + 2];
      bounds.Y.Max = box_rec_buff[offset + 3];
      bounds.Z.Min = box_rec_buff[offset + 4];
      bounds.Z.Max = box_rec_buff[offset + 5];

      m_rank_map[dom_id] = rank_map[i];
    }

    delete[] dom_send_buff;
    delete[] dom_rec_buff;
    delete[] box_send_buff;
    delete[] box_rec_buff;
    delete[] dom_offsets;
    delete[] box_offsets;
    delete[] dom_counts;
    delete[] box_counts;
    delete[] rank_map;
#endif

    //Get the global bounds.
    for (auto &it : bm)
        globalBounds.Include(it.second);

    // Placeholder for more complex representatoin like a bvh that needs to
    // be constructed
  }

  std::map<int, vtkm::Bounds> bm; // map<dom_id, bounds>
  std::map<int, int> m_rank_map;  // map<dom_id,rank>
  vtkm::Bounds globalBounds;
protected:

};

inline std::ostream &operator<<(std::ostream &os, const vtkh::BoundsMap &bm)
  {
//    os<<"BoundsMap: "<<bm.globalBounds<<" d->r "<<bm.m_rank_map<<" d->b "<<bm.bm;
    return os;
  }



} // namespace vtkh
#endif
