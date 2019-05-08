#ifndef VTK_H_SPATIAL_QUERY_DETAIL_HPP
#define VTK_H_SPATIAL_QUERY_DETAIL_HPP

#include <vtkm/worklet/WorkletMapField.h>

namespace vtkh
{
namespace detail
{

template <typename Device>
class AABBDataExec
{
public:
  using FloatHandle = vtkm::cont::ArrayHandle<vtkm::Float32>;
  using FloatPortal = typename FloatHandle::ExecutionTypes<Device>::PortalConst;

  FloatPortal Xmins;
  FloatPortal Ymins;
  FloatPortal Zmins;

  FloatPortal Xmaxs;
  FloatPortal Ymaxs;
  FloatPortal Zmaxs;

public:
  AABBDataExec() = default;

  AABBDataExec(const FloatHandle& xmins,
               const FloatHandle& ymins,
               const FloatHandle& zmins,
               const FloatHandle& xmaxs,
               const FloatHandle& ymaxs,
               const FloatHandle& zmaxs)
    : Xmins(xmins.PrepareForInput(Device())),
      Ymins(ymins.PrepareForInput(Device())),
      Zmins(zmins.PrepareForInput(Device())),
      Xmaxs(xmaxs.PrepareForInput(Device())),
      Ymaxs(ymaxs.PrepareForInput(Device())),
      Zmaxs(zmaxs.PrepareForInput(Device()))
  {
  }

  void print_aabb(const int index) const
  {
    printf("AABB [%f, %f. %f] - [%f, %f, %f]\n", Xmins.Get(index),
                                                 Ymins.Get(index),
                                                 Zmins.Get(index),
                                                 Xmaxs.Get(index),
                                                 Ymaxs.Get(index),
                                                 Zmaxs.Get(index));
  }

  template <typename Precision>
  VTKM_EXEC inline void IntersectAABB(const vtkm::Id &aabb_index,
                                      const vtkm::Vec<Precision, 3>& originDir,
                                      const vtkm::Vec<Precision, 3>& invDir,
                                      const Precision& closestDistance,
                                      const Precision& minDistance,
                                      vtkm::Vec<Precision,2> &min_max) const
  {
      const vtkm::Float32 xmin = Xmins.Get(aabb_index);
      const vtkm::Float32 ymin = Ymins.Get(aabb_index);
      const vtkm::Float32 zmin = Zmins.Get(aabb_index);

      const vtkm::Float32 xmax = Xmaxs.Get(aabb_index);
      const vtkm::Float32 ymax = Ymaxs.Get(aabb_index);
      const vtkm::Float32 zmax = Zmaxs.Get(aabb_index);
      //printf("AABB [%f, %f. %f] - [%f, %f, %f]\n", xmin, ymin, zmin, xmax, ymax, zmax);
      // do intersect
      const Precision xmin0 = xmin * invDir[0] - originDir[0];
      const Precision ymin0 = ymin * invDir[1] - originDir[1];
      const Precision zmin0 = zmin * invDir[2] - originDir[2];
      const Precision xmax0 = xmax * invDir[0] - originDir[0];
      const Precision ymax0 = ymax * invDir[1] - originDir[1];
      const Precision zmax0 = zmax * invDir[2] - originDir[2];

      min_max[0] = vtkm::Max(
        vtkm::Max(vtkm::Max(vtkm::Min(ymin0, ymax0), vtkm::Min(xmin0, xmax0)), vtkm::Min(zmin0, zmax0)),
        minDistance);
      min_max[1] = vtkm::Min(
        vtkm::Min(vtkm::Min(vtkm::Max(ymin0, ymax0), vtkm::Max(xmin0, xmax0)), vtkm::Max(zmin0, zmax0)),
        closestDistance);

  }

  template <typename Precision, typename LeafPortalType>
  VTKM_EXEC inline void ClosestInterval(const vtkm::Int32& currentNode,
                                        const vtkm::Vec<Precision, 3>& originDir,
                                        const vtkm::Vec<Precision, 3>& invDir,
                                        vtkm::Vec<Precision,2> &min_max,
                                        Precision& closestDistance,
                                        LeafPortalType &leafs,
                                        const Precision& minDistance) const
  {
    const vtkm::Id count = leafs.Get(currentNode);
    for (vtkm::Id i = 1; i <= count; ++i)
    {
      const vtkm::Id aabb_index = leafs.Get(currentNode + i);
      //std::cout<<"********\n";
      //std::cout<<"* Current min "<<minDistance<<" closest "<<closestDistance<<"\n";
      //std::cout<<"* index  "<<aabb_index<<"\n";
      //std::cout<<"* "; print_aabb(aabb_index);
      Precision distance = -1.;
      vtkm::Vec<Precision,2> mm;
      IntersectAABB(aabb_index, originDir, invDir,closestDistance, minDistance, mm);
      if(mm[1] > mm[0])
      {
        //std::cout<<"* Valid intersection \n";
        distance = mm[0];
      }
      //std::cout<<"* "<<mm[0]<<" "<<mm[1]<<"\n";
      if (distance != -1. && distance < closestDistance && distance >= minDistance)
      {
        //std::cout<<"* Closest so far\n";
        closestDistance = distance;
        min_max = mm;
      }
      //std::cout<<"*********\n\n";
    } // for
  }

  template <typename Precision, typename LeafPortalType>
  VTKM_EXEC inline int InsideInterval(const vtkm::Int32& currentNode,
                                      const vtkm::Vec<Precision, 3>& originDir,
                                      const vtkm::Vec<Precision, 3>& invDir,
                                      const vtkm::Vec<Precision,2> &min_max,
                                      LeafPortalType &leafs) const
  {
    int count = 0;
    const vtkm::Id ncount = leafs.Get(currentNode);
    for (vtkm::Id i = 1; i <= ncount; ++i)
    {
      const vtkm::Id aabb_index = leafs.Get(currentNode + i);
      //std::cout<<"min "<<minDistance<<" closest "<<closestDistance<<"\n";
      //std::cout<<"inde  "<<aabb_index<<"\n";
      Precision distance = -1.;
      vtkm::Vec<Precision,2> mm;
      IntersectAABB(aabb_index, originDir, invDir, min_max[1], min_max[0], mm);
      //std::cout<<"Inside interval "<<min_max[0]<<" - "<<min_max[1]<<" aabb "<<mm[0]<<"- "<<mm[1]<<"\n";
      bool valid = mm[1] > mm[0];
      if(valid && mm[1] >= mm[0] && mm[0] >= min_max[0] && mm[0] < min_max[1])
      {
        count++;
      }
    } // for
    return count;
  }

  template <typename Precision, typename LeafPortalType, typename OutputType>
  VTKM_EXEC inline void InsideInterval(const vtkm::Int32& currentNode,
                                      const vtkm::Vec<Precision, 3>& originDir,
                                      const vtkm::Vec<Precision, 3>& invDir,
                                      const vtkm::Vec<Precision,2> &min_max,
                                      LeafPortalType &leafs,
                                      vtkm::Id &index,
                                      OutputType &output) const
  {
    const vtkm::Id ncount = leafs.Get(currentNode);
    for (vtkm::Id i = 1; i <= ncount; ++i)
    {
      const vtkm::Id aabb_index = leafs.Get(currentNode + i);
      //std::cout<<"min "<<minDistance<<" closest "<<closestDistance<<"\n";
      //std::cout<<"inde  "<<aabb_index<<"\n";
      Precision distance = -1.;
      vtkm::Vec<Precision,2> mm;
      IntersectAABB(aabb_index, originDir, invDir,min_max[1], min_max[0], mm);
      bool valid = mm[1] > mm[0];
      if(valid && mm[1] >= mm[0] && mm[0] >= min_max[0] && mm[0] < min_max[1])
      {
        output.Set(index, aabb_index);
        index++;
      }
    } // for
    //return count;
  }
};

class AABBData : public vtkm::cont::ExecutionObjectBase
{
protected:
  using FloatHandle = vtkm::cont::ArrayHandle<vtkm::Float32>;

  FloatHandle Xmins;
  FloatHandle Ymins;
  FloatHandle Zmins;

  FloatHandle Xmaxs;
  FloatHandle Ymaxs;
  FloatHandle Zmaxs;

public:
  AABBData(const FloatHandle& xmins,
           const FloatHandle& ymins,
           const FloatHandle& zmins,
           const FloatHandle& xmaxs,
           const FloatHandle& ymaxs,
           const FloatHandle& zmaxs)
    : Xmins(xmins),
      Ymins(ymins),
      Zmins(zmins),
      Xmaxs(xmaxs),
      Ymaxs(ymaxs),
      Zmaxs(zmaxs)
  {}

  template <typename Device>
  VTKM_CONT AABBDataExec<Device> PrepareForExecution(Device) const
  {
    return AABBDataExec<Device>(Xmins, Ymins, Zmins, Xmaxs, Ymaxs, Zmaxs);
  }
};


VTKM_EXEC
inline vtkm::Float32 rcp(vtkm::Float32 f) { return 1.0f / f; }

VTKM_EXEC
inline vtkm::Float32 rcp_safe(vtkm::Float32 f)
{
  return rcp((vtkm::Abs(f) < 1e-8f) ? 1e-8f : f);
}

VTKM_EXEC
inline vtkm::Float64 rcp(vtkm::Float64 f) { return 1.0 / f; }

VTKM_EXEC
inline vtkm::Float64 rcp_safe(vtkm::Float64 f)
{
  return rcp((vtkm::Abs(f) < 1e-8f) ? 1e-8f : f);
}

class ClosestCandidates : public vtkm::worklet::WorkletMapField
{
private:

public:
  VTKM_CONT
  ClosestCandidates() {}
  using ControlSignature = void(FieldIn,
                                FieldIn,
                                FieldIn,
                                FieldIn,
                                FieldOut,
                                WholeArrayIn,
                                WholeArrayIn,
                                ExecObject);

  using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8);


  template <typename Precision,
            typename BVHInnerType,
            typename BVHLeafType,
            typename AABBType>
  VTKM_EXEC void operator()(const vtkm::Vec<Precision, 3>& dir,
                            const vtkm::Vec<Precision, 3>& origin,
                            const Precision& minDistance,
                            const Precision& maxDistance,
                            vtkm::Vec<vtkm::Float32,2>& min_max,
                            const BVHInnerType& bvh,
                            const BVHLeafType& leafs,
                            AABBType &aabbs) const
  {
    vtkm::Id hitIndex = -1;
    min_max[0] = vtkm::Infinity32();
    min_max[1] = vtkm::NegativeInfinity32();

    Precision closestDistance = maxDistance;

    const vtkm::Vec<Precision, 3> invDir(rcp_safe(dir[0]),
                                         rcp_safe(dir[1]),
                                         rcp_safe(dir[2]));

    const vtkm::Vec<Precision, 3> originDir = origin * invDir;

    vtkm::Int32 currentNode;
    vtkm::Int32 todo[64];
    vtkm::Int32 stackptr = 0;
    constexpr vtkm::Int32 barrier = (vtkm::Int32)2000000000;
    currentNode = 0;

    todo[stackptr] = barrier;

    while (currentNode != barrier)
    {
      if (currentNode > -1)
      {

        bool hitLeftChild, hitRightChild;

        const vtkm::Vec<vtkm::Float32, 4> first4 = bvh.Get(currentNode);
        const vtkm::Vec<vtkm::Float32, 4> second4 = bvh.Get(currentNode + 1);
        const vtkm::Vec<vtkm::Float32, 4> third4 = bvh.Get(currentNode + 2);

        const Precision xmin0 = first4[0] * invDir[0] - originDir[0];
        const Precision ymin0 = first4[1] * invDir[1] - originDir[1];
        const Precision zmin0 = first4[2] * invDir[2] - originDir[2];
        const Precision xmax0 = first4[3] * invDir[0] - originDir[0];
        const Precision ymax0 = second4[0] * invDir[1] - originDir[1];
        const Precision zmax0 = second4[1] * invDir[2] - originDir[2];

        const Precision min0 = vtkm::Max(
          vtkm::Max(vtkm::Max(vtkm::Min(ymin0, ymax0), vtkm::Min(xmin0, xmax0)), vtkm::Min(zmin0, zmax0)),
          minDistance);
        const Precision max0 = vtkm::Min(
          vtkm::Min(vtkm::Min(vtkm::Max(ymin0, ymax0), vtkm::Max(xmin0, xmax0)), vtkm::Max(zmin0, zmax0)),
          closestDistance);

        hitLeftChild = (max0 >= min0);

        const Precision xmin1 = second4[2] * invDir[0] - originDir[0];
        const Precision ymin1 = second4[3] * invDir[1] - originDir[1];
        const Precision zmin1 = third4[0] * invDir[2] - originDir[2];
        const Precision xmax1 = third4[1] * invDir[0] - originDir[0];
        const Precision ymax1 = third4[2] * invDir[1] - originDir[1];
        const Precision zmax1 = third4[3] * invDir[2] - originDir[2];

        const Precision min1 = vtkm::Max(
          vtkm::Max(vtkm::Max(vtkm::Min(ymin1, ymax1), vtkm::Min(xmin1, xmax1)), vtkm::Min(zmin1, zmax1)),
          minDistance);
        const Precision max1 = vtkm::Min(
          vtkm::Min(vtkm::Min(vtkm::Max(ymin1, ymax1), vtkm::Max(xmin1, xmax1)), vtkm::Max(zmin1, zmax1)),
          closestDistance);
        hitRightChild = (max1 >= min1);

        if (!hitLeftChild && !hitRightChild)
        {
          currentNode = todo[stackptr];
          stackptr--;
        }
        else
        {
          vtkm::Vec<vtkm::Float32, 4> children = bvh.Get(currentNode + 3);
          vtkm::Int32 leftChild;
          memcpy(&leftChild, &children[0], 4);
          vtkm::Int32 rightChild;
          memcpy(&rightChild, &children[1], 4);
          currentNode = (hitLeftChild) ? leftChild : rightChild;
          if (hitLeftChild && hitRightChild)
          {
            if (min1 >  min0)
            {
              currentNode = rightChild;
              stackptr++;
              todo[stackptr] = leftChild;
            }
            else
            {
              stackptr++;
              todo[stackptr] = rightChild;
            }
          }
        }
      } // if inner node

      if (currentNode < 0 && currentNode != barrier) //check register usage
      {
        currentNode = -currentNode - 1; //swap the neg address

        aabbs.ClosestInterval(currentNode,
                              originDir,
                              invDir,
                              min_max,
                              closestDistance,
                              leafs,
                              minDistance);

        currentNode = todo[stackptr];
        stackptr--;
      } // if leaf node

    } //while
    //std::cout<<"Min "<<min_max[0]<<"\n";
    //std::cout<<"Max "<<min_max[1]<<"\n";
    //if (hitIndex != -1)
    //  distance = closestDistance;
  } // ()
};

class CountCandidates : public vtkm::worklet::WorkletMapField
{
private:

public:
  VTKM_CONT
  CountCandidates() {}
  using ControlSignature = void(FieldIn,
                                FieldIn,
                                FieldIn,
                                FieldIn,
                                FieldOut,
                                WholeArrayIn,
                                WholeArrayIn,
                                ExecObject);

  using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8);


  template <typename Precision,
            typename BVHType,
            typename LeafType,
            typename AABBType>
  VTKM_EXEC void operator()(const vtkm::Vec<Precision, 3>& dir,
                            const vtkm::Vec<Precision, 3>& origin,
                            const vtkm::Vec<Precision,2>& min_max,
                            const vtkm::UInt8& status,
                            vtkm::Id& count,
                            const BVHType& bvh,
                            const LeafType& leafs,
                            const AABBType &aabbs) const
  {
    count = 0;

    if(status != RAY_ACTIVE)
    {
      return;
    }

    Precision closestDistance = min_max[1];
    Precision minDistance = min_max[0];

    const vtkm::Vec<Precision, 3> invDir(rcp_safe(dir[0]),
                                         rcp_safe(dir[1]),
                                         rcp_safe(dir[2]));

    const vtkm::Vec<Precision, 3> originDir = origin * invDir;

    vtkm::Int32 currentNode;
    vtkm::Int32 todo[64];
    vtkm::Int32 stackptr = 0;
    constexpr vtkm::Int32 barrier = (vtkm::Int32)2000000000;
    currentNode = 0;

    todo[stackptr] = barrier;

    while (currentNode != barrier)
    {
      if (currentNode > -1)
      {

        bool hitLeftChild, hitRightChild;

        const vtkm::Vec<vtkm::Float32, 4> first4 = bvh.Get(currentNode);
        const vtkm::Vec<vtkm::Float32, 4> second4 = bvh.Get(currentNode + 1);
        const vtkm::Vec<vtkm::Float32, 4> third4 = bvh.Get(currentNode + 2);

        const Precision xmin0 = first4[0] * invDir[0] - originDir[0];
        const Precision ymin0 = first4[1] * invDir[1] - originDir[1];
        const Precision zmin0 = first4[2] * invDir[2] - originDir[2];
        const Precision xmax0 = first4[3] * invDir[0] - originDir[0];
        const Precision ymax0 = second4[0] * invDir[1] - originDir[1];
        const Precision zmax0 = second4[1] * invDir[2] - originDir[2];

        const Precision min0 = vtkm::Max(
          vtkm::Max(vtkm::Max(vtkm::Min(ymin0, ymax0), vtkm::Min(xmin0, xmax0)), vtkm::Min(zmin0, zmax0)),
          minDistance);
        const Precision max0 = vtkm::Min(
          vtkm::Min(vtkm::Min(vtkm::Max(ymin0, ymax0), vtkm::Max(xmin0, xmax0)), vtkm::Max(zmin0, zmax0)),
          closestDistance);

        hitLeftChild = (max0 >= min0);

        const Precision xmin1 = second4[2] * invDir[0] - originDir[0];
        const Precision ymin1 = second4[3] * invDir[1] - originDir[1];
        const Precision zmin1 = third4[0] * invDir[2] - originDir[2];
        const Precision xmax1 = third4[1] * invDir[0] - originDir[0];
        const Precision ymax1 = third4[2] * invDir[1] - originDir[1];
        const Precision zmax1 = third4[3] * invDir[2] - originDir[2];

        const Precision min1 = vtkm::Max(
          vtkm::Max(vtkm::Max(vtkm::Min(ymin1, ymax1), vtkm::Min(xmin1, xmax1)), vtkm::Min(zmin1, zmax1)),
          minDistance);
        const Precision max1 = vtkm::Min(
          vtkm::Min(vtkm::Min(vtkm::Max(ymin1, ymax1), vtkm::Max(xmin1, xmax1)), vtkm::Max(zmin1, zmax1)),
          closestDistance);
        hitRightChild = (max1 >= min1);

        if (!hitLeftChild && !hitRightChild)
        {
          currentNode = todo[stackptr];
          stackptr--;
        }
        else
        {
          vtkm::Vec<vtkm::Float32, 4> children = bvh.Get(currentNode + 3);
          vtkm::Int32 leftChild;
          memcpy(&leftChild, &children[0], 4);
          vtkm::Int32 rightChild;
          memcpy(&rightChild, &children[1], 4);
          currentNode = (hitLeftChild) ? leftChild : rightChild;
          if (hitLeftChild && hitRightChild)
          {
            if (min1 >  min0)
            {
              currentNode = rightChild;
              stackptr++;
              todo[stackptr] = leftChild;
            }
            else
            {
              stackptr++;
              todo[stackptr] = rightChild;
            }
          }
        }
      } // if inner node

      if (currentNode < 0 && currentNode != barrier) //check register usage
      {
        currentNode = -currentNode - 1; //swap the neg address
        count += aabbs.InsideInterval(currentNode,
                                      originDir,
                                      invDir,
                                      min_max,
                                      leafs);

        currentNode = todo[stackptr];
        stackptr--;
      } // if leaf node

    } //while
    //std::cout<<"count "<<count<<"\n";
    //if (hitIndex != -1)
    //  distance = closestDistance;
  } // ()
};

class GetCandidates : public vtkm::worklet::WorkletMapField
{
private:

public:
  VTKM_CONT
  GetCandidates() {}
  using ControlSignature = void(FieldIn,
                                FieldIn,
                                FieldIn,
                                FieldIn,
                                WholeArrayIn,
                                WholeArrayIn,
                                ExecObject,
                                FieldIn,
                                WholeArrayOut);

  using ExecutionSignature = void(_1, _2, _3, _4, _5, _6, _7, _8, _9);


  template <typename Precision,
            typename BVHType,
            typename LeafType,
            typename AABBType,
            typename OutputType>
  VTKM_EXEC void operator()(const vtkm::Vec<Precision, 3>& dir,
                            const vtkm::Vec<Precision, 3>& origin,
                            const vtkm::Vec<Precision,2>& min_max,
                            const vtkm::Id& count,
                            const BVHType& bvh,
                            const LeafType& leafs,
                            const AABBType &aabbs,
                            const vtkm::Id& offset,
                            OutputType &output) const
  {
    if(count == 0) return;
    vtkm::Id index = offset;
    Precision closestDistance = min_max[1];
    Precision minDistance = min_max[0];

    const vtkm::Vec<Precision, 3> invDir(rcp_safe(dir[0]),
                                         rcp_safe(dir[1]),
                                         rcp_safe(dir[2]));

    const vtkm::Vec<Precision, 3> originDir = origin * invDir;

    vtkm::Int32 currentNode;
    vtkm::Int32 todo[64];
    vtkm::Int32 stackptr = 0;
    constexpr vtkm::Int32 barrier = (vtkm::Int32)2000000000;
    currentNode = 0;

    todo[stackptr] = barrier;

    while (currentNode != barrier)
    {
      if (currentNode > -1)
      {

        bool hitLeftChild, hitRightChild;

        const vtkm::Vec<vtkm::Float32, 4> first4 = bvh.Get(currentNode);
        const vtkm::Vec<vtkm::Float32, 4> second4 = bvh.Get(currentNode + 1);
        const vtkm::Vec<vtkm::Float32, 4> third4 = bvh.Get(currentNode + 2);

        const Precision xmin0 = first4[0] * invDir[0] - originDir[0];
        const Precision ymin0 = first4[1] * invDir[1] - originDir[1];
        const Precision zmin0 = first4[2] * invDir[2] - originDir[2];
        const Precision xmax0 = first4[3] * invDir[0] - originDir[0];
        const Precision ymax0 = second4[0] * invDir[1] - originDir[1];
        const Precision zmax0 = second4[1] * invDir[2] - originDir[2];

        const Precision min0 = vtkm::Max(
          vtkm::Max(vtkm::Max(vtkm::Min(ymin0, ymax0), vtkm::Min(xmin0, xmax0)), vtkm::Min(zmin0, zmax0)),
          minDistance);
        const Precision max0 = vtkm::Min(
          vtkm::Min(vtkm::Min(vtkm::Max(ymin0, ymax0), vtkm::Max(xmin0, xmax0)), vtkm::Max(zmin0, zmax0)),
          closestDistance);

        hitLeftChild = (max0 >= min0);

        const Precision xmin1 = second4[2] * invDir[0] - originDir[0];
        const Precision ymin1 = second4[3] * invDir[1] - originDir[1];
        const Precision zmin1 = third4[0] * invDir[2] - originDir[2];
        const Precision xmax1 = third4[1] * invDir[0] - originDir[0];
        const Precision ymax1 = third4[2] * invDir[1] - originDir[1];
        const Precision zmax1 = third4[3] * invDir[2] - originDir[2];

        const Precision min1 = vtkm::Max(
          vtkm::Max(vtkm::Max(vtkm::Min(ymin1, ymax1), vtkm::Min(xmin1, xmax1)), vtkm::Min(zmin1, zmax1)),
          minDistance);
        const Precision max1 = vtkm::Min(
          vtkm::Min(vtkm::Min(vtkm::Max(ymin1, ymax1), vtkm::Max(xmin1, xmax1)), vtkm::Max(zmin1, zmax1)),
          closestDistance);
        hitRightChild = (max1 >= min1);

        if (!hitLeftChild && !hitRightChild)
        {
          currentNode = todo[stackptr];
          stackptr--;
        }
        else
        {
          vtkm::Vec<vtkm::Float32, 4> children = bvh.Get(currentNode + 3);
          vtkm::Int32 leftChild;
          memcpy(&leftChild, &children[0], 4);
          vtkm::Int32 rightChild;
          memcpy(&rightChild, &children[1], 4);
          currentNode = (hitLeftChild) ? leftChild : rightChild;
          if (hitLeftChild && hitRightChild)
          {
            if (min1 >  min0)
            {
              currentNode = rightChild;
              stackptr++;
              todo[stackptr] = leftChild;
            }
            else
            {
              stackptr++;
              todo[stackptr] = rightChild;
            }
          }
        }
      } // if inner node

      if (currentNode < 0 && currentNode != barrier) //check register usage
      {
        currentNode = -currentNode - 1; //swap the neg address
        aabbs.InsideInterval(currentNode,
                             originDir,
                             invDir,
                             min_max,
                             leafs,
                             index,
                             output);

        currentNode = todo[stackptr];
        stackptr--;
      } // if leaf node

    } //while
    //std::cout<<"count "<<count<<"\n";
    //if (hitIndex != -1)
    //  distance = closestDistance;
  } // ()
};

}} // namespace vtkh::detail
#endif
