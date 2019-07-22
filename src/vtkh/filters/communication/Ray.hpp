#ifndef VTK_H_RAY_HPP
#define VTK_H_RAY_HPP

#include <iostream>
#include <vector>
#include <vtkm/Types.h>

#include <vtkh/filters/communication/MemStream.h>
namespace vtkh
{

class Ray
{
public:
  vtkm::Vec<vtkm::Float32,3> m_origin;
  vtkm::Vec<vtkm::Float32,3> m_dir;
  vtkm::Vec<vtkm::Float32,3> m_color;
  vtkm::Vec<vtkm::Float32,3> m_throughput;
  vtkm::UInt8                m_status;
  vtkm::Int32                m_depth;
  vtkm::Float32              m_distance;
  vtkm::Float32              m_max_distance;
  vtkm::Int32                m_pixel_id;
  vtkm::Int32                m_dest_dom;

  friend std::ostream &operator<<(std::ostream &os, const vtkh::Ray &r)
  {
    os<<"Ray("<<r.m_pixel_id<<"):\n";
    os<<"  Origin     "<<r.m_origin[0]<<", "<<r.m_origin[1]<<", "<<r.m_origin[0]<<"\n";
    os<<"  Dir        "<<r.m_dir[0]<<", "<<r.m_dir[1]<<", "<<r.m_dir[2]<<"\n";
    os<<"  Color      "<<r.m_color[0]<<", "<<r.m_color[1]<<", "<<r.m_color[2]<<"\n";
    os<<"  Throughput "<<r.m_throughput[0]<<", "<<r.m_throughput[1]<<", "<<r.m_throughput[2]<<"\n";
    return os;
  }
};

template<>
struct Serialization<vtkh::Ray>
{
  static void write(MemStream &memstream, const vtkh::Ray &data)
  {
    vtkh::write(memstream, data.m_origin);
    vtkh::write(memstream, data.m_dir);
    vtkh::write(memstream, data.m_color);
    vtkh::write(memstream, data.m_throughput);
    vtkh::write(memstream, data.m_distance);
    vtkh::write(memstream, data.m_max_distance);
    vtkh::write(memstream, data.m_pixel_id);
    vtkh::write(memstream, data.m_depth);
    vtkh::write(memstream, data.m_dest_dom);
  }//

  static void read(MemStream &memstream, vtkh::Ray &data)
  {
    vtkh::read(memstream, data.m_origin);
    vtkh::read(memstream, data.m_dir);
    vtkh::read(memstream, data.m_color);
    vtkh::read(memstream, data.m_throughput);
    vtkh::read(memstream, data.m_distance);
    vtkh::read(memstream, data.m_max_distance);
    vtkh::read(memstream, data.m_pixel_id);
    vtkh::read(memstream, data.m_depth);
    vtkh::read(memstream, data.m_dest_dom);
  }
};
} //namespace vtkh
#endif
