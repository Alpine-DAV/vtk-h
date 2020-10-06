#ifndef VTK_H_PARTICLE_HPP
#define VTK_H_PARTICLE_HPP

#include <iostream>
#include <vector>
#include <vtkm/Particle.h>
#include <vtkm/Types.h>

#include <vtkh/filters/communication/MemStream.h>
#include <vtkh/utils/StreamUtil.hpp>

namespace vtkh
{

class Particle
{
public:
    Particle() : p(), blockIds() {}
    Particle(const std::vector<float> &c, int id) : p(vtkm::Vec3f(c[0],c[1],c[2]), id), blockIds() {}
    Particle(const std::vector<double> &c, int id) : p(vtkm::Vec3f(c[0],c[1],c[2]), id), blockIds() {}
    Particle(const double *c, int id) : p(vtkm::Vec3f(c[0],c[1],c[2]), id), blockIds() {}
    Particle(const float *c, int id) : p(vtkm::Vec3f(c[0],c[1],c[2]), id), blockIds() {}
    Particle(const vtkm::Vec<double,3> &c, int id) : p(c, id), blockIds() {}
    Particle(const vtkm::Vec<float,3> &c, int id) : p(c, id), blockIds() {}
    Particle(const Particle &part) : p(part.p), blockIds(part.blockIds) {}

    vtkm::Particle p;
    std::vector<int> blockIds;

    friend std::ostream &operator<<(std::ostream &os, const vtkh::Particle p)
    {
        os<<"(";
        os<<"P_"<<p.p.ID<<" ";
        os<<"["<<p.p.Pos[0]<<" "<<p.p.Pos[1]<<" "<<p.p.Pos[2]<<"] #s "<<p.p.NumSteps<<" ";
        os<<"{ ";
        if (p.p.Status.CheckOk()) os<<"OK ";
        if (p.p.Status.CheckTerminate()) os<<"TERM ";
        if (p.p.Status.CheckSpatialBounds()) os<<"OOB ";
        if (!p.p.Status.CheckTookAnySteps()) os<<"NO_STEPS ";
        os<<"}";
        os<<" bid = "<<p.blockIds;
        os<<")";
        return os;
    }
};

template<>
struct Serialization<vtkh::Particle>
{
  static void write(MemStream &memstream, const vtkh::Particle &data)
  {
    vtkh::write(memstream, data.p.Pos[0]);
    vtkh::write(memstream, data.p.Pos[1]);
    vtkh::write(memstream, data.p.Pos[2]);
    vtkh::write(memstream, data.p.ID);
    vtkh::write(memstream, data.p.Status);
    vtkh::write(memstream, data.p.NumSteps);
    vtkh::write(memstream, data.p.Time);
    vtkh::write(memstream, data.blockIds);
  }

  static void read(MemStream &memstream, vtkh::Particle &data)
  {
    vtkh::read(memstream, data.p.Pos[0]);
    vtkh::read(memstream, data.p.Pos[1]);
    vtkh::read(memstream, data.p.Pos[2]);
    vtkh::read(memstream, data.p.ID);
    vtkh::read(memstream, data.p.Status);
    vtkh::read(memstream, data.p.NumSteps);
    vtkh::read(memstream, data.p.Time);
    vtkh::read(memstream, data.blockIds);
  }
};
} //namespace vtkh
#endif
