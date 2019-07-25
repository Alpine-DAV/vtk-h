#ifndef VTK_H_PARTICLE_HPP
#define VTK_H_PARTICLE_HPP

#include <iostream>
#include <vector>
#include <vtkm/Types.h>

#include <vtkh/filters/communication/MemStream.h>
#include <vtkh/utils/StreamUtil.hpp>

namespace vtkh
{

class Particle
{
public:
    enum Status {ACTIVE, TERMINATE, OUTOFBOUNDS, WRONG_DOMAIN};

    Particle() : coords(), id(-1), nSteps(0), status(ACTIVE), blockIds() {}
    Particle(const std::vector<float> &c, int _id) : coords(c[0],c[1],c[2]), id(_id), nSteps(0), status(ACTIVE), blockIds() {}
    Particle(const std::vector<double> &c, int _id) : coords(c[0],c[1],c[2]), id(_id), nSteps(0), status(ACTIVE), blockIds() {}
    Particle(const double *c, int _id) : coords(c[0],c[1],c[2]), id(_id), nSteps(0), status(ACTIVE), blockIds() {}
    Particle(const float *c, int _id) : coords(c[0],c[1],c[2]), id(_id), nSteps(0), status(ACTIVE), blockIds() {}
    Particle(const vtkm::Vec<double,3> &c, int _id) : coords(c), id(_id), nSteps(0), status(ACTIVE), blockIds() {}
    Particle(const vtkm::Vec<float,3> &c, int _id) : coords(c[0],c[1],c[2]), id(_id), nSteps(0), status(ACTIVE), blockIds() {}
    Particle(const Particle &p) : coords(p.coords), nSteps(p.nSteps), id(p.id), status(p.status), blockIds(p.blockIds) {}

    vtkm::Vec<double,3> coords;
    int id, nSteps;
    std::vector<int> blockIds;
    Status status;

    friend std::ostream &operator<<(std::ostream &os, const vtkh::Particle p)
    {
        os<<"(";
        os<<"P_"<<p.id<<" ";
        os<<"["<<p.coords[0]<<" "<<p.coords[1]<<" "<<p.coords[2]<<"] #s "<<p.nSteps<<" ";
        if (p.status == Particle::ACTIVE) os<<"ACTIVE";
        else if (p.status == Particle::TERMINATE) os<<"TERM";
        else if (p.status == Particle::OUTOFBOUNDS) os<<"OOB";
        else if (p.status == Particle::WRONG_DOMAIN) os<<"WRONG_DOMAIN";
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
    vtkh::write(memstream, data.coords[0]);
    vtkh::write(memstream, data.coords[1]);
    vtkh::write(memstream, data.coords[2]);
    vtkh::write(memstream, data.id);
    vtkh::write(memstream, data.nSteps);
    vtkh::write(memstream, data.status);
    vtkh::write(memstream, data.blockIds);
  }

  static void read(MemStream &memstream, vtkh::Particle &data)
  {
    vtkh::read(memstream, data.coords[0]);
    vtkh::read(memstream, data.coords[1]);
    vtkh::read(memstream, data.coords[2]);
    vtkh::read(memstream, data.id);
    vtkh::read(memstream, data.nSteps);
    vtkh::read(memstream, data.status);
    vtkh::read(memstream, data.blockIds);
  }
};
} //namespace vtkh
#endif
