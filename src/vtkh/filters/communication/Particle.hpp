#ifndef VTK_H_PARTICLE_HPP
#define VTK_H_PARTICLE_HPP

#include <iostream>
#include <vector>
#include <vtkm/Types.h>

#include <vtkh/filters/communication/MemStream.h>
namespace vtkh
{

class Particle
{
public:
    enum Status {ACTIVE, TERMINATE, OUTOFBOUNDS};

    Particle() : coords(), id(-1), nSteps(0), status(ACTIVE), blockId(-1) {}
    Particle(const std::vector<float> &c, int _id) : coords(c[0],c[1],c[2]), id(_id), nSteps(0), status(ACTIVE), blockId(-1) {}
    Particle(const std::vector<double> &c, int _id) : coords(c[0],c[1],c[2]), id(_id), nSteps(0), status(ACTIVE), blockId(-1) {}
    Particle(const double *c, int _id) : coords(c[0],c[1],c[2]), id(_id), nSteps(0), status(ACTIVE), blockId(-1) {}
    Particle(const float *c, int _id) : coords(c[0],c[1],c[2]), id(_id), nSteps(0), status(ACTIVE), blockId(-1) {}
    Particle(const vtkm::Vec<double,3> &c, int _id) : coords(c), id(_id), nSteps(0), status(ACTIVE), blockId(-1) {}
    Particle(const vtkm::Vec<float,3> &c, int _id) : coords(c[0],c[1],c[2]), id(_id), nSteps(0), status(ACTIVE), blockId(-1) {}
    Particle(const Particle &p) : coords(p.coords), nSteps(p.nSteps), id(p.id), status(p.status), blockId(p.blockId) {}

    vtkm::Vec<double,3> coords;
    int id, nSteps;
    int blockId;
    Status status;
};

inline std::ostream &operator<<(std::ostream &os, const Particle &p)
{
    os<<"P_"<<p.id<<": ["<<p.coords[0]<<" "<<p.coords[1]<<" "<<p.coords[2]<<"] #s "<<p.nSteps<<" ";
    if (p.status == Particle::ACTIVE) os<<"ACTIVE";
    else if (p.status == Particle::TERMINATE) os<<"TERM";
    else if (p.status == Particle::OUTOFBOUNDS) os<<"OOB";
    os<<" bid = "<<p.blockId;
    return os;
}

template<>
struct Serialization<Particle>
{
  static void write(MemStream &memstream, const Particle &data)
  {
    vtkh::write(memstream, data.coords[0]);
    vtkh::write(memstream, data.coords[1]);
    vtkh::write(memstream, data.coords[2]);
    vtkh::write(memstream, data.id);
    vtkh::write(memstream, data.nSteps);
    vtkh::write(memstream, data.status);
    vtkh::write(memstream, data.blockId);
  }

  static void read(MemStream &memstream, Particle &data)
  {
    vtkh::read(memstream, data.coords[0]);
    vtkh::read(memstream, data.coords[1]);
    vtkh::read(memstream, data.coords[2]);
    vtkh::read(memstream, data.id);
    vtkh::read(memstream, data.nSteps);
    vtkh::read(memstream, data.status);
    vtkh::read(memstream, data.blockId);
  }
};
} //namespace vtkh
#endif
