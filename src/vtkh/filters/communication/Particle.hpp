#ifndef VTK_H_PARTICLE_HPP
#define VTK_H_PARTICLE_HPP

#include <iostream>
#include <vector>
#include <vtkm/Types.h>

#include <vtkh/filters/communication/MemStream.h>

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
    memstream.write(data.coords[0]);
    memstream.write(data.coords[1]);
    memstream.write(data.coords[2]);
    memstream.write(data.id);
    memstream.write(data.nSteps);
    memstream.write(data.status);
    memstream.write(data.blockId);
  }

  static void read(MemStream &memstream, Particle &data)
  {
    memstream.read(data.coords[0]);
    memstream.read(data.coords[1]);
    memstream.read(data.coords[2]);
    memstream.read(data.id);
    memstream.read(data.nSteps);
    memstream.read(data.status);
    memstream.read(data.blockId);
  }
};

#endif
