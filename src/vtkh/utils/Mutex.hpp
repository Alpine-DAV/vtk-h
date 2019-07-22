#ifndef VTK_H_MUTEX_HPP
#define VTK_H_MUTEX_HPP
#include <memory>

namespace vtkh
{

//Mutex class for both openmp and std::mutex
class Mutex
{

public:
  Mutex();
  ~Mutex();
  void Lock();
  void Unlock();
private:
  struct InternalsType;
  std::shared_ptr<InternalsType> m_internals;
};

} //namespace vtkh

#endif //VTK_H_MUTEX_HPP
