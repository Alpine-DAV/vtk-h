#ifndef VTK_H_MUTEX_HPP
#define VTK_H_MUTEX_HPP

#ifdef ENABLE_OPENMP
#include <omp.h>
#else
#include <thread>
#include <mutex>
#endif

namespace vtkh
{

//Mutex class for both openmp and std::mutex
class Mutex
{

//openMP version
#ifdef ENABLE_OPENMP
public:
  Mutex() { omp_init_lock(&lock); }
  ~Mutex() { omp_destroy_lock(&lock); }
  void Lock() {omp_set_lock(&lock);}
  void Unlock() {omp_unset_lock(&lock);}
private:
  omp_lock_t lock;

//std::mutex version
#else
public:
  Mutex() {}
  ~Mutex() {}
  void Lock() {lock.lock();}
  void Unlock() {lock.unlock();}
private:
  std::mutex lock;
};

#endif


} //namespace vtkh

#endif //VTK_H_MUTEX_HPP
