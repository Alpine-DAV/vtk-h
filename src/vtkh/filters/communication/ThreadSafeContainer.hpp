#ifndef VTK_H_THREAD_SAFE_CONTAINER_HPP
#define VTK_H_THREAD_SAFE_CONTAINER_HPP

#include <algorithm>

namespace vtkh
{

template <typename T,
          template <typename, typename> class Container,
          typename Lock>
class ThreadSafeContainer
{
public:
 ThreadSafeContainer()
 {
 }

 template <template <typename, typename> class C,
           typename Allocator=std::allocator<T>>
 ThreadSafeContainer(const C<T,Allocator> &c)
 {
     Assign(c);
 }

 ThreadSafeContainer(const ThreadSafeContainer<T, Container, Lock> &c)
 {
     if (this == &c)
         throw "ERROR: attempting to assign thread identical thread safe containers.";

     Container<T, std::allocator<T>> tmp;
     c.Get(tmp);
     Assign(tmp);
 }

 ~ThreadSafeContainer()
 {
     Clear();
 }

 bool Empty()
 {
     lock.set();
     bool val = data.empty();
     lock.unset();
     return val;
 }
 size_t Size()
 {
     lock.set();
     size_t sz = data.size();
     lock.unset();
     return sz;
 }
 void Clear()
 {
     lock.set();
     data.clear();
     lock.unset();
 }

 //Add/set elements
 template <template <typename, typename> class C,
           typename Allocator=std::allocator<T>>
 void Insert(const C<T, Allocator> &c)
 {
     if (c.empty())
         return;

     lock.set();
     data.insert(data.end(), c.begin(), c.end());
     lock.unset();
 }

 template <template <typename, typename> class C,
           typename Allocator=std::allocator<T>>
 void Assign(const C<T, Allocator> &c)
 {
     lock.set();
     data.clear();
     data.insert(data.end(), c.begin(), c.end());
     lock.unset();
 }

 template <template <typename, typename> class C,
           typename Allocator=std::allocator<T>>
 Container<T, std::allocator<T>>& operator=(const C<T, Allocator> &c)
 {
     lock.set();
     data.clear();
     data.insert(data.end(), c.begin(), c.end());
     lock.unset();
     return *this;
 }

 //Get elements
 template <template <typename, typename> class C,
           typename Allocator=std::allocator<T>>
 bool Get(C<T, Allocator> &c)
 {
     lock.set();
     c.insert(c.end(), data.begin(), data.end());
     data.clear();
     lock.unset();

     return !c.empty();
 }

 //Get elements
 template <template <typename, typename> class C,
           typename Allocator=std::allocator<T>>
 bool Get(C<T, Allocator> &c, std::size_t N)
 {
     lock.set();
     std::size_t n = std::min(N, data.size());
     c.insert(c.end(), data.begin(), data.begin()+n);
     data.erase(data.begin(), data.begin()+n);
     lock.unset();

     return !c.empty();
 }

 template <template <typename, typename> class C,
           typename Allocator=std::allocator<T>>
 void Put(C<T, Allocator> &c)
 {
     if (c.empty())
         return;

     lock.set();
     data.insert(data.end(), c.begin(), c.end());
     lock.unset();
 }

 friend std::ostream &
 operator<<(std::ostream &os, ThreadSafeContainer<T, Container, Lock> &c)
 {
    c.lock.set();
    os<<"ts_(("<<c.data<<"))";
    c.lock.unset();
    return os;
 }

protected:
 Container<T, std::allocator<T>> data;
 Lock lock;
};

};

#endif //VTK_H_THREAD_SAFE_CONTAINER_HPP
