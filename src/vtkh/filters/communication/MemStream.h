#ifndef VTKH_MEM_STREAM_H
#define VTKH_MEM_STREAM_H

#include <iostream>
#include <string>
#include <vector>
#include <list>
#include <cstring>

#include <vtkh/vtkh_exports.h>

namespace vtkh
{

class VTKH_API MemStream
{
  public:
    enum Mode
    {
        READ = 0,
        WRITE
    };

    MemStream(size_t sz0= 32);
    MemStream(size_t sz, const unsigned char *buff);
    MemStream(const MemStream &s);
    ~MemStream();

    void rewind() {_pos = 0;}
    size_t pos() const {return _pos;}
    void setPos(size_t p);
    size_t len() const { return _len; }
    size_t capacity() const { return _maxLen; }
    unsigned char *data() const { return _data; }

    // General read/write routines.
    template <typename T> void io(Mode mode, T *pt, size_t num) {return (mode == READ ? read(pt,num) : write(pt,num));}
    template <typename T> void io(Mode mode, T& t) {size_t s=1; io( mode, &t, s );}
    template <typename T> void io(Mode mode, std::vector<T> &v)  {return (mode == READ ? read(v) : write(v));}
    template <typename T> void io(Mode mode, std::list<T> &l)  {return (mode == READ ? read(l) : write(l));}

    //Read from buffer.

    void read_binary(unsigned char *data, const size_t &size);
    //void read(std::string &str);

    //Write to buffer.
    void write_binary(const unsigned char *data, size_t size);

    void SaveFile( const char *filename );
    void LoadFile( const char *filename );
    void ClearMemStream();

  private:
    // data members
    unsigned char *_data;
    size_t _len, _maxLen, _pos;

    void CheckSize(size_t sz);

    friend std::ostream& operator<<(std::ostream &out, const MemStream &m)
    {
        out<<" MemStream(p= "<<m.pos()<<", l= "<<m.len()<<"["<<m.capacity()<<"]) data=[";
        /*
        for (size_t i=0; i < m.len(); i++)
            out<<(int)(m._data[i])<<" ";
        */
        out<<"]";
        return out;
    }
};

inline void MemStream::read_binary(unsigned char *data, const size_t &size)
{
    size_t nBytes = sizeof(unsigned char)*size;
    memcpy(data, &_data[_pos], nBytes);
    _pos += nBytes;
}

inline void MemStream::write_binary(const unsigned char *data, size_t size)
{
    size_t nBytes = sizeof(unsigned char)*size;
    CheckSize(nBytes);
    memcpy(&_data[_pos], data, nBytes );
    _pos += nBytes;

    if (_pos > _len)
        _len = _pos;
}

inline void MemStream::setPos(size_t p)
{
    _pos = p;
    if (_pos > len())
        throw "MemStream::setPos failed";
}

template<typename T>
struct Serialization
{
#if (defined(__clang__) && !defined(__ppc64__)) || (defined(__GNUC__) && __GNUC__ >= 5)
    static_assert(std::is_trivially_copyable<T>::value, "Default serialization works only for trivially copyable types");
#endif
  static void write(MemStream &memstream, const T &data)
  {
    memstream.write_binary((const unsigned char*) &data, sizeof(T));
  }
  static void read(MemStream &memstream, T &data)
  {
    memstream.read_binary((unsigned char*) &data, sizeof(T));
  }
};

template<typename T>
static void write(MemStream &memstream, const T &data)
{
  Serialization<T>::write(memstream, data);
}

template<typename T>
static void read(MemStream &memstream, T &data)
{
  Serialization<T>::read(memstream, data);
}

template<class T>
struct Serialization<std::vector<T>>
{
  static void write(MemStream &memstream, const std::vector<T> &data)
  {
    const size_t sz = data.size();
    vtkh::write(memstream, sz);
    for (size_t i = 0; i < sz; i++)
        vtkh::write(memstream, data[i]);
  }

  static void read(MemStream &memstream, std::vector<T> &data)
  {
    size_t sz;
    vtkh::read(memstream, sz);
    data.resize(sz);
    for ( size_t i = 0; i < sz; i++ )
        vtkh::read(memstream, data[i]);
  }
};

template<class T>
struct Serialization<std::list<T>>
{
  static void write(MemStream &memstream, const std::list<T> &data)
  {
    vtkh::write(memstream, data.size());
    typename std::list<T>::const_iterator it;
    for (it = data.begin(); it != data.end(); it++)
        vtkh::write(memstream, *it);
  }

  static void read(MemStream &memstream, std::list<T> &data)
  {
    size_t sz;
    vtkh::read(memstream, sz);
    for (size_t i = 0; i < sz; i++)
    {
        T v;
        vtkh::read(memstream, v);
        data.push_back(v);
    }
  }
};

//template<>
//struct Serialization<std::string>
//{
//  static void write(MemStream &memstream, const std::string &data)
//  {
//    size_t sz = data.size();
//    memstream.write(sz);
//    memstream.write(data.data(), sz);
//  }
//
//  static void read(MemStream &memstream, std::string &data)
//  {
//    size_t sz;
//    memstream.read(sz);
//    data.resize(sz);
//    memstream.read(&data[0], sz);
//  }
//};

} // namespace vtkh
#endif
