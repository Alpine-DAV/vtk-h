#include "MemStream.h"

using namespace std;

namespace vtkh
{

MemStream::MemStream(size_t sz0)
{
    _pos = 0; _len = 0;
    _maxLen = _len;
    _data = NULL;
    CheckSize(sz0);
}

MemStream::MemStream(size_t sz, const unsigned char *buff)
{
    _pos = 0;
    _len = sz;
    _maxLen = _len;

    _data = new unsigned char[_len];
    memcpy(_data, buff, _len);
}

MemStream::MemStream(const MemStream &s)
{
    _pos = 0;
    _len = s.len();
    _maxLen = _len;
    _data = new unsigned char[_len];
    memcpy(_data, s.data(), _len);
}

MemStream::~MemStream()
{
    ClearMemStream();
}

void
MemStream::ClearMemStream()
{
    if (_data)
    {
        delete [] _data;
        _data = NULL;
    }
    _pos = 0;
    _len = 0;
    _maxLen = 0;
}

void
MemStream::CheckSize(size_t sz)
{
    size_t reqLen = _pos+sz;

    if (reqLen > _maxLen)
    {
        size_t newLen = 2*_maxLen; // double current size.
        if (newLen < reqLen)
            newLen = reqLen;

        unsigned char *newData = new unsigned char[newLen];

        if (_data)
        {
            memcpy(newData, _data, _len); // copy existing data to new buffer.
            delete [] _data;
        }
        _data = newData;
        _maxLen = newLen;
    }
}

void
MemStream::SaveFile( const char *filename )
{
    FILE *fp = fopen( filename, "wb" );

    if( fp )
    {
        fwrite( &_len, sizeof(_len), 1, fp );
        fwrite( _data, sizeof(_data[0]), _len, fp );

        fflush( fp );
        fclose( fp );
    }
}

void
MemStream::LoadFile( const char *filename )
{
    FILE *fp = fopen( filename, "rb" );

    if( fp )
    {
        int res = 0;
        ClearMemStream();

        res = fread( &_len, sizeof(_len), 1, fp );
        if (res != sizeof(_len))
        {
            cerr << "Bad read of MemStream from " << filename << endl;
        }

        CheckSize( _len );
        res = fread( _data, sizeof(_data[0]), _len, fp );
        if ((size_t)res != sizeof(_data[0])*_len)
        {
            cerr << "Bad read of MemStream from " << filename << endl;
        }
        fclose( fp );
    }
}

} // namespace vtkh
