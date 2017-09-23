#include "WaveletCompressor.hpp"

#include <vtkh/Error.hpp>
#include <vtkm/worklet/WaveletCompressor.h>

namespace vtkh
{

struct WaveletCompressor::InternalsType
{

};

WaveletCompressor::WaveletCompressor()
  : m_internals(new InternalsType)
{

}

WaveletCompressor::~WaveletCompressor()
{

}

void
WaveletCompressor::PreExecute()
{
  int topo_dims;
}

void
WaveletCompressor::DoExecute()
{

}

void
WaveletCompressor::PostExecute()
{

}

} // namespace vtkh
