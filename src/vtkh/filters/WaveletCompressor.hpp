#ifndef VTK_H_WAVELET_COMPRESSOR_HPP
#define VTK_H_WAVELET_COMPRESSOR_HPP

#include <vtkh/vtkh.hpp>
#include <vtkh/filters/Filter.hpp>
#include <vtkh/DataSet.hpp>

#include <memory>

namespace vtkh
{

class WaveletCompressor : public Filter
{
public:
  WaveletCompressor(); 
  virtual ~WaveletCompressor(); 

protected:
  void PreExecute() override;
  void PostExecute() override;
  void DoExecute() override;

  struct InternalsType;
  std::shared_ptr<InternalsType> m_internals;
};

} //namespace vtkh
#endif
