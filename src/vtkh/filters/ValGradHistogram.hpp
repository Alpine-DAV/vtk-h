#ifndef VTK_H_VALGRADHISTOGRAM_HPP
#define VTK_H_VALGRADHISTOGRAM_HPP

#include <vtkh/vtkh.hpp>
#include <vtkh/vtkh_exports.h>
#include <vtkh/filters/Filter.hpp>
#include <vtkh/DataSet.hpp>
#include <vtkh/filters/Histogram.hpp>

#include <vector>
#include <iostream>

namespace vtkh
{

class VTKH_API ValGradHistogram
{
public:
  ValGradHistogram();
  virtual ~ValGradHistogram();

  struct SparseHistogramResult
  {
    std::vector<vtkm::cont::ArrayHandle<vtkm::Id>> m_binIds;
    vtkm::cont::ArrayHandle<vtkm::Id> m_freqs;
    int m_nob;
    
    vtkm::Range m_range_val, m_range_grad;
    vtkm::Float64 m_bin_delta_val, m_bin_delta_grad;
    void Print(std::ostream &out);
  };

  struct NonSparseHistogramResult
  {
    //this is a flattened 2d histogram
    vtkm::cont::ArrayHandle<vtkm::Id> m_freqs;
    int m_nob;
    
    vtkm::Range m_range_val, m_range_grad;
    vtkm::Float64 m_bin_delta_val, m_bin_delta_grad;
    void Print(std::ostream &out);
    Histogram::HistogramResult getValDist();
    Histogram::HistogramResult getGradDist();
    Histogram::HistogramResult getConditionalGradDist(int val_bin_id);
    vtkm::Id totalCount();
  };

  NonSparseHistogramResult Run(vtkh::DataSet &data_set, const std::string &field_name, vtkh::DataSet &grad_mag_ds);
  void SetNumBins(const int num_bins);
protected:
  int m_num_bins;
};

} //namespace vtkh
#endif
