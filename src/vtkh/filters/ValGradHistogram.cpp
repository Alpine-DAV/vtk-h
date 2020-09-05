#include <vtkh/filters/ValGradHistogram.hpp>
#include <vtkh/filters/Histogram.hpp>

#include <vtkh/Error.hpp>
#include <vtkh/Logger.hpp>
#include <vtkh/utils/vtkm_array_utils.hpp>

#include <vtkm/worklet/FieldHistogram.h>
#include <vtkm/worklet/NDimsHistogram.h>

#ifdef VTKH_PARALLEL
#include <mpi.h>
#endif

namespace vtkh
{

namespace valgrad_detail
{

ValGradHistogram::NonSparseHistogramResult
merge_histograms(std::vector<ValGradHistogram::SparseHistogramResult> &histograms)
{
  const int size = histograms.size();
  if(size < 1)
  {
    throw Error("2d histogram size is 0");
  }
  ValGradHistogram::NonSparseHistogramResult res;

  //initialize the frequencies
  std::vector<vtkm::Id> freq_vec(histograms[0].m_nob*histograms[0].m_nob, 0);

  for(int h=0; h<size; ++h)
  {
    vtkm::Id nonSparseBins = histograms[h].m_binIds[0].ReadPortal().GetNumberOfValues();
    for (int i = 0; i < nonSparseBins; i++)
    {
      vtkm::Id idx0 = histograms[h].m_binIds[0].ReadPortal().Get(i);
      vtkm::Id idx1 = histograms[h].m_binIds[1].ReadPortal().Get(i);
      vtkm::Id f = histograms[h].m_freqs.ReadPortal().Get(i);

      freq_vec[idx0*histograms[0].m_nob + idx1] += f;
    }
  }

  res.m_freqs = vtkm::cont::make_ArrayHandle(freq_vec, vtkm::CopyFlag::On);

  res.m_nob = histograms[0].m_nob;
  res.m_range_val = histograms[0].m_range_val;
  res.m_range_grad = histograms[0].m_range_grad;
  res.m_bin_delta_val = histograms[0].m_bin_delta_val;
  res.m_bin_delta_grad = histograms[0].m_bin_delta_grad;  

  return res;
}

template<typename T>
void reduce(T *array, int size);

template<>
void reduce<vtkm::Int32>(vtkm::Int32 *array, int size)
{
#ifdef VTKH_PARALLEL
  MPI_Comm mpi_comm = MPI_Comm_f2c(vtkh::GetMPICommHandle());
  MPI_Allreduce(MPI_IN_PLACE,array,size, MPI_INT,MPI_SUM,mpi_comm);
#else
  (void) array;
  (void) size;
#endif
}

template<>
void reduce<vtkm::Int64>(vtkm::Int64 *array, int size)
{
#ifdef VTKH_PARALLEL
  MPI_Comm mpi_comm = MPI_Comm_f2c(vtkh::GetMPICommHandle());
  MPI_Allreduce(MPI_IN_PLACE,array,size, MPI_LONG_LONG,MPI_SUM,mpi_comm);
#else
  (void) array;
  (void) size;
#endif
}


} // namespace valgrad_detail

ValGradHistogram::ValGradHistogram()
  : m_num_bins(256)
{

}

ValGradHistogram::~ValGradHistogram()
{

}

void
ValGradHistogram::SetNumBins(const int num_bins)
{
  m_num_bins = num_bins;
}

ValGradHistogram::NonSparseHistogramResult
ValGradHistogram::Run(vtkh::DataSet &data_set, const std::string &field_name, vtkh::DataSet &grad_mag_ds)
{
  VTKH_DATA_OPEN("histogram");
  VTKH_DATA_ADD("device", GetCurrentDevice());
  VTKH_DATA_ADD("bins", m_num_bins);
  VTKH_DATA_ADD("input_cells", data_set.GetNumberOfCells());
  VTKH_DATA_ADD("input_domains", data_set.GetNumberOfDomains());

  if(!data_set.GlobalFieldExists(field_name))
  {
    throw Error("Histogram: field '"+field_name+"' does not exist");
  }

  vtkm::Range val_range;
  vtkm::cont::ArrayHandle<vtkm::Range> val_ranges = data_set.GetGlobalRange(field_name);
  if(val_ranges.GetNumberOfValues() != 1)
  {
    throw Error("Histogram: field must have a single component");
  }
  val_range = val_ranges.ReadPortal().Get(0);

  vtkm::Range grad_range;
  vtkm::cont::ArrayHandle<vtkm::Range> grad_ranges = grad_mag_ds.GetGlobalRange("mag");
  if(grad_ranges.GetNumberOfValues() != 1)
  {
    throw Error("Grad Mag Histogram: field must have a single component");
  }
  grad_range = grad_ranges.ReadPortal().Get(0);

  // std::cout << "[new]scalar val global range: " << val_range << std::endl;
  // std::cout << "[new]gradient mag global range: " << grad_range << std::endl;

  const int num_domains = data_set.GetNumberOfDomains(); 


  std::vector<ValGradHistogram::SparseHistogramResult> local_sparse_histograms;
  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id val_domain_id, grad_domain_id;
    vtkm::cont::DataSet val_dom, grad_dom;
    
    data_set.GetDomain(i, val_dom, val_domain_id);
    grad_mag_ds.GetDomain(i, grad_dom, grad_domain_id);

    if(!val_dom.HasField(field_name)) continue;
    vtkm::cont::Field val_field = val_dom.GetField(field_name);
    
    if(!grad_dom.HasField("mag")) continue;
    vtkm::cont::Field grad_field = grad_dom.GetField("mag");
    
    

    if(val_field.GetNumberOfValues() != grad_field.GetNumberOfValues() )
    {
      throw Error("value and gradient magintude should have same number of datapoints");
    }

    vtkm::worklet::NDimsHistogram ndHistogram;
    // Set the number of data points
    ndHistogram.SetNumOfDataPoints(val_field.GetNumberOfValues());
    vtkm::Float64 delta_val;
    ndHistogram.AddField(val_field.GetData(), m_num_bins, val_range, delta_val, true);
    vtkm::Float64 delta_grad;
    ndHistogram.AddField(grad_field.GetData(), m_num_bins, grad_range, delta_grad, true);


    std::vector<vtkm::cont::ArrayHandle<vtkm::Id>> binIds;
    vtkm::cont::ArrayHandle<vtkm::Id> freqs;
    ndHistogram.Run(binIds, freqs);

    ValGradHistogram::SparseHistogramResult dom_2d_hist;
    dom_2d_hist.m_binIds = binIds;
    dom_2d_hist.m_freqs = freqs;
    dom_2d_hist.m_nob = m_num_bins;
    dom_2d_hist.m_range_val = val_range;
    dom_2d_hist.m_range_grad = grad_range;
    dom_2d_hist.m_bin_delta_val = delta_val;
    dom_2d_hist.m_bin_delta_grad = delta_grad;

    //dom_2d_hist.Print(std::cout);
    local_sparse_histograms.push_back(dom_2d_hist);



  }

  ValGradHistogram::NonSparseHistogramResult local = valgrad_detail::merge_histograms(local_sparse_histograms);

  
  vtkm::Id * bin_ptr = GetVTKMPointer(local.m_freqs);
  valgrad_detail::reduce(bin_ptr, m_num_bins*m_num_bins);

  VTKH_DATA_CLOSE();
  return local;
}

void
ValGradHistogram::SparseHistogramResult::Print(std::ostream &out)
{
  vtkm::Id nonSparseBins = m_binIds[0].ReadPortal().GetNumberOfValues();

  vtkm::Id f_sum = 0;
  for (int i = 0; i < nonSparseBins; i++)
  {
    vtkm::Id idx0 = m_binIds[0].ReadPortal().Get(i);
    vtkm::Id idx1 = m_binIds[1].ReadPortal().Get(i);
    vtkm::Id f = m_freqs.ReadPortal().Get(i);
    f_sum += f;
    out << "[" << idx0 << ", " << idx1 << "] --> " << f << "\n";
  }
  out << "Total num of points in 2d histogram = " << f_sum << "\n"; 
}

void
ValGradHistogram::NonSparseHistogramResult::Print(std::ostream &out)
{

  auto freqPortal = m_freqs.ReadPortal();  
  
  vtkm::Id f_sum = 0;
  for(int i=0; i<m_nob; ++i)
  {
    for(int j=0; j<m_nob; ++j)
    {
      out << freqPortal.Get(i*m_nob + j) << "  ";
      f_sum += freqPortal.Get(i*m_nob + j);
    }
    out << "\n";
  }

  out << "Total num of points in 2d histogram = " << f_sum << "\n"; 
}

vtkm::Id
ValGradHistogram::NonSparseHistogramResult::totalCount()
{
  auto freqPortal = m_freqs.ReadPortal();  
  vtkm::Id f_sum = 0;
  for(int i=0; i<m_nob*m_nob; ++i)
  {
   f_sum += freqPortal.Get(i);
  }
  return f_sum;
}

Histogram::HistogramResult
ValGradHistogram::NonSparseHistogramResult::getValDist()
{
  auto freqPortal = m_freqs.ReadPortal();
  std::vector<vtkm::Id> marginal_freqs(m_nob, 0);
  for(int i=0; i<m_nob; ++i)
  {
    for(int j=0; j<m_nob; ++j)
    {
      marginal_freqs[i] += freqPortal.Get(i*m_nob + j);
    }
  }

  Histogram::HistogramResult marginal_hist;
  marginal_hist.m_bins = vtkm::cont::make_ArrayHandle(marginal_freqs, vtkm::CopyFlag::On); 
  marginal_hist.m_range = m_range_val;
  marginal_hist.m_bin_delta = m_bin_delta_val;

  return marginal_hist;
}

Histogram::HistogramResult
ValGradHistogram::NonSparseHistogramResult::getGradDist()
{
  auto freqPortal = m_freqs.ReadPortal();
  std::vector<vtkm::Id> marginal_freqs(m_nob, 0);
  for(int i=0; i<m_nob; ++i)
  {
    for(int j=0; j<m_nob; ++j)
    {
      marginal_freqs[j] += freqPortal.Get(i*m_nob + j);
    }
  }

  Histogram::HistogramResult marginal_hist;
  marginal_hist.m_bins = vtkm::cont::make_ArrayHandle(marginal_freqs, vtkm::CopyFlag::On); 
  marginal_hist.m_range = m_range_grad;
  marginal_hist.m_bin_delta = m_bin_delta_grad;

  return marginal_hist;
}

Histogram::HistogramResult
ValGradHistogram::NonSparseHistogramResult::getConditionalGradDist(int val_bin_id)
{
  if(val_bin_id < 0 || val_bin_id >= m_nob)
  {
    throw Error("Conditional Grad Histogram: value bin id out of range");
  }

  auto freqPortal = m_freqs.ReadPortal();
  std::vector<vtkm::Id> conditional_freqs(m_nob, 0);
  for(int j=0; j<m_nob; ++j)
  {
    conditional_freqs[j] = freqPortal.Get(val_bin_id*m_nob + j);
  }

  Histogram::HistogramResult conditional_hist;
  conditional_hist.m_bins = vtkm::cont::make_ArrayHandle(conditional_freqs, vtkm::CopyFlag::On); 
  conditional_hist.m_range = m_range_grad;
  conditional_hist.m_bin_delta = m_bin_delta_grad;

  return conditional_hist;
}

} //  namespace vtkh
