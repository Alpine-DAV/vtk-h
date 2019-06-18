#include <vtkm/worklet/FieldHistogram.h>
#include <vtkh/filters/Histogram.hpp>
#include <vtkh/utils/vtkm_array_utils.hpp>

#ifdef VTKH_PARALLEL
#include <mpi.h>
#endif
namespace vtkh
{

namespace detail
{

class InitBins : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  InitBins()
  {}

  typedef void ControlSignature(FieldOut);
  typedef void ExecutionSignature(_1);

  template<typename T>
  VTKM_EXEC
  void operator()(T &value) const
  {
    value = 0;
  }
}; //class InitBins

class AccBins : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  AccBins()
  {}

  typedef void ControlSignature(FieldIn, FieldInOut);
  typedef void ExecutionSignature(_1, _2);

  template<typename T>
  VTKM_EXEC
  void operator()(const T &value, T &acc) const
  {
    acc += value;
  }
}; //class InitBins

} // namespace detail

Histogram::Histogram()
  : m_num_bins(256)
{

}

Histogram::~Histogram()
{

}

void
Histogram::SetField(const std::string &field_name)
{
  m_field_name = field_name;
}

void Histogram::PreExecute()
{
  Filter::PreExecute();
  Filter::CheckForRequiredField(m_field_name);
}

void Histogram::PostExecute()
{
  Filter::PostExecute();
}
struct HistRes
{
  vtkm::cont::ArrayHandle<vtkm::Id> bins;
  vtkm::Int32 num_bins;
  vtkm::Float64 min;
  vtkm::Float64 max;
  vtkm::Float64 delta;
};

struct RunHist
{
  HistRes *result;
  RunHist(HistRes *res)
    : result(res)
  {

  }

  template<typename T, typename S>
  VTKM_CONT void operator()(const vtkm::cont::ArrayHandle<T,S> &array) const
  {
    vtkm::worklet::FieldHistogram hist;
    T delta;
    T min = static_cast<T>(result->min);
    T max = static_cast<T>(result->max);
    hist.Run(array,
             result->num_bins,
             min,
             max,
             delta,
             result->bins);
    result->delta = static_cast<vtkm::Float64>(delta);
  }
};

void Histogram::DoExecute()
{
  this->m_output = new DataSet();
  const int num_domains = this->m_input->GetNumberOfDomains();

  vtkm::Range scalar_range = m_input->GetGlobalRange(m_field_name).GetPortalControl().Get(0);
  vtkm::Float64 delta = -1.f;
  vtkm::cont::ArrayHandle<vtkm::Int32> local_bins;
  local_bins.Allocate(m_num_bins);
  vtkm::worklet::DispatcherMapField<detail::InitBins>()
    .Invoke(local_bins);

  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);
    HistRes histres;
    histres.num_bins = m_num_bins;
    histres.min = scalar_range.Min;
    histres.max = scalar_range.Max;

    dom.GetField(m_field_name).
        GetData().
        ResetTypes(vtkm::TypeListTagFieldScalar()).CastAndCall(RunHist(&histres));

    vtkm::worklet::DispatcherMapField<detail::AccBins>()
      .Invoke(histres.bins, local_bins);

  }
#ifdef VTKH_PARALLEL
  MPI_Comm mpi_comm = MPI_Comm_f2c(vtkh::GetMPICommHandle());

  vtkm::cont::ArrayHandle<vtkm::Int32> global_bins;
  global_bins.Allocate(m_num_bins);

  int *local_ptr = GetVTKMPointer(local_bins);
  int *global_ptr = GetVTKMPointer(global_bins);
  MPI_Allreduce( local_ptr, global_ptr, m_num_bins, MPI_INT, MPI_SUM, mpi_comm);

  vtkm::cont::Field field("histogram",
                          vtkm::cont::Field::Association::WHOLE_MESH,
                          global_bins);
  vtkm::cont::DataSet res;
  res.AddField(field);
  m_output->AddDomain(res, vtkh::GetMPIRank());

#else
  vtkm::cont::Field field("histogram",
                          vtkm::cont::Field::Association::WHOLE_MESH,
                          local_bins);
  vtkm::cont::DataSet res;
  res.AddField(field);
  m_output->AddDomain(res, 0);
#endif
}

std::string
Histogram::GetName() const
{
  return "vtkh::Histogram";
}

} //  namespace vtkh
