#include <vtkh/filters/GhostStripper.hpp>
#include <vtkh/utils/vtkm_dataset_info.hpp>
#include <vtkh/vtkm_filters/vtkmThreshold.hpp>
#include <vtkh/vtkm_filters/vtkmCleanGrid.hpp>
#include <vtkh/vtkm_filters/vtkmExtractStructured.hpp>

#include <vtkm/worklet/DispatcherMapField.h>
#include <vtkm/cont/Algorithm.h>

#include <limits>

namespace vtkh
{

namespace detail
{

template<int DIMS>
VTKM_EXEC_CONT
vtkm::Vec<vtkm::Id,3> get_logical(const vtkm::Id &index, const vtkm::Vec<vtkm::Id,3> &cell_dims);

template<>
VTKM_EXEC_CONT
vtkm::Vec<vtkm::Id,3> get_logical<3>(const vtkm::Id &index, const vtkm::Vec<vtkm::Id,3> &cell_dims)
{
  vtkm::Vec<vtkm::Id,3> res(0,0,0);
  res[0] = index % cell_dims[0];
  res[1] = (index / (cell_dims[0])) % (cell_dims[1]);
  res[2] = index / ((cell_dims[0]) * (cell_dims[1]));
  return res;
}

template<>
VTKM_EXEC_CONT
vtkm::Vec<vtkm::Id,3> get_logical<2>(const vtkm::Id &index, const vtkm::Vec<vtkm::Id,3> &cell_dims)
{
  vtkm::Vec<vtkm::Id,3> res(0,0,0);
  res[0] = index % cell_dims[0];
  res[1] = index / cell_dims[0];
  return res;
}

template<>
VTKM_EXEC_CONT
vtkm::Vec<vtkm::Id,3> get_logical<1>(const vtkm::Id &index, const vtkm::Vec<vtkm::Id,3> &cell_dims)
{
  vtkm::Vec<vtkm::Id,3> res(0,0,0);
  res[0] = index;
  return res;
}

template<int DIMS>
class RealMinMax : public vtkm::worklet::WorkletMapField
{
protected:
  vtkm::Vec<vtkm::Id,3> m_cell_dims;
  vtkm::Int32 m_min_value;
  vtkm::Int32 m_max_value;
public:
  VTKM_CONT
  RealMinMax(vtkm::Vec<vtkm::Id,3> cell_dims, vtkm::Int32 min_value, vtkm::Int32 max_value)
    : m_cell_dims(cell_dims),
      m_min_value(min_value),
      m_max_value(max_value)
  {
  }

  typedef void ControlSignature(FieldIn, AtomicArrayInOut);
  typedef void ExecutionSignature(_1, WorkIndex, _2);

  template<typename Atomic>
  VTKM_EXEC void Max(Atomic &boom,
                     const vtkm::Int32 &val,
                     const vtkm::Id &index) const
  {
    vtkm::Int32 old = -1;
    do
    {
      old = boom.CompareAndSwap(index, val, old);
    }
    while (old < val);
  }

  template<typename Atomic>
  VTKM_EXEC void Min(Atomic &boom,
                     const vtkm::Int32 &val,
                     const vtkm::Id &index) const
  {
    vtkm::Int32 old = 1000000000;
    do
    {
      old = boom.CompareAndSwap(index, val, old);
    }
    while (old > val);
  }

  template<typename T, typename AtomicType>
  VTKM_EXEC
  void operator()(const T &value, const vtkm::Id &index, AtomicType &boom) const
  {
    // we are finding the logical min max of valid zones
    if( value < m_min_value || value > m_max_value) return;

    vtkm::Vec<vtkm::Id,3> logical = get_logical<DIMS>(index, m_cell_dims);

    Min(boom, logical[0], 0);
    Min(boom, logical[1], 1);
    Min(boom, logical[2], 2);

    Max(boom, logical[0], 3);
    Max(boom, logical[1], 4);
    Max(boom, logical[2], 5);

  }
}; //class TheRealMinMax

template<int DIMS>
class Validate : public vtkm::worklet::WorkletMapField
{
protected:
  vtkm::Vec<vtkm::Id,3> m_cell_dims;
  vtkm::Int32 m_min_value;
  vtkm::Int32 m_max_value;
  vtkm::Vec<vtkm::Id,3> m_valid_min;
  vtkm::Vec<vtkm::Id,3> m_valid_max;
public:
  VTKM_CONT
  Validate(vtkm::Vec<vtkm::Id,3> cell_dims,
           vtkm::Int32 min_value,
           vtkm::Int32 max_value,
           vtkm::Vec<vtkm::Id,3> valid_min,
           vtkm::Vec<vtkm::Id,3> valid_max)
    : m_cell_dims(cell_dims),
      m_min_value(min_value),
      m_max_value(max_value),
      m_valid_min(valid_min),
      m_valid_max(valid_max)
  {
  }

  typedef void ControlSignature(FieldIn, FieldOut);
  typedef void ExecutionSignature(_1, WorkIndex, _2);

  template<typename T>
  VTKM_EXEC
  void operator()(const T &value, const vtkm::Id &index, vtkm::UInt8 &valid) const
  {
    valid = 0; // this is a valid zone
    // we are validating if non-valid cells fall completely outside
    // the min max range of valid cells
    if( value >= m_min_value || value <= m_max_value) return;

    vtkm::Vec<vtkm::Id,3> logical = get_logical<DIMS>(index, m_cell_dims);
    for(vtkm::Int32 i = 0; i < DIMS; ++i)
    {
      if(logical[i] >= m_valid_min[i] || logical[i] <= m_valid_max[i])
      {
        valid = 1;
      }
    }

  }
}; //class Validate

template<int DIMS>
bool CanStrip(vtkm::cont::Field  &ghost_field,
              const vtkm::Int32 min_value,
              const vtkm::Int32 max_value,
              vtkm::Vec<vtkm::Id,3> &min,
              vtkm::Vec<vtkm::Id,3> &max,
              vtkm::Vec<vtkm::Id,3> cell_dims,
              vtkm::Id size)
{
  vtkm::cont::ArrayHandle<vtkm::Id> minmax;
  minmax.Allocate(6);
  minmax.GetPortalControl().Set(0,std::numeric_limits<int>::max());
  minmax.GetPortalControl().Set(1,std::numeric_limits<int>::max());
  minmax.GetPortalControl().Set(2,std::numeric_limits<int>::max());
  minmax.GetPortalControl().Set(3,std::numeric_limits<int>::min());
  minmax.GetPortalControl().Set(4,std::numeric_limits<int>::min());
  minmax.GetPortalControl().Set(5,std::numeric_limits<int>::min());

  vtkm::worklet::DispatcherMapField<RealMinMax<3>>(RealMinMax<3>(cell_dims, min_value, max_value))
     .Invoke(ghost_field.GetData().ResetTypes(vtkm::TypeListTagScalarAll()), minmax);

  vtkm::Vec<vtkm::Id, 3> valid_min, valid_max;
  valid_min[0] = minmax.GetPortalConstControl().Get(0);
  valid_min[1] = minmax.GetPortalConstControl().Get(1);
  valid_min[2] = minmax.GetPortalConstControl().Get(2);

  valid_max[0] = minmax.GetPortalConstControl().Get(3);
  valid_max[1] = minmax.GetPortalConstControl().Get(4);
  valid_max[2] = minmax.GetPortalConstControl().Get(5);

  vtkm::cont::ArrayHandle<vtkm::UInt8> valid_flags;
  valid_flags.Allocate(size);

  min = valid_min;
  max = valid_max;

  vtkm::worklet::DispatcherMapField<Validate<DIMS>>(Validate<DIMS>(cell_dims,
                                                             min_value,
                                                             max_value,
                                                             valid_min,
                                                             valid_max))
     .Invoke(ghost_field.GetData().ResetTypes(vtkm::TypeListTagScalarAll()), valid_flags);

  vtkm::UInt8 res = vtkm::cont::Algorithm::Reduce(valid_flags, vtkm::UInt8(0), vtkm::Maximum());
  return res == 0;
}

bool StructuredStrip(vtkm::cont::DataSet &dataset,
                     vtkm::cont::Field   &ghost_field,
                     const vtkm::Int32 min_value,
                     const vtkm::Int32 max_value,
                     vtkm::Vec<vtkm::Id,3> &min,
                     vtkm::Vec<vtkm::Id,3> &max)
{
  vtkm::cont::DynamicCellSet cell_set = dataset.GetCellSet();
  int dims[3];
  VTKMDataSetInfo::GetPointDims(cell_set, dims);
  vtkm::Vec<vtkm::Id,3> cell_dims(0,0,0);


  bool can_strip = false;
  vtkm::Id size = 0;
  if(cell_set.IsSameType(vtkm::cont::CellSetStructured<1>()))
  {
    cell_dims[0] = dims[0] - 1;
    size = cell_dims[0];

    can_strip = CanStrip<1>(ghost_field,
                            min_value,
                            max_value,
                            min,
                            max,
                            cell_dims,
                            size);
  }
  else if(cell_set.IsSameType(vtkm::cont::CellSetStructured<2>()))
  {
    cell_dims[0] = dims[0] - 1;
    cell_dims[1] = dims[1] - 1;
    size = cell_dims[0] * cell_dims[1];

    can_strip = CanStrip<2>(ghost_field,
                            min_value,
                            max_value,
                            min,
                            max,
                            cell_dims,
                            size);
  }
  else if(cell_set.IsSameType(vtkm::cont::CellSetStructured<3>()))
  {
    cell_dims[0] = dims[0] - 1;
    cell_dims[1] = dims[1] - 1;
    cell_dims[2] = dims[2] - 1;
    size = cell_dims[0] * cell_dims[1] * cell_dims[2];

    can_strip = CanStrip<3>(ghost_field,
                            min_value,
                            max_value,
                            min,
                            max,
                            cell_dims,
                            size);
  }

  return can_strip;
}

} // namespace detail

GhostStripper::GhostStripper()
  : m_min_value(0),  // default to real zones only
    m_max_value(0)   // 0 = real, 1 = valid ghost, 2 = garbage ghost
{

}

GhostStripper::~GhostStripper()
{

}

void
GhostStripper::SetField(const std::string &field_name)
{
  m_field_name = field_name;
}

void
GhostStripper::SetMinValue(const vtkm::Int32 min_value)
{
  m_min_value = min_value;
}

void
GhostStripper::SetMaxValue(const vtkm::Int32 max_value)
{
  m_max_value = max_value;
}

void GhostStripper::PreExecute()
{
  Filter::PreExecute();
  Filter::CheckForRequiredField(m_field_name);
}

void GhostStripper::PostExecute()
{
  Filter::PostExecute();
}

void GhostStripper::DoExecute()
{
  this->m_output = new DataSet();

  const int num_domains = this->m_input->GetNumberOfDomains();

  for(int i = 0; i < num_domains; ++i)
  {

    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);

    if(!dom.HasField(m_field_name))
    {
      continue;
    }

    vtkm::cont::Field field = dom.GetField(m_field_name);

    int topo_dims = 0;
    bool do_threshold = true;

    if(VTKMDataSetInfo::IsStructured(dom, topo_dims))
    {
      vtkm::Vec<vtkm::Id,3> min, max;
      bool can_strip = detail::StructuredStrip(dom, field, m_min_value, m_max_value, min, max);
      if(can_strip)
      {
        do_threshold = false;
        vtkm::RangeId3 range(min[0],max[0]+2, min[1], max[1]+2, min[2], max[2]+2);
        vtkm::Id3 sample(1, 1, 1);

        vtkh::vtkmExtractStructured extract;
        auto output = extract.Run(dom,
                                  range,
                                  sample,
                                  this->GetFieldSelection());

        m_output->AddDomain(output, domain_id);
      }

    }

    if(do_threshold)
    {
      vtkmThreshold thresholder;
      auto tout = thresholder.Run(dom,
                                  m_field_name,
                                  m_min_value,
                                  m_max_value,
                                  this->GetFieldSelection());

      vtkh::vtkmCleanGrid cleaner;
      auto clout = cleaner.Run(dom, this->GetFieldSelection());
      m_output->AddDomain(clout, domain_id);
    }

  }

}

std::string
GhostStripper::GetName() const
{
  return "vtkh::GhostStripper";
}

} //  namespace vtkh
