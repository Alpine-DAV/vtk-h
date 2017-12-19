
#include <vtkh/filters/Slice.hpp>
#include <vtkh/Error.hpp>
#include <vtkh/filters/MarchingCubes.hpp>

#include <vtkm/VectorAnalysis.h>
#include <vtkm/cont/TryExecute.h>
#include <vtkm/worklet/WorkletMapField.h>

namespace vtkh
{

namespace detail
{
  
class SliceField : public vtkm::worklet::WorkletMapField
{
protected:
  vtkm::Vec<vtkm::Float32,3> m_point;
  vtkm::Vec<vtkm::Float32,3> m_normal;
public:
  VTKM_CONT
  SliceField(vtkm::Vec<vtkm::Float32,3> point, vtkm::Vec<vtkm::Float32,3> normal)
    : m_point(point),
      m_normal(normal)
  {
    vtkm::Normalize(m_normal);
  }

  typedef void ControlSignature(FieldIn<>, FieldOut<>);
  typedef void ExecutionSignature(_1, _2);
  
  template<typename T>
  VTKM_EXEC
  void operator()(const vtkm::Vec<T,3> &point, vtkm::Float32& distance) const
  {
    vtkm::Vec<vtkm::Float32,3> f_point(point[0], point[1], point[2]);
    distance = vtkm::dot(m_point - f_point, m_normal);
  }
}; //class SliceField 

struct SliceCaller 
{

  template <typename Device>
  VTKM_CONT bool operator()(Device, 
                            const vtkm::cont::CoordinateSystem &coords,
                            vtkm::cont::ArrayHandle<vtkm::Float32> &output,
                            vtkm::Vec<vtkm::Float32, 3> point,
                            vtkm::Vec<vtkm::Float32, 3> normal) const
  {
    VTKM_IS_DEVICE_ADAPTER_TAG(Device);
    vtkm::worklet::DispatcherMapField<SliceField, Device>(SliceField(point, normal))
      .Invoke(coords.GetData(), output);
    return true;
  }
};
}

Slice::Slice()
  : m_point(20.f, 20.1f, 20.f),
    m_normal(.0f, 0.0f, 1.f)
{

}

Slice::~Slice()
{

}

void 
Slice::SetPlane(vtkm::Vec<vtkm::Float32,3> point, vtkm::Vec<vtkm::Float32,3> normal)
{
  m_point = point;
  m_normal = normal;
}

void
Slice::PreExecute()
{

}

void
Slice::DoExecute()
{
  const std::string fname = "slice_field";
  const int num_domains = this->m_input->GetNumberOfDomains(); 
  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet &dom = this->m_input->GetDomain(i);

    vtkm::cont::ArrayHandle<vtkm::Float32> slice_field;
    vtkm::cont::TryExecute(detail::SliceCaller(), dom.GetCoordinateSystem(), slice_field, m_point, m_normal);
     
    dom.AddField(vtkm::cont::Field(fname,
                                    vtkm::cont::Field::ASSOC_POINTS,
                                    slice_field));
  }

  vtkh::MarchingCubes marcher;
  marcher.SetInput(this->m_input);
  marcher.SetIsoValue(0.);
  marcher.SetField(fname);
  marcher.Update();
  this->m_output = marcher.GetOutput();
}

void
Slice::PostExecute()
{

}

} // namespace vtkh
