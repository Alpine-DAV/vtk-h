#include "Log.hpp"
#include <vtkh/Error.hpp>

#include <vtkm/Math.h>
#include <vtkm/worklet/WorkletMapField.h>

namespace vtkh 
{

namespace detail
{
class LogField : public vtkm::worklet::WorkletMapField
{
public:
  VTKM_CONT
  LogField()
  {}

  typedef void ControlSignature(FieldIn<Scalar>, FieldOut<>);
  typedef void ExecutionSignature(_1, _2);
  
  template<typename T>
  VTKM_EXEC
  void operator()(const T &value, vtkm::Float32& log_value) const
  {
    vtkm::Float32 f_value = static_cast<vtkm::Float32>(value);
    log_value = vtkm::Log(f_value);
  }
}; //class SliceField 

} // namespace detail

Log::Log()
{
}

Log::~Log()
{

}

void 
Log::SetField(const std::string &field_name)
{
  m_field_name = field_name;
}

void 
Log::SetResultField(const std::string &field_name)
{
  m_result_name = field_name;
}


std::string
Log::GetField() const
{
  return m_field_name;
}

std::string
Log::GetResultField() const
{
  return m_result_name;
}

void Log::PreExecute() 
{
  Filter::PreExecute();
  if(m_result_name== "")
  {
    m_result_name= "log(" + m_field_name + ")";
  }
   
  vtkm::Range scalar_range = m_input->GetGlobalRange(m_field_name).GetPortalControl().Get(0);
  if(scalar_range.Min <= 0.f)
  {
    std::stringstream msg;
    msg<<"Log : error cannot perform log on field with negative values ";
    msg<<scalar_range;
    throw Error(msg.str());
  }
}

void Log::PostExecute()
{
  Filter::PostExecute();
}

void Log::DoExecute()
{
  
  this->m_output = new DataSet();
  // shallow copy input data set and bump internal ref counts
  *m_output = *m_input;

  const int num_domains = this->m_input->GetNumberOfDomains(); 

  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::cont::DataSet &dom =  this->m_output->GetDomain(i);

    if(!dom.HasField(m_field_name))
    {
      continue;
    }

    vtkm::cont::Field::Association in_assoc = dom.GetField(m_field_name).GetAssociation(); 
    bool is_cell_assoc = in_assoc == vtkm::cont::Field::Association::CELL_SET; 
    bool is_point_assoc = in_assoc == vtkm::cont::Field::Association::POINTS; 

    if(!is_cell_assoc && !is_point_assoc)
    {
      throw Error("Log: input field must be zonal or nodal");
    }


    vtkm::cont::ArrayHandle<vtkm::Float32> log_field;
    vtkm::cont::Field in_field = dom.GetField(m_field_name);
    vtkm::worklet::DispatcherMapField<detail::LogField>()
      .Invoke(in_field, log_field);
  
    if(is_cell_assoc)
    {
      vtkm::cont::Field out_field(m_result_name,
                                  in_assoc,
                                  in_field.GetAssocCellSet(),
                                  log_field);
      dom.AddField(out_field);
    }
    else
    {
      vtkm::cont::Field out_field(m_result_name,
                                  in_assoc,
                                  log_field);
      dom.AddField(out_field);
    }
  }
}

std::string
Log::GetName() const
{
  return "vtkh::Log";
}

} //  namespace vtkh
