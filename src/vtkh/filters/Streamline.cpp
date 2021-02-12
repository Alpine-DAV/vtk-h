#include <iostream>
#include <vtkh/vtkm_filters/vtkmStreamline.hpp>
#include <vtkh/filters/Streamline.hpp>
#include <vtkh/vtkh.hpp>
#include <vtkh/Error.hpp>

namespace vtkh
{

Streamline::Streamline()
{
}

Streamline::~Streamline()
{

}

void Streamline::PreExecute()
{
  Filter::PreExecute();
  Filter::CheckForRequiredField(m_field_name);
}

void Streamline::PostExecute()
{
  Filter::PostExecute();
}

void Streamline::DoExecute()
{
  vtkmStreamline streamlineFilter;

  this->m_output = new DataSet();
  const int num_domains = this->m_input->GetNumberOfDomains();

  vtkm::cont::PartitionedDataSet inputs;
  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);
    if(dom.HasField(m_field_name))
    {
      using vectorField_d = vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float64, 3>>;
      using vectorField_f = vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float32, 3>>;
      auto field = dom.GetField(m_field_name).GetData();
      if(!field.IsType<vectorField_d>() && !field.IsType<vectorField_f>())
      {
        throw Error("Vector field type does not match <vtkm::Vec<vtkm::Float32,3>> or <vtkm::Vec<vtkm::Float64,3>>");
      }
    }
    else
    {
      throw Error("Domain does not contain specified vector field for Streamline analysis.");
    }

    inputs.AppendPartition(dom);
  }

  auto out = streamlineFilter.Run(inputs, m_field_name, m_step_size, m_num_steps, m_seeds);

  for (vtkm::Id i = 0; i < out.GetNumberOfPartitions(); i++)
  {
    this->m_output->AddDomain(out.GetPartition(i), i);
  }
}

} //  namespace vtkh
