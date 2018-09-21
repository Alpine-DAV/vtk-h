#include <iostream>
#include <vtkm/worklet/ParticleAdvection.h>
#include <vtkh/filters/ParticleAdvection.hpp>
#include <vtkh/vtkh.hpp>
#include <vtkh/Error.hpp>

namespace vtkh
{

ParticleAdvection::ParticleAdvection()
{

}

ParticleAdvection::~ParticleAdvection()
{

}

void
ParticleAdvection::SetField(const std::string &field_name)
{
  m_field_name = field_name;
}

void
ParticleAdvection::SetStepSize(const double &step_size)
{
  m_step_size = step_size;
}

void
ParticleAdvection::SetWriteFrequency(const int &write_frequency)
{
  m_write_frequency = write_frequency;
}

void
ParticleAdvection::SetCustomSeedResolution(const int &cust_res)
{
	m_cust_res = cust_res;
}

void
ParticleAdvection::SetSeedResolutionInX(const int &x_res)
{
	m_x_res = x_res;
}

void
ParticleAdvection::SetSeedResolutionInY(const int &y_res)
{
	m_y_res = y_res;
}
void
ParticleAdvection::SetSeedResolutionInZ(const int &z_res)
{
	m_z_res = z_res;
}


void ParticleAdvection::PreExecute()
{
  Filter::PreExecute();
}

void ParticleAdvection::PostExecute()
{
  Filter::PostExecute();
}

void ParticleAdvection::DoExecute()
{
  vtkm::worklet::ParticleAdvection particleAdvectionWorklet;

  this->m_output = new DataSet();
  const int num_domains = this->m_input->GetNumberOfDomains();

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
      if(!field.IsSameType(vectorField_d()) && !field.IsSameType(vectorField_f()))
      {
        throw Error("Vector field type does not match <vtkm::Vec<vtkm::Float32,3>> or <vtkm::Vec<vtkm::Float64,3>>");
      }
    }
    else
    {
      throw Error("Domain does not contain specified vector field for ParticleAdvection analysis.");
    }

    vtkm::cont::DataSet res;
    m_output->AddDomain(res, domain_id);
  }
}

std::string
ParticleAdvection::GetName() const
{
  return "vtkh::ParticleAdvection";
}

} //  namespace vtkh
