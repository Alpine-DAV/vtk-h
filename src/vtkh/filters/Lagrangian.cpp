#include <iostream>
#include <vtkm/filter/Lagrangian.h>
#include <vtkh/filters/Lagrangian.hpp>
#include <vtkh/vtkh.hpp>

namespace vtkh 
{

Lagrangian::Lagrangian()
{

}

Lagrangian::~Lagrangian()
{

}

void 
Lagrangian::SetField(const std::string &field_name)
{
  m_field_name = field_name;
}

void 
Lagrangian::SetStepSize(const double &step_size)
{
  m_step_size = step_size;
}

void 
Lagrangian::SetWriteFrequency(const int &write_frequency)
{
  m_write_frequency = write_frequency;
}

void 
Lagrangian::SetCustomSeedResolution(const int &cust_res)
{
	m_cust_res = cust_res;
}

void 
Lagrangian::SetSeedResolutionInX(const int &x_res)
{
	m_x_res = x_res;
}

void 
Lagrangian::SetSeedResolutionInY(const int &y_res)
{
	m_y_res = y_res;
}
void 
Lagrangian::SetSeedResolutionInZ(const int &z_res)
{
	m_z_res = z_res;
}


void Lagrangian::PreExecute() 
{
  Filter::PreExecute();
}

void Lagrangian::PostExecute()
{
  Filter::PostExecute();
}

void Lagrangian::DoExecute()
{
	vtkm::filter::Lagrangian lagrangianFilter;
  lagrangianFilter.SetStepSize(m_step_size);
	lagrangianFilter.SetWriteFrequency(m_write_frequency);
	lagrangianFilter.SetRank(vtkh::GetMPIRank());
	lagrangianFilter.SetActiveField(m_field_name);
	lagrangianFilter.SetCustomSeedResolution(m_cust_res);
	lagrangianFilter.SetSeedResolutionInX(m_x_res);
	lagrangianFilter.SetSeedResolutionInY(m_y_res);
	lagrangianFilter.SetSeedResolutionInZ(m_z_res);
	
	this->m_output = new DataSet();
  const int num_domains = this->m_input->GetNumberOfDomains();
  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);
		
		if(!dom.HasField(m_field_name))
		{
			// Cloverleaf3D has a velocity field stored as velocity_x, velocity_y, velocity_z. Composite these vectors.
			if(dom.HasField("velocity_x") && dom.HasField("velocity_y") && dom.HasField("velocity_z"))
			{
				vtkm::cont::Field x_field = dom.GetField("velocity_x");
				vtkm::cont::Field y_field = dom.GetField("velocity_y");
				vtkm::cont::Field z_field = dom.GetField("velocity_z");

				vtkm::cont::DynamicArrayHandle xIn = x_field.GetData();	
				vtkm::cont::DynamicArrayHandle yIn = y_field.GetData();	
				vtkm::cont::DynamicArrayHandle zIn = z_field.GetData();	
			
				vtkm::cont::ArrayHandle<vtkm::Float64> xData;
				vtkm::cont::ArrayHandle<vtkm::Float64> yData;
				vtkm::cont::ArrayHandle<vtkm::Float64> zData;

				xIn.CopyTo(xData);
				yIn.CopyTo(yData);
				zIn.CopyTo(zData);
			
				vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float64,3>> velocityField;
				auto composite = vtkm::cont::make_ArrayHandleCompositeVector(xData, yData, zData);
				vtkm::cont::ArrayCopy(composite, velocityField);
				
				vtkm::cont::Field velocity(m_field_name, vtkm::cont::Field::Association::POINTS, velocityField);
				dom.AddField(velocity);
			}
      else if(dom.HasField("vector_data"))
      {
				vtkm::cont::Field vel_field = dom.GetField("vector_data");
				vtkm::cont::DynamicArrayHandle velIn = vel_field.GetData();	
				vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float64,3>> velocityField;
				velIn.CopyTo(velocityField);
				vtkm::cont::Field velocity(m_field_name, vtkm::cont::Field::Association::POINTS, velocityField);
				dom.AddField(velocity);
      }
			else
			{
				std::cout << "Lagrangian filter not called: velocity field names do not match." << std::endl;
				continue;
			}
		}
		std::cout << "Lagrangian filter call on rank: " << vtkh::GetMPIRank() << std::endl;
		vtkm::cont::DataSet extractedBasis = lagrangianFilter.Execute(dom);
		std::cout << "Returned from Lagrangian filter call on rank: " << vtkh::GetMPIRank() << std::endl;
    m_output->AddDomain(extractedBasis, domain_id);
  }
}

std::string
Lagrangian::GetName() const
{
  return "vtkh::Lagrangian";
}

} //  namespace vtkh
