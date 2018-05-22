#include <vtkm/filter/Lagrangian.h>
#include <vtkh/filters/Lagrangian.hpp>

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
	vtkm::Float32 stepSize = 0.01;
  lagrangianFilter.SetStepSize(stepSize);
	lagrangianFilter.SetWriteFrequency(10);
	lagrangianFilter.SetActiveField("velocity");

	this->m_output = new DataSet();
  const int num_domains = this->m_input->GetNumberOfDomains();
  
  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);
		
		if(!dom.HasField("velocity"))
		{
			// Cloverleaf3D has vectors stored in the form xvec, yvec, zvec. Composite these vectors.
			if(dom.HasField("xvec") && dom.HasField("yvec") && dom.HasField("zvec"))
			{
				vtkm::cont::Field x_field = dom.GetField("xvec");
				vtkm::cont::Field y_field = dom.GetField("yvec");
				vtkm::cont::Field z_field = dom.GetField("zvec");

				vtkm::cont::DynamicArrayHandle xIn = x_field.GetData();	
				vtkm::cont::DynamicArrayHandle yIn = x_field.GetData();	
				vtkm::cont::DynamicArrayHandle zIn = x_field.GetData();	
			
				vtkm::cont::ArrayHandle<vtkm::Float64> xData;
				vtkm::cont::ArrayHandle<vtkm::Float64> yData;
				vtkm::cont::ArrayHandle<vtkm::Float64> zData;

				xIn.CopyTo(xData);
				yIn.CopyTo(yData);
				zIn.CopyTo(zData);
			
				vtkm::cont::ArrayHandle<vtkm::Vec<vtkm::Float64,3>> velocityField;
				auto composite = vtkm::cont::make_ArrayHandleCompositeVector(xData, 0, yData, 0, zData, 0);
				vtkm::cont::ArrayCopy(composite, velocityField);

				vtkm::cont::Field velocity("velocity", vtkm::cont::Field::ASSOC_POINTS, velocityField);
				dom.AddField(velocity);
			}
			else
			{
				continue;
			}
		}
		vtkm::cont::DataSet extractedBasis = lagrangianFilter.Execute(dom);
    m_output->AddDomain(extractedBasis, domain_id);
  }
}

std::string
Lagrangian::GetName() const
{
  return "vtkh::Lagrangian";
}

} //  namespace vtkh
