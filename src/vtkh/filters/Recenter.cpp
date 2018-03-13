#include <vtkh/filters/Recenter.hpp>
#include <vtkm/filter/PointAverage.h>
#include <vtkm/filter/CellAverage.h>

namespace vtkh 
{

Recenter::Recenter()
{

}

Recenter::~Recenter()
{

}

void 
Recenter::SetField(const std::string &field_name)
{
  m_field_name = field_name;
}

void Recenter::PreExecute() 
{
  Filter::PreExecute();
  assert(m_field_name != "");
}

void Recenter::PostExecute()
{
  Filter::PostExecute();
}

void Recenter::DoExecute()
{
  this->m_output = new DataSet();
  const int num_domains = this->m_input->GetNumberOfDomains();
  
  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);
   
    vtkm::cont::DataSet out_data, temp;
    // Since there is no way to remove a field from a dataset
    // we have to iterate over the data set to create a shallow
    // copy of everything else
    
    const vtkm::Id num_fields = dom.GetNumberOfFields();

    for(vtkm::Id f = 0; f < num_fields; ++f)
    {
      vtkm::cont::Field field = dom.GetField(f);
      if(field.GetName() != m_field_name)
      {
        out_data.AddField(field); 
      }
      else
      {
        temp.AddField(field);
      }
    }
    
    const vtkm::Id num_coords = dom.GetNumberOfCoordinateSystems();

    for(vtkm::Id f = 0; f < num_coords; ++f)
    {
      vtkm::cont::CoordinateSystem coords= dom.GetCoordinateSystem(f);
      out_data.AddCoordinateSystem(coords); 
      temp.AddCoordinateSystem(coords);
    }
    
    const vtkm::Id num_cellsets= dom.GetNumberOfCellSets();

    for(vtkm::Id f = 0; f < num_cellsets; ++f)
    {
      vtkm::cont::DynamicCellSet cellset = dom.GetCellSet(f);
      out_data.AddCellSet(cellset); 
      temp.AddCellSet(cellset);
    }

    if(temp.HasField(m_field_name))
    {
      bool isCellAssoc = temp.GetField(m_field_name).GetAssociation() == 
                         vtkm::cont::Field::ASSOC_CELL_SET; 

      vtkm::cont::DataSet dataset;
      std::string out_name = m_field_name + "_out";

      if(isCellAssoc)
      {
        vtkm::filter::PointAverage avg;
        avg.SetOutputFieldName(out_name);
        avg.SetActiveField(m_field_name);
        dataset = avg.Execute(dom);
      }
      else
      {
        vtkm::filter::CellAverage avg;
        avg.SetOutputFieldName(out_name);
        avg.SetActiveField(m_field_name);
        dataset = avg.Execute(dom);
      }

      vtkm::cont::Field recentered_field;
      recentered_field = vtkm::cont::Field(m_field_name, 
                                           dataset.GetField(out_name).GetAssociation(),
                                           dataset.GetField(out_name).GetData());
      out_data.AddField(recentered_field);
      
    }
    
    m_output->AddDomain(out_data, domain_id);
  }
}

std::string
Recenter::GetName() const
{
  return "vtkh::Recenter";
}

} //  namespace vtkh
