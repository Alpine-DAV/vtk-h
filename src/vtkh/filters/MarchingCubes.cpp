#include <vtkh/filters/MarchingCubes.hpp>
#include <vtkh/filters/CleanGrid.hpp>
#include <vtkh/filters/Recenter.hpp>
#include <vtkh/vtkm_filters/vtkmMarchingCubes.hpp>

namespace vtkh
{

MarchingCubes::MarchingCubes()
 : m_levels(10)
{

}

MarchingCubes::~MarchingCubes()
{

}

void
MarchingCubes::SetIsoValue(const double &iso_value)
{
  m_iso_values.clear();
  m_iso_values.push_back(iso_value);
  m_levels = -1;
}

void
MarchingCubes::SetLevels(const int &levels)
{
  m_iso_values.clear();
  m_levels = levels;
}

void
MarchingCubes::SetIsoValues(const double *iso_values, const int &num_values)
{
  assert(num_values > 0);
  m_iso_values.clear();
  for(int i = 0; i < num_values; ++i)
  {
    m_iso_values.push_back(iso_values[i]);
  }
  m_levels = -1;
}

void
MarchingCubes::SetField(const std::string &field_name)
{
  m_field_name = field_name;
}

void MarchingCubes::PreExecute()
{
  Filter::PreExecute();
  Filter::CheckForRequiredField(m_field_name);

  if(m_levels != -1)
  {
    vtkm::Range scalar_range = m_input->GetGlobalRange(m_field_name).GetPortalControl().Get(0);
    float length = scalar_range.Length();
    float step = length / (m_levels + 1.f);

    m_iso_values.clear();
    for(int i = 1; i <= m_levels; ++i)
    {
      float iso = scalar_range.Min + float(i) * step;
      m_iso_values.push_back(iso);
    }
  }

  assert(m_iso_values.size() > 0);
}

void MarchingCubes::PostExecute()
{
  Filter::PostExecute();
}

void MarchingCubes::DoExecute()
{
  DataSet temp_data;
  vtkh::DataSet *old_input = this->m_input;


  // make sure we have a node-centered field
  bool valid_field = false;
  bool is_cell_assoc = m_input->GetFieldAssociation(m_field_name, valid_field) ==
                       vtkm::cont::Field::Association::CELL_SET;
  bool delete_input = false;
  if(valid_field && is_cell_assoc)
  {
    Recenter recenter;
    recenter.SetInput(m_input);
    recenter.SetField(m_field_name);
    recenter.SetResultAssoc(vtkm::cont::Field::Association::POINTS);
    recenter.Update();
    m_input = recenter.GetOutput();
    delete_input = true;
  }

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

    vtkh::vtkmMarchingCubes marcher;

    auto dataset = marcher.Run(dom,
                               m_field_name,
                               m_iso_values,
                               this->GetFieldSelection());

    temp_data.AddDomain(dataset, domain_id);

  }

  CleanGrid cleaner;
  cleaner.SetInput(&temp_data);
  cleaner.Update();
  this->m_output = cleaner.GetOutput();

  if(delete_input)
  {
    delete m_input;
    this->m_input = old_input;
  }
}

std::string
MarchingCubes::GetName() const
{
  return "vtkh::MarchingCubes";
}

} //  namespace vtkh
