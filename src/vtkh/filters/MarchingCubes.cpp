#include <vtkh/filters/MarchingCubes.hpp>
#include <vtkh/filters/CleanGrid.hpp>
#include <vtkm/filter/MarchingCubes.h>

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
  if(m_levels != -1) {
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
  assert(m_field_name != "");
  if(m_map_fields.size() == 0)
  {
    this->MapAllFields(); 
  }
}

void MarchingCubes::PostExecute()
{

}

bool 
MarchingCubes::ContainsIsoValues(vtkm::cont::DataSet &dom)
{
  vtkm::cont::Field field = dom.GetField(m_field_name);
  vtkm::cont::ArrayHandle<vtkm::Range> ranges = field.GetRange();
  assert(ranges.GetNumberOfValues() == 1);
  vtkm::Range range = ranges.GetPortalControl().Get(0);
  for(size_t i = 0; i < m_iso_values.size(); ++i)
  {
    double val = m_iso_values[i];
    if(range.Contains(val))
    {
      return true;
    }
  }

  return false;

}

void MarchingCubes::DoExecute()
{
  this->m_output = new DataSet();
  vtkm::filter::MarchingCubes marcher;

  marcher.SetIsoValues(m_iso_values);
  marcher.SetMergeDuplicatePoints(true); 
  const int num_domains = this->m_input->GetNumberOfDomains(); 
  int valid = 0;
  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);
    bool valid_domain = ContainsIsoValues(dom);
    if(!valid_domain)
    {
      // vtkm does not like it if we ask it to contour
      // values that do not exist in the field, so
      // we have to check.
      continue;
    }
    valid++;
    vtkm::filter::Result res = marcher.Execute(dom, m_field_name);
    for(size_t f = 0; f < m_map_fields.size(); ++f)
    {
      marcher.MapFieldOntoOutput(res, dom.GetField(m_map_fields[f]));
    }
    m_output->AddDomain(res.GetDataSet(), domain_id);
  }

  this->PropagateMetadata();
  
}

} //  namespace vtkh
