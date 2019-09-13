#include "Clip.hpp"

#include <vtkh/filters/CleanGrid.hpp>
#include <vtkh/vtkm_filters/vtkmClip.hpp>

namespace vtkh
{
struct Clip::InternalsType
{
  vtkm::cont::ImplicitFunctionHandle m_func;
  InternalsType()
  {}
};

Clip::Clip()
  : m_internals(new InternalsType),
    m_invert(false)
{

}

Clip::~Clip()
{

}

void
Clip::SetInvertClip(bool invert)
{
  m_invert = invert;
}

void
Clip::SetBoxClip(const vtkm::Bounds &clipping_bounds)
{
   auto box =  vtkm::cont::make_ImplicitFunctionHandle(
                 vtkm::Box({ clipping_bounds.X.Min,
                             clipping_bounds.Y.Min,
                             clipping_bounds.Z.Min},
                           { clipping_bounds.X.Max,
                             clipping_bounds.Y.Max,
                             clipping_bounds.Z.Max}));


  m_internals->m_func = box;
}

void
Clip::SetSphereClip(const double center[3], const double radius)
{
  vtkm::Vec<vtkm::FloatDefault,3> vec_center;
  vec_center[0] = center[0];
  vec_center[1] = center[1];
  vec_center[2] = center[2];
  vtkm::FloatDefault r = radius;

  auto sphere = vtkm::cont::make_ImplicitFunctionHandle(vtkm::Sphere(vec_center, r));
  m_internals->m_func = sphere;
}

void
Clip::SetPlaneClip(const double origin[3], const double normal[3])
{
  vtkm::Vec<vtkm::FloatDefault,3> vec_origin;
  vec_origin[0] = origin[0];
  vec_origin[1] = origin[1];
  vec_origin[2] = origin[2];

  vtkm::Vec<vtkm::FloatDefault,3> vec_normal;
  vec_normal[0] = normal[0];
  vec_normal[1] = normal[1];
  vec_normal[2] = normal[2];

  auto plane = vtkm::cont::make_ImplicitFunctionHandle(vtkm::Plane(vec_origin, vec_normal));
  m_internals->m_func = plane;
}

void Clip::PreExecute()
{
  Filter::PreExecute();
}

void Clip::PostExecute()
{
  Filter::PostExecute();
}

void Clip::DoExecute()
{

  DataSet data_set;

  const int num_domains = this->m_input->GetNumberOfDomains();
  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::Id domain_id;
    vtkm::cont::DataSet dom;
    this->m_input->GetDomain(i, dom, domain_id);

    vtkh::vtkmClip clipper;

    auto dataset = clipper.Run(dom,
                               m_internals->m_func,
                               m_invert,
                               this->GetFieldSelection());

    data_set.AddDomain(dataset, domain_id);
  }

  CleanGrid cleaner;
  cleaner.SetInput(&data_set);
  cleaner.Update();
  this->m_output = cleaner.GetOutput();
}

std::string
Clip::GetName() const
{
  return "vtkh::Clip";
}

} //  namespace vtkh
