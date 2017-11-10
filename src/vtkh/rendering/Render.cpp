#include "Render.hpp"
#include <vtkm/rendering/MapperRayTracer.h>
#include <vtkm/rendering/View2D.h>
#include <vtkm/rendering/View3D.h>

namespace vtkh 
{

Render::Render()
  : m_color_table("cool2warm"),
    m_has_color_table(false)
{
}

Render::~Render()
{
}

Render::vtkmCanvasPtr 
Render::GetDomainCanvas(const vtkm::Id &domain_id)
{
  vtkm::Id dom = -1;
  for(size_t i = 0; i < m_domain_ids.size(); ++i)
  {
    if(m_domain_ids[i] == domain_id)
    {
      dom = i;
      break;
    }
  }

  if(dom == -1)
  {
    std::stringstream ss;
    ss<<"Render: canvas with domain id "<< domain_id <<" not found ";
    throw Error(ss.str());
  }
  return m_canvases[dom];
}

Render::vtkmCanvasPtr 
Render::GetCanvas(const vtkm::Id index)
{
  assert(index >= 0 && index < m_canvases.size());
  return m_canvases[index];
}

vtkm::Bounds
Render::GetSceneBounds() const
{
  return m_scene_bounds;
}

void
Render::SetSceneBounds(const vtkm::Bounds &bounds) 
{
  m_scene_bounds = bounds;
}

vtkm::Range
Render::GetScalarRange() const
{
  return m_scalar_range;
}

void
Render::SetScalarRange(const vtkm::Range &range) 
{
  m_scalar_range = range;
}

void 
Render::AddCanvas(vtkmCanvasPtr canvas, vtkm::Id domain_id)
{
  m_canvases.push_back(canvas);
  m_domain_ids.push_back(domain_id);
}

int 
Render::GetNumberOfCanvases() const
{
  return static_cast<int>(m_canvases.size());
}

bool 
Render::HasCanvas(const vtkm::Id &domain_id) const 
{
  vtkm::Id dom = -1;
  for(size_t i = 0; i < m_domain_ids.size(); ++i)
  {
    if(m_domain_ids[i] == domain_id)
    {
      dom = i;
      break;
    }
  }

  return dom != -1;
}

const vtkm::rendering::Camera& 
Render::GetCamera() const
{ 
  return m_camera;
}

void 
Render::SetCamera(const vtkm::rendering::Camera &camera)
{ 
   m_camera = camera;
}

void 
Render::SetImageName(const std::string &name)
{
  m_image_name = name;
}

std::string 
Render::GetImageName() const 
{
  return m_image_name;
}

void 
Render::SetColorTable(const vtkm::rendering::ColorTable &color_table)
{
  m_color_table = color_table;
  m_has_color_table = true;
}

bool
Render::HasColorTable() const
{
  return m_has_color_table;
}

vtkm::rendering::ColorTable
Render::GetColorTable() const
{
  return m_color_table;
}

void 
Render::RenderWorldAnnotations()
{
  int size = m_canvases.size(); 
  if(size < 1) return;

  bool is3d = m_camera.GetMode() == vtkm::rendering::Camera::MODE_3D;
  if( !is3d ) return; // no world annotations for 2d
  // create a a dummy data set so we can create a vtkm::view 
  // to handle the annotations for us
  vtkm::rendering::MapperRayTracer dummy_mapper;
  vtkm::Vec<vtkm::Float32,3> origin(m_scene_bounds.X.Min,
                                    m_scene_bounds.Y.Min,
                                    m_scene_bounds.Z.Min);

  vtkm::Vec<vtkm::Float32,3> spacing(m_scene_bounds.X.Max - m_scene_bounds.X.Min,
                                     m_scene_bounds.Y.Max - m_scene_bounds.Y.Min,
                                     m_scene_bounds.Z.Max - m_scene_bounds.Z.Min);

  vtkm::Id3 dims(2,2,2);
  vtkm::cont::CoordinateSystem coords("coords", dims, origin, spacing);
  vtkm::cont::DynamicCellSet cells;
  vtkm::cont::ArrayHandle<vtkm::Float32> range;
  range.Allocate(2);

  auto portal = range.GetPortalControl();
  portal.Set(0, m_scalar_range.Min);
  portal.Set(1, m_scalar_range.Max);
  vtkm::cont::Field field("field_name", vtkm::cont::Field::ASSOC_POINTS, range);

  for(int i = 0; i < size; ++i)
  {
    m_canvases[i]->SetViewToWorldSpace(m_camera, true);
    vtkm::rendering::Scene scene;
    vtkm::rendering::Actor actor(cells, coords, field, m_color_table);
    scene.AddActor(actor);
    vtkm::rendering::View3D view(scene, dummy_mapper, *m_canvases[i], m_camera, m_canvases[i]->GetBackgroundColor());
    view.RenderWorldAnnotations();
  }
    
}

void 
Render::RenderScreenAnnotations()
{
  int size = m_canvases.size(); 
  if(size < 1) return;

  bool is3d = m_camera.GetMode() == vtkm::rendering::Camera::MODE_3D;
  // create a a dummy data set so we can create a vtkm::view 
  // to handle the annotations for us
  vtkm::rendering::MapperRayTracer dummy_mapper;
  vtkm::Vec<vtkm::Float32,3> origin(m_scene_bounds.X.Min,
                                    m_scene_bounds.Y.Min,
                                    m_scene_bounds.Z.Min);

  vtkm::Vec<vtkm::Float32,3> spacing(m_scene_bounds.X.Max - m_scene_bounds.X.Min,
                                     m_scene_bounds.Y.Max - m_scene_bounds.Y.Min,
                                     m_scene_bounds.Z.Max - m_scene_bounds.Z.Min);

  vtkm::Id3 dims(2,2,2);
  vtkm::cont::CoordinateSystem coords("coords", dims, origin, spacing);
  vtkm::cont::DynamicCellSet cells;
  vtkm::cont::ArrayHandle<vtkm::Float32> range;
  range.Allocate(2);

  auto portal = range.GetPortalControl();
  portal.Set(0, m_scalar_range.Min);
  portal.Set(1, m_scalar_range.Max);
  vtkm::cont::Field field("field_name", vtkm::cont::Field::ASSOC_POINTS, range);

  for(int i = 0; i < size; ++i)
  {
    m_canvases[i]->SetViewToScreenSpace(m_camera, false);
    vtkm::rendering::Scene scene;
    vtkm::rendering::Actor actor(cells, coords, field, m_color_table);
    scene.AddActor(actor);
    if(is3d)
    {
      vtkm::rendering::View3D view(scene, 
                                   dummy_mapper, 
                                   *m_canvases[i], 
                                   m_camera, 
                                   m_canvases[i]->GetBackgroundColor());
      view.RenderScreenAnnotations();
    }
    else
    {
      vtkm::rendering::View2D view(scene, 
                                   dummy_mapper, 
                                   *m_canvases[i], 
                                   m_camera, 
                                   m_canvases[i]->GetBackgroundColor());
      view.RenderScreenAnnotations();
    }
  }
}

} // namespace vtkh
