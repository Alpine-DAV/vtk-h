#include "Render.hpp"
#include <vtkh/rendering/Annotator.hpp>
#include <vtkh/utils/PNGEncoder.hpp>
#include <vtkh/utils/vtkm_array_utils.hpp>
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

  for(int i = 0; i < size; ++i)
  {

    Annotator annotator(*m_canvases[i], m_camera, m_scene_bounds);
    annotator.RenderWorldAnnotations();
  }
    
}

void 
Render::RenderScreenAnnotations()
{
  int size = m_canvases.size(); 
  if(size < 1) return;
  
  Annotator annotator(*m_canvases[0], m_camera, m_scene_bounds);
  annotator.RenderScreenAnnotations("field_name", m_scalar_range, m_color_table);
}

void
Render::Save()
{
  // After rendering and compositing 
  // Rank 0 domain 0 contains the complete image.
  int size = m_canvases.size(); 
  if(size < 1) return;
#ifdef PARALLEL
  if(vtkh::GetMPIRank() != 0) return;
#endif
  
  float* color_buffer = &GetVTKMPointer(m_canvases[0]->GetColorBuffer())[0][0]; 
  int height = m_canvases[0]->GetHeight(); 
  int width = m_canvases[0]->GetWidth(); 
  PNGEncoder encoder;
  encoder.Encode(color_buffer, width, height);
  encoder.Save(m_image_name + ".png");
}

} // namespace vtkh
