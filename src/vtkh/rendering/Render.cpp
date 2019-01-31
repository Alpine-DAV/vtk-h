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
  : m_width(1024),
    m_height(1024),
    m_render_annotations(true),
    m_render_background(true),
    m_shading(true)
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

  if(m_canvases[dom] == nullptr)
  {
    m_canvases[dom] = this->CreateCanvas();
  }

  return m_canvases[dom];
}

Render::vtkmCanvasPtr
Render::GetCanvas(const vtkm::Id index)
{
  assert(index >= 0 && index < m_canvases.size());

  if(m_canvases[index] == nullptr)
  {
    m_canvases[index] = this->CreateCanvas();
  }
  return m_canvases[index];
}

vtkm::Bounds
Render::GetSceneBounds() const
{
  return m_scene_bounds;
}

void
Render::DoRenderAnnotations(bool on)
{
  m_render_annotations = on;
}

void
Render::DoRenderBackground(bool on)
{
  m_render_background = on;
}

vtkm::Int32
Render::GetWidth() const
{
  return m_width;
}

vtkm::Int32
Render::GetHeight() const
{
  return m_height;
}

void
Render::SetWidth(const vtkm::Int32 width)
{
  m_width = width;
}

void
Render::SetShadingOn(bool on)
{
  m_shading = on;
}

bool
Render::GetShadingOn() const
{
  return m_shading;
}

void
Render::SetHeight(const vtkm::Int32 height)
{
  m_height = height;
}

void
Render::SetSceneBounds(const vtkm::Bounds &bounds)
{
  m_scene_bounds = bounds;
}

void
Render::AddDomain(vtkm::Id domain_id)
{
  m_canvases.push_back(nullptr);
  m_domain_ids.push_back(domain_id);
}

int
Render::GetNumberOfCanvases() const
{
  return static_cast<int>(m_canvases.size());
}

void
Render::ClearCanvases()
{
  int num = static_cast<int>(m_canvases.size());
  for(int i = 0; i < num; ++i)
  {
    if(m_canvases[i] != nullptr)
    {
      m_canvases[i] = nullptr;
    }
  }
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

void
Render::SetBackgroundColor(float bg_color[4])
{
  m_bg_color.Components[0] = bg_color[0];
  m_bg_color.Components[1] = bg_color[1];
  m_bg_color.Components[2] = bg_color[2];
  m_bg_color.Components[3] = bg_color[3];
}

void
Render::SetForegroundColor(float fg_color[4])
{
  m_fg_color.Components[0] = fg_color[0];
  m_fg_color.Components[1] = fg_color[1];
  m_fg_color.Components[2] = fg_color[2];
  m_fg_color.Components[3] = fg_color[3];
}

std::string
Render::GetImageName() const
{
  return m_image_name;
}

vtkm::rendering::Color
Render::GetBackgroundColor() const
{
  return m_bg_color;
}

void
Render::RenderWorldAnnotations()
{
  if(!m_render_annotations) return;
  int size = m_canvases.size();
  if(size < 1) return;

#ifdef VTKH_PARALLEL
  if(vtkh::GetMPIRank() != 0) return;
#endif

  Annotator annotator(*m_canvases[0], m_camera, m_scene_bounds);
  annotator.RenderWorldAnnotations();

}

void
Render::RenderScreenAnnotations(const std::vector<std::string> &field_names,
                                const std::vector<vtkm::Range> &ranges,
                                const std::vector<vtkm::cont::ColorTable> &colors)
{
  if(!m_render_annotations) return;
  int size = m_canvases.size();
  if(size < 1) return;

  if(m_render_background) m_canvases[0]->BlendBackground();
  Annotator annotator(*m_canvases[0], m_camera, m_scene_bounds);
  annotator.RenderScreenAnnotations(field_names, ranges, colors);
}

void
Render::RenderBackground()
{
  int size = m_canvases.size();
  if(size < 1) return;
  if(m_render_background) m_canvases[0]->BlendBackground();
}

Render::vtkmCanvasPtr
Render::CreateCanvas()
{
  Render::vtkmCanvasPtr canvas = std::make_shared<vtkm::rendering::CanvasRayTracer>(m_width, m_height);
  canvas->SetBackgroundColor(m_bg_color);
  canvas->SetForegroundColor(m_fg_color);
  canvas->Clear();
  return canvas;
}

void
Render::Save()
{
  // After rendering and compositing
  // Rank 0 domain 0 contains the complete image.
  int size = m_canvases.size();
  if(size < 1) return;
#ifdef VTKH_PARALLEL
  if(vtkh::GetMPIRank() != 0) return;
#endif
  float* color_buffer = &GetVTKMPointer(m_canvases[0]->GetColorBuffer())[0][0];
  int height = m_canvases[0]->GetHeight();
  int width = m_canvases[0]->GetWidth();
  PNGEncoder encoder;
  encoder.Encode(color_buffer, width, height);
  encoder.Save(m_image_name + ".png");
}

vtkh::Render
MakeRender(int width,
           int height,
           vtkm::Bounds scene_bounds,
           const std::vector<vtkm::Id> &domain_ids,
           const std::string &image_name,
           float bg_color[4],
           float fg_color[4])
{
  vtkh::Render render;
  vtkm::rendering::Camera camera;
  camera.ResetToBounds(scene_bounds);
  render.SetSceneBounds(scene_bounds);
  render.SetWidth(width);
  render.SetHeight(height);
  //
  // detect a 2d data set
  //
  camera.SetModeTo3D();

  bool is_2d = scene_bounds.Z.Min == 0. && scene_bounds.Z.Max == 0.;

  if(is_2d)
  {
    camera.SetModeTo2D();
    render.SetShadingOn(false);
  }

  render.SetCamera(camera);
  render.SetImageName(image_name);
  render.SetBackgroundColor(bg_color);
  render.SetForegroundColor(fg_color);

  const size_t num_domains = domain_ids.size();

  for(size_t i = 0; i < num_domains; ++i)
  {
    //auto canvas = RendererType::GetNewCanvas(1, 1);
    //canvas->Clear();
    //canvas->SetBackgroundColor(render.GetBackgroundColor());
    render.AddDomain(domain_ids[i]);
  }

  if(num_domains == 0)
  {
    //auto canvas = RendererType::GetNewCanvas(width, height);
    //canvas->SetBackgroundColor(render.GetBackgroundColor());
    //canvas->Clear();
    vtkm::Id empty_data_set_id = -1;
    render.AddDomain(empty_data_set_id);
  }
  return render;
}

vtkh::Render
MakeRender(int width,
           int height,
           vtkm::rendering::Camera camera,
           vtkh::DataSet &data_set,
           const std::string &image_name,
           float bg_color[4],
           float fg_color[4])
{
  vtkh::Render render;
  render.SetCamera(camera);
  render.SetImageName(image_name);
  vtkm::Bounds bounds = data_set.GetGlobalBounds();
  render.SetSceneBounds(bounds);
  render.SetWidth(width);
  render.SetHeight(height);
  //
  // detect a 2d data set
  //
  bool is_2d = bounds.Z.Min == 0. && bounds.Z.Max == 0.;

  if(is_2d)
  {
    camera.SetModeTo2D();
    render.SetShadingOn(false);
  }

  render.SetBackgroundColor(bg_color);
  render.SetForegroundColor(fg_color);

  int num_domains = static_cast<int>(data_set.GetNumberOfDomains());
  for(int i = 0; i < num_domains; ++i)
  {
    vtkm::cont::DataSet ds;
    vtkm::Id domain_id;
    data_set.GetDomain(i, ds, domain_id);
    //auto canvas = RendererType::GetNewCanvas(width, height);

    //canvas->SetBackgroundColor(render.GetBackgroundColor());
    //canvas->Clear();

    render.AddDomain(domain_id);
  }

  if(num_domains == 0)
  {
    //auto canvas = RendererType::GetNewCanvas(1,1);
    //canvas->SetBackgroundColor(render.GetBackgroundColor());
    //canvas->Clear();
    vtkm::Id empty_data_set_id = -1;
    render.AddDomain(empty_data_set_id);
  }

  return render;
}
} // namespace vtkh
