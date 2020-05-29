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
    m_shading(true),
    m_canvas(m_width, m_height)
{
}

Render::~Render()
{
}

Render::vtkmCanvas&
Render::GetCanvas()
{
  return m_canvas;
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
  if(width == m_width) return;
  m_width = width;
  m_canvas.ResizeBuffers(m_width, m_height);
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
  if(height == m_height) return;
  m_height = height;
  m_canvas.ResizeBuffers(m_width, m_height);
}

void
Render::SetSceneBounds(const vtkm::Bounds &bounds)
{
  m_scene_bounds = bounds;
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
#ifdef VTKH_PARALLEL
  if(vtkh::GetMPIRank() != 0) return;
#endif

  Annotator annotator(m_canvas, m_camera, m_scene_bounds);
  annotator.RenderWorldAnnotations();

}

void
Render::RenderScreenAnnotations(const std::vector<std::string> &field_names,
                                const std::vector<vtkm::Range> &ranges,
                                const std::vector<vtkm::cont::ColorTable> &colors)
{
  if(!m_render_annotations) return;

  if(m_render_background) m_canvas.BlendBackground();
  Annotator annotator(m_canvas, m_camera, m_scene_bounds);
  annotator.RenderScreenAnnotations(field_names, ranges, colors);
}

void
Render::RenderBackground()
{
  if(m_render_background) m_canvas.BlendBackground();
}

Render::vtkmCanvas
Render::CreateCanvas()
{
  Render::vtkmCanvas canvas(m_width, m_height);
  canvas.SetBackgroundColor(m_bg_color);
  canvas.SetForegroundColor(m_fg_color);
  canvas.Clear();
  return canvas;
}

void
Render::Save()
{
  // After rendering and compositing
  // Rank 0 contains the complete image.
#ifdef VTKH_PARALLEL
  if(vtkh::GetMPIRank() != 0) return;
#endif
  float* color_buffer = &GetVTKMPointer(m_canvas.GetColorBuffer())[0][0];
  int height = m_canvas.GetHeight();
  int width = m_canvas.GetWidth();
  PNGEncoder encoder;
  encoder.Encode(color_buffer, width, height);
  encoder.Save(m_image_name + ".png");
}

vtkh::Render
MakeRender(int width,
           int height,
           vtkm::Bounds scene_bounds,
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

  return render;
}

vtkh::Render
MakeRender(int width,
           int height,
           vtkm::Bounds scene_bounds,
           vtkm::rendering::Camera camera,
           const std::string &image_name,
           float bg_color[4],
           float fg_color[4])
{
  vtkh::Render render;
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

  return render;
}
} // namespace vtkh
