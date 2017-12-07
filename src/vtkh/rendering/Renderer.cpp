#include "Renderer.hpp"
#include "Image.hpp"
#include "compositing/Compositor.hpp"

#include <vtkh/utils/vtkm_array_utils.hpp>
#include <vtkh/utils/vtkm_dataset_info.hpp>
#include <vtkh/utils/PNGEncoder.hpp>
#include <vtkm/rendering/raytracing/Logger.h>
#ifdef PARALLEL
#include "compositing/DIYCompositor.hpp"
#endif

#include <assert.h>

namespace vtkh {

Renderer::Renderer()
  : m_do_composite(true),
    m_color_table("cool2warm"),
    m_field_index(0)
{
  m_compositor  = NULL; 
#ifdef PARALLEL
  m_compositor  = new DIYCompositor(); 
#else
  m_compositor  = new Compositor(); 
#endif

}

Renderer::~Renderer()
{
  delete m_compositor;
}

void 
Renderer::SetField(const std::string field_name)
{
  m_field_name = field_name; 
}

void 
Renderer::SetDoComposite(bool do_composite)
{
  m_do_composite = do_composite;
}

void
Renderer::AddRender(vtkh::Render &render)
{
  m_renders.push_back(render); 
}

std::vector<vtkh::Render>
Renderer::GetRenders()
{
  return m_renders; 
}

void
Renderer::SetRenders(const std::vector<vtkh::Render> &renders)
{
  m_renders = renders; 
}

int
Renderer::GetNumberOfRenders() const
{
  return static_cast<int>(m_renders.size());
}

void
Renderer::ClearRenders()
{
  m_renders.clear(); 
}

void Renderer::SetColorTable(const vtkm::rendering::ColorTable &color_table)
{
  m_color_table = color_table;
}

vtkm::rendering::ColorTable Renderer::GetColorTable() const
{
  return m_color_table;
}

void 
Renderer::Composite(const int &num_images)
{

  m_compositor->SetCompositeMode(Compositor::Z_BUFFER_SURFACE);
  for(int i = 0; i < num_images; ++i)
  {
    const int num_canvases = m_renders[i].GetNumberOfCanvases();

    for(int dom = 0; dom < num_canvases; ++dom)
    {
      float* color_buffer = &GetVTKMPointer(m_renders[i].GetCanvas(dom)->GetColorBuffer())[0][0]; 
      float* depth_buffer = GetVTKMPointer(m_renders[i].GetCanvas(dom)->GetDepthBuffer()); 

      int height = m_renders[i].GetCanvas(dom)->GetHeight();
      int width = m_renders[i].GetCanvas(dom)->GetWidth();

      m_compositor->AddImage(color_buffer,
                             depth_buffer,
                             width,
                             height);
    } //for dom

    Image result = m_compositor->Composite();

    float bg_color[4];
    vtkm::rendering::Color color = m_renders[i].GetCanvas(0)->GetBackgroundColor();
    bg_color[0] = color.Components[0];
    bg_color[1] = color.Components[1];
    bg_color[2] = color.Components[2];
    bg_color[3] = color.Components[3];
    result.CompositeBackground(bg_color);
#ifdef PARALLEL
    if(vtkh::GetMPIRank() == 0)
    {
      ImageToCanvas(result, *m_renders[i].GetCanvas(0), true); 
    }
#else
    ImageToCanvas(result, *m_renders[i].GetCanvas(0), true); 
#endif
    m_compositor->ClearImages();
  } // for image
}

void 
Renderer::PreExecute() 
{
  vtkm::cont::ArrayHandle<vtkm::Range> ranges = m_input->GetGlobalRange(m_field_name);
  int num_components = ranges.GetPortalControl().GetNumberOfValues();
  //
  // current vtkm renderers only supports single component scalar fields
  //
  assert(num_components == 1);
  m_range = ranges.GetPortalControl().Get(0);
  m_bounds = m_input->GetGlobalBounds();
  int total_renders = static_cast<int>(m_renders.size());
  for(int i = 0; i < total_renders; ++i)
  {
    m_renders[i].SetScalarRange(m_range);
    m_renders[i].RenderWorldAnnotations();
  }

}

void 
Renderer::Update() 
{
  PreExecute();
  DoExecute();
  PostExecute();
}

void 
Renderer::PostExecute() 
{
  int total_renders = static_cast<int>(m_renders.size());
  if(m_do_composite)
  {
    this->Composite(total_renders);
  }
}

void 
Renderer::DoExecute() 
{
  if(m_mapper.get() == 0)
  {
    std::string msg = "Renderer Error: no renderer was set by sub-class"; 
    throw Error(msg);
  }

  int total_renders = static_cast<int>(m_renders.size());
  int num_domains = static_cast<int>(m_input->GetNumberOfDomains());
  for(int dom = 0; dom < num_domains; ++dom)
  {
    vtkm::cont::DataSet data_set; 
    vtkm::Id domain_id;
    m_input->GetDomain(dom, data_set, domain_id);
    const vtkm::cont::DynamicCellSet &cellset = data_set.GetCellSet();
    const vtkm::cont::Field &field = data_set.GetField(m_field_name);
    const vtkm::cont::CoordinateSystem &coords = data_set.GetCoordinateSystem();
    if(cellset.GetNumberOfCells() == 0) continue;

    for(int i = 0; i < total_renders; ++i)
    {
      // paint
      if(m_renders[i].HasColorTable())
      {
        m_mapper->SetActiveColorTable(m_renders[i].GetColorTable());
      }
      else
      {
        m_mapper->SetActiveColorTable(m_color_table);
      }
      
      vtkmCanvasPtr p_canvas = m_renders[i].GetDomainCanvas(domain_id);
      const vtkmCamera &camera = m_renders[i].GetCamera(); 
      m_mapper->SetCanvas(&(*p_canvas));
      m_mapper->RenderCells(cellset,
                            coords,
                            field,
                            m_color_table,
                            camera,
                            m_range);

    }
  }


}

void 
Renderer::ImageToCanvas(Image &image, vtkm::rendering::Canvas &canvas, bool get_depth) 
{
  const int width = canvas.GetWidth(); 
  const int height = canvas.GetHeight(); 
  const int size = width * height;
  const int color_size = size * 4;
  float* color_buffer = &GetVTKMPointer(canvas.GetColorBuffer())[0][0]; 
  float one_over_255 = 1.f / 255.f;
#ifdef VTKH_USE_OPENMP
  #pragma omp parallel for 
#endif
  for(int i = 0; i < color_size; ++i)
  {
    color_buffer[i] = static_cast<float>(image.m_pixels[i]) * one_over_255;
  }

  float* depth_buffer = GetVTKMPointer(canvas.GetDepthBuffer()); 
  if(get_depth) memcpy(depth_buffer, &image.m_depths[0], sizeof(float) * size);
}

std::vector<Render> 
Renderer::GetRenders() const
{
  return m_renders;
}

vtkh::DataSet *
Renderer::GetInput() 
{
  return m_input;
}

} // namespace vtkh
