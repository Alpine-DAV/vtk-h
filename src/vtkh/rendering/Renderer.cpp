#include "Renderer.hpp"
#include "compositing/Compositor.hpp"

#include <vtkh/Logger.hpp>
#include <vtkh/utils/vtkm_array_utils.hpp>
#include <vtkh/utils/vtkm_dataset_info.hpp>
#include <vtkh/utils/PNGEncoder.hpp>
#include <vtkm/rendering/raytracing/Logger.h>

#include <assert.h>
#include <chrono>

namespace vtkh {

Renderer::Renderer()
  : m_do_composite(false),
    m_color_table("Cool to Warm"),
    m_field_index(0),
    m_has_color_table(true)
{
  m_compositor  = new Compositor();
}

Renderer::~Renderer()
{
  delete m_compositor;
}

void
Renderer::SetShadingOn(bool on)

{
  // do nothing by default;
}

void
Renderer::SetField(const std::string field_name)
{
  m_field_name = field_name;
}

std::string
Renderer::GetFieldName() const
{
  return m_field_name;
}

bool
Renderer::GetHasColorTable() const
{
  return m_has_color_table;
}

double 
Renderer::GetLastRenderTime() const
{
  return m_render_times.back();
}

int 
Renderer::GetMpiRank() const
{
  return vtkh::GetMPIRank();
}


std::vector<double>
Renderer::GetRenderTimes() 
{
  return m_render_times;
}

std::vector<std::vector<unsigned char> >
Renderer::GetColorBuffers() 
{
  return m_color_buffers;
}

std::vector<std::vector<float> >
Renderer::GetDepthBuffers() 
{
  return m_depth_buffers;
}

std::vector<float>
Renderer::GetDepths() 
{
  return m_depths;
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

void Renderer::SetColorTable(const vtkm::cont::ColorTable &color_table)
{
  m_color_table = color_table;
}

vtkm::cont::ColorTable Renderer::GetColorTable() const
{
  return m_color_table;
}

bool Renderer::HasContribution(const vtkm::Range &plot_scalar_range,
                               const vtkm::cont::DataSet &dom,
                               const vtkm::Float64 threshold)
{
  int num_alphas = m_color_table.GetNumberOfPointsAlpha();
  // i don't know what the defualt behavior is for only one
  // alpha point, so if we have 0 or 1 just say that this
  // domain contributes
  if(num_alphas < 2)
  {
    return true;
  }

  if(!dom.HasField(m_field_name))
  {
    return false;
  }

  const vtkm::cont::Field &field = dom.GetField(m_field_name);
  vtkm::cont::ArrayHandle<vtkm::Range> sub_range;
  sub_range = field.GetRange();

  vtkm::Range field_range = sub_range.GetPortalControl().Get(0);

  vtkm::Float64 min_value = plot_scalar_range.Min;
  vtkm::Float64 max_value = plot_scalar_range.Max;
  vtkm::Float64 length = min_value == max_value ? 1.0 : max_value - min_value;

  vtkm::Float64 domain_min = field_range.Min;
  vtkm::Float64 domain_max = field_range.Max;
  // normalize to color table positions
  domain_min = (domain_min - min_value) / length;
  domain_max = (domain_max - min_value) / length;

  vtkm::Float64 max_alpha = 0;

  for(int i = 0; i < num_alphas-1; ++i)
  {
    vtkm::Vec<vtkm::Float64,4> point0, point1;
    bool valid = m_color_table.GetPointAlpha(i, point0);
    valid = m_color_table.GetPointAlpha(i+1, point1);
    // alpha points are location, alpha, mid point, sharpness
    if(point0[0] <= domain_max && domain_min <= point1[0])
    {
      max_alpha = std::max(max_alpha, point0[1]);
      max_alpha = std::max(max_alpha, point1[1]);
    }
  }

  return max_alpha > threshold;

}


void
Renderer::Composite(const int &num_images)
{
  VTKH_DATA_OPEN("Composite");
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

#ifdef VTKH_PARALLEL
    if(vtkh::GetMPIRank() == 0)
    {
      ImageToCanvas(result, *m_renders[i].GetCanvas(0), true);
    }
#else
    ImageToCanvas(result, *m_renders[i].GetCanvas(0), true);
#endif
    m_compositor->ClearImages();
  } // for image
  VTKH_DATA_CLOSE();
}

void
Renderer::PreExecute()
{
  bool range_set = m_range.IsNonEmpty();
  // TODO: work around the global call
  // bool field_exists = m_input->GlobalFieldExists(m_field_name);
  bool field_exists = m_input->FieldExists(m_field_name); 
  if(!range_set && field_exists)
  {
    // we have not been given a range, so ask the data set
    // TODO: work around the global call -> validate
    // vtkm::cont::ArrayHandle<vtkm::Range> ranges = m_input->GetGlobalRange(m_field_name);
    vtkm::cont::ArrayHandle<vtkm::Range> ranges = m_input->GetRange(m_field_name);
    int num_components = ranges.GetPortalControl().GetNumberOfValues();
    //
    // current vtkm renderers only supports single component scalar fields
    //
    assert(num_components == 1);
    if(num_components != 1)
    {
      std::stringstream msg;
      msg<<"Renderer '"<<this->GetName()<<"' cannot render a field with ";
      msg<<"'"<<num_components<<"' components. Field must be a scalar field.";
      throw Error(msg.str());
    }

    vtkm::Range global_range = ranges.GetPortalControl().Get(0);
    // a min or max may be been set by the user, check to see
    if(m_range.Min == vtkm::Infinity64())
    {
      m_range.Min = global_range.Min;
    }

    if(m_range.Max == vtkm::NegativeInfinity64())
    {
      m_range.Max = global_range.Max;
    }
  }

  // TODO: work around the global call -> validate
  m_bounds = m_input->GetBounds(); // m_input->GetGlobalBounds();

  // std::cout << "*** bounds: " << m_bounds.X.Min << "  " 
  //           << m_bounds.X.Max << " / " 
  //           << m_bounds.Y.Min << " " 
  //           << m_bounds.Y.Max << " / "
  //           << m_bounds.Z.Min << " " 
  //           << m_bounds.Z.Max << std::endl;

  m_bounds.X.Min =  0;
  m_bounds.X.Max = 10;
  m_bounds.Y.Min =  0;
  m_bounds.Y.Max = 10; 
  m_bounds.Z.Min =  0;
  m_bounds.Z.Max = 10;
}

void
Renderer::Update()
{
  VTKH_DATA_OPEN(this->GetName());
#ifdef VTKH_ENABLE_LOGGING
  long long int in_cells = this->m_input->GetNumberOfCells();
  VTKH_DATA_ADD("input_cells", in_cells);
#endif
  // PreExecute();
  DoExecute();
  PostExecute();
  VTKH_DATA_CLOSE();
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
  int skipped_dom = 0;
  int total_renders = static_cast<int>(m_renders.size());
  int num_domains = static_cast<int>(m_input->GetNumberOfDomains());
  for(int dom = 0; dom < num_domains; ++dom)
  {
    vtkm::cont::DataSet data_set;
    vtkm::Id domain_id;
    m_input->GetDomain(dom, data_set, domain_id);

    if(!data_set.HasField(m_field_name))
    {
      continue;
    }

    // TODO: m_range ->  m_scalar_range ; threshold == 0.01 ?
    if (!HasContribution(m_range, data_set, vtkm::Float64(0.01)))
    {
      int rank = 0;
#ifdef VTKH_PARALLEL
      rank = vtkh::GetMPIRank();
#endif
      std::cout << "---Skip block on rank " << rank << std::endl;
      ++skipped_dom;
      AddRenderTime(0.0);
      continue;
    }

    const vtkm::cont::DynamicCellSet &cellset = data_set.GetCellSet();
    const vtkm::cont::Field &field = data_set.GetField(m_field_name);
    const vtkm::cont::CoordinateSystem &coords = data_set.GetCoordinateSystem();
    if(cellset.GetNumberOfCells() == 0) continue;

    log_global_time("begin rendering", vtkh::GetMPIRank());
    for(int i = 0; i < total_renders; ++i)
    {
      if(m_renders[i].GetShadingOn())
      {
        this->SetShadingOn(true);
      }
      else
      {
        this->SetShadingOn(false);
      }

      m_mapper->SetActiveColorTable(m_color_table);

      auto t1 = std::chrono::high_resolution_clock::now();
      
      vtkmCanvasPtr p_canvas = m_renders[i].GetDomainCanvas(domain_id);
      const vtkmCamera &camera = m_renders[i].GetCamera();
      m_mapper->SetCanvas(&(*p_canvas));
      m_mapper->RenderCells(cellset,
                            coords,
                            field,
                            m_color_table,
                            camera,
                            m_range);
      
      auto t2 = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1).count();
      AddRenderTime(duration);
    }
    log_global_time("end rendering", vtkh::GetMPIRank());
  }
  
  if (skipped_dom == num_domains)
    m_skipped = true;
  else
    m_skipped = false;
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

void Renderer::AddRenderTime(double t)
{
  m_render_times.push_back(t);
  if (m_render_times.size() > 10000) // only keep the last 10k frame times for now
    m_render_times.erase(m_render_times.begin());
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

vtkm::Range
Renderer::GetRange() const
{
  return m_range;
}

void
Renderer::SetRange(const vtkm::Range &range)
{
  m_range = range;
}

} // namespace vtkh
