#include "ScalarRenderer.hpp"
#include <vtkh/compositing/PayloadCompositor.hpp>

#include <vtkh/Logger.hpp>
#include <vtkh/utils/vtkm_array_utils.hpp>
#include <vtkh/utils/vtkm_dataset_info.hpp>
#include <vtkh/utils/PNGEncoder.hpp>
#include <vtkm/rendering/raytracing/Logger.h>
#include <vtkm/rendering/ScalarRenderer.h>
#include <vtkm/io/writer/VTKDataSetWriter.h>

#include <assert.h>

namespace vtkh {

ScalarRenderer::ScalarRenderer()
  : m_width(1024),
    m_height(1024)
{
}

ScalarRenderer::~ScalarRenderer()
{
}

std::string
ScalarRenderer::GetName() const
{
  return "vtkh::ScalarRenderer";
}

void
ScalarRenderer::AddCamera(vtkmCamera &camera)
{
  m_cameras.push_back(camera);
}

void
ScalarRenderer::SetCameras(const std::vector<vtkmCamera> &cameras)
{
  m_cameras = cameras;
}

void
ScalarRenderer::ClearCameras()
{
  m_cameras.clear();
}

void
ScalarRenderer::Composite(const int &num_images)
{
//  VTKH_DATA_OPEN("Composite");
//  for(int i = 0; i < num_images; ++i)
//  {
//    const int num_canvases = m_renders[i].GetNumberOfCanvases();
//
//    for(int dom = 0; dom < num_canvases; ++dom)
//    {
//      float* color_buffer = &GetVTKMPointer(m_renders[i].GetCanvas(dom)->GetColorBuffer())[0][0];
//      float* depth_buffer = GetVTKMPointer(m_renders[i].GetCanvas(dom)->GetDepthBuffer());
//
//      int height = m_renders[i].GetCanvas(dom)->GetHeight();
//      int width = m_renders[i].GetCanvas(dom)->GetWidth();
//
//      m_compositor->AddImage(color_buffer,
//                             depth_buffer,
//                             width,
//                             height);
//    } //for dom
//
//    Image result = m_compositor->Composite();
//
//#ifdef VTKH_PARALLEL
//    if(vtkh::GetMPIRank() == 0)
//    {
//      ImageToCanvas(result, *m_renders[i].GetCanvas(0), true);
//    }
//#else
//    ImageToCanvas(result, *m_renders[i].GetCanvas(0), true);
//#endif
//    m_compositor->ClearImages();
//  } // for image
//  VTKH_DATA_CLOSE();
}

void
ScalarRenderer::PreExecute()
{
}

void
ScalarRenderer::Update()
{
  VTKH_DATA_OPEN(this->GetName());
#ifdef VTKH_ENABLE_LOGGING
  long long int in_cells = this->m_input->GetNumberOfCells();
  VTKH_DATA_ADD("input_cells", in_cells);
#endif
  PreExecute();
  DoExecute();
  PostExecute();
  VTKH_DATA_CLOSE();
}

void
ScalarRenderer::PostExecute()
{
  int total_cameras = static_cast<int>(m_cameras.size());
  this->Composite(total_cameras);
}

void
ScalarRenderer::DoExecute()
{

  int total_cameras = static_cast<int>(m_cameras.size());
  int num_domains = static_cast<int>(m_input->GetNumberOfDomains());

  //
  // There external faces + bvh construction happens
  // when we set the input for the renderer, which
  // we don't want to repeat for every camera. Also,
  // We could be processing AMR patches, numbering
  // in the 1000s, and with 100 images * 1000s amr
  // patches we could blow memory. We will set the input
  // once and composite after every image (todo: batch images
  // in groups of X).
  //
  std::vector<vtkm::rendering::ScalarRenderer> renderers;
  std::vector<vtkm::Id> cell_counts;
  renderers.resize(num_domains);
  cell_counts.resize(num_domains);
  for(int dom = 0; dom < num_domains; ++dom)
  {
    vtkm::cont::DataSet data_set;
    vtkm::Id domain_id;
    m_input->GetDomain(dom, data_set, domain_id);
    renderers[dom].SetInput(data_set);
    renderers[dom].SetWidth(m_width);
    renderers[dom].SetHeight(m_height);

    // all the data sets better be the same
    cell_counts.push_back(data_set.GetCellSet().GetNumberOfCells());
  }

  // basic sanity checking
  int min_p = std::numeric_limits<int>::max();
  int max_p = std::numeric_limits<int>::min();
  bool do_once = true;

  for(int i = 0; i < total_cameras; ++i)
  {
    const vtkmCamera &camera = m_cameras[i];
    std::vector<std::string> field_names;
    PayloadCompositor compositor;

    for(int dom = 0; dom < num_domains; ++dom)
    {
      Result res = renderers[dom].Render(camera);

      //vtkm::cont::DataSet dset = res.ToDataSet();
      //std::stringstream ss;
      //ss<<"part"<<dom<<".vtk";
      //vtkm::io::writer::VTKDataSetWriter writer(ss.str());
      //writer.WriteDataSet(dset);

      field_names = res.ScalarNames;
      PayloadImage *pimage = Convert(res);
      min_p = std::min(min_p, pimage->m_payload_bytes);
      max_p = std::max(max_p, pimage->m_payload_bytes);
      compositor.AddImage(*pimage);
      delete pimage;
    }

    if(min_p != max_p)
    {
      std::cout<<"VERY BAD\n";
    }

    PayloadImage final_image = compositor.Composite();
    if(vtkh::GetMPIRank() == 0)
    {
      Result final_result = Convert(final_image, field_names);
      vtkm::cont::DataSet dset = final_result.ToDataSet();
      vtkm::io::writer::VTKDataSetWriter writer("test.vtk");
      writer.WriteDataSet(dset);

    }
    // WE BETTER BE SURE THAT payload sizes are equal and in the same
    // order !!!
    if(do_once)
    {
      //if(min_p
    }
  }


}

ScalarRenderer::Result ScalarRenderer::Convert(PayloadImage &image, std::vector<std::string> &names)
{
  Result result;
  result.ScalarNames = names;
  const int num_fields = names.size();

  const int dx  = image.m_bounds.X.Max - image.m_bounds.X.Min + 1;
  const int dy  = image.m_bounds.Y.Max - image.m_bounds.Y.Min + 1;
  const int size = dx * dy;

  result.Width = dx;
  result.Height = dy;

  std::vector<float*> buffers;
  for(int i = 0; i < num_fields; ++i)
  {
    vtkm::cont::ArrayHandle<vtkm::Float32> array;
    array.Allocate(size);
    result.Scalars.push_back(array);
    float* buffer = GetVTKMPointer(result.Scalars[i]);
    buffers.push_back(buffer);
  }

  const unsigned char *loads = &image.m_payloads[0];
  const size_t payload_size = image.m_payload_bytes;

  for(size_t x = 0; x < size; ++x)
  {
    for(int i = 0; i < num_fields; ++i)
    {
      const size_t offset = x * payload_size + i * sizeof(float);
      memcpy(&buffers[i][x], loads + offset, sizeof(float));
    }
  }

  //
  result.Depths.Allocate(size);
  float* dbuffer = GetVTKMPointer(result.Depths);
  memcpy(dbuffer, &image.m_depths[0], sizeof(float) * size);

  return result;
}

PayloadImage * ScalarRenderer::Convert(Result &result)
{
  const int num_fields = result.Scalars.size();
  const int payload_size = num_fields * sizeof(float);
  vtkm::Bounds bounds;
  bounds.X.Min = 1;
  bounds.Y.Min = 1;
  bounds.X.Max = result.Width;
  bounds.Y.Max = result.Height;
  const size_t size = result.Width * result.Height;

  PayloadImage *image = new PayloadImage(bounds, payload_size);
  unsigned char *loads = &image->m_payloads[0];

  float* dbuffer = GetVTKMPointer(result.Depths);
  std::cout<<"Size "<<size<<"\n";
  std::cout<<"DSize "<<result.Depths.GetNumberOfValues()<<"\n";
  std::cout<<"mdSize "<<image->m_depths.size()<<"\n";
  memcpy(&image->m_depths[0], dbuffer, sizeof(float) * size);
  // copy scalars into payload
  std::vector<float*> buffers;
  for(int i = 0; i < num_fields; ++i)
  {
    vtkm::cont::ArrayHandle<vtkm::Float32> scalar = result.Scalars[i];
    float* buffer = GetVTKMPointer(scalar);
    buffers.push_back(buffer);
  }
#ifdef VTKH_USE_OPENMP
    #pragma omp parallel for
#endif
  for(size_t x = 0; x < size; ++x)
  {
    for(int i = 0; i < num_fields; ++i)
    {
      const size_t offset = x * payload_size + i * sizeof(float);
      memcpy(loads + offset, &buffers[i][x], sizeof(float));
    }
  }
  return image;
}

//void
//ScalarRenderer::ImageToCanvas(Image &image, vtkm::rendering::Canvas &canvas, bool get_depth)
//{
//  const int width = canvas.GetWidth();
//  const int height = canvas.GetHeight();
//  const int size = width * height;
//  const int color_size = size * 4;
//  float* color_buffer = &GetVTKMPointer(canvas.GetColorBuffer())[0][0];
//  float one_over_255 = 1.f / 255.f;
//#ifdef VTKH_USE_OPENMP
//  #pragma omp parallel for
//#endif
//  for(int i = 0; i < color_size; ++i)
//  {
//    color_buffer[i] = static_cast<float>(image.m_pixels[i]) * one_over_255;
//  }
//
//  float* depth_buffer = GetVTKMPointer(canvas.GetDepthBuffer());
//  if(get_depth) memcpy(depth_buffer, &image.m_depths[0], sizeof(float) * size);
//}

vtkh::DataSet *
ScalarRenderer::GetInput()
{
  return m_input;
}

} // namespace vtkh
