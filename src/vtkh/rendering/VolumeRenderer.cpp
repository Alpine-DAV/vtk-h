#include "VolumeRenderer.hpp"

#include <vtkh/utils/vtkm_array_utils.hpp>
#include <vtkh/compositing/Compositor.hpp>
#include <vtkh/Logger.hpp>

#include <vtkm/rendering/CanvasRayTracer.h>

#include <memory>

#ifdef VTKH_PARALLEL
#include <mpi.h>
#endif


#define VTKH_OPACITY_CORRECTION 10.f

namespace vtkh {

namespace detail
{
  struct VisOrdering
  {
    int m_rank;
    int m_domain_index;
    int m_order;
    float m_minz;
  };

  struct DepthOrder
  {
    inline bool operator()(const VisOrdering &lhs, const VisOrdering &rhs)
    {
      return lhs.m_minz < rhs.m_minz;
    }
  };

  struct RankOrder
  {
    inline bool operator()(const VisOrdering &lhs, const VisOrdering &rhs)
    {
      if(lhs.m_rank < rhs.m_rank)
      {
        return true;
      }
      else if(lhs.m_rank == rhs.m_rank)
      {
        return lhs.m_domain_index < rhs.m_domain_index;
      }
      return false;
    }
  };
} //  namespace detail

VolumeRenderer::VolumeRenderer()
{
  typedef vtkm::rendering::MapperVolume TracerType;
  m_tracer = std::make_shared<TracerType>();
  this->m_mapper = m_tracer;
  m_tracer->SetCompositeBackground(false);
  //
  // add some default opacity to the color table
  //
  m_uncorrected_color_table.AddPointAlpha(0.0f, .02);
  m_uncorrected_color_table.AddPointAlpha(.0f, .5);
  m_num_samples = 100.f;
  CorrectOpacity();
}

VolumeRenderer::~VolumeRenderer()
{
}

void
VolumeRenderer::Update()
{
  VTKH_DATA_OPEN(this->GetName());
#ifdef VTKH_ENABLE_LOGGING
  VTKH_DATA_ADD("device", GetCurrentDevice());
  long long int in_cells = this->m_input->GetNumberOfCells();
  VTKH_DATA_ADD("input_cells", in_cells);
  VTKH_DATA_ADD("input_domains", this->m_input->GetNumberOfDomains());
  int in_topo_dims;
  bool in_structured = this->m_input->IsStructured(in_topo_dims);
  if(in_structured)
  {
    VTKH_DATA_ADD("in_topology", "structured");
  }
  else
  {
    VTKH_DATA_ADD("in_topology", "unstructured");
  }
#endif

  PreExecute();
  DoExecute();
  PostExecute();

  VTKH_DATA_CLOSE();
}

void VolumeRenderer::SetColorTable(const vtkm::cont::ColorTable &color_table)
{
  m_uncorrected_color_table = color_table;
  CorrectOpacity();
}

void VolumeRenderer::CorrectOpacity()
{
  const float correction_scalar = VTKH_OPACITY_CORRECTION;
  float samples = m_num_samples;

  float ratio = correction_scalar / samples;
  vtkm::cont::ColorTable corrected;
  corrected = m_uncorrected_color_table;
  int num_points = corrected.GetNumberOfPointsAlpha();
  for(int i = 0; i < num_points; i++)
  {
    vtkm::Vec<vtkm::Float64,4> point;
    corrected.GetPointAlpha(i,point);
    point[1] = 1. - vtkm::Pow((1. - point[1]), double(ratio));
    corrected.UpdatePointAlpha(i,point);
  }

  this->m_color_table = corrected;
}

void
VolumeRenderer::DoExecute()
{
  std::cout<<"DO EXECUTE\n";
  if(m_input->OneDomainPerRank())
  {
    RenderOneDomainPerRank();
  }
  else
  {
    RenderMultipleDomainsPerRank();
  }
}

void
VolumeRenderer::RenderOneDomainPerRank()
{
  std::cout<<"Doing it\n";
  if(m_mapper.get() == 0)
  {
    std::string msg = "Renderer Error: no renderer was set by sub-class";
    throw Error(msg);
  }

  int total_renders = static_cast<int>(m_renders.size());
  int num_domains = static_cast<int>(m_input->GetNumberOfDomains());
  if(num_domains > 1)
  {
    std::cout<<"VERY BAD\n";
  }
  for(int dom = 0; dom < num_domains; ++dom)
  {
    std::cout<<"DOM "<<dom<<"\n";
    vtkm::cont::DataSet data_set;
    vtkm::Id domain_id;
    m_input->GetDomain(0, data_set, domain_id);

    if(!data_set.HasField(m_field_name))
    {
      continue;
    }

    const vtkm::cont::DynamicCellSet &cellset = data_set.GetCellSet();
    const vtkm::cont::Field &field = data_set.GetField(m_field_name);
    const vtkm::cont::CoordinateSystem &coords = data_set.GetCoordinateSystem();

    if(cellset.GetNumberOfCells() == 0) continue;

    for(int i = 0; i < total_renders; ++i)
    {
      m_mapper->SetActiveColorTable(m_color_table);

      Render::vtkmCanvas &canvas = m_renders[i].GetCanvas();
      const vtkmCamera &camera = m_renders[i].GetCamera();
      m_mapper->SetCanvas(&canvas);
      m_mapper->RenderCells(cellset,
                            coords,
                            field,
                            m_color_table,
                            camera,
                            m_range);
      canvas.SaveAs("bananas.pnm");
    }
  }
}

void
VolumeRenderer::RenderMultipleDomainsPerRank()
{
  std::cout<<"NOT IMPLEMENTED\n";
}

void
VolumeRenderer::PreExecute()
{
  Renderer::PreExecute();

  vtkm::Vec<vtkm::Float32,3> extent;
  extent[0] = static_cast<vtkm::Float32>(this->m_bounds.X.Length());
  extent[1] = static_cast<vtkm::Float32>(this->m_bounds.Y.Length());
  extent[2] = static_cast<vtkm::Float32>(this->m_bounds.Z.Length());
  vtkm::Float32 dist = vtkm::Magnitude(extent) / m_num_samples;
  m_tracer->SetSampleDistance(dist);
}

void
VolumeRenderer::PostExecute()
{
  int total_renders = static_cast<int>(m_renders.size());
  if(m_do_composite)
  {
    this->Composite(total_renders);
  }
}

void
VolumeRenderer::SetNumberOfSamples(const int num_samples)
{
  if(num_samples < 1)
  {
    throw Error("Volume rendering samples must be greater than 0");
  }
  m_num_samples = num_samples;
  CorrectOpacity();
}

Renderer::vtkmCanvasPtr
VolumeRenderer::GetNewCanvas(int width, int height)
{
  return std::make_shared<vtkm::rendering::CanvasRayTracer>(width, height);
}

float
VolumeRenderer::FindMinDepth(const vtkm::rendering::Camera &camera,
                                 const vtkm::Bounds &bounds) const
{

  vtkm::Vec<vtkm::Float64,3> center = bounds.Center();
  vtkm::Vec<vtkm::Float64,3> fcenter;
  fcenter[0] = static_cast<vtkm::Float32>(center[0]);
  fcenter[1] = static_cast<vtkm::Float32>(center[1]);
  fcenter[2] = static_cast<vtkm::Float32>(center[2]);
  vtkm::Vec<vtkm::Float32,3> pos = camera.GetPosition();
  vtkm::Float32 dist = vtkm::Magnitude(fcenter - pos);
  return dist;
}

void
VolumeRenderer::Composite(const int &num_images)
{
  std::cout<<"Compositing\n";
  const int num_domains = static_cast<int>(m_input->GetNumberOfDomains());

  m_compositor->SetCompositeMode(Compositor::VIS_ORDER_BLEND);

  FindVisibilityOrdering();

  for(int i = 0; i < num_images; ++i)
  {
    float* color_buffer =
      &GetVTKMPointer(m_renders[i].GetCanvas().GetColorBuffer())[0][0];
    float* depth_buffer =
      GetVTKMPointer(m_renders[i].GetCanvas().GetDepthBuffer());
    int height = m_renders[i].GetCanvas().GetHeight();
    int width = m_renders[i].GetCanvas().GetWidth();

    m_compositor->AddImage(color_buffer,
                           depth_buffer,
                           width,
                           height,
                           m_visibility_orders[i][0]);

    Image result = m_compositor->Composite();
    const std::string image_name = m_renders[i].GetImageName() + ".png";
#ifdef VTKH_PARALLEL
    if(vtkh::GetMPIRank() == 0)
    {
#endif
      ImageToCanvas(result, m_renders[i].GetCanvas(), true);
#ifdef VTKH_PARALLEL
    }
#endif
    m_compositor->ClearImages();
  } // for image
}

void
VolumeRenderer::DepthSort(int num_domains,
                          std::vector<float> &min_depths,
                          std::vector<int> &local_vis_order)
{
  if(min_depths.size() != num_domains)
  {
    throw Error("min depths size does not equal the number of domains");
  }
  if(local_vis_order.size() != num_domains)
  {
    throw Error("local vis order not equal to number of domains");
  }
#ifdef VTKH_PARALLEL
  int root = 0;
  MPI_Comm comm = MPI_Comm_f2c(vtkh::GetMPICommHandle());
  int num_ranks = vtkh::GetMPISize();
  int rank = vtkh::GetMPIRank();
  int *domain_counts = NULL;
  int *domain_offsets = NULL;
  int *vis_order = NULL;
  float *depths = NULL;

  if(rank == root)
  {
    domain_counts = new int[num_ranks];
    domain_offsets = new int[num_ranks];
  }

  MPI_Gather(&num_domains,
             1,
             MPI_INT,
             domain_counts,
             1,
             MPI_INT,
            root,
             comm);

  int depths_size = 0;
  if(rank == root)
  {
    //scan for dispacements
    domain_offsets[0] = 0;
    for(int i = 1; i < num_ranks; ++i)
    {
      domain_offsets[i] = domain_offsets[i - 1] + domain_counts[i - 1];
    }

    for(int i = 0; i < num_ranks; ++i)
    {
      depths_size += domain_counts[i];
    }

    depths = new float[depths_size];

  }

  MPI_Gatherv(&min_depths[0],
              num_domains,
              MPI_FLOAT,
              depths,
              domain_counts,
              domain_offsets,
              MPI_FLOAT,
              root,
              comm);

  if(rank == root)
  {
    std::vector<detail::VisOrdering> order;
    order.resize(depths_size);

    for(int i = 0; i < num_ranks; ++i)
    {
      for(int c = 0; c < domain_counts[i]; ++c)
      {
        int index = domain_offsets[i] + c;
        order[index].m_rank = i;
        order[index].m_domain_index = c;
        order[index].m_minz = depths[index];
      }
    }

    std::sort(order.begin(), order.end(), detail::DepthOrder());

    for(int i = 0; i < depths_size; ++i)
    {
      order[i].m_order = i;
    }

    std::sort(order.begin(), order.end(), detail::RankOrder());

    vis_order = new int[depths_size];
    for(int i = 0; i < depths_size; ++i)
    {
      vis_order[i] = order[i].m_order;
    }
  }

  MPI_Scatterv(vis_order,
               domain_counts,
               domain_offsets,
               MPI_INT,
               &local_vis_order[0],
               num_domains,
               MPI_INT,
               root,
               comm);

  if(rank == root)
  {
    delete[] domain_counts;
    delete[] domain_offsets;
    delete[] vis_order;
    delete[] depths;
  }
#else

  std::vector<detail::VisOrdering> order;
  order.resize(num_domains);

  for(int i = 0; i < num_domains; ++i)
  {
      order[i].m_rank = 0;
      order[i].m_domain_index = i;
      order[i].m_minz = min_depths[i];
  }
  std::sort(order.begin(), order.end(), detail::DepthOrder());

  for(int i = 0; i < num_domains; ++i)
  {
    order[i].m_order = i;
  }

  std::sort(order.begin(), order.end(), detail::RankOrder());

  for(int i = 0; i < num_domains; ++i)
  {
    local_vis_order[i] = order[i].m_order;
  }
#endif
}

void
VolumeRenderer::FindVisibilityOrdering()
{
  const int num_domains = static_cast<int>(m_input->GetNumberOfDomains());
  const int num_cameras = static_cast<int>(m_renders.size());
  m_visibility_orders.resize(num_cameras);

  for(int i = 0; i < num_cameras; ++i)
  {
    m_visibility_orders[i].resize(num_domains);
  }

  //
  // In order for parallel volume rendering to composite correctly,
  // we nee to establish a visibility ordering to pass to IceT.
  // We will transform the data extents into camera space and
  // take the minimum z value. Then sort them while keeping
  // track of rank, then pass the list in.
  //
  std::vector<float> min_depths;
  min_depths.resize(num_domains);

  for(int i = 0; i < num_cameras; ++i)
  {
    const vtkm::rendering::Camera &camera = m_renders[i].GetCamera();
    for(int dom = 0; dom < num_domains; ++dom)
    {
      vtkm::Bounds bounds = this->m_input->GetDomainBounds(dom);
      min_depths[dom] = FindMinDepth(camera, bounds);
    }

    DepthSort(num_domains, min_depths, m_visibility_orders[i]);

  } // for each camera
}

std::string
VolumeRenderer::GetName() const
{
  return "vtkh::VolumeRenderer";
}

} // namespace vtkh
