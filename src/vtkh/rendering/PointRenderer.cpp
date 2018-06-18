#include "PointRenderer.hpp"

#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/MapperPoint.h>
#include <memory>

namespace vtkh {
  
PointRenderer::PointRenderer()
  : m_use_nodes(true),
    m_radius_set(false),
    m_use_variable_radius(false),
    m_base_radius(0.5f),
    m_delta_radius(0.5f)
{
  typedef vtkm::rendering::MapperPoint TracerType;
  auto mapper = std::make_shared<TracerType>();
  mapper->SetCompositeBackground(false);
  this->m_mapper = mapper;
}

PointRenderer::~PointRenderer()
{
}

Renderer::vtkmCanvasPtr 
PointRenderer::GetNewCanvas(int width, int height)
{
  return std::make_shared<vtkm::rendering::CanvasRayTracer>(width, height);
}

std::string
PointRenderer::GetName() const
{
  return "vtkh::PointRenderer";
}

void 
PointRenderer::UseCells()
{
  m_use_nodes = false;
}

void 
PointRenderer::UseNodes()
{
  m_use_nodes = true;
}

void 
PointRenderer::UseVariableRadius(bool useVariableRadius)
{
  m_use_variable_radius = useVariableRadius;
}

void 
PointRenderer::SetBaseRadius(vtkm::Float32 radius)
{
  m_base_radius = radius;
  m_radius_set = true;
}

void 
PointRenderer::SetRadiusDelta(vtkm::Float32 delta)
{
  m_delta_radius = delta;
}

void
PointRenderer::PreExecute()
{
  Renderer::PreExecute();

  typedef vtkm::rendering::MapperPoint MapperType;
  std::shared_ptr<MapperType> mesh_mapper = 
    std::dynamic_pointer_cast<MapperType>(this->m_mapper);
  
  if(m_use_nodes)
  {
    mesh_mapper->UseNodes();
  }
  else
  {
    mesh_mapper->UseCells();
  }

  // allow for athe default mapper radius
  if(m_radius_set)
  {
    mesh_mapper->SetRadius(m_base_radius);
  }

  mesh_mapper->UseVariableRadius(m_use_variable_radius);
  mesh_mapper->SetRadiusDelta(m_delta_radius);

}

} // namespace vtkh
