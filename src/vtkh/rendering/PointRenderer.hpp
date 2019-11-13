#ifndef VTK_H_RENDERER_POINTS_HPP
#define VTK_H_RENDERER_POINTS_HPP

#include <vtkh/rendering/Renderer.hpp>

namespace vtkh {

#include<vtkh/vtkh_exports.h>
class VTKH_API PointRenderer : public Renderer
{
public:
  PointRenderer();
  virtual ~PointRenderer();
  std::string GetName() const override;
  static Renderer::vtkmCanvasPtr GetNewCanvas(int width = 1024, int height = 1024);
  void PreExecute() override;
  void UseCells();
  void UseNodes();
  void UseVariableRadius(bool useVariableRadius);
  void SetBaseRadius(vtkm::Float32 radius);
  void SetRadiusDelta(vtkm::Float32 delta);
private:
  bool m_use_nodes;
  bool m_radius_set;
  bool m_use_variable_radius;
  vtkm::Float32 m_base_radius;
  vtkm::Float32 m_delta_radius;

};

} // namespace vtkh
#endif
