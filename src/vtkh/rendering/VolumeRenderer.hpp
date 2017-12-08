#ifndef VTK_H_RENDERER_VOLUME_HPP
#define VTK_H_RENDERER_VOLUME_HPP

#include <vtkh/rendering/Renderer.hpp>
#include <vtkm/rendering/MapperVolume.h>

namespace vtkh {

class VolumeRenderer : public Renderer
{
public:
  VolumeRenderer();
  virtual ~VolumeRenderer();
  void SetNumberOfSamples(const int num_samples);
  static Renderer::vtkmCanvasPtr GetNewCanvas(int width = 1024, int height = 1024);

  void Update() override;
protected:  
  virtual void Composite(const int &num_images) override;
  virtual void PreExecute() override;
  virtual void PostExecute() override;

  std::vector<std::vector<int>> m_visibility_orders;
  void FindVisibilityOrdering();
  void DepthSort(int num_domains, 
                 std::vector<float> &min_depths,
                 std::vector<int> &local_vis_order);
  float FindMinDepth(const vtkm::rendering::Camera &camera, 
                     const vtkm::Bounds &bounds) const;
  
  std::shared_ptr<vtkm::rendering::MapperVolume> m_tracer;
  int m_num_samples;
};

} // namespace vtkh
#endif
