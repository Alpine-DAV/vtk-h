#ifndef VTK_H_RENDER_HPP
#define VTK_H_RENDER_HPP

#include <vector>
#include <vtkh/DataSet.hpp>
#include <vtkh/Error.hpp>

#include <vtkm/rendering/Camera.h>
#include <vtkm/rendering/CanvasRayTracer.h>
#include <vtkm/rendering/Mapper.h>

namespace vtkh {
//
// A Render contains the information needed to create a single image.
// There are 'n' canvases that matches the number of domains in the 
// data set. It is possible to chain multiple plots together that 
// are rendering separate data, i.e. the result of different data
// transformations, to handle this we keep track of the domain ids
// that each canvas is associated with.
//

class Render
{
public: 
  typedef std::shared_ptr<vtkm::rendering::CanvasRayTracer> vtkmCanvasPtr; 

  Render();
  ~Render();
  vtkmCanvasPtr                   GetDomainCanvas(const vtkm::Id &domain_id);
  vtkmCanvasPtr                   GetCanvas(const vtkm::Id index);
  int                             GetNumberOfCanvases() const;
  const vtkm::rendering::Camera&  GetCamera() const;
  std::string                     GetImageName() const;
  vtkm::Bounds                    GetSceneBounds() const;
  vtkm::Int32                     GetHeight() const;
  vtkm::Int32                     GetWidth() const;
  vtkm::rendering::Color          GetBackgroundColor() const;

  void                            SetWidth(const vtkm::Int32 width);
  void                            SetHeight(const vtkm::Int32 height);
  void                            SetSceneBounds(const vtkm::Bounds &bounds);
  void                            SetCamera(const vtkm::rendering::Camera &camera);
  void                            SetImageName(const std::string &name);
  void                            SetBackgroundColor(float bg_color[4]);
  void                            ClearCanvases();
  bool                            HasCanvas(const vtkm::Id &domain_id) const;
  void                            AddDomain(vtkm::Id domain_id);
  void                            RenderWorldAnnotations();
  void                            RenderScreenAnnotations(const std::vector<std::string> &field_names,
                                                          const std::vector<vtkm::Range> &ranges,
                                                          const std::vector<vtkm::rendering::ColorTable> &colors);
  void                            Save();
protected:
  std::vector<vtkmCanvasPtr>   m_canvases;
  std::vector<vtkm::Id>        m_domain_ids;
  vtkm::rendering::Camera      m_camera; 
  std::string                  m_image_name;
  vtkm::Bounds                 m_scene_bounds;
  vtkm::Int32                  m_width;
  vtkm::Int32                  m_height;
  vtkm::rendering::Color       m_bg_color;
  vtkmCanvasPtr                CreateCanvas();
}; 

static float vtkh_default_bg_color[4] = {0.f, 0.f, 0.f, 1.f};

//template<typename RendererType>
vtkh::Render 
MakeRender(int width,
           int height, 
           vtkm::Bounds scene_bounds,
           const std::vector<vtkm::Id> &domain_ids,
           const std::string &image_name,
           float bg_color[4] = vtkh_default_bg_color);

vtkh::Render 
MakeRender(int width,
           int height, 
           vtkm::rendering::Camera camera,
           vtkh::DataSet &data_set,
           const std::string &image_name,
           float bg_color[4] = vtkh_default_bg_color);

} // namespace vtkh
#endif
