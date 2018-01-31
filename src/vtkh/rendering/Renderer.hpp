#ifndef VTK_H_RENDERER_HPP
#define VTK_H_RENDERER_HPP

#include <vector>
#include <vtkh/Error.hpp>
#include <vtkh/filters/Filter.hpp>
#include <vtkh/rendering/Render.hpp>
#include <vtkh/rendering/Image.hpp>

#include <vtkm/rendering/Camera.h>
#include <vtkm/rendering/Canvas.h>
#include <vtkm/rendering/Mapper.h>

namespace vtkh {

class Compositor;

class Renderer : public Filter
{
public:
  typedef std::shared_ptr<vtkm::rendering::Canvas> vtkmCanvasPtr; 
  typedef std::shared_ptr<vtkm::rendering::Mapper> vtkmMapperPtr; 
  typedef vtkm::rendering::Camera vtkmCamera; 

  Renderer(); 
  virtual ~Renderer(); 
  virtual void Update();

  void AddRender(vtkh::Render &render); 
  void ClearRenders(); 

  void SetField(const std::string field_name);
  void SetColorTable(const vtkm::rendering::ColorTable &color_table);
  void SetDoComposite(bool do_composite);
  void SetRenders(const std::vector<Render> &renders);
  void SetRange(const vtkm::Range &range);

  vtkm::rendering::ColorTable GetColorTable() const;
  std::string                 GetFieldName() const;
  int                         GetNumberOfRenders() const; 
  std::vector<Render>         GetRenders() const;
  vtkh::DataSet              *GetInput();
  vtkm::Range                 GetRange() const;
  bool                        GetHasColorTable() const;
protected:
  
  // image related data with cinema support
  std::vector<vtkh::Render>                m_renders;
  int                                      m_field_index;
  Compositor                              *m_compositor;
  std::string                              m_field_name;
  bool                                     m_do_composite;   
  vtkmMapperPtr                            m_mapper;
  vtkm::Bounds                             m_bounds;
  vtkm::Range                              m_range;
  vtkm::rendering::ColorTable              m_color_table;
  vtkm::rendering::ColorTable              m_corrected_color_table;
  bool                                     m_has_color_table;  
  // methods
  virtual void PreExecute() override;
  virtual void PostExecute() override;
  virtual void DoExecute() override;

  virtual void Composite(const int &num_images);
  void ImageToCanvas(Image &image, vtkm::rendering::Canvas &canvas, bool get_depth);
};

} // namespace vtkh
#endif
