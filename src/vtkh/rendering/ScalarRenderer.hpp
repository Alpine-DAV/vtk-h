#ifndef VTK_H_SCALAR_RENDERER_HPP
#define VTK_H_SCALAR_RENDERER_HPP

#include <vector>
#include <vtkh/vtkh_exports.h>
#include <vtkh/Error.hpp>
#include <vtkh/filters/Filter.hpp>
#include <vtkh/rendering/Render.hpp>
#include <vtkh/compositing/PayloadImage.hpp>

#include <vtkm/rendering/Camera.h>
#include <vtkm/rendering/ScalarRenderer.h>

namespace vtkh {

class VTKH_API ScalarRenderer : public Filter
{
public:
  typedef vtkm::rendering::Camera vtkmCamera;
  using Result = vtkm::rendering::ScalarRenderer::Result;

  ScalarRenderer();
  virtual ~ScalarRenderer();
  virtual void Update();
  virtual std::string GetName() const override;

  void AddCamera(vtkmCamera &camera);
  void ClearCameras();

  void SetCameras(const std::vector<vtkmCamera> &cameras);

  int GetNumberOfCameras() const;
  vtkh::DataSet *GetInput();
protected:

  int m_width;
  int m_height;
  // image related data with cinema support
  std::vector<vtkmCamera>    m_cameras;
  // methods
  virtual void PreExecute() override;
  virtual void PostExecute() override;
  virtual void DoExecute() override;

  virtual void Composite(const int &num_images);
  PayloadImage * Convert(Result &result);
  ScalarRenderer::Result Convert(PayloadImage &image, std::vector<std::string> &names);
  //void ImageToDataSet(Image &image, vtkm::rendering::Canvas &canvas, bool get_depth);
};

} // namespace vtkh
#endif
