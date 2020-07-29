#ifndef VTK_H_RENDERER_HPP
#define VTK_H_RENDERER_HPP

#include <vector>
#include <memory>
#include <chrono>
#include <iomanip>

#include <vtkh/vtkh_exports.h>
#include <vtkh/Error.hpp>
#include <vtkh/filters/Filter.hpp>
#include <vtkh/rendering/Render.hpp>
#include <vtkh/rendering/Image.hpp>

#include <vtkm/rendering/Camera.h>
#include <vtkm/rendering/Canvas.h>
#include <vtkm/rendering/Mapper.h>

namespace vtkh {

class Compositor;

class VTKH_API Renderer : public Filter
{
public:
  typedef std::shared_ptr<vtkm::rendering::Canvas> vtkmCanvasPtr;
  typedef std::shared_ptr<vtkm::rendering::Mapper> vtkmMapperPtr;
  typedef vtkm::rendering::Camera vtkmCamera;

  Renderer();
  virtual ~Renderer();
  virtual void SetShadingOn(bool on);
  virtual void Update();

  void AddRender(vtkh::Render &render);
  void ClearRenders();

  void SetField(const std::string field_name);
  virtual void SetColorTable(const vtkm::cont::ColorTable &color_table);
  void SetDoComposite(bool do_composite);
  void SetRenders(const std::vector<Render> &renders);
  void SetRange(const vtkm::Range &range);

  bool HasContribution(const vtkm::Range &plot_scalar_range,
                       const vtkm::cont::DataSet &dom,
                       const vtkm::Float64 threshold);

  vtkm::cont::ColorTable      GetColorTable() const;
  std::string                 GetFieldName() const;
  int                         GetNumberOfRenders() const;
  std::vector<Render>         GetRenders() const;
  vtkh::DataSet              *GetInput();
  vtkm::Range                 GetRange() const;
  bool                        GetHasColorTable() const;
  double                      GetLastRenderTime() const;
  std::vector<double>         GetRenderTimes();
  int                         GetMpiRank() const;
  std::vector<std::vector<unsigned char> > GetColorBuffers();
  std::vector<std::vector<float> >         GetDepthBuffers();
  std::vector<float>          GetDepths();

static std::string get_timing_file_name(const int value, const int precision, 
                                        const std::string &prefix,
                                        const std::string &path = "timings")
{
    std::ostringstream oss;
    oss << path;
    oss << "/";
    oss << prefix;
    oss << "_";
    oss << std::setw(precision) << std::setfill('0') << value;
    oss << ".txt";
    return oss.str();
}

static void log_global_time(const std::string &description,
                            const int rank)
{
    auto now = std::chrono::system_clock::now();
    auto duration = now.time_since_epoch();
    auto millis = std::chrono::duration_cast<std::chrono::milliseconds>(duration).count();

    std::ofstream out(get_timing_file_name(rank, 5, "global"), std::ios_base::app);
    out << description << " : " << millis << std::endl;
    out.close();
}

protected:

  // image related data with cinema support
  std::vector<vtkh::Render>                m_renders;
  int                                      m_field_index;
  Compositor                              *m_compositor;
  std::string                              m_field_name;
  bool                                     m_do_composite = false;
  vtkmMapperPtr                            m_mapper;
  vtkm::Bounds                             m_bounds;
  vtkm::Range                              m_range;
  vtkm::cont::ColorTable                   m_color_table;
  bool                                     m_has_color_table;
  std::vector<double>                      m_render_times;  // render time in milliseconds
  std::vector<std::vector<unsigned char>>  m_color_buffers;
  std::vector<std::vector<float> >         m_depth_buffers;
  std::vector<float>                       m_depths;
  bool                                     m_skipped = false;

  // methods
  virtual void PreExecute() override;
  virtual void PostExecute() override;
  virtual void DoExecute() override;

  virtual void Composite(const int &num_images);
  void ImageToCanvas(Image &image, vtkm::rendering::Canvas &canvas, bool get_depth);
  void AddRenderTime(double t);
};

} // namespace vtkh
#endif
