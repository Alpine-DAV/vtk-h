#include <vtkh/rendering/Scene.hpp>
#include <vtkh/rendering/MeshRenderer.hpp>
#include <vtkh/rendering/VolumeRenderer.hpp>


namespace vtkh
{

Scene::Scene()
  : m_has_volume(false),
    m_batch_size(10)
{

}

Scene::~Scene()
{

}

void
Scene::SetRenderBatchSize(int batch_size)
{
  assert(batch_size > 0);
  m_batch_size = batch_size;
}

int
Scene::GetRenderBatchSize() const
{
  return m_batch_size;
}

void
Scene::AddRender(vtkh::Render &render)
{
  m_renders.push_back(render);
}

void
Scene::SetRenders(const std::vector<vtkh::Render> &renders)
{
  m_renders = renders;
}

bool
Scene::IsMesh(vtkh::Renderer *renderer)
{
  bool is_mesh = false;

  if(dynamic_cast<vtkh::MeshRenderer*>(renderer) != nullptr)
  {
    is_mesh = true;
  }
  return is_mesh;
}

bool
Scene::IsVolume(vtkh::Renderer *renderer)
{
  bool is_volume = false;

  if(dynamic_cast<vtkh::VolumeRenderer*>(renderer) != nullptr)
  {
    is_volume = true;
  }
  return is_volume;
}

void
Scene::AddRenderer(vtkh::Renderer *renderer)
{
  bool is_volume = IsVolume(renderer);
  bool is_mesh = IsMesh(renderer);

  if(is_volume)
  {
    if(m_has_volume)
    {
      throw Error("Scenes only support a single volume plot");
    }

    m_has_volume = true;
    // make sure that the volume render is last
    m_renderers.push_back(renderer);
  }
  else if(is_mesh)
  {
    // make sure that the mesh plot is last
    // and before the volume pl0t
    if(m_has_volume)
    {
      if(m_renderers.size() == 1)
      {
        m_renderers.push_front(renderer);
      }
      else
      {
        auto it = m_renderers.end();
        it--;
        it--;
        m_renderers.insert(it,renderer);
      }
    }
    else
    {
      m_renderers.push_back(renderer);
    }
  }
  else
  {
    m_renderers.push_front(renderer);
  }
}

void
Scene::Render()
{

  std::vector<vtkm::Range> ranges;
  std::vector<std::string> field_names;
  std::vector<vtkm::cont::ColorTable> color_tables;

  bool do_once = true;

  //
  // We are going to render images in batches. With databases
  // like Cinema, we could be rendering hundres of images. Keeping
  // all the canvases around can hog memory so we will conserve it.
  // For example, if we rendered 360 images at 1024^2, all the canvases
  // would consume 7GB of space. Not good on the GPU, where resources
  // are limited.
  //
  const int render_size = m_renders.size();
  int batch_start = 0;
  while(batch_start < render_size)
  {
    int batch_end = std::min(m_batch_size + batch_start, render_size);
    auto begin = m_renders.begin() + batch_start;
    auto end = m_renders.begin() + batch_end;

    std::vector<vtkh::Render> current_batch(begin, end);
    const int plot_size = m_renderers.size();
    auto renderer = m_renderers.begin();

    for(int i = 0; i < plot_size; ++i)
    {
      if(i == plot_size - 1)
      {
        (*renderer)->SetDoComposite(true);
      }
      else
      {
        (*renderer)->SetDoComposite(false);
      }

      (*renderer)->SetRenders(current_batch);
      (*renderer)->Update();

      // we only need to get the ranges and color tables once
      if(do_once)
      {
        if((*renderer)->GetHasColorTable())
        {
          ranges.push_back((*renderer)->GetRange());
          field_names.push_back((*renderer)->GetFieldName());
          color_tables.push_back((*renderer)->GetColorTable());
        }
        do_once = false;
      }

      current_batch  = (*renderer)->GetRenders();
      (*renderer)->ClearRenders();

      renderer++;
    }

    // render screen annotations last and save
    for(int i = 0; i < current_batch.size(); ++i)
    {
      current_batch[i].RenderWorldAnnotations();
      current_batch[i].RenderScreenAnnotations(field_names, ranges, color_tables);
      current_batch[i].RenderBackground();
      current_batch[i].Save();
      // free buffers
      m_renders[batch_start + i].ClearCanvases();
    }

    batch_start = batch_end;
  } // while
}

void
Scene::Save()
{

}

} // namespace vtkh
