#include <vtkh/rendering/Scene.hpp>
#include <vtkh/rendering/MeshRenderer.hpp>
#include <vtkh/rendering/VolumeRenderer.hpp>


namespace vtkh 
{

Scene::Scene()
  : m_has_volume(false)
{

}

Scene::~Scene()
{

}

void 
Scene::AddRender(vtkh::Render &render)
{
  m_renders.push_back(render);
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
  const int render_size = m_renders.size();

  std::vector<vtkm::Range> ranges; 
  std::vector<std::string> field_names; 
  std::vector<vtkm::rendering::ColorTable> color_tables; 

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

    (*renderer)->SetRenders(m_renders);
    (*renderer)->Update();

    if((*renderer)->GetHasColorTable())
    {
      std::cout<<"has color table\n";
      ranges.push_back((*renderer)->GetRange());
      field_names.push_back((*renderer)->GetFieldName());
      color_tables.push_back((*renderer)->GetColorTable());
    }
    else std::cout<<"NO has color table\n";

    m_renders = (*renderer)->GetRenders();
    renderer++;
  }
  
  // render screen annotations last and save
  for(int i = 0; i < render_size; ++i)
  {
    m_renders[i].RenderWorldAnnotations();
    m_renders[i].RenderScreenAnnotations(field_names, ranges, color_tables);
    m_renders[i].Save();
  }
}

void 
Scene::Save()
{

}

} // namespace vtkh
