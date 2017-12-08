#include <vtkh/rendering/Scene.hpp>
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

void 
Scene::AddRenderer(vtkh::Renderer *renderer)
{
  bool is_volume = false;

  if(dynamic_cast<vtkh::VolumeRenderer*>(renderer) != nullptr)
  {
    is_volume = true;
  }

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
  else
  {
    m_renderers.push_front(renderer);
  }
}

void 
Scene::Render()
{
  const int render_size = m_renders.size();

  // Always render world annotations first
  for(int i = 0; i < render_size; ++i)
  {
    m_renders[i].RenderWorldAnnotations();
  }
  
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

    m_renders = (*renderer)->GetRenders();
    renderer++;
  }
  
  // render screen annotations last and save
  for(int i = 0; i < render_size; ++i)
  {
    m_renders[i].RenderScreenAnnotations();
    m_renders[i].Save();
  }
}

void 
Scene::Save()
{

}

} // namespace vtkh
